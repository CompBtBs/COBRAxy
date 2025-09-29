"""
Utilities for generating and manipulating COBRA models and related metadata.

This module includes helpers to:
- extract rules, reactions, bounds, objective coefficients, and compartments
- build a COBRA model from a tabular file
- set objective and medium from dataframes
- validate a model and convert gene identifiers
- translate model GPRs using mapping tables
"""
import os
import cobra
import pandas as pd
import re
import logging
from typing import Optional, Tuple, Union, List, Dict, Set
from collections import defaultdict
import utils.rule_parsing  as rulesUtils
import utils.reaction_parsing as reactionUtils
from cobra import Model as cobraModel, Reaction, Metabolite
import sys


############################ check_methods ####################################
def gene_type(l :str, name :str) -> str:
    """
    Determine the type of gene ID.

    Args:
        l (str): The gene identifier to check.
        name (str): The name of the dataset, used in error messages.

    Returns:
        str: The type of gene ID ('hugo_id', 'ensembl_gene_id', 'symbol', or 'entrez_id').

    Raises:
        sys.exit: If the gene ID type is not supported, the execution is aborted.
    """
    if check_hgnc(l):
        return 'hugo_id'
    elif check_ensembl(l):
        return 'ensembl_gene_id'
    elif check_symbol(l):
        return 'symbol'
    elif check_entrez(l):
        return 'entrez_id'
    else:
        sys.exit('Execution aborted:\n' +
                 'gene ID type in ' + name + ' not supported. Supported ID'+
                 'types are: HUGO ID, Ensemble ID, HUGO symbol, Entrez ID\n')

def check_hgnc(l :str) -> bool:
    """
    Check if a gene identifier follows the HGNC format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the HGNC format, False otherwise.
    """
    if len(l) > 5:
        if (l.upper()).startswith('HGNC:'):
            return l[5:].isdigit()
        else:
            return False
    else:
        return False

def check_ensembl(l :str) -> bool:
    """
    Check if a gene identifier follows the Ensembl format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the Ensembl format, False otherwise.
    """
    return l.upper().startswith('ENS')
 

def check_symbol(l :str) -> bool:
    """
    Check if a gene identifier follows the symbol format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the symbol format, False otherwise.
    """
    if len(l) > 0:
        if l[0].isalpha() and l[1:].isalnum():
            return True
        else:
            return False
    else:
        return False

def check_entrez(l :str) -> bool:
    """
    Check if a gene identifier follows the Entrez ID format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the Entrez ID format, False otherwise.
    """ 
    if len(l) > 0:
        return l.isdigit()
    else: 
        return False

################################- DATA GENERATION -################################
ReactionId = str
def generate_rules(model: cobraModel, *, asParsed = True) -> Union[Dict[ReactionId, rulesUtils.OpList], Dict[ReactionId, str]]:
    """
    Generate a dictionary mapping reaction IDs to GPR rules from the model.

    Args:
        model: COBRA model to derive data from.
        asParsed: If True, parse rules into a nested list structure; otherwise keep raw strings.

    Returns:
        Dict[ReactionId, rulesUtils.OpList]: Parsed rules by reaction ID.
        Dict[ReactionId, str]: Raw rules by reaction ID.
    """
    _ruleGetter   =  lambda reaction : reaction.gene_reaction_rule
    ruleExtractor = (lambda reaction :
        rulesUtils.parseRuleToNestedList(_ruleGetter(reaction))) if asParsed else _ruleGetter

    return {
        reaction.id : ruleExtractor(reaction)
        for reaction in model.reactions
        if reaction.gene_reaction_rule }

def generate_reactions(model :cobraModel, *, asParsed = True) -> Dict[ReactionId, str]:
    """
    Generate a dictionary mapping reaction IDs to reaction formulas from the model.

    Args:
        model: COBRA model to derive data from.
        asParsed: If True, convert formulas into a parsed representation; otherwise keep raw strings.

    Returns:
        Dict[ReactionId, str]: Reactions by reaction ID (parsed if requested).
    """

    unparsedReactions = {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction 
    }

    if not asParsed: return unparsedReactions
    
    return reactionUtils.create_reaction_dict(unparsedReactions)

def get_medium(model:cobraModel) -> pd.DataFrame:
    """
    Extract the uptake reactions representing the model medium.

    Returns a DataFrame with a single column 'reaction' listing exchange reactions
    with negative lower bound and no positive stoichiometric coefficients (uptake only).
    """
    trueMedium=[]
    for r in model.reactions:
        positiveCoeff=0
        for m in r.metabolites:
            if r.get_coefficient(m.id)>0:
                positiveCoeff=1;
        if (positiveCoeff==0 and r.lower_bound<0):
            trueMedium.append(r.id)

    df_medium = pd.DataFrame()
    df_medium["reaction"] = trueMedium
    return df_medium

def extract_objective_coefficients(model: cobraModel) -> pd.DataFrame:
    """
    Extract objective coefficients for each reaction.

    Args:
        model: COBRA model

    Returns:
        pd.DataFrame with columns: ReactionID, ObjectiveCoefficient
    """
    coeffs = []
    # model.objective.expression is a linear expression
    objective_expr = model.objective.expression.as_coefficients_dict()
    
    for reaction in model.reactions:
        coeff = objective_expr.get(reaction.forward_variable, 0.0)
        coeffs.append({
            "ReactionID": reaction.id,
            "ObjectiveCoefficient": coeff
        })
    
    return pd.DataFrame(coeffs)

def generate_bounds(model:cobraModel) -> pd.DataFrame:
    """
    Build a DataFrame of lower/upper bounds for all reactions.

    Returns:
        pd.DataFrame indexed by reaction IDs with columns ['lower_bound', 'upper_bound'].
    """

    rxns = []
    for reaction in model.reactions:
        rxns.append(reaction.id)

    bounds = pd.DataFrame(columns = ["lower_bound", "upper_bound"], index=rxns)

    for reaction in model.reactions:
        bounds.loc[reaction.id] = [reaction.lower_bound, reaction.upper_bound]
    return bounds



def generate_compartments(model: cobraModel) -> pd.DataFrame:
    """
    Generates a DataFrame containing compartment information for each reaction.
    Creates columns for each compartment position (Compartment_1, Compartment_2, etc.)
    
    Args:
        model: the COBRA model to extract compartment data from.
        
    Returns:
        pd.DataFrame: DataFrame with ReactionID and compartment columns
    """
    pathway_data = []

    # First pass: determine the maximum number of pathways any reaction has
    max_pathways = 0
    reaction_pathways = {}

    for reaction in model.reactions:
        # Get unique pathways from all metabolites in the reaction
        if type(reaction.annotation['pathways']) == list:
            reaction_pathways[reaction.id] = reaction.annotation['pathways']
            max_pathways = max(max_pathways, len(reaction.annotation['pathways']))
        else:
            reaction_pathways[reaction.id] = [reaction.annotation['pathways']]

    # Create column names for pathways
    pathway_columns = [f"Pathway_{i+1}" for i in range(max_pathways)]

    # Second pass: create the data
    for reaction_id, pathways in reaction_pathways.items():
        row = {"ReactionID": reaction_id}
        
        # Fill pathway columns
        for i in range(max_pathways):
            col_name = pathway_columns[i]
            if i < len(pathways):
                row[col_name] = pathways[i]
            else:
                row[col_name] = None  # or "" if you prefer empty strings

        pathway_data.append(row)

    return pd.DataFrame(pathway_data)



def build_cobra_model_from_csv(csv_path: str, model_id: str = "new_model") -> cobraModel:
    """
    Build a COBRApy model from a tabular file with reaction data.

    Args:
        csv_path: Path to the tab-separated file.
        model_id: ID for the newly created model.

    Returns:
        cobra.Model: The constructed COBRApy model.
    """
    
    df = pd.read_csv(csv_path, sep='\t')
    
    model = cobraModel(model_id)
    
    metabolites_dict = {}
    compartments_dict = {}
    
    print(f"Building model from {len(df)} reactions...")
    
    for idx, row in df.iterrows():
        reaction_formula = str(row['Formula']).strip()
        if not reaction_formula or reaction_formula == 'nan':
            continue
            
        metabolites = extract_metabolites_from_reaction(reaction_formula)
        
        for met_id in metabolites:
            compartment = extract_compartment_from_metabolite(met_id)
            
            if compartment not in compartments_dict:
                compartments_dict[compartment] = compartment
            
            if met_id not in metabolites_dict:
                metabolites_dict[met_id] = Metabolite(
                    id=met_id,
                    compartment=compartment,
                    name=met_id.replace(f"_{compartment}", "").replace("__", "_")
                )
    
    model.compartments = compartments_dict
    
    model.add_metabolites(list(metabolites_dict.values()))
    
    print(f"Added {len(metabolites_dict)} metabolites and {len(compartments_dict)} compartments")
    
    reactions_added = 0
    reactions_skipped = 0
    
    for idx, row in df.iterrows():

        reaction_id = str(row['ReactionID']).strip()
        reaction_formula = str(row['Formula']).strip()
        
        if not reaction_formula or reaction_formula == 'nan':
            raise ValueError(f"Missing reaction formula for {reaction_id}")
        
        reaction = Reaction(reaction_id)
        reaction.name = reaction_id
        
        reaction.lower_bound = float(row['lower_bound']) if pd.notna(row['lower_bound']) else -1000.0
        reaction.upper_bound = float(row['upper_bound']) if pd.notna(row['upper_bound']) else 1000.0
        
        if pd.notna(row['GPR']) and str(row['GPR']).strip():
            reaction.gene_reaction_rule = str(row['GPR']).strip()
        
        try:
            parse_reaction_formula(reaction, reaction_formula, metabolites_dict)
        except Exception as e:
            print(f"Error parsing reaction {reaction_id}: {e}")
            reactions_skipped += 1
            continue
        
        model.add_reactions([reaction])
        reactions_added += 1
            
    
    print(f"Added {reactions_added} reactions, skipped {reactions_skipped} reactions")
    
    # set objective function
    set_objective_from_csv(model, df, obj_col="ObjectiveCoefficient")

    set_medium_from_data(model, df)
    
    print(f"Model completed: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    return model


# Estrae tutti gli ID metaboliti nella formula (gestisce prefissi numerici + underscore)
def extract_metabolites_from_reaction(reaction_formula: str) -> Set[str]:
    """
    Extract metabolite IDs from a reaction formula.
    Robust pattern: tokens ending with _<compartment> (e.g., _c, _m, _e),
    allowing leading digits/underscores.
    """
    metabolites = set()
    # optional coefficient followed by a token ending with _<letters>
    pattern = r'(?:\d+(?:\.\d+)?\s+)?([A-Za-z0-9_]+_[a-z]+)'
    matches = re.findall(pattern, reaction_formula)
    metabolites.update(matches)
    return metabolites


def extract_compartment_from_metabolite(metabolite_id: str) -> str:
    """Extract the compartment from a metabolite ID."""
    if '_' in metabolite_id:
        return metabolite_id.split('_')[-1]
    return 'c'  # default cytoplasm


def parse_reaction_formula(reaction: Reaction, formula: str, metabolites_dict: Dict[str, Metabolite]):
    """Parse a reaction formula and set metabolites with their coefficients."""

    if '<=>' in formula:
        left, right = formula.split('<=>')
        reversible = True
    elif '<--' in formula:
        left, right = formula.split('<--')
        reversible = False
    elif '-->' in formula:
        left, right = formula.split('-->')
        reversible = False
    elif '<-' in formula:
        left, right = formula.split('<-')
        reversible = False
    else:
        raise ValueError(f"Unrecognized reaction format: {formula}")
    
    reactants = parse_metabolites_side(left.strip())
    products = parse_metabolites_side(right.strip())
    
    metabolites_to_add = {}
    
    for met_id, coeff in reactants.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = -coeff
    
    for met_id, coeff in products.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = coeff
    
    reaction.add_metabolites(metabolites_to_add)


def parse_metabolites_side(side_str: str) -> Dict[str, float]:
    """Parse one side of a reaction and extract metabolites with coefficients."""
    metabolites = {}
    if not side_str or side_str.strip() == '':
        return metabolites

    terms = side_str.split('+')
    for term in terms:
        term = term.strip()
        if not term:
            continue

        # optional coefficient + id ending with _<compartment>
        match = re.match(r'(?:(\d+\.?\d*)\s+)?([A-Za-z0-9_]+_[a-z]+)', term)
        if match:
            coeff_str, met_id = match.groups()
            coeff = float(coeff_str) if coeff_str else 1.0
            metabolites[met_id] = coeff

    return metabolites



def set_objective_from_csv(model: cobra.Model, df: pd.DataFrame, obj_col: str = "ObjectiveCoefficient"):
    """
    Sets the model's objective function based on a column of coefficients in the CSV.
    Can be any reaction(s), not necessarily biomass.
    """
    obj_dict = {}
    
    for idx, row in df.iterrows():
        reaction_id = str(row['ReactionID']).strip()
        coeff = float(row[obj_col]) if pd.notna(row[obj_col]) else 0.0
        if coeff != 0:
            if reaction_id in model.reactions:
                obj_dict[model.reactions.get_by_id(reaction_id)] = coeff
            else:
                print(f"Warning: reaction {reaction_id} not found in model, skipping for objective.")

    if not obj_dict:
        raise ValueError("No reactions found with non-zero objective coefficient.")

    model.objective = obj_dict
    print(f"Objective set with {len(obj_dict)} reactions.")




def set_medium_from_data(model: cobraModel, df: pd.DataFrame):
    """Set the medium based on the 'InMedium' column in the dataframe."""
    medium_reactions = df[df['InMedium'] == True]['ReactionID'].tolist()
    
    medium_dict = {}
    for rxn_id in medium_reactions:
        if rxn_id in [r.id for r in model.reactions]:
            reaction = model.reactions.get_by_id(rxn_id)
            if reaction.lower_bound < 0:  # Solo reazioni di uptake
                medium_dict[rxn_id] = abs(reaction.lower_bound)
    
    if medium_dict:
        model.medium = medium_dict
        print(f"Medium set with {len(medium_dict)} components")


def validate_model(model: cobraModel) -> Dict[str, any]:
    """Validate the model and return basic statistics."""
    validation = {
        'num_reactions': len(model.reactions),
        'num_metabolites': len(model.metabolites),
        'num_genes': len(model.genes),
        'num_compartments': len(model.compartments),
        'objective': str(model.objective),
        'medium_size': len(model.medium),
        'reversible_reactions': len([r for r in model.reactions if r.reversibility]),
        'exchange_reactions': len([r for r in model.reactions if r.id.startswith('EX_')]),
    }
    
    try:
        # Growth test
        solution = model.optimize()
        validation['growth_rate'] = solution.objective_value
        validation['status'] = solution.status
    except Exception as e:
        validation['growth_rate'] = None
        validation['status'] = f"Error: {e}"
    
    return validation

def convert_genes(model, annotation):
    """Rename genes using a selected annotation key in gene.notes; returns a model copy."""
    from cobra.manipulation import rename_genes
    model2=model.copy()
    try:
        dict_genes={gene.id:gene.notes[annotation]  for gene in model2.genes}
    except:
        print("No annotation in gene dict!")
        return -1
    rename_genes(model2,dict_genes)

    return model2

# ---------- Utility helpers ----------
def _normalize_colname(col: str) -> str:
    return col.strip().lower().replace(' ', '_')

def _choose_columns(mapping_df: 'pd.DataFrame') -> Dict[str, str]:
    """
    Find useful columns and return a dict {ensg: colname1, hgnc_id: colname2, ...}.
    Raise ValueError if no suitable mapping is found.
    """
    cols = { _normalize_colname(c): c for c in mapping_df.columns }
    chosen = {}
    # candidate names for each category
    candidates = {
        'ensg': ['ensg', 'ensembl_gene_id', 'ensembl'],
        'hgnc_id': ['hgnc_id', 'hgnc', 'hgnc:'],
        'hgnc_symbol': ['hgnc_symbol', 'hgnc symbol', 'symbol'],
        'entrez_id': ['entrez', 'entrez_id', 'entrezgene'],
        'gene_number': ['gene_number']
    }
    for key, names in candidates.items():
        for n in names:
            if n in cols:
                chosen[key] = cols[n]
                break
    return chosen

def _validate_target_uniqueness(mapping_df: 'pd.DataFrame',
                                source_col: str,
                                target_col: str,
                                model_source_genes: Optional[Set[str]] = None,
                                logger: Optional[logging.Logger] = None) -> None:
    """
        Check that, within the filtered mapping_df, each target maps to at most one source.
        Log examples if duplicates are found.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if mapping_df.empty:
        logger.warning("Mapping dataframe is empty for the requested source genes; skipping uniqueness validation.")
        return

    # normalize temporary columns for grouping (without altering the original df)
    tmp = mapping_df[[source_col, target_col]].copy()
    tmp['_src_norm'] = tmp[source_col].astype(str).map(_normalize_gene_id)
    tmp['_tgt_norm'] = tmp[target_col].astype(str).str.strip()

    # optionally filter to the set of model source genes
    if model_source_genes is not None:
        tmp = tmp[tmp['_src_norm'].isin(model_source_genes)]

    if tmp.empty:
        logger.warning("After filtering to model source genes, mapping table is empty — nothing to validate.")
        return

    # build reverse mapping: target -> set(sources)
    grouped = tmp.groupby('_tgt_norm')['_src_norm'].agg(lambda s: set(s.dropna()))
    # find targets with more than one source
    problematic = {t: sorted(list(s)) for t, s in grouped.items() if len(s) > 1}

    if problematic:
    # prepare warning message with examples (limited subset)
        sample_items = list(problematic.items())
        msg_lines = ["Mapping validation failed: some target IDs are associated with multiple source IDs."]
        for tgt, sources in sample_items:
            msg_lines.append(f"  - target '{tgt}' <- sources: {', '.join(sources)}")
        full_msg = "\n".join(msg_lines)
    # log warning
        logger.warning(full_msg)

    # if everything is fine
    logger.info("Mapping validation passed: no target ID is associated with multiple source IDs (within filtered set).")


def _normalize_gene_id(g: str) -> str:
    """Normalize a gene ID for use as a key (removes prefixes like 'HGNC:' and strips)."""
    if g is None:
        return ""
    g = str(g).strip()
    # remove common prefixes
    g = re.sub(r'^(HGNC:)', '', g, flags=re.IGNORECASE)
    g = re.sub(r'^(ENSG:)', '', g, flags=re.IGNORECASE)
    return g

def _simplify_boolean_expression(expr: str) -> str:
    """
    Simplify a boolean expression by removing duplicates and redundancies.
    Handles expressions with 'and' and 'or'.
    """
    if not expr or not expr.strip():
        return expr
    
    # normalize operators
    expr = expr.replace(' AND ', ' and ').replace(' OR ', ' or ')
    
    # recursive helper to process expressions
    def process_expression(s: str) -> str:
        s = s.strip()
        if not s:
            return s
            
    # handle parentheses
        while '(' in s:
            # find the innermost parentheses
            start = -1
            for i, c in enumerate(s):
                if c == '(':
                    start = i
                elif c == ')' and start != -1:
                    # process inner content
                    inner = s[start+1:i]
                    processed_inner = process_expression(inner)
                    s = s[:start] + processed_inner + s[i+1:]
                    break
            else:
                break
        
    # split by 'or' at top level
        or_parts = []
        current_part = ""
        paren_count = 0
        
        tokens = s.split()
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token == 'or' and paren_count == 0:
                if current_part.strip():
                    or_parts.append(current_part.strip())
                current_part = ""
            else:
                if token.count('(') > token.count(')'):
                    paren_count += token.count('(') - token.count(')')
                elif token.count(')') > token.count('('):
                    paren_count -= token.count(')') - token.count('(')
                current_part += token + " "
            i += 1
        
        if current_part.strip():
            or_parts.append(current_part.strip())
        
    # process each OR part
        processed_or_parts = []
        for or_part in or_parts:
            # split by 'and' within each OR part
            and_parts = []
            current_and = ""
            paren_count = 0
            
            and_tokens = or_part.split()
            j = 0
            while j < len(and_tokens):
                token = and_tokens[j]
                if token == 'and' and paren_count == 0:
                    if current_and.strip():
                        and_parts.append(current_and.strip())
                    current_and = ""
                else:
                    if token.count('(') > token.count(')'):
                        paren_count += token.count('(') - token.count(')')
                    elif token.count(')') > token.count('('):
                        paren_count -= token.count(')') - token.count('(')
                    current_and += token + " "
                j += 1
            
            if current_and.strip():
                and_parts.append(current_and.strip())
            
            # deduplicate AND parts
            unique_and_parts = list(dict.fromkeys(and_parts))  # mantiene l'ordine
            
            if len(unique_and_parts) == 1:
                processed_or_parts.append(unique_and_parts[0])
            elif len(unique_and_parts) > 1:
                processed_or_parts.append(" and ".join(unique_and_parts))
        
    # deduplicate OR parts
        unique_or_parts = list(dict.fromkeys(processed_or_parts))
        
        if len(unique_or_parts) == 1:
            return unique_or_parts[0]
        elif len(unique_or_parts) > 1:
            return " or ".join(unique_or_parts)
        else:
            return ""
    
    try:
        return process_expression(expr)
    except Exception:
    # if simplification fails, return the original expression
        return expr

# ---------- Main public function ----------
def translate_model_genes(model: 'cobra.Model',
                         mapping_df: 'pd.DataFrame',
                         target_nomenclature: str,
                         source_nomenclature: str = 'hgnc_id',
                         allow_many_to_one: bool = False,
                         logger: Optional[logging.Logger] = None) -> 'cobra.Model':
    """
    Translate model genes from source_nomenclature to target_nomenclature using a mapping table.
    mapping_df should contain columns enabling mapping (e.g., ensg, hgnc_id, hgnc_symbol, entrez).

    Args:
        model: COBRA model to translate.
        mapping_df: DataFrame containing the mapping information.
        target_nomenclature: Desired target key (e.g., 'hgnc_symbol').
        source_nomenclature: Current source key in the model (default 'hgnc_id').
        allow_many_to_one: If True, allow many-to-one mappings and handle duplicates in GPRs.
        logger: Optional logger.
    """
    if logger is None:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        logger = logging.getLogger(__name__)

    logger.info(f"Translating genes from '{source_nomenclature}' to '{target_nomenclature}'")

    # normalize column names and choose relevant columns
    chosen = _choose_columns(mapping_df)
    if not chosen:
        raise ValueError("Could not detect useful columns in mapping_df. Expected at least one of: ensg, hgnc_id, hgnc_symbol, entrez.")

    # map source/target to actual dataframe column names (allow user-specified source/target keys)
    # normalize input args
    src_key = source_nomenclature.strip().lower()
    tgt_key = target_nomenclature.strip().lower()

    # try to find the actual column names for requested keys
    col_for_src = None
    col_for_tgt = None
    # first, try exact match
    for k, actual in chosen.items():
        if k == src_key:
            col_for_src = actual
        if k == tgt_key:
            col_for_tgt = actual

    # if not found, try mapping common names
    if col_for_src is None:
        possible_src_names = {k: v for k, v in chosen.items()}
        # try to match by contained substring
        for k, actual in possible_src_names.items():
            if src_key in k:
                col_for_src = actual
                break

    if col_for_tgt is None:
        for k, actual in chosen.items():
            if tgt_key in k:
                col_for_tgt = actual
                break

    if col_for_src is None:
        raise ValueError(f"Source column for '{source_nomenclature}' not found in mapping dataframe.")
    if col_for_tgt is None:
        raise ValueError(f"Target column for '{target_nomenclature}' not found in mapping dataframe.")

    model_source_genes = { _normalize_gene_id(g.id) for g in model.genes }
    logger.info(f"Filtering mapping to {len(model_source_genes)} source genes present in model (normalized).")

    tmp_map = mapping_df[[col_for_src, col_for_tgt]].dropna().copy()
    tmp_map[col_for_src + "_norm"] = tmp_map[col_for_src].astype(str).map(_normalize_gene_id)

    filtered_map = tmp_map[tmp_map[col_for_src + "_norm"].isin(model_source_genes)].copy()

    if filtered_map.empty:
        logger.warning("No mapping rows correspond to source genes present in the model after filtering. Proceeding with empty mapping (no translation will occur).")

    if not allow_many_to_one:
        _validate_target_uniqueness(filtered_map, col_for_src, col_for_tgt, model_source_genes=model_source_genes, logger=logger)

    # Crea il mapping
    gene_mapping = _create_gene_mapping(filtered_map, col_for_src, col_for_tgt, logger)

    # copy model
    model_copy = model.copy()

    # statistics
    stats = {'translated': 0, 'one_to_one': 0, 'one_to_many': 0, 'not_found': 0, 'simplified_gprs': 0}
    unmapped = []
    multi = []

    original_genes = {g.id for g in model_copy.genes}
    logger.info(f"Original genes count: {len(original_genes)}")

    # translate GPRs
    for rxn in model_copy.reactions:
        gpr = rxn.gene_reaction_rule
        if gpr and gpr.strip():
            new_gpr = _translate_gpr(gpr, gene_mapping, stats, unmapped, multi, logger)
            if new_gpr != gpr:
                simplified_gpr = _simplify_boolean_expression(new_gpr)
                if simplified_gpr != new_gpr:
                    stats['simplified_gprs'] += 1
                    logger.debug(f"Simplified GPR for {rxn.id}: '{new_gpr}' -> '{simplified_gpr}'")
                rxn.gene_reaction_rule = simplified_gpr
                logger.debug(f"Reaction {rxn.id}: '{gpr}' -> '{simplified_gpr}'")

    # update model genes based on new GPRs
    _update_model_genes(model_copy, logger)

    # final logging
    _log_translation_statistics(stats, unmapped, multi, original_genes, model_copy.genes, logger)

    logger.info("Translation finished")
    return model_copy


# ---------- helper functions ----------
def _create_gene_mapping(mapping_df, source_col: str, target_col: str, logger: logging.Logger) -> Dict[str, List[str]]:
    """
    Build mapping dict: source_id -> list of target_ids
    Normalizes IDs (removes prefixes like 'HGNC:' etc).
    """
    df = mapping_df[[source_col, target_col]].dropna().copy()
    # normalize to string
    df[source_col] = df[source_col].astype(str).map(_normalize_gene_id)
    df[target_col] = df[target_col].astype(str).str.strip()

    df = df.drop_duplicates()

    logger.info(f"Creating mapping from {len(df)} rows")

    mapping = defaultdict(list)
    for _, row in df.iterrows():
        s = row[source_col]
        t = row[target_col]
        if t not in mapping[s]:
            mapping[s].append(t)

    # stats
    one_to_one = sum(1 for v in mapping.values() if len(v) == 1)
    one_to_many = sum(1 for v in mapping.values() if len(v) > 1)
    logger.info(f"Mapping: {len(mapping)} source keys, {one_to_one} 1:1, {one_to_many} 1:many")
    return dict(mapping)


def _translate_gpr(gpr_string: str,
                   gene_mapping: Dict[str, List[str]],
                   stats: Dict[str, int],
                   unmapped_genes: List[str],
                   multi_mapping_genes: List[Tuple[str, List[str]]],
                   logger: logging.Logger) -> str:
    """
    Translate genes inside a GPR string using gene_mapping.
    Returns new GPR string.
    """
    # Generic token pattern: letters, digits, :, _, -, ., (captures HGNC:1234, ENSG000..., symbols)
    token_pattern = r'\b[A-Za-z0-9:_.-]+\b'
    tokens = re.findall(token_pattern, gpr_string)

    logical = {'and', 'or', 'AND', 'OR', '(', ')'}
    tokens = [t for t in tokens if t not in logical]

    new_gpr = gpr_string

    for token in sorted(set(tokens), key=lambda x: -len(x)):  # longer tokens first to avoid partial replacement
        norm = _normalize_gene_id(token)
        if norm in gene_mapping:
            targets = gene_mapping[norm]
            stats['translated'] += 1
            if len(targets) == 1:
                stats['one_to_one'] += 1
                replacement = targets[0]
            else:
                stats['one_to_many'] += 1
                multi_mapping_genes.append((token, targets))
                replacement = "(" + " or ".join(targets) + ")"

            pattern = r'\b' + re.escape(token) + r'\b'
            new_gpr = re.sub(pattern, replacement, new_gpr)
        else:
            stats['not_found'] += 1
            if token not in unmapped_genes:
                unmapped_genes.append(token)
            logger.debug(f"Token not found in mapping (left as-is): {token}")

    return new_gpr


def _update_model_genes(model: 'cobra.Model', logger: logging.Logger):
    """
    Rebuild model.genes from gene_reaction_rule content.
    Removes genes not referenced and adds missing ones.
    """
    # collect genes in GPRs
    gene_pattern = r'\b[A-Za-z0-9:_.-]+\b'
    logical = {'and', 'or', 'AND', 'OR', '(', ')'}
    genes_in_gpr: Set[str] = set()

    for rxn in model.reactions:
        gpr = rxn.gene_reaction_rule
        if gpr and gpr.strip():
            toks = re.findall(gene_pattern, gpr)
            toks = [t for t in toks if t not in logical]
            # normalize IDs consistent with mapping normalization
            toks = [_normalize_gene_id(t) for t in toks]
            genes_in_gpr.update(toks)

    # existing gene ids
    existing = {g.id for g in model.genes}

    # remove obsolete genes
    to_remove = [gid for gid in existing if gid not in genes_in_gpr]
    removed = 0
    for gid in to_remove:
        try:
            gene_obj = model.genes.get_by_id(gid)
            model.genes.remove(gene_obj)
            removed += 1
        except Exception:
            # safe-ignore
            pass

    # add new genes
    added = 0
    for gid in genes_in_gpr:
        if gid not in existing:
            new_gene = cobra.Gene(gid)
            try:
                model.genes.add(new_gene)
            except Exception:
                # fallback: if model.genes doesn't support add, try append or model.add_genes
                try:
                    model.genes.append(new_gene)
                except Exception:
                    try:
                        model.add_genes([new_gene])
                    except Exception:
                        logger.warning(f"Could not add gene object for {gid}")
            added += 1

    logger.info(f"Model genes updated: removed {removed}, added {added}")


def _log_translation_statistics(stats: Dict[str, int],
                               unmapped_genes: List[str],
                               multi_mapping_genes: List[Tuple[str, List[str]]],
                               original_genes: Set[str],
                               final_genes,
                               logger: logging.Logger):
    logger.info("=== TRANSLATION STATISTICS ===")
    logger.info(f"Translated: {stats.get('translated', 0)} (1:1 = {stats.get('one_to_one', 0)}, 1:many = {stats.get('one_to_many', 0)})")
    logger.info(f"Not found tokens: {stats.get('not_found', 0)}")
    logger.info(f"Simplified GPRs: {stats.get('simplified_gprs', 0)}")

    final_ids = {g.id for g in final_genes}
    logger.info(f"Genes in model: {len(original_genes)} -> {len(final_ids)}")

    if unmapped_genes:
        logger.warning(f"Unmapped tokens ({len(unmapped_genes)}): {', '.join(unmapped_genes[:20])}{(' ...' if len(unmapped_genes)>20 else '')}")
    if multi_mapping_genes:
        logger.info(f"Multi-mapping examples ({len(multi_mapping_genes)}):")
        for orig, targets in multi_mapping_genes[:10]:
            logger.info(f"  {orig} -> {', '.join(targets)}")