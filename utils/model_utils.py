import os
import csv
import cobra
import pickle
import argparse
import pandas as pd
import re
import logging
from typing import Optional, Tuple, Union, List, Dict, Set
from collections import defaultdict
import utils.general_utils as utils
import utils.rule_parsing  as rulesUtils
import utils.reaction_parsing as reactionUtils
from cobra import Model as cobraModel, Reaction, Metabolite

################################- DATA GENERATION -################################
ReactionId = str
def generate_rules(model: cobraModel, *, asParsed = True) -> Union[Dict[ReactionId, rulesUtils.OpList], Dict[ReactionId, str]]:
    """
    Generates a dictionary mapping reaction ids to rules from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the rules to an optimized runtime format, otherwise leaves them as strings.

    Returns:
        Dict[ReactionId, rulesUtils.OpList] : the generated dictionary of parsed rules.
        Dict[ReactionId, str] : the generated dictionary of raw rules.
    """
    # Is the below approach convoluted? yes
    # Ok but is it inefficient? probably
    # Ok but at least I don't have to repeat the check at every rule (I'm clinically insane)
    _ruleGetter   =  lambda reaction : reaction.gene_reaction_rule
    ruleExtractor = (lambda reaction :
        rulesUtils.parseRuleToNestedList(_ruleGetter(reaction))) if asParsed else _ruleGetter

    return {
        reaction.id : ruleExtractor(reaction)
        for reaction in model.reactions
        if reaction.gene_reaction_rule }

def generate_reactions(model :cobraModel, *, asParsed = True) -> Dict[ReactionId, str]:
    """
    Generates a dictionary mapping reaction ids to reaction formulas from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the reactions to an optimized runtime format, otherwise leaves them as they are.

    Returns:
        Dict[ReactionId, str] : the generated dictionary.
    """

    unparsedReactions = {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction 
    }

    if not asParsed: return unparsedReactions
    
    return reactionUtils.create_reaction_dict(unparsedReactions)

def get_medium(model:cobraModel) -> pd.DataFrame:
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
    Estrae i coefficienti della funzione obiettivo per ciascuna reazione del modello.
    
    Args:
        model : cobra.Model
        
    Returns:
        pd.DataFrame con colonne: ReactionID, ObjectiveCoefficient
    """
    coeffs = []
    # model.objective.expression è un'espressione lineare
    objective_expr = model.objective.expression.as_coefficients_dict()
    
    for reaction in model.reactions:
        coeff = objective_expr.get(reaction.forward_variable, 0.0)
        coeffs.append({
            "ReactionID": reaction.id,
            "ObjectiveCoefficient": coeff
        })
    
    return pd.DataFrame(coeffs)

def generate_bounds(model:cobraModel) -> pd.DataFrame:

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
    Costruisce un modello COBRApy a partire da un file CSV con i dati delle reazioni.
    
    Args:
        csv_path: Path al file CSV (separato da tab)
        model_id: ID del modello da creare
        
    Returns:
        cobra.Model: Il modello COBRApy costruito
    """
    
    # Leggi i dati dal CSV
    df = pd.read_csv(csv_path, sep='\t')
    
    # Crea il modello vuoto
    model = cobraModel(model_id)
    
    # Dict per tenere traccia di metaboliti e compartimenti
    metabolites_dict = {}
    compartments_dict = {}
    
    print(f"Costruendo modello da {len(df)} reazioni...")
    
    # Prima passata: estrai metaboliti e compartimenti dalle formule delle reazioni
    for idx, row in df.iterrows():
        reaction_formula = str(row['Reaction']).strip()
        if not reaction_formula or reaction_formula == 'nan':
            continue
            
        # Estrai metaboliti dalla formula della reazione
        metabolites = extract_metabolites_from_reaction(reaction_formula)
        
        for met_id in metabolites:
            compartment = extract_compartment_from_metabolite(met_id)
            
            # Aggiungi compartimento se non esiste
            if compartment not in compartments_dict:
                compartments_dict[compartment] = compartment
            
            # Aggiungi metabolita se non esiste
            if met_id not in metabolites_dict:
                metabolites_dict[met_id] = Metabolite(
                    id=met_id,
                    compartment=compartment,
                    name=met_id.replace(f"_{compartment}", "").replace("__", "_")
                )
    
    # Aggiungi compartimenti al modello
    model.compartments = compartments_dict
    
    # Aggiungi metaboliti al modello  
    model.add_metabolites(list(metabolites_dict.values()))
    
    print(f"Aggiunti {len(metabolites_dict)} metaboliti e {len(compartments_dict)} compartimenti")
    
    # Seconda passata: aggiungi le reazioni
    reactions_added = 0
    reactions_skipped = 0
    
    for idx, row in df.iterrows():

        reaction_id = str(row['ReactionID']).strip()
        reaction_formula = str(row['Formula']).strip()
        
        # Salta reazioni senza formula
        if not reaction_formula or reaction_formula == 'nan':
            raise ValueError(f"Formula della reazione mancante {reaction_id}")
        
        # Crea la reazione
        reaction = Reaction(reaction_id)
        reaction.name = reaction_id
        
        # Imposta bounds
        reaction.lower_bound = float(row['lower_bound']) if pd.notna(row['lower_bound']) else -1000.0
        reaction.upper_bound = float(row['upper_bound']) if pd.notna(row['upper_bound']) else 1000.0
        
        # Aggiungi gene rule se presente
        if pd.notna(row['GPR']) and str(row['GPR']).strip():
            reaction.gene_reaction_rule = str(row['GPR']).strip()
        
        # Parse della formula della reazione
        try:
            parse_reaction_formula(reaction, reaction_formula, metabolites_dict)
        except Exception as e:
            print(f"Errore nel parsing della reazione {reaction_id}: {e}")
            reactions_skipped += 1
            continue
        
        # Aggiungi la reazione al modello
        model.add_reactions([reaction])
        reactions_added += 1
            
    
    print(f"Aggiunte {reactions_added} reazioni, saltate {reactions_skipped} reazioni")
    
    # set objective function
    set_objective_from_csv(model, df, obj_col="ObjectiveCoefficient")

    # Imposta il medium
    set_medium_from_data(model, df)
    
    print(f"Modello completato: {len(model.reactions)} reazioni, {len(model.metabolites)} metaboliti")
    
    return model


# Estrae tutti gli ID metaboliti nella formula (gestisce prefissi numerici + underscore)
def extract_metabolites_from_reaction(reaction_formula: str) -> Set[str]:
    """
    Estrae gli ID dei metaboliti da una formula di reazione.
    Pattern robusto: cattura token che terminano con _<compartimento> (es. _c, _m, _e)
    e permette che comincino con cifre o underscore.
    """
    metabolites = set()
    # coefficiente opzionale seguito da un token che termina con _<letters>
    pattern = r'(?:\d+(?:\.\d+)?\s+)?([A-Za-z0-9_]+_[a-z]+)'
    matches = re.findall(pattern, reaction_formula)
    metabolites.update(matches)
    return metabolites


def extract_compartment_from_metabolite(metabolite_id: str) -> str:
    """
    Estrae il compartimento dall'ID del metabolita.
    """
    # Il compartimento è solitamente l'ultima lettera dopo l'underscore
    if '_' in metabolite_id:
        return metabolite_id.split('_')[-1]
    return 'c'  # default cytoplasm


def parse_reaction_formula(reaction: Reaction, formula: str, metabolites_dict: Dict[str, Metabolite]):
    """
    Parsa una formula di reazione e imposta i metaboliti con i loro coefficienti.
    """

    if reaction.id == 'EX_thbpt_e':
        print(reaction.id)
        print(formula)
    # Dividi in parte sinistra e destra
    if '<=>' in formula:
        left, right = formula.split('<=>')
        reversible = True
    elif '<--' in formula:
        left, right = formula.split('<--')
        reversible = False
        left, right = left, right
    elif '-->' in formula:
        left, right = formula.split('-->')
        reversible = False
    elif '<-' in formula:
        left, right = formula.split('<-')
        reversible = False
        left, right = left, right
    else:
        raise ValueError(f"Formato reazione non riconosciuto: {formula}")
    
    # Parse dei metaboliti e coefficienti
    reactants = parse_metabolites_side(left.strip())
    products = parse_metabolites_side(right.strip())
    
    # Aggiungi metaboliti alla reazione
    metabolites_to_add = {}
    
    # Reagenti (coefficienti negativi)
    for met_id, coeff in reactants.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = -coeff
    
    # Prodotti (coefficienti positivi)
    for met_id, coeff in products.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = coeff
    
    reaction.add_metabolites(metabolites_to_add)


def parse_metabolites_side(side_str: str) -> Dict[str, float]:
    """
    Parsa un lato della reazione per estrarre metaboliti e coefficienti.
    """
    metabolites = {}
    if not side_str or side_str.strip() == '':
        return metabolites

    terms = side_str.split('+')
    for term in terms:
        term = term.strip()
        if not term:
            continue

        # pattern allineato: coefficiente opzionale + id che termina con _<compartimento>
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
            # Assicuriamoci che la reazione esista nel modello
            if reaction_id in model.reactions:
                obj_dict[model.reactions.get_by_id(reaction_id)] = coeff
            else:
                print(f"Warning: reaction {reaction_id} not found in model, skipping for objective.")

    if not obj_dict:
        raise ValueError("No reactions found with non-zero objective coefficient.")

    # Imposta la funzione obiettivo
    model.objective = obj_dict
    print(f"Objective set with {len(obj_dict)} reactions.")




def set_medium_from_data(model: cobraModel, df: pd.DataFrame):
    """
    Imposta il medium basato sulla colonna InMedium.
    """
    medium_reactions = df[df['InMedium'] == True]['ReactionID'].tolist()
    
    medium_dict = {}
    for rxn_id in medium_reactions:
        if rxn_id in [r.id for r in model.reactions]:
            reaction = model.reactions.get_by_id(rxn_id)
            if reaction.lower_bound < 0:  # Solo reazioni di uptake
                medium_dict[rxn_id] = abs(reaction.lower_bound)
    
    if medium_dict:
        model.medium = medium_dict
        print(f"Medium impostato con {len(medium_dict)} componenti")


def validate_model(model: cobraModel) -> Dict[str, any]:
    """
    Valida il modello e fornisce statistiche di base.
    """
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
        # Test di crescita
        solution = model.optimize()
        validation['growth_rate'] = solution.objective_value
        validation['status'] = solution.status
    except Exception as e:
        validation['growth_rate'] = None
        validation['status'] = f"Error: {e}"
    
    return validation

def convert_genes(model,annotation):
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
    Cerca colonne utili e ritorna dict {ensg: colname1, hgnc_id: colname2, ...}
    Lancia ValueError se non trova almeno un mapping utile.
    """
    cols = { _normalize_colname(c): c for c in mapping_df.columns }
    chosen = {}
    # possibili nomi per ciascuna categoria
    candidates = {
        'ensg': ['ensg', 'ensembl_gene_id', 'ensembl'],
        'hgnc_id': ['hgnc_id', 'hgnc', 'hgnc:'],
        'hgnc_symbol': ['hgnc_symbol', 'hgnc symbol', 'symbol'],
        'entrez_id': ['entrez', 'entrez_id', 'entrezgene']
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
    Verifica che, nel mapping_df (eventualmente già filtrato sui source di interesse),
    ogni target sia associato ad al massimo un source. Se trova target associati a
    >1 source solleva ValueError mostrando esempi.

    - mapping_df: DataFrame con colonne source_col, target_col
    - model_source_genes: se fornito, è un set di source normalizzati che stiamo traducendo
      (se None, si usa tutto mapping_df)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if mapping_df.empty:
        logger.warning("Mapping dataframe is empty for the requested source genes; skipping uniqueness validation.")
        return

    # normalizza le colonne temporanee per gruppi (senza modificare il df originale)
    tmp = mapping_df[[source_col, target_col]].copy()
    tmp['_src_norm'] = tmp[source_col].astype(str).map(_normalize_gene_id)
    tmp['_tgt_norm'] = tmp[target_col].astype(str).str.strip()

    # se è passato un insieme di geni modello, filtra qui (già fatto nella chiamata, ma doppio-check ok)
    if model_source_genes is not None:
        tmp = tmp[tmp['_src_norm'].isin(model_source_genes)]

    if tmp.empty:
        logger.warning("After filtering to model source genes, mapping table is empty — nothing to validate.")
        return

    # costruisci il reverse mapping target -> set(sources)
    grouped = tmp.groupby('_tgt_norm')['_src_norm'].agg(lambda s: set(s.dropna()))
    # trova target con più di 1 source
    problematic = {t: sorted(list(s)) for t, s in grouped.items() if len(s) > 1}

    if problematic:
        # prepara messaggio di errore con esempi (fino a 20)
        sample_items = list(problematic.items())[:20]
        msg_lines = ["Mapping validation failed: some target IDs are associated with multiple source IDs."]
        for tgt, sources in sample_items:
            msg_lines.append(f"  - target '{tgt}' <- sources: {', '.join(sources)}")
        if len(problematic) > len(sample_items):
            msg_lines.append(f"  ... and {len(problematic) - len(sample_items)} more cases.")
        full_msg = "\n".join(msg_lines)
        # loggare e sollevare errore
        logger.error(full_msg)
        raise ValueError(full_msg)

    # se tutto ok
    logger.info("Mapping validation passed: no target ID is associated with multiple source IDs (within filtered set).")


def _normalize_gene_id(g: str) -> str:
    """Rendi consistente un gene id per l'uso come chiave (rimuove prefissi come 'HGNC:' e strip)."""
    if g is None:
        return ""
    g = str(g).strip()
    # remove common prefixes
    g = re.sub(r'^(HGNC:)', '', g, flags=re.IGNORECASE)
    g = re.sub(r'^(ENSG:)', '', g, flags=re.IGNORECASE)
    return g

# ---------- Main public function ----------
def translate_model_genes(model: 'cobra.Model',
                         mapping_df: 'pd.DataFrame',
                         target_nomenclature: str,
                         source_nomenclature: str = 'hgnc_id',
                         logger: Optional[logging.Logger] = None) -> 'cobra.Model':
    """
    Translate model genes from source_nomenclature to target_nomenclature.
    mapping_df should contain at least columns that allow the mapping
    (e.g. ensg, hgnc_id, hgnc_symbol, entrez).
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
    # support synonyms: user may pass "ensg" or "ENSG" etc.
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
        # fallback: if user passed 'hgnc_id' but chosen has only 'hgnc_symbol', it's not useful
        # we require at least the source column to exist
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

    # Se non ci sono righe rilevanti, avvisa (possono non esserci mapping per i geni presenti)
    if filtered_map.empty:
        logger.warning("No mapping rows correspond to source genes present in the model after filtering. Proceeding with empty mapping (no translation will occur).")

    # --- VALIDAZIONE: nessun target deve essere mappato da piu' di un source (nell'insieme filtrato) ---
    # Se vuoi la verifica su tutto il dataframe (non solo sui geni del modello), passa model_source_genes=None.
    _validate_target_uniqueness(filtered_map, col_for_src, col_for_tgt, model_source_genes=model_source_genes, logger=logger)

    # Ora crea il mapping solo sul sottoinsieme filtrato (piu' efficiente)
    # ATTENZIONE: _create_gene_mapping si aspetta i nomi originali delle colonne
    # quindi passiamo filtered_map con le colonne rimappate (senza la col_for_src + "_norm")
    gene_mapping = _create_gene_mapping(filtered_map, col_for_src, col_for_tgt, logger)

    # copy model
    model_copy = model.copy()

    # statistics
    stats = {'translated': 0, 'one_to_one': 0, 'one_to_many': 0, 'not_found': 0}
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
                rxn.gene_reaction_rule = new_gpr
                logger.debug(f"Reaction {rxn.id}: '{gpr}' -> '{new_gpr}'")

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

    final_ids = {g.id for g in final_genes}
    logger.info(f"Genes in model: {len(original_genes)} -> {len(final_ids)}")

    if unmapped_genes:
        logger.warning(f"Unmapped tokens ({len(unmapped_genes)}): {', '.join(unmapped_genes[:20])}{(' ...' if len(unmapped_genes)>20 else '')}")
    if multi_mapping_genes:
        logger.info(f"Multi-mapping examples ({len(multi_mapping_genes)}):")
        for orig, targets in multi_mapping_genes[:10]:
            logger.info(f"  {orig} -> {', '.join(targets)}")