import os
import csv
import cobra
import pickle
import argparse
import pandas as pd
import re
from typing import Optional, Tuple, Union, List, Dict, Set
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
        reaction_formula = str(row['Reaction']).strip()
        
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
        if pd.notna(row['Rule']) and str(row['Rule']).strip():
            reaction.gene_reaction_rule = str(row['Rule']).strip()
        
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
    
    # Imposta l'obiettivo di biomassa
    set_biomass_objective(model)
    
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



def set_biomass_objective(model: cobraModel):
    """
    Imposta la reazione di biomassa come obiettivo.
    """
    biomass_reactions = [r for r in model.reactions if 'biomass' in r.id.lower()]
    
    if biomass_reactions:
        model.objective = biomass_reactions[0].id
        print(f"Obiettivo impostato su: {biomass_reactions[0].id}")
    else:
        print("Nessuna reazione di biomassa trovata")


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
