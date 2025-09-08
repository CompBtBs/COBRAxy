import argparse
import utils.general_utils as utils
from typing import Optional, Dict, Set, List, Tuple
import os
import numpy as np
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite
import re
import sys
import csv
from joblib import Parallel, delayed, cpu_count

# , medium

################################# process args ###############################
def process_args(args :List[str] = None) -> argparse.Namespace:
    """
    Processes command-line arguments.

    Args:
        args (list): List of command-line arguments.

    Returns:
        Namespace: An object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s')
    
    
    parser.add_argument("-mo", "--model_upload", type = str,
        help = "path to input file with custom rules, if provided")
    
    parser.add_argument("-meo", "--medium", type = str,
        help = "path to input file with custom medium, if provided")

    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")
    
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    
    parser.add_argument('-ir', '--input_ras',
                        type=str,
                        required = False,
                        help = 'input ras')
    
    parser.add_argument('-rn', '--name',
                type=str,
                help = 'ras class names')
    
    parser.add_argument('-rs', '--ras_selector',
                        required = True,
                        type=utils.Bool("using_RAS"),
                        help = 'ras selector')

    parser.add_argument('-cc', '--cell_class',
                    type = str,
                    help = 'output of cell class')
    parser.add_argument(
        '-idop', '--output_path', 
        type = str,
        default='ras_to_bounds/',
        help = 'output path for maps')
    
    
    ARGS = parser.parse_args(args)
    return ARGS

########################### warning ###########################################
def warning(s :str) -> None:
    """
    Log a warning message to an output log file and print it to the console.

    Args:
        s (str): The warning message to be logged and printed.
    
    Returns:
      None
    """
    with open(ARGS.out_log, 'a') as log:
        log.write(s + "\n\n")
    print(s)

############################ dataset input ####################################
def read_dataset(data :str, name :str) -> pd.DataFrame:
    """
    Read a dataset from a CSV file and return it as a pandas DataFrame.

    Args:
        data (str): Path to the CSV file containing the dataset.
        name (str): Name of the dataset, used in error messages.

    Returns:
        pandas.DataFrame: DataFrame containing the dataset.

    Raises:
        pd.errors.EmptyDataError: If the CSV file is empty.
        sys.exit: If the CSV file has the wrong format, the execution is aborted.
    """
    try:
        dataset = pd.read_csv(data, sep = '\t', header = 0, engine='python')
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    return dataset


def apply_ras_bounds(bounds, ras_row):
    """
    Adjust the bounds of reactions in the model based on RAS values.

    Args:
        bounds (pd.DataFrame): Model bounds.
        ras_row (pd.Series): A row from a RAS DataFrame containing scaling factors for reaction bounds.
    Returns:
        new_bounds (pd.DataFrame): integrated bounds.
    """
    new_bounds = bounds.copy()
    for reaction in ras_row.index:
        scaling_factor = ras_row[reaction]
        if not np.isnan(scaling_factor):
            lower_bound=bounds.loc[reaction, "lower_bound"]
            upper_bound=bounds.loc[reaction, "upper_bound"]
            valMax=float((upper_bound)*scaling_factor)
            valMin=float((lower_bound)*scaling_factor)
            if upper_bound!=0 and lower_bound==0:
                new_bounds.loc[reaction, "upper_bound"] = valMax
            if upper_bound==0 and lower_bound!=0:
                new_bounds.loc[reaction, "lower_bound"] = valMin
            if upper_bound!=0 and lower_bound!=0:
                new_bounds.loc[reaction, "lower_bound"] = valMin
                new_bounds.loc[reaction, "upper_bound"] = valMax
    return new_bounds

def process_ras_cell(cellName, ras_row, model, rxns_ids, output_folder):
    """
    Process a single RAS cell, apply bounds, and save the bounds to a CSV file.

    Args:
        cellName (str): The name of the RAS cell (used for naming the output file).
        ras_row (pd.Series): A row from a RAS DataFrame containing scaling factors for reaction bounds.
        model (cobra.Model): The metabolic model to be modified.
        rxns_ids (list of str): List of reaction IDs to which the scaling factors will be applied.
        output_folder (str): Folder path where the output CSV file will be saved.
    
    Returns:
        None
    """
    bounds = pd.DataFrame([(rxn.lower_bound, rxn.upper_bound) for rxn in model.reactions], index=rxns_ids, columns=["lower_bound", "upper_bound"])
    new_bounds = apply_ras_bounds(bounds, ras_row)
    new_bounds.to_csv(output_folder + cellName + ".csv", sep='\t', index=True)
    pass

def generate_bounds(model: cobra.Model, ras=None, output_folder='output/') -> pd.DataFrame:
    """
    Generate reaction bounds for a metabolic model based on medium conditions and optional RAS adjustments.
    
    Args:
        model (cobra.Model): The metabolic model for which bounds will be generated.
        medium (dict): A dictionary where keys are reaction IDs and values are the medium conditions.
        ras (pd.DataFrame, optional): RAS pandas dataframe. Defaults to None.
        output_folder (str, optional): Folder path where output CSV files will be saved. Defaults to 'output/'.

    Returns:
        pd.DataFrame: DataFrame containing the bounds of reactions in the model.
    """
    rxns_ids = [rxn.id for rxn in model.reactions]            
            
    # Perform Flux Variability Analysis (FVA) on this medium
    df_FVA = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0, processes=1).round(8)
    
    # Set FVA bounds
    for reaction in rxns_ids:
        model.reactions.get_by_id(reaction).lower_bound = float(df_FVA.loc[reaction, "minimum"])
        model.reactions.get_by_id(reaction).upper_bound = float(df_FVA.loc[reaction, "maximum"])

    if ras is not None:
        Parallel(n_jobs=cpu_count())(delayed(process_ras_cell)(cellName, ras_row, model, rxns_ids, output_folder) for cellName, ras_row in ras.iterrows())
    else:
        bounds = pd.DataFrame([(rxn.lower_bound, rxn.upper_bound) for rxn in model.reactions], index=rxns_ids, columns=["lower_bound", "upper_bound"])
        newBounds = apply_ras_bounds(bounds, pd.Series([1]*len(rxns_ids), index=rxns_ids))
        newBounds.to_csv(output_folder + "bounds.csv", sep='\t', index=True)
    pass

# TODO: VALUTARE QUALI DI QUESTE FUNZIONI METTERE IN UTILS.PY
def build_cobra_model_from_csv(csv_path: str, model_id: str = "ENGRO2_custom") -> cobra.Model:
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
    model = Model(model_id)
    
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
        try:
            reaction_id = str(row['ReactionID']).strip()
            if reaction_id == 'EX_thbpt_e':
              print('qui')
              print(reaction_id)
              print(str(row['Reaction']).strip())
              print('qui')
            reaction_formula = str(row['Reaction']).strip()
            
            # Salta reazioni senza formula
            if not reaction_formula or reaction_formula == 'nan':
                reactions_skipped += 1
                continue
            
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
            
        except Exception as e:
            print(f"Errore nell'aggiungere la reazione {reaction_id}: {e}")
            reactions_skipped += 1
            continue
    
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



def set_biomass_objective(model: Model):
    """
    Imposta la reazione di biomassa come obiettivo.
    """
    biomass_reactions = [r for r in model.reactions if 'biomass' in r.id.lower()]
    
    if biomass_reactions:
        model.objective = biomass_reactions[0].id
        print(f"Obiettivo impostato su: {biomass_reactions[0].id}")
    else:
        print("Nessuna reazione di biomassa trovata")


def set_medium_from_data(model: Model, df: pd.DataFrame):
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


def validate_model(model: Model) -> Dict[str, any]:
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


############################# main ###########################################
def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    if not os.path.exists('ras_to_bounds'):
        os.makedirs('ras_to_bounds')


    global ARGS
    ARGS = process_args(args)

    if(ARGS.ras_selector == True):
        ras_file_list = ARGS.input_ras.split(",")
        ras_file_names = ARGS.name.split(",")
        if len(ras_file_names) != len(set(ras_file_names)):
            error_message = "Duplicated file names in the uploaded RAS matrices."
            warning(error_message)
            raise ValueError(error_message)
            pass
        ras_class_names = []
        for file in ras_file_names:
            ras_class_names.append(file.rsplit(".", 1)[0])
        ras_list = []
        class_assignments = pd.DataFrame(columns=["Patient_ID", "Class"])
        for ras_matrix, ras_class_name in zip(ras_file_list, ras_class_names):
            ras = read_dataset(ras_matrix, "ras dataset")
            ras.replace("None", None, inplace=True)
            ras.set_index("Reactions", drop=True, inplace=True)
            ras = ras.T
            ras = ras.astype(float)
            if(len(ras_file_list)>1):
                #append class name to patient id (dataframe index)
                ras.index = [f"{idx}_{ras_class_name}" for idx in ras.index]
            else:
                ras.index = [f"{idx}" for idx in ras.index]
            ras_list.append(ras)
            for patient_id in ras.index:
                class_assignments.loc[class_assignments.shape[0]] = [patient_id, ras_class_name]
        
        
        # Concatenate all ras DataFrames into a single DataFrame
        ras_combined = pd.concat(ras_list, axis=0)
        # Normalize the RAS values by max RAS
        ras_combined = ras_combined.div(ras_combined.max(axis=0))
        ras_combined.dropna(axis=1, how='all', inplace=True)


    
    #model_type :utils.Model = ARGS.model_selector
    #if model_type is utils.Model.Custom:
    #    model = model_type.getCOBRAmodel(customPath = utils.FilePath.fromStrPath(ARGS.model), customExtension = utils.FilePath.fromStrPath(ARGS.model_name).ext)
    #else:
    #    model = model_type.getCOBRAmodel(toolDir=ARGS.tool_dir)

    # TODO LOAD MODEL FROM UPLOAD

    model = build_cobra_model_from_csv(ARGS.model_upload)

    validation = validate_model(model)

    print("\n=== VALIDAZIONE MODELLO ===")
    for key, value in validation.items():
        print(f"{key}: {value}")

    #if(ARGS.medium_selector == "Custom"):
    #    medium = read_dataset(ARGS.medium, "medium dataset")
    #    medium.set_index(medium.columns[0], inplace=True)
    #    medium = medium.astype(float)
    #    medium = medium[medium.columns[0]].to_dict()
    #else:
    #    df_mediums = pd.read_csv(ARGS.tool_dir + "/local/medium/medium.csv", index_col = 0)
    #    ARGS.medium_selector = ARGS.medium_selector.replace("_", " ")
    #    medium = df_mediums[[ARGS.medium_selector]]
    #    medium = medium[ARGS.medium_selector].to_dict()

    if(ARGS.ras_selector == True):
        generate_bounds(model, ras = ras_combined, output_folder=ARGS.output_path)
        class_assignments.to_csv(ARGS.cell_class, sep = '\t', index = False)
    else:
        generate_bounds(model, output_folder=ARGS.output_path)

    pass
        
##############################################################################
if __name__ == "__main__":
    main()