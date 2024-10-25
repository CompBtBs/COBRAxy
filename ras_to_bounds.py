import argparse
import utils.general_utils as utils
from typing import Optional, List
import os
import numpy as np
import pandas as pd
import cobra
import sys
import csv
from joblib import Parallel, delayed, cpu_count

################################# process args ###############################
def process_args(args :List[str]) -> argparse.Namespace:
    """
    Processes command-line arguments.

    Args:
        args (list): List of command-line arguments.

    Returns:
        Namespace: An object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s')
    
    parser.add_argument(
        '-ms', '--model_selector', 
        type = utils.Model, default = utils.Model.ENGRO2, choices = [utils.Model.ENGRO2, utils.Model.Custom],
        help = 'chose which type of model you want use')
    
    parser.add_argument("-mo", "--model", type = str,
        help = "path to input file with custom rules, if provided")
    
    parser.add_argument("-mn", "--model_name", type = str, help = "custom mode name")

    parser.add_argument(
        '-mes', '--medium_selector', 
        default = "allOpen",
        help = 'chose which type of medium you want use')
    
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
    
    
    ARGS = parser.parse_args()
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


def apply_ras_bounds(model, ras_row):
    """
    Adjust the bounds of reactions in the model based on RAS values.

    Args:
        model (cobra.Model): The metabolic model to be modified.
        ras_row (pd.Series): A row from a RAS DataFrame containing scaling factors for reaction bounds.
    Returns:
        None
    """
    for reaction in ras_row.index:
        scaling_factor = ras_row[reaction]
        lower_bound=model.reactions.get_by_id(reaction).lower_bound
        upper_bound=model.reactions.get_by_id(reaction).upper_bound
        valMax=float((upper_bound)*scaling_factor)
        valMin=float((lower_bound)*scaling_factor)
        if upper_bound!=0 and lower_bound==0:
            model.reactions.get_by_id(reaction).upper_bound=valMax
        if upper_bound==0 and lower_bound!=0:
            model.reactions.get_by_id(reaction).lower_bound=valMin
        if upper_bound!=0 and lower_bound!=0:
            model.reactions.get_by_id(reaction).lower_bound=valMin
            model.reactions.get_by_id(reaction).upper_bound=valMax
    pass

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
    model_new = model.copy()
    apply_ras_bounds(model_new, ras_row)
    bounds = pd.DataFrame([(rxn.lower_bound, rxn.upper_bound) for rxn in model_new.reactions], index=rxns_ids, columns=["lower_bound", "upper_bound"])
    bounds.to_csv(output_folder + cellName + ".csv", sep='\t', index=True)
    pass

def generate_bounds(model: cobra.Model, medium: dict, ras=None, output_folder='output/') -> pd.DataFrame:
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

    # Set all reactions to zero in the medium
    for rxn_id, _ in model.medium.items():
        model.reactions.get_by_id(rxn_id).lower_bound = float(0.0)
    
    # Set medium conditions
    for reaction, value in medium.items():
        if value is not None:
            model.reactions.get_by_id(reaction).lower_bound = -float(value)
            
            
    # Perform Flux Variability Analysis (FVA) on this medium
    df_FVA = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0, processes=1).round(8)
    
    # Set FVA bounds
    for reaction in rxns_ids:
        model.reactions.get_by_id(reaction).lower_bound = float(df_FVA.loc[reaction, "minimum"])
        model.reactions.get_by_id(reaction).upper_bound = float(df_FVA.loc[reaction, "maximum"])

    if ras is not None:
        Parallel(n_jobs=cpu_count())(delayed(process_ras_cell)(cellName, ras_row, model, rxns_ids, output_folder) for cellName, ras_row in ras.iterrows())
    else:
        model_new = model.copy()
        apply_ras_bounds(model_new, pd.Series([1]*len(rxns_ids), index=rxns_ids))
        bounds = pd.DataFrame([(rxn.lower_bound, rxn.upper_bound) for rxn in model_new.reactions], index=rxns_ids, columns=["lower_bound", "upper_bound"])
        bounds.to_csv(output_folder + "bounds.csv", sep='\t', index=True)
    pass



############################# main ###########################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    if not os.path.exists('ras_to_bounds'):
        os.makedirs('ras_to_bounds')


    global ARGS
    ARGS = process_args(sys.argv)

    ARGS.output_folder = 'ras_to_bounds/'

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
            ras_class_names.append(file.split(".")[0])
        ras_list = []
        class_assignments = pd.DataFrame(columns=["Patient_ID", "Class"])
        for ras_matrix, ras_class_name in zip(ras_file_list, ras_class_names):
            ras = read_dataset(ras_matrix, "ras dataset")
            ras.replace("None", None, inplace=True)
            ras.set_index("Reactions", drop=True, inplace=True)
            ras = ras.T
            ras = ras.astype(float)
            #append class name to patient id (dataframe index)
            ras.index = [f"{idx}_{ras_class_name}" for idx in ras.index]
            ras_list.append(ras)
            for patient_id in ras.index:
                class_assignments.loc[class_assignments.shape[0]] = [patient_id, ras_class_name]
        
        
        # Concatenate all ras DataFrames into a single DataFrame
        ras_combined = pd.concat(ras_list, axis=0)
        # Normalize the RAS values by max RAS
        ras_combined = ras_combined.div(ras_combined.max(axis=0))
        ras_combined.dropna(axis=1, how='all', inplace=True)


    
    model_type :utils.Model = ARGS.model_selector
    if model_type is utils.Model.Custom:
        model = model_type.getCOBRAmodel(customPath = utils.FilePath.fromStrPath(ARGS.model), customExtension = utils.FilePath.fromStrPath(ARGS.model_name).ext)
    else:
        model = model_type.getCOBRAmodel(toolDir=ARGS.tool_dir)

    if(ARGS.medium_selector == "Custom"):
        medium = read_dataset(ARGS.medium, "medium dataset")
        medium.set_index(medium.columns[0], inplace=True)
        medium = medium.astype(float)
        medium = medium[medium.columns[0]].to_dict()
    else:
        df_mediums = pd.read_csv(ARGS.tool_dir + "/local/medium/medium.csv", index_col = 0)
        ARGS.medium_selector = ARGS.medium_selector.replace("_", " ")
        medium = df_mediums[[ARGS.medium_selector]]
        medium = medium[ARGS.medium_selector].to_dict()

    if(ARGS.ras_selector == True):
        generate_bounds(model, medium, ras = ras_combined, output_folder=ARGS.output_folder)
        class_assignments.to_csv(ARGS.cell_class, sep = '\t', index = False)
    else:
        generate_bounds(model, medium, output_folder=ARGS.output_folder)

    pass
        
##############################################################################
if __name__ == "__main__":
    main()