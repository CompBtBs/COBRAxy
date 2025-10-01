"""
Apply RAS-based scaling to reaction bounds and optionally save updated models.

Workflow:
- Read one or more RAS matrices (patients/samples x reactions)
- Normalize and merge them, optionally adding class suffixes to sample IDs
- Build a COBRA model from a tabular CSV
- Run FVA to initialize bounds, then scale per-sample based on RAS values
- Save bounds per sample and optionally export updated models in chosen formats
"""
import argparse
import utils.general_utils as utils
from typing import Optional, Dict, Set, List, Tuple, Union
import os
import numpy as np
import pandas as pd
import cobra
from cobra import Model
import sys
from joblib import Parallel, delayed, cpu_count
import utils.model_utils as modelUtils

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

    parser.add_argument('-cc', '--cell_class',
                    type = str,
                    help = 'output of cell class')
    parser.add_argument(
        '-idop', '--output_path', 
        type = str,
        default='ras_to_bounds/',
        help = 'output path for maps')
    
    parser.add_argument('-sm', '--save_models',
                    type=utils.Bool("save_models"),
                    default=False,
                    help = 'whether to save models with applied bounds')
    
    parser.add_argument('-smp', '--save_models_path',
                        type = str,
                        default='saved_models/',
                        help = 'output path for saved models')
    
    parser.add_argument('-smf', '--save_models_format',
                        type = str,
                        default='csv',
                        help = 'format for saved models (csv, xml, json, mat, yaml, tabular)')

    
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
    if ARGS.out_log:
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


def save_model(model, filename, output_folder, file_format='csv'):
    """
    Save a COBRA model to file in the specified format.
    
    Args:
        model (cobra.Model): The model to save.
        filename (str): Base filename (without extension).
        output_folder (str): Output directory.
        file_format (str): File format ('xml', 'json', 'mat', 'yaml', 'tabular', 'csv').
    
    Returns:
        None
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    try:
        if file_format == 'tabular' or file_format == 'csv':
            # Special handling for tabular format using utils functions
            filepath = os.path.join(output_folder, f"{filename}.csv")
            
            rules = modelUtils.generate_rules(model, asParsed = False)
            reactions = modelUtils.generate_reactions(model, asParsed = False)
            bounds = modelUtils.generate_bounds(model)
            medium = modelUtils.get_medium(model)
            
            compartments = modelUtils.generate_compartments(model)

            df_rules = pd.DataFrame(list(rules.items()), columns = ["ReactionID", "Rule"])
            df_reactions = pd.DataFrame(list(reactions.items()), columns = ["ReactionID", "Reaction"])
            df_bounds = bounds.reset_index().rename(columns = {"index": "ReactionID"})
            df_medium = medium.rename(columns = {"reaction": "ReactionID"})
            df_medium["InMedium"] = True

            merged = df_reactions.merge(df_rules, on = "ReactionID", how = "outer")
            merged = merged.merge(df_bounds, on = "ReactionID", how = "outer")
            # Add compartments only if they exist
            if compartments is not None:
                merged = merged.merge(compartments, on = "ReactionID", how = "outer")
            
            merged = merged.merge(df_medium, on = "ReactionID", how = "left")
            merged["InMedium"] = merged["InMedium"].fillna(False)
            merged = merged.sort_values(by = "InMedium", ascending = False)
            
            merged.to_csv(filepath, sep="\t", index=False)
            
        else:
            # Standard COBRA formats
            filepath = os.path.join(output_folder, f"{filename}.{file_format}")
            
            if file_format == 'xml':
                cobra.io.write_sbml_model(model, filepath)
            elif file_format == 'json':
                cobra.io.save_json_model(model, filepath)
            elif file_format == 'mat':
                cobra.io.save_matlab_model(model, filepath)
            elif file_format == 'yaml':
                cobra.io.save_yaml_model(model, filepath)
            else:
                raise ValueError(f"Unsupported format: {file_format}")
        
        print(f"Model saved: {filepath}")
        
    except Exception as e:
        warning(f"Error saving model {filename}: {str(e)}")

def apply_bounds_to_model(model, bounds):
    """
    Apply bounds from a DataFrame to a COBRA model.
    
    Args:
        model (cobra.Model): The metabolic model to modify.
        bounds (pd.DataFrame): DataFrame with reaction bounds.
    
    Returns:
        cobra.Model: Modified model with new bounds.
    """
    model_copy = model.copy()
    for reaction_id in bounds.index:
        try:
            reaction = model_copy.reactions.get_by_id(reaction_id)
            reaction.lower_bound = bounds.loc[reaction_id, "lower_bound"]
            reaction.upper_bound = bounds.loc[reaction_id, "upper_bound"]
        except KeyError:
            # Reaction not found in model, skip
            continue
    return model_copy

def process_ras_cell(cellName, ras_row, model, rxns_ids, output_folder, save_models=False, save_models_path='saved_models/', save_models_format='csv'):
    """
    Process a single RAS cell, apply bounds, and save the bounds to a CSV file.

    Args:
        cellName (str): The name of the RAS cell (used for naming the output file).
        ras_row (pd.Series): A row from a RAS DataFrame containing scaling factors for reaction bounds.
        model (cobra.Model): The metabolic model to be modified.
        rxns_ids (list of str): List of reaction IDs to which the scaling factors will be applied.
        output_folder (str): Folder path where the output CSV file will be saved.
        save_models (bool): Whether to save models with applied bounds.
        save_models_path (str): Path where to save models.
        save_models_format (str): Format for saved models.
    
    Returns:
        None
    """
    bounds = pd.DataFrame([(rxn.lower_bound, rxn.upper_bound) for rxn in model.reactions], index=rxns_ids, columns=["lower_bound", "upper_bound"])
    new_bounds = apply_ras_bounds(bounds, ras_row)
    new_bounds.to_csv(output_folder + cellName + ".csv", sep='\t', index=True)
    
    # Save model if requested
    if save_models:
        modified_model = apply_bounds_to_model(model, new_bounds)
        save_model(modified_model, cellName, save_models_path, save_models_format)
    
    return

def generate_bounds_model(model: cobra.Model, ras=None, output_folder='output/', save_models=False, save_models_path='saved_models/', save_models_format='csv') -> pd.DataFrame:
    """
    Generate reaction bounds for a metabolic model based on medium conditions and optional RAS adjustments.
    
    Args:
        model (cobra.Model): The metabolic model for which bounds will be generated.
        ras (pd.DataFrame, optional): RAS pandas dataframe. Defaults to None.
        output_folder (str, optional): Folder path where output CSV files will be saved. Defaults to 'output/'.
        save_models (bool): Whether to save models with applied bounds.
        save_models_path (str): Path where to save models.
        save_models_format (str): Format for saved models.

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
        Parallel(n_jobs=cpu_count())(delayed(process_ras_cell)(
            cellName, ras_row, model, rxns_ids, output_folder, 
            save_models, save_models_path, save_models_format
        ) for cellName, ras_row in ras.iterrows())
    else:
        raise ValueError("RAS DataFrame is None. Cannot generate bounds without RAS data.")
    return

############################# main ###########################################
def main(args:List[str] = None) -> None:
    """
    Initialize and execute RAS-to-bounds pipeline based on the frontend input arguments.

    Returns:
        None
    """
    if not os.path.exists('ras_to_bounds'):
        os.makedirs('ras_to_bounds')

    global ARGS
    ARGS = process_args(args)


    ras_file_list = ARGS.input_ras.split(",")
    ras_file_names = ARGS.name.split(",")
    if len(ras_file_names) != len(set(ras_file_names)):
        error_message = "Duplicated file names in the uploaded RAS matrices."
        warning(error_message)
        raise ValueError(error_message)
        
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
            # Append class name to patient id (DataFrame index)
            ras.index = [f"{idx}_{ras_class_name}" for idx in ras.index]
        else:
            ras.index = [f"{idx}" for idx in ras.index]
        ras_list.append(ras)
        for patient_id in ras.index:
            class_assignments.loc[class_assignments.shape[0]] = [patient_id, ras_class_name]
    
        
    # Concatenate all RAS DataFrames into a single DataFrame
        ras_combined = pd.concat(ras_list, axis=0)
    # Normalize RAS values column-wise by max RAS
        ras_combined = ras_combined.div(ras_combined.max(axis=0))
        ras_combined.dropna(axis=1, how='all', inplace=True)

    model = modelUtils.build_cobra_model_from_csv(ARGS.model_upload)

    validation = modelUtils.validate_model(model)

    print("\n=== MODEL VALIDATION ===")
    for key, value in validation.items():
        print(f"{key}: {value}")


    generate_bounds_model(model, ras=ras_combined, output_folder=ARGS.output_path,
                    save_models=ARGS.save_models, save_models_path=ARGS.save_models_path,
                    save_models_format=ARGS.save_models_format)
    class_assignments.to_csv(ARGS.cell_class, sep='\t', index=False)


    return
        
##############################################################################
if __name__ == "__main__":
    main()