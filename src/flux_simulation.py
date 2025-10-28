"""
Flux sampling and analysis utilities for COBRA models.

This script supports two modes:
- Mode 1 (model_and_bounds=True): load a base model and apply bounds from
    separate files before sampling.
- Mode 2 (model_and_bounds=False): load complete models and sample directly.

Sampling algorithms supported: OPTGP and CBS. Outputs include flux samples
and optional analyses (pFBA, FVA, sensitivity), saved as tabular files.
"""

import argparse
from typing import List
import os
import pandas as pd
import numpy as np
import cobra
from joblib import Parallel, delayed, cpu_count
from cobra.sampling import OptGPSampler
import sys

try:
    from .utils import general_utils as utils
    from .utils import CBS_backend
    from .utils import model_utils
except:
    import utils.general_utils as utils
    import utils.CBS_backend as CBS_backend
    import utils.model_utils as model_utils


################################# process args ###############################
def process_args(args: List[str] = None) -> argparse.Namespace:
    """
    Processes command-line arguments.
    
    Args:
        args (list): List of command-line arguments.
    
    Returns:
        Namespace: An object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='process some value\'s')
    
    parser.add_argument("-mo", "--model_upload", type=str,
        help="path to input file with custom rules, if provided")
    
    parser.add_argument("-mab", "--model_and_bounds", type=str,
        choices=['True', 'False'],
        required=True,
        help="upload mode: True for model+bounds, False for complete models")
    
    parser.add_argument("-ens", "--sampling_enabled", type=str,
        choices=['true', 'false'],
        required=True,
        help="enable sampling: 'true' for sampling, 'false' for no sampling")
    
    parser.add_argument('-ol', '--out_log',
                        help="Output log")
    
    parser.add_argument('-td', '--tool_dir',
                        type=str,
                        default=os.path.dirname(os.path.abspath(__file__)),
                        help='your tool directory (default: auto-detected package location)')
    
    parser.add_argument('-inf', '--input_file',
                        required=True,
                        type=str,
                        help='path to file containing list of input bounds files or complete model files (one per line)')
    
    parser.add_argument('-nif', '--name_file',
                        required=True,
                        type=str,
                        help='path to file containing list of cell names (one per line)')
    
    parser.add_argument('-a', '--algorithm',
                        type=str,
                        choices=['OPTGP', 'CBS'],
                        required=True,
                        help='choose sampling algorithm')
    
    parser.add_argument('-th', '--thinning',
                        type=int,
                        default=100,
                        required=True,
                        help='choose thinning')
    
    parser.add_argument('-ns', '--n_samples',
                        type=int,
                        required=True,
                        help='choose how many samples (set to 0 for optimization only)')
    
    parser.add_argument('-sd', '--seed',
                        type=int,
                        required=True,
                        help='seed for random number generation')
    
    parser.add_argument('-nb', '--n_batches',
                        type=int,
                        required=True,
                        help='choose how many batches')
    
    parser.add_argument('-opt', '--perc_opt',
                        type=float,
                        default=0.9,
                        required=False,
                        help='choose the fraction of optimality for FVA (0-1)')
    
    parser.add_argument('-ot', '--output_type',
                        type=str,
                        required=True,
                        help='output type for sampling results')
    
    parser.add_argument('-ota', '--output_type_analysis',
                        type=str,
                        required=False,
                        help='output type analysis (optimization methods)')

    parser.add_argument('-idop', '--output_path',
                        type=str,
                        default='flux_simulation/',
                        help = 'output path for fluxes')
    
    parser.add_argument('-otm', '--out_mean',
                    type = str,
                    required=False,
                    help = 'output of mean of fluxes')
    
    parser.add_argument('-otmd', '--out_median',
                    type = str,
                    required=False,
                    help = 'output of median of fluxes')

    parser.add_argument('-otq', '--out_quantiles',
                    type = str,
                    required=False,
                    help = 'output of quantiles of fluxes')
    
    parser.add_argument('-otfva', '--out_fva',
                    type = str, 
                    required=False,
                    help = 'output of FVA results')
    parser.add_argument('-otp', '--out_pfba',
                    type = str,
                    required=False,
                    help = 'output of pFBA results')
    parser.add_argument('-ots', '--out_sensitivity',
                    type = str,
                    required=False,
                    help = 'output of sensitivity results')
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


def write_to_file(dataset: pd.DataFrame, path: str, keep_index:bool=False, name:str=None)->None:
    """
    Write a DataFrame to a TSV file under path with a given base name.

    Args:
        dataset: The DataFrame to write.
        name: Base file name (without extension). If None, 'path' is treated as the full file path.
        path: Directory path where the file will be saved.
        keep_index: Whether to keep the DataFrame index in the file.

    Returns:
        None
    """
    dataset.index.name = 'Reactions'
    if name:
        dataset.to_csv(os.path.join(path, name + ".csv"), sep = '\t', index = keep_index)
    else:
        dataset.to_csv(path, sep = '\t', index = keep_index)

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
        dataset = pd.read_csv(data, sep = '\t', header = 0, index_col=0, engine='python')
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    return dataset



def OPTGP_sampler(model: cobra.Model, model_name: str, n_samples: int = 1000, thinning: int = 100, n_batches: int = 1, seed: int = 0) -> None:
    """
    Samples from the OPTGP (Optimal Global Perturbation) algorithm and saves the results to CSV files.
    
    Args:
        model (cobra.Model): The COBRA model to sample from.
        model_name (str): The name of the model, used in naming output files.
        n_samples (int, optional): Number of samples per batch. Default is 1000.
        thinning (int, optional): Thinning parameter for the sampler. Default is 100.
        n_batches (int, optional): Number of batches to run. Default is 1.
        seed (int, optional): Random seed for reproducibility. Default is 0.
        
    Returns:
        None
    """
    import numpy as np
    
    # Get reaction IDs for consistent column ordering
    reaction_ids = [rxn.id for rxn in model.reactions]
    
    # Sample and save each batch as numpy file
    for i in range(n_batches):
        optgp = OptGPSampler(model, thinning, seed)
        samples = optgp.sample(n_samples)
        
        # Save as numpy array (more memory efficient)
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_OPTGP.npy"
        np.save(batch_filename, samples.to_numpy())
        
        seed += 1
    
    # Merge all batches into a single DataFrame
    all_samples = []
    
    for i in range(n_batches):
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_OPTGP.npy"
        batch_data = np.load(batch_filename, allow_pickle=True)
        all_samples.append(batch_data)
    
    # Concatenate all batches
    samplesTotal_array = np.vstack(all_samples)
    
    # Convert back to DataFrame with proper column names
    samplesTotal = pd.DataFrame(samplesTotal_array, columns=reaction_ids)
    
    # Save the final merged result as CSV
    write_to_file(samplesTotal.T, ARGS.output_path, True, name=model_name)
    
    # Clean up temporary numpy files
    for i in range(n_batches):
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_OPTGP.npy"
        if os.path.exists(batch_filename):
            os.remove(batch_filename)


def CBS_sampler(model: cobra.Model, model_name: str, n_samples: int = 1000, n_batches: int = 1, seed: int = 0) -> None:
    """
    Samples using the CBS (Constraint-based Sampling) algorithm and saves the results to CSV files.
    
    Args:
        model (cobra.Model): The COBRA model to sample from.
        model_name (str): The name of the model, used in naming output files.
        n_samples (int, optional): Number of samples per batch. Default is 1000.
        n_batches (int, optional): Number of batches to run. Default is 1.
        seed (int, optional): Random seed for reproducibility. Default is 0.
        
    Returns:
        None
    """
    import numpy as np
    
    # Get reaction IDs for consistent column ordering
    reaction_ids = [reaction.id for reaction in model.reactions]
    
    # Perform FVA analysis once for all batches
    df_FVA = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0).round(6)
    
    # Generate random objective functions for all samples across all batches
    df_coefficients = CBS_backend.randomObjectiveFunction(model, n_samples * n_batches, df_FVA, seed=seed)
    
    # Sample and save each batch as numpy file
    for i in range(n_batches):
        samples = pd.DataFrame(columns=reaction_ids, index=range(n_samples))
        
        try:
            CBS_backend.randomObjectiveFunctionSampling(
                model, 
                n_samples, 
                df_coefficients.iloc[:, i * n_samples:(i + 1) * n_samples], 
                samples
            )
        except Exception as e:
            utils.logWarning(
                f"Warning: GLPK solver has failed for {model_name}. Trying with COBRA interface. Error: {str(e)}",
                ARGS.out_log
            )
            CBS_backend.randomObjectiveFunctionSampling_cobrapy(
                model, 
                n_samples, 
                df_coefficients.iloc[:, i * n_samples:(i + 1) * n_samples],
                samples
            )
        
        # Save as numpy array (more memory efficient)
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_CBS.npy"
        utils.logWarning(batch_filename, ARGS.out_log)
        np.save(batch_filename, samples.to_numpy())
    
    # Merge all batches into a single DataFrame
    all_samples = []
    
    for i in range(n_batches):
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_CBS.npy"
        batch_data = np.load(batch_filename, allow_pickle=True)
        all_samples.append(batch_data)
    
    # Concatenate all batches
    samplesTotal_array = np.vstack(all_samples)
    
    # Convert back to DataFrame with proper column namesq
    samplesTotal = pd.DataFrame(samplesTotal_array, columns=reaction_ids)
    
    # Save the final merged result as CSV
    write_to_file(samplesTotal.T, ARGS.output_path, True, name=model_name)
    
    # Clean up temporary numpy files
    for i in range(n_batches):
        batch_filename = f"{ARGS.output_path}/{model_name}_{i}_CBS.npy"
        if os.path.exists(batch_filename):
            os.remove(batch_filename)



def model_sampler_with_bounds(model_path: str, bounds_path: str, cell_name: str) -> List[pd.DataFrame]:
    """
    MODE 1: Loads model from file, applies bounds from separate bounds file and performs sampling.

    Args:
        model_path (str): Path to the tabular model file.
        bounds_path (str): Path to the CSV file containing the bounds dataset.
        cell_name (str): Name of the cell, used to generate filenames for output.

    Returns:
        List[pd.DataFrame]: A list of DataFrames containing statistics and analysis results.
    """
    
    model_input = model_utils.build_cobra_model_from_csv(model_path)

    validation = model_utils.validate_model(model_input)

    print("\n=== MODEL VALIDATION ===")
    for key, value in validation.items():
        print(f"{key}: {value}")

    model_input.solver.configuration.verbosity = 1
    
    bounds_df = read_dataset(bounds_path, "bounds dataset")
    
    # Apply bounds to model
    for rxn_index, row in bounds_df.iterrows():
        try:
            model_input.reactions.get_by_id(rxn_index).lower_bound = row.lower_bound
            model_input.reactions.get_by_id(rxn_index).upper_bound = row.upper_bound
        except KeyError:
            warning(f"Warning: Reaction {rxn_index} not found in model. Skipping.")

    return perform_sampling_and_analysis(None, cell_name, model_input=model_input)


def perform_sampling_and_analysis(model_path: str, cell_name: str, model_input: cobra.Model=None) -> List[pd.DataFrame]:
    """
    Common function to perform sampling and analysis on a prepared model.

    Args:
        model_path (str): Path to the tabular model file. model with bounds applied.
        cell_name (str): Name of the cell, used to generate filenames for output.

    Returns:
        List[pd.DataFrame]: A list of DataFrames containing statistics and analysis results.
    """

    returnList = []
    
    if model_input is None:
        model_input = model_utils.build_cobra_model_from_csv(model_path)

    if ARGS.sampling_enabled == "true":
    
        if ARGS.algorithm == 'OPTGP':
            OPTGP_sampler(model_input, cell_name, ARGS.n_samples, ARGS.thinning, ARGS.n_batches, ARGS.seed)
        elif ARGS.algorithm == 'CBS':
            CBS_sampler(model_input, cell_name, ARGS.n_samples, ARGS.n_batches, ARGS.seed)

        df_mean, df_median, df_quantiles = fluxes_statistics(cell_name, ARGS.output_types)

        if("fluxes" not in ARGS.output_types):
            os.remove(ARGS.output_path + "/" + cell_name + '.csv')

        returnList = [df_mean, df_median, df_quantiles]

    df_pFBA, df_FVA, df_sensitivity = fluxes_analysis(model_input, cell_name, ARGS.output_type_analysis)

    if("pFBA" in ARGS.output_type_analysis):
        returnList.append(df_pFBA)
    if("FVA" in ARGS.output_type_analysis):
        returnList.append(df_FVA)
    if("sensitivity" in ARGS.output_type_analysis):
        returnList.append(df_sensitivity)

    return returnList

def fluxes_statistics(model_name: str,  output_types:List)-> List[pd.DataFrame]:
    """
    Computes statistics (mean, median, quantiles) for the fluxes.

    Args:
        model_name (str): Name of the model, used in filename for input.
        output_types (List[str]): Types of statistics to compute (mean, median, quantiles).

    Returns:
        List[pd.DataFrame]: List of DataFrames containing mean, median, and quantiles statistics.
    """

    df_mean = pd.DataFrame()
    df_median= pd.DataFrame()
    df_quantiles= pd.DataFrame()

    df_samples = pd.read_csv(ARGS.output_path + "/"  +  model_name + '.csv', sep = '\t', index_col = 0).T
    df_samples = df_samples.round(8)

    for output_type in output_types:
        if(output_type == "mean"):
            df_mean = df_samples.mean()
            df_mean = df_mean.to_frame().T
            df_mean = df_mean.reset_index(drop=True)
            df_mean.index = [model_name]
        elif(output_type == "median"):
            df_median = df_samples.median()
            df_median = df_median.to_frame().T
            df_median = df_median.reset_index(drop=True)
            df_median.index = [model_name]
        elif(output_type == "quantiles"):
            newRow = []
            cols = []
            for rxn in df_samples.columns:
                quantiles = df_samples[rxn].quantile([0.25, 0.50, 0.75])
                newRow.append(quantiles[0.25])
                cols.append(rxn + "_q1")
                newRow.append(quantiles[0.5])
                cols.append(rxn + "_q2")
                newRow.append(quantiles[0.75])
                cols.append(rxn + "_q3")
            df_quantiles = pd.DataFrame(columns=cols)
            df_quantiles.loc[0] = newRow
            df_quantiles = df_quantiles.reset_index(drop=True)
            df_quantiles.index = [model_name]
    
    return df_mean, df_median, df_quantiles

def fluxes_analysis(model:cobra.Model,  model_name:str, output_types:List)-> List[pd.DataFrame]:
    """
    Performs flux analysis including pFBA, FVA, and sensitivity analysis. The objective function
    is assumed to be already set in the model.

    Args:
        model (cobra.Model): The COBRA model to analyze.
        model_name (str): Name of the model, used in filenames for output.
        output_types (List[str]): Types of analysis to perform (pFBA, FVA, sensitivity).

    Returns:
        List[pd.DataFrame]: List of DataFrames containing pFBA, FVA, and sensitivity analysis results.
    """

    df_pFBA = pd.DataFrame()
    df_FVA= pd.DataFrame()
    df_sensitivity= pd.DataFrame()

    for output_type in output_types:
        if(output_type == "pFBA"):
            solution = cobra.flux_analysis.pfba(model)
            fluxes = solution.fluxes
            df_pFBA.loc[0,[rxn.id for rxn in model.reactions]] = fluxes.tolist()
            df_pFBA = df_pFBA.reset_index(drop=True)
            df_pFBA.index = [model_name]
            df_pFBA = df_pFBA.astype(float).round(6)
        elif(output_type == "FVA"):
            fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=ARGS.perc_opt, processes=1).round(8)
            columns = []
            for rxn in fva.index.to_list():
                columns.append(rxn + "_min")
                columns.append(rxn + "_max")
            df_FVA= pd.DataFrame(columns = columns)
            for index_rxn, row in fva.iterrows():
                df_FVA.loc[0, index_rxn+ "_min"] = fva.loc[index_rxn, "minimum"]
                df_FVA.loc[0, index_rxn+ "_max"] = fva.loc[index_rxn, "maximum"]
            df_FVA = df_FVA.reset_index(drop=True)
            df_FVA.index = [model_name]
            df_FVA = df_FVA.astype(float).round(6)
        elif(output_type == "sensitivity"):
            solution_original = model.optimize().objective_value
            reactions = model.reactions
            single = cobra.flux_analysis.single_reaction_deletion(model)
            newRow = []
            df_sensitivity = pd.DataFrame(columns = [rxn.id for rxn in reactions], index = [model_name])
            for rxn in reactions:
                newRow.append(single.knockout[rxn.id].growth.values[0]/solution_original)
            df_sensitivity.loc[model_name] = newRow
            df_sensitivity = df_sensitivity.astype(float).round(6)
    return df_pFBA, df_FVA, df_sensitivity

############################# main ###########################################
def main(args: List[str] = None) -> None:
    """
    Initialize and run sampling/analysis based on the frontend input arguments.

    Returns:
        None
    """

    num_processors = max(1, cpu_count() - 1)

    global ARGS
    ARGS = process_args(args)

    if not os.path.exists('flux_simulation'):
        os.makedirs('flux_simulation')

    # --- Read input files and names from the provided file paths ---
    with open(ARGS.input_file, 'r') as f:
        ARGS.input_files = [line.strip() for line in f if line.strip()]
    
    with open(ARGS.name_file, 'r') as f:
        ARGS.file_names = [line.strip() for line in f if line.strip()]
    
    # output types (required) -> list
    ARGS.output_types = ARGS.output_type.split(",") if ARGS.output_type else []
    # optional analysis output types -> list or empty
    ARGS.output_type_analysis = ARGS.output_type_analysis.split(",") if ARGS.output_type_analysis else []

    # Determine if sampling should be performed
    if ARGS.sampling_enabled == "true":
        perform_sampling = True
    else:
        perform_sampling = False

    print("=== INPUT FILES ===")
    print(f"{ARGS.input_files}")
    print(f"{ARGS.file_names}")
    print(f"{ARGS.output_type}")
    print(f"{ARGS.output_types}")
    print(f"{ARGS.output_type_analysis}")
    print(f"Sampling enabled: {perform_sampling} (n_samples: {ARGS.n_samples})")
    
    if ARGS.model_and_bounds == "True":
        # MODE 1: Model + bounds (separate files)
        print("=== MODE 1: Model + Bounds (separate files) ===")
        
        # Load base model
        if not ARGS.model_upload:
            sys.exit("Error: model_upload is required for Mode 1")

        # Process each bounds file with the base model
        results = Parallel(n_jobs=num_processors)(
            delayed(model_sampler_with_bounds)(ARGS.model_upload, bounds_file, cell_name) 
            for bounds_file, cell_name in zip(ARGS.input_files, ARGS.file_names)
        )

    else:
        # MODE 2: Multiple complete models
        print("=== MODE 2: Multiple complete models ===")
        
        # Process each complete model file
        results = Parallel(n_jobs=num_processors)(
            delayed(perform_sampling_and_analysis)(model_file, cell_name) 
            for model_file, cell_name in zip(ARGS.input_files, ARGS.file_names)
        )

    # Handle sampling outputs (only if sampling was performed)
    if perform_sampling:
        print("=== PROCESSING SAMPLING RESULTS ===")
        
        all_mean = pd.concat([result[0] for result in results], ignore_index=False)
        all_median = pd.concat([result[1] for result in results], ignore_index=False)
        all_quantiles = pd.concat([result[2] for result in results], ignore_index=False)

        if "mean" in ARGS.output_types:
            all_mean = all_mean.fillna(0.0)
            all_mean = all_mean.sort_index()
            write_to_file(all_mean.T, ARGS.out_mean, True)

        if "median" in ARGS.output_types:
            all_median = all_median.fillna(0.0)
            all_median = all_median.sort_index()
            write_to_file(all_median.T, ARGS.out_median, True)
        
        if "quantiles" in ARGS.output_types:
            all_quantiles = all_quantiles.fillna(0.0)
            all_quantiles = all_quantiles.sort_index()
            write_to_file(all_quantiles.T, ARGS.out_quantiles, True)
    else:
        print("=== SAMPLING SKIPPED (n_samples = 0 or sampling disabled) ===")

    # Handle optimization analysis outputs (always available)
    print("=== PROCESSING OPTIMIZATION RESULTS ===")
    
    # Determine the starting index for optimization results
    # If sampling was performed, optimization results start at index 3
    # If no sampling, optimization results start at index 0
    index_result = 3 if perform_sampling else 0
    
    if "pFBA" in ARGS.output_type_analysis:
        all_pFBA = pd.concat([result[index_result] for result in results], ignore_index=False)
        all_pFBA = all_pFBA.sort_index()
        write_to_file(all_pFBA.T, ARGS.out_pfba, True)
        index_result += 1
        
    if "FVA" in ARGS.output_type_analysis:
        all_FVA = pd.concat([result[index_result] for result in results], ignore_index=False)
        all_FVA = all_FVA.sort_index()
        write_to_file(all_FVA.T, ARGS.out_fva, True)
        index_result += 1
        
    if "sensitivity" in ARGS.output_type_analysis:
        all_sensitivity = pd.concat([result[index_result] for result in results], ignore_index=False)
        all_sensitivity = all_sensitivity.sort_index()
        write_to_file(all_sensitivity.T, ARGS.out_sensitivity, True)

    return
        
##############################################################################
if __name__ == "__main__":
    main()