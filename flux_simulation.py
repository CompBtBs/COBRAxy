import argparse
import utils.general_utils as utils
from typing import Optional, List
import os
import numpy as np
import pandas as pd
import cobra
import utils.CBS_backend as CBS_backend
from joblib import Parallel, delayed, cpu_count
from cobra.sampling import OptGPSampler
import sys

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

    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")
    
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    
    parser.add_argument('-in', '--input',
                        required = True,
                        type=str,
                        help = 'inputs bounds')
    
    parser.add_argument('-ni', '--names',
                        required = True,
                        type=str,
                        help = 'cell names')
 
    parser.add_argument(
        '-ms', '--model_selector', 
        type = utils.Model, default = utils.Model.ENGRO2, choices = [utils.Model.ENGRO2, utils.Model.Custom],
        help = 'chose which type of model you want use')
    
    parser.add_argument("-mo", "--model", type = str)
    
    parser.add_argument("-mn", "--model_name", type = str, help = "custom mode name")
    
    parser.add_argument('-a', '--algorithm',
                        type = str,
                        choices = ['OPTGP', 'CBS'],
                        required = True,
                        help = 'choose sampling algorithm')
    
    parser.add_argument('-th', '--thinning', 
                        type = int,
                        default= 100,
                        required=False,
                        help = 'choose thinning')
    
    parser.add_argument('-ns', '--n_samples', 
                        type = int,
                        required = True,
                        help = 'choose how many samples')
    
    parser.add_argument('-sd', '--seed', 
                        type = int,
                        required = True,
                        help = 'seed')
    
    parser.add_argument('-nb', '--n_batches', 
                        type = int,
                        required = True,
                        help = 'choose how many batches')
    
    parser.add_argument('-ot', '--output_type', 
                        type = str,
                        required = True,
                        help = 'output type')
    
    parser.add_argument('-ota', '--output_type_analysis', 
                        type = str,
                        required = False,
                        help = 'output type analysis')
    
    parser.add_argument('-idop', '--output_path', 
                        type = str,
                        default='flux_simulation',
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


def write_to_file(dataset: pd.DataFrame, name: str, keep_index:bool=False)->None:
    dataset.index.name = 'Reactions'
    dataset.to_csv(ARGS.output_path + "/" + name + ".csv", sep = '\t', index = keep_index)

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



def OPTGP_sampler(model:cobra.Model, model_name:str, n_samples:int=1000, thinning:int=100, n_batches:int=1, seed:int=0)-> None:
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

    for i in range(0, n_batches):
        optgp = OptGPSampler(model, thinning, seed)
        samples = optgp.sample(n_samples)
        samples.to_csv(ARGS.output_path + "/" +  model_name + '_'+ str(i)+'_OPTGP.csv', index=False)
        seed+=1
    samplesTotal = pd.DataFrame()
    for i in range(0, n_batches):
        samples_batch = pd.read_csv(ARGS.output_path + "/"  +  model_name + '_'+ str(i)+'_OPTGP.csv')
        samplesTotal = pd.concat([samplesTotal, samples_batch], ignore_index = True)

    write_to_file(samplesTotal.T, model_name, True)

    for i in range(0, n_batches):
        os.remove(ARGS.output_path + "/" +   model_name + '_'+ str(i)+'_OPTGP.csv')
    pass


def CBS_sampler(model:cobra.Model, model_name:str, n_samples:int=1000, n_batches:int=1, seed:int=0)-> None:
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

    df_FVA = cobra.flux_analysis.flux_variability_analysis(model,fraction_of_optimum=0).round(6)
    
    df_coefficients = CBS_backend.randomObjectiveFunction(model, n_samples*n_batches, df_FVA, seed=seed)

    for i in range(0, n_batches):
        samples = pd.DataFrame(columns =[reaction.id for reaction in model.reactions], index = range(n_samples))
        try:
            CBS_backend.randomObjectiveFunctionSampling(model, n_samples, df_coefficients.iloc[:,i*n_samples:(i+1)*n_samples], samples)
        except Exception as e:
            utils.logWarning(
            "Warning: GLPK solver has failed for " + model_name + ". Trying with COBRA interface. Error:" + str(e),
            ARGS.out_log)
            CBS_backend.randomObjectiveFunctionSampling_cobrapy(model, n_samples, df_coefficients.iloc[:,i*n_samples:(i+1)*n_samples], 
                                                    samples)
        utils.logWarning(ARGS.output_path + "/" +  model_name + '_'+ str(i)+'_CBS.csv', ARGS.out_log)
        samples.to_csv(ARGS.output_path + "/" +  model_name + '_'+ str(i)+'_CBS.csv', index=False)

    samplesTotal = pd.DataFrame()
    for i in range(0, n_batches):
        samples_batch = pd.read_csv(ARGS.output_path + "/"  +  model_name + '_'+ str(i)+'_CBS.csv')
        samplesTotal = pd.concat([samplesTotal, samples_batch], ignore_index = True)

    write_to_file(samplesTotal.T, model_name, True)

    for i in range(0, n_batches):
        os.remove(ARGS.output_path + "/" + model_name + '_'+ str(i)+'_CBS.csv')
    pass


def model_sampler(model_input_original:cobra.Model, bounds_path:str, cell_name:str)-> List[pd.DataFrame]:
    """
    Prepares the model with bounds from the dataset and performs sampling and analysis based on the selected algorithm.

    Args:
        model_input_original (cobra.Model): The original COBRA model.
        bounds_path (str): Path to the CSV file containing the bounds dataset.
        cell_name (str): Name of the cell, used to generate filenames for output.

    Returns:
        List[pd.DataFrame]: A list of DataFrames containing statistics and analysis results.
    """

    model_input = model_input_original.copy()
    bounds_df = read_dataset(bounds_path, "bounds dataset")
    for rxn_index, row in bounds_df.iterrows():
        model_input.reactions.get_by_id(rxn_index).lower_bound = row.lower_bound
        model_input.reactions.get_by_id(rxn_index).upper_bound = row.upper_bound
    
    
    if ARGS.algorithm == 'OPTGP':
        OPTGP_sampler(model_input, cell_name, ARGS.n_samples, ARGS.thinning, ARGS.n_batches, ARGS.seed)

    elif ARGS.algorithm == 'CBS':
        CBS_sampler(model_input,  cell_name, ARGS.n_samples, ARGS.n_batches, ARGS.seed)

    df_mean, df_median, df_quantiles = fluxes_statistics(cell_name, ARGS.output_types)

    if("fluxes" not in ARGS.output_types):
        os.remove(ARGS.output_path + "/"  +  cell_name + '.csv')

    returnList = []
    returnList.append(df_mean)
    returnList.append(df_median)
    returnList.append(df_quantiles)

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
    Performs flux analysis including pFBA, FVA, and sensitivity analysis.

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
            model.objective = "Biomass"
            solution = cobra.flux_analysis.pfba(model)
            fluxes = solution.fluxes
            df_pFBA.loc[0,[rxn._id for rxn in model.reactions]] = fluxes.tolist()
            df_pFBA = df_pFBA.reset_index(drop=True)
            df_pFBA.index = [model_name]
            df_pFBA = df_pFBA.astype(float).round(6)
        elif(output_type == "FVA"):
            fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0, processes=1).round(8)
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
            model.objective = "Biomass"
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
def main(args :List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """

    num_processors = cpu_count()

    global ARGS
    ARGS = process_args(args)

    if not os.path.exists(ARGS.output_path):
        os.makedirs(ARGS.output_path)
    
    model_type :utils.Model = ARGS.model_selector
    if model_type is utils.Model.Custom:
        model = model_type.getCOBRAmodel(customPath = utils.FilePath.fromStrPath(ARGS.model), customExtension = utils.FilePath.fromStrPath(ARGS.model_name).ext)
    else:
        model = model_type.getCOBRAmodel(toolDir=ARGS.tool_dir)

    #Set solver verbosity to 1 to see warning and error messages only.
    model.solver.configuration.verbosity = 1
    
    ARGS.bounds = ARGS.input.split(",")
    ARGS.bounds_name = ARGS.names.split(",")
    ARGS.output_types = ARGS.output_type.split(",")
    ARGS.output_type_analysis = ARGS.output_type_analysis.split(",")


    results = Parallel(n_jobs=num_processors)(delayed(model_sampler)(model, bounds_path, cell_name) for bounds_path, cell_name in zip(ARGS.bounds, ARGS.bounds_name))

    all_mean = pd.concat([result[0] for result in results], ignore_index=False)
    all_median = pd.concat([result[1] for result in results], ignore_index=False)
    all_quantiles = pd.concat([result[2] for result in results], ignore_index=False)

    if("mean" in ARGS.output_types):
        all_mean = all_mean.fillna(0.0)
        all_mean = all_mean.sort_index()
        write_to_file(all_mean.T, "mean", True)

    if("median" in ARGS.output_types):
        all_median = all_median.fillna(0.0)
        all_median = all_median.sort_index()
        write_to_file(all_median.T, "median", True)
    
    if("quantiles" in ARGS.output_types):
        all_quantiles = all_quantiles.fillna(0.0)
        all_quantiles = all_quantiles.sort_index()
        write_to_file(all_quantiles.T, "quantiles", True)

    index_result = 3
    if("pFBA" in ARGS.output_type_analysis):
        all_pFBA = pd.concat([result[index_result] for result in results], ignore_index=False)
        all_pFBA = all_pFBA.sort_index()
        write_to_file(all_pFBA.T, "pFBA", True)
        index_result+=1
    if("FVA" in ARGS.output_type_analysis):
        all_FVA= pd.concat([result[index_result] for result in results], ignore_index=False)
        all_FVA = all_FVA.sort_index()
        write_to_file(all_FVA.T, "FVA", True)
        index_result+=1
    if("sensitivity" in ARGS.output_type_analysis):
        all_sensitivity = pd.concat([result[index_result] for result in results], ignore_index=False)
        all_sensitivity = all_sensitivity.sort_index()
        write_to_file(all_sensitivity.T, "sensitivity", True)

    pass
        
##############################################################################
if __name__ == "__main__":
    main()