import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple
import sys
import math
import argparse
import pickle as pk

########################## argparse ##########################################
def process_args(args :List[str]) -> argparse.Namespace:
    """
    Processes command-line arguments.

    Args:
        args (list): List of command-line arguments.

    Returns:
        Namespace: An object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(usage = '%(prog)s [options]',
                                     description = 'process some value\'s'+
                                     ' abundances and reactions to create RPS scores.')
    parser.add_argument('-rc', '--reaction_choice', 
                        type = str,
                        default = 'default',
                        choices = ['default','custom'], 
                        help = 'chose which type of reaction dataset you want use')
    parser.add_argument('-cm', '--custom',
                        type = str,
                        help='your dataset if you want custom reactions')
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")    
    parser.add_argument('-id', '--input',
                        type = str,
                        help = 'input dataset')
    parser.add_argument('-rp', '--rps_output',
                        type = str,
                        required = True,
                        help = 'ras output')
    
    args = parser.parse_args()
    return args



############################ read_dataset ####################################
def read_dataset(data :str, name :str) -> pd.DataFrame:
    """
    Tries to read the dataset from its path (data) as a tsv and turns it into a DataFrame.

    Args:
        data : filepath of a dataset (from frontend input params or literals upon calling)
        name : name associated with the dataset (from frontend input params or literals upon calling)

    Returns:
        pd.DataFrame : dataset in a runtime operable shape
    
    Raises:
        sys.exit : if there's no data (pd.errors.EmptyDataError) or if the dataset has less than 2 columns
    """
    try:
        dataset = pd.read_csv(data, sep = '\t', header = 0, engine='python')
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    return dataset


############################ dataset name #####################################
def name_dataset(name_data :str, count :int) -> str:
    """
    Produces a unique name for a dataset based on what was provided by the user. The default name for any dataset is "Dataset", thus if the user didn't change it this function appends f"_{count}" to make it unique.

    Args:
        name_data : name associated with the dataset (from frontend input params)
        count : counter from 1 to make these names unique (external)

    Returns:
        str : the name made unique
    """
    if str(name_data) == 'Dataset':
        return str(name_data) + '_' + str(count)
    else:
        return str(name_data)
    

############################ get_abund_data ####################################
def get_abund_data(dataset: pd.DataFrame, cell_line_index:int) -> Optional[pd.Series]:
    """
    Extracts abundance data and turns it into a series for a specific cell line from the dataset.

    Args:
        dataset (pandas.DataFrame): The DataFrame containing abundance data for all cell lines and metabolites.
        cell_line_index (int): The index of the cell line of interest in the dataset.

    Returns:
        pd.Series or None: A series containing abundance values for the specified cell line.
                           The name of the series is the name of the cell line.
                           Returns None if the cell index is invalid.
    """
    if cell_line_index < 0 or cell_line_index >= len(dataset.index):
        print(f"Errore: L'indice della cell line '{cell_line_index}' non è valido.")
        return None

    cell_line_name = dataset.iloc[cell_line_index, 0] 
    abundances_series = dataset.iloc[cell_line_index, 1:]  
    abundances_series.name = cell_line_name

    return abundances_series 


############################ update_metabolite_names ####################################
def update_metabolite_names(abundances_dict: Dict[str, float], syn_dict: Dict[str, List[str]]) -> Dict[str, float]:
    """
    Update metabolite names in the abundance series based on synonyms provided in the synonym dictionary.

    Args:
        abundances_series (pandas.Series): A series containing metabolite names as index and abundance values as values.
        syn_dict (dict): A dictionary where keys are general metabolite names and values are lists of possible synonyms.

    Returns:
        dict: An updated series where metabolite names have been replaced with their corresponding general names
              according to the synonym dictionary. If a metabolite name doesn't have a synonym it is deleted, whereas 
              if it has already the general name, it remains unchanged.
    """
    updated_abundances = abundances_dict.copy()
    for metabolite, abundance in list(abundances_dict.items()):
        for key, synonyms in syn_dict.items():
            if metabolite in synonyms or metabolite in key:
                updated_abundances[key] = updated_abundances.pop(metabolite)
                break
        else:
            del updated_abundances[metabolite]

    return updated_abundances


############################ check_missing_metab ####################################
def check_missing_metab(reactions: Dict[str, Dict[str, int]], updated_abundances: Dict[str, float]) -> Tuple[Dict[str, float], List[str]]:
    """
    Check for missing metabolites in the abundances dictionary compared to the reactions dictionary and update abundances accordingly.

    Parameters:
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        updated_abundances (dict): A dictionary representing abundances where keys are metabolite names and values are their corresponding abundances.

    Returns:
        tuple[dict, list]: A tuple containing:
            - A dictionary representing updated abundances where metabolite names have been added or set to 1 if missing in the original abundances dictionary.
            - A list of metabolite names that were missing in the original abundances dictionary and thus their aboundances were set to 1.

    Side effects:
        The input updated_abundances is modified in place      
    """
    missing_list=[]
    for reaction in reactions.values():
        for metabolite in reaction.keys():
          if metabolite not in updated_abundances:
            updated_abundances[metabolite]= 1
            missing_list.append(metabolite)

    return updated_abundances, missing_list


############################ calculate_rps ####################################
def calculate_rps(reactions: Dict[str, Dict[str, int]], abundances: Dict[str, float], black_list: List[str], missing_list: List[str]) -> Dict[str, float]:
    """
    Calculate the Reaction Propensity scores (RPS) based on the availability of reaction substrates, for (ideally) each input model reaction and for each sample.
    The score is computed as the product of the concentrations of the reacting substances, with each concentration raised to a power equal to its stoichiometric coefficient
    for each reaction using the provided coefficient and abundance values.

    Parameters:
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        abundances (dict): A dictionary representing metabolite abundances where keys are metabolite names and values are their corresponding abundances.
        black_list (list): A list containing metabolite names that should be excluded from the RPS calculation.
        missing_list (list): A list containing metabolite names that were missing in the original abundances dictionary and thus their values were set to 1.

    Returns:
        dict: A dictionary containing Reaction Propensity Scores (RPS) where keys are reaction names and values are the corresponding RPS scores.
    """
    rps_scores = {}
 
    for reaction_name, substrates in reactions.items():
        total_contribution = 0
        metab_significant = False
        for metabolite, stoichiometry in substrates.items():
            temp = 1 if math.isnan(abundances[metabolite]) else abundances[metabolite]
            if metabolite not in black_list and metabolite not in missing_list:
              metab_significant = True
            total_contribution *= temp ** stoichiometry
        
        rps_scores[reaction_name] = total_contribution if metab_significant else math.nan
    
    return rps_scores


############################ rps_for_cell_lines ####################################
def rps_for_cell_lines(dataframe: pd.DataFrame, reactions: Dict[str, Dict[str, int]], black_list: List[str], syn_dict: Dict[str, List[str]], file: str, flag: bool) -> None:
    """
    Calculate Reaction Propensity Scores (RPS) for each cell line represented in the dataframe and creates an output file.

    Parameters:
        dataframe (pandas.DataFrame): A DataFrame containing metabolite abundances for different cell lines.
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        black_list (list): A list containing metabolite names that should be excluded from the RPS calculation.
        syn_dict (dict): A dictionary where keys are general metabolite names and values are lists of possible synonyms.
        file (str): Path to the output RPS file.
        flag(bool): True if default reaction dict is being used, False otherwise.

    Returns:
        None
    """
    rps_scores=[]
    
    for series in dataframe.iterrows():
        updated_abundances = update_metabolite_names(series.to_dict(), syn_dict) if flag else series.to_dict()
        abundances, missing_list = check_missing_metab(reactions, updated_abundances)
        rps_scores.append(calculate_rps(reactions, abundances, black_list, missing_list))
     

    text_file = open(file, "w")
    text_file.write(rps_scores)
    text_file.close()


############################ main ####################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """

#TODO: sistemare i nomi dei path in modo che corrispondano a dei file esistenti

    args = process_args(sys.argv)

    with open(args.tool_dir + '/local/black_list.p', 'rb') as bl:
        black_list = pk.load(bl)

    with open(args.tool_dir + '/local/syn_dict.p', 'rb') as sd:
        syn_dict = pk.load(sd)


    flag = True
    if args.reactions_selector == 'default':        
        reactions = pk.load(open(args.tool_dir + '/local/default_reactions.p', 'rb'))
    elif args.rules_selector == 'Custom':
        #custom_reactions = parse_custom_rules(args.custom)   
        flag = False
    
    name = "RPS Dataset"
    dataset = read_dataset(args.input, "dataset")
        
    if args.rules_selector != 'custom':
        rps_for_cell_lines(dataset, reactions, black_list, syn_dict, name, flag)
    #elif args.rules_selector == 'Custom':
        #rps_for_cell_lines(dataset, custom_reactions, black_list, syn_dict, name, flag)
 
    print('Execution succeded')
    return None


##############################################################################
if __name__ == "__main__":
    main()