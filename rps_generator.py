import math
import argparse

import numpy  as np
import pickle as pk
import pandas as pd

from typing import Optional, List, Dict

import utils.general_utils as utils
import utils.reaction_parsing as reactionUtils

########################## argparse ##########################################
ARGS :argparse.Namespace
def process_args(args:List[str] = None) -> argparse.Namespace:
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
                        required = True,
                        help = 'input dataset')
    parser.add_argument('-rp', '--rps_output',
                        type = str,
                        required = True,
                        help = 'rps output')
    
    args = parser.parse_args(args)
    return args

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
    Extracts abundance data and turns it into a series for a specific cell line from the dataset, which rows are
    metabolites and columns are cell lines.

    Args:
        dataset (pandas.DataFrame): The DataFrame containing abundance data for all cell lines and metabolites.
        cell_line_index (int): The index of the cell line of interest in the dataset.

    Returns:
        pd.Series or None: A series containing abundance values for the specified cell line.
                           The name of the series is the name of the cell line.
                           Returns None if the cell index is invalid.
    """
    if cell_line_index < 0 or cell_line_index >= len(dataset.index):
        print(f"Errore: This cell line index: '{cell_line_index}' is not valid.")
        return None

    cell_line_name = dataset.columns[cell_line_index]
    abundances_series = dataset[cell_line_name][1:]

    return abundances_series

############################ clean_metabolite_name ####################################
def clean_metabolite_name(name :str) -> str:
    """
    Removes some characters from a metabolite's name, provided as input, and makes it lowercase in order to simplify
    the search of a match in the dictionary of synonyms.

    Args:
        name : the metabolite's name, as given in the dataset.
    
    Returns:
        str : a new string with the cleaned name.
    """
    return "".join(ch for ch in name if ch not in ",;-_'([{ }])").lower()

############################ get_metabolite_id ####################################
def get_metabolite_id(name :str, syn_dict :Dict[str, List[str]]) -> str:
    """
    Looks through a dictionary of synonyms to find a match for a given metabolite's name.

    Args:
        name : the metabolite's name, as given in the dataset.
        syn_dict : the dictionary of synonyms, using unique identifiers as keys and lists of clean synonyms as values.
    
    Returns:
        str : the internal :str unique identifier of that metabolite, used in all other parts of the model in use.
        An empty string is returned if a match isn't found.
    """
    name = clean_metabolite_name(name)
    for id, synonyms in syn_dict.items():
        if name in synonyms: return id
    
    return ""

############################ check_missing_metab ####################################
def check_missing_metab(reactions: Dict[str, Dict[str, int]], dataset_by_rows: Dict[str, List[float]], cell_lines_amt :int) -> List[str]:
    """
    Check for missing metabolites in the abundances dictionary compared to the reactions dictionary and update abundances accordingly.

    Parameters:
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        dataset_by_rows (dict): A dictionary representing abundances where keys are metabolite names and values are their corresponding abundances for all cell lines.
        cell_lines_amt : amount of cell lines, needed to add a new list of abundances for missing metabolites.

    Returns:
        list[str] : list of metabolite names that were missing in the original abundances dictionary and thus their aboundances were set to 1.

    Side effects:
        dataset_by_rows : mut
    """
    missing_list = []
    for reaction in reactions.values():
        for metabolite in reaction.keys():
          if metabolite not in dataset_by_rows:
            dataset_by_rows[metabolite] = [1] * cell_lines_amt
            missing_list.append(metabolite)

    return missing_list

############################ calculate_rps ####################################
def calculate_rps(reactions: Dict[str, Dict[str, int]], abundances: Dict[str, float], black_list: List[str], missing_list: List[str], substrateFreqTable: Dict[str, int]) -> Dict[str, float]:
    """
    Calculate the Reaction Propensity scores (RPS) based on the availability of reaction substrates, for (ideally) each input model reaction and for each sample.
    The score is computed as the product of the concentrations of the reacting substances, with each concentration raised to a power equal to its stoichiometric coefficient
    for each reaction using the provided coefficient and abundance values. The value is then normalized, based on how frequent the metabolite is in the selected model's reactions,
    and log-transformed.

    Parameters:
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        abundances (dict): A dictionary representing metabolite abundances where keys are metabolite names and values are their corresponding abundances.
        black_list (list): A list containing metabolite names that should be excluded from the RPS calculation.
        missing_list (list): A list containing metabolite names that were missing in the original abundances dictionary and thus their values were set to 1.
        substrateFreqTable (dict): A dictionary where each metabolite name (key) is associated with how many times it shows up in the model's reactions (value).
        
    Returns:
        dict: A dictionary containing Reaction Propensity Scores (RPS) where keys are reaction names and values are the corresponding RPS scores.
    """
    rps_scores = {}

    for reaction_name, substrates in reactions.items():
        total_contribution = 0
        metab_significant  = False
        for metabolite, stoichiometry in substrates.items():
            abundance = 1 if math.isnan(abundances[metabolite]) else abundances[metabolite]
            if metabolite not in black_list and metabolite not in missing_list:
                metab_significant = True
            
            total_contribution += math.log((abundance + np.finfo(float).eps) / substrateFreqTable[metabolite]) * stoichiometry
        
        rps_scores[reaction_name] = total_contribution if metab_significant else math.nan
    
    return rps_scores

############################ rps_for_cell_lines ####################################
def rps_for_cell_lines(dataset: List[List[str]], reactions: Dict[str, Dict[str, int]], black_list: List[str], syn_dict: Dict[str, List[str]], substrateFreqTable: Dict[str, int]) -> None:
    """
    Calculate Reaction Propensity Scores (RPS) for each cell line represented in the dataframe and creates an output file.

    Parameters:
        dataset : the dataset's data, by rows
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        black_list (list): A list containing metabolite names that should be excluded from the RPS calculation.
        syn_dict (dict): A dictionary where keys are general metabolite names and values are lists of possible synonyms.
        substrateFreqTable (dict): A dictionary where each metabolite name (key) is associated with how many times it shows up in the model's reactions (value).

    Returns:
        None
    """
    cell_lines = dataset[0][1:]
    abundances_dict = {}

    translationIsApplied = ARGS.reaction_choice == "default"
    for row in dataset[1:]:
        id = get_metabolite_id(row[0], syn_dict) if translationIsApplied else row[0]
        if id: abundances_dict[id] = list(map(utils.Float(), row[1:]))
    
    missing_list = check_missing_metab(reactions, abundances_dict, len((cell_lines)))
    
    rps_scores :Dict[Dict[str, float]] = {}
    for pos, cell_line_name in enumerate(cell_lines):
        abundances = { metab : abundances[pos] for metab, abundances in abundances_dict.items() }
        rps_scores[cell_line_name] = calculate_rps(reactions, abundances, black_list, missing_list, substrateFreqTable)
    
    df = pd.DataFrame.from_dict(rps_scores)
    
    df.index.name = 'Reactions'
    df.to_csv(ARGS.rps_output, sep='\t', na_rep='None', index=True)

############################ main ####################################
def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    global ARGS
    ARGS = process_args(args)

    # TODO:use utils functions vvv
    with open(ARGS.tool_dir + '/local/pickle files/black_list.pickle', 'rb') as bl:
        black_list = pk.load(bl)

    with open(ARGS.tool_dir + '/local/pickle files/synonyms.pickle', 'rb') as sd:
        syn_dict = pk.load(sd)

    dataset = utils.readCsv(utils.FilePath.fromStrPath(ARGS.input), '\t', skipHeader = False)

    if ARGS.reaction_choice == 'default':
        reactions = pk.load(open(ARGS.tool_dir + '/local/pickle files/reactions.pickle', 'rb'))
        substrateFreqTable = pk.load(open(ARGS.tool_dir + '/local/pickle files/substrate_frequencies.pickle', 'rb'))
    
    elif ARGS.reaction_choice == 'custom':
        reactions = reactionUtils.parse_custom_reactions(ARGS.custom)
        substrateFreqTable = {}
        for _, substrates in reactions.items():
            for substrateName, _ in substrates.items():
                if substrateName not in substrateFreqTable: substrateFreqTable[substrateName] = 0
                substrateFreqTable[substrateName] += 1

    rps_for_cell_lines(dataset, reactions, black_list, syn_dict, substrateFreqTable)
    print('Execution succeded')

##############################################################################
if __name__ == "__main__": main()
