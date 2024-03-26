import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple
import sys
import math


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
                           Non-numeric or missing (NaN) values are converted to zero.
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
            total_contribution += temp * stoichiometry
        
        rps_scores[reaction_name] = total_contribution if metab_significant else math.nan
    
    return rps_scores


############################ rps_for_cell_lines ####################################
def rps_for_cell_lines(dataframe: pd.DataFrame, reactions: Dict[str, dict[str, int]], black_list: List[str], syn_dict: Dict[str, List[str]]) -> List[Dict[str, float]]:
    """
    Calculate Reaction Propensity Scores (RPS) for each cell line represented in the dataframe.

    Parameters:
        dataframe (pandas.DataFrame): A DataFrame containing metabolite abundances for different cell lines.
        reactions (dict): A dictionary representing reactions where keys are reaction names and values are dictionaries containing metabolite names as keys and stoichiometric coefficients as values.
        black_list (list): A list containing metabolite names that should be excluded from the RPS calculation.
        syn_dict (dict): A dictionary where keys are general metabolite names and values are lists of possible synonyms.

    Returns:
        list[dict]: A list of dictionaries containing Reaction Propensity Scores (RPS) for each cell line. Each dictionary contains reaction names as keys and their corresponding RPS scores as values.
    
    """
    rps_scores=[]
    for series in dataframe.rows():
     updated_abundances = update_metabolite_names(series.to_dict(), syn_dict)
     abundances, missing_list = check_missing_metab(reactions, updated_abundances)
     rps_scores.append(calculate_rps(reactions, abundances, black_list, missing_list))
    return rps_scores


############################ main ####################################