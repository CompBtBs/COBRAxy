import numpy as np
import pandas as pd
from typing import Union, Optional, List, Dict, Tuple, TypeVar
import sys


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

    # Converti i valori non numerici o mancanti (NaN) a zero
    abundances_series = pd.to_numeric(abundances_series, errors='coerce').fillna(0)
    abundances_series.name = cell_line_name

    return abundances_series 


############################ update_metabolite_names ####################################
def update_metabolite_names(abundances_series: pd.Series, syn_dict: Dict[str, List[str]]) -> Optional[pd.Series]:
    """
    Update metabolite names in the abundance series based on synonyms provided in the synonym dictionary.

    Args:
        abundances_series (pandas.Series): A series containing metabolite names as index and abundance values as values.
        syn_dict (dict): A dictionary where keys are general metabolite names and values are lists of possible synonyms.

    Returns:
        pd.Series or None : An updated series where metabolite names have been replaced with their corresponding general names
                            according to the synonym dictionary. If a metabolite name doesn't have a synonym or it has already the general name, it remains unchanged.
                            Returns None if abundances_series is None.
    """
    if abundances_series is None:
        print("Error: abundances_series cannot be None.")
        return None

    updated_abundances = abundances_series.copy()

    for metabolite, abundance in abundances_series.items():
        for key, synonyms in syn_dict.items():
            if metabolite in synonyms:
                updated_abundances = updated_abundances.rename(index={metabolite: key})
                break
        else:
            updated_abundances = updated_abundances.drop(index=metabolite)

    return updated_abundances


############################ sort_reactions ####################################
def sort_reactions(reactions: Dict[str, Dict[str, int]], updated_abundances: pd.Series) -> Optional[Dict[str, Dict[str, int]]]:
    """
    Sort reactions based on the order of metabolites in the updated abundance series.

    Args:
        reactions (dict): A dictionary containing reactions data.
                          Keys are metabolite names and values are dictionaries representing stoichiometric coefficients.
        updated_abundances (pandas.Series): A pandas Series containing updated metabolite abundances.
                                             Metabolite names in the series will determine the order of reactions.

    Returns:
        dict or None: An updated dictionary of reactions that are sorted based on the order in which metabolites are in the series.
                      If the metabolite in a reaction is not found in the updated abundance series, the reaction is not included in the updated dictionary.
                      Returns None if updated_abundances is None

    """
    if updated_abundances is None:
       print("Error: abundances_series cannot be None.")
       return None

    ordered_reactions = {}

    for reaction_id in updated_abundances.index:
        if reaction_id in reactions:
            ordered_reactions[reaction_id] = reactions[reaction_id]

    return ordered_reactions


############################ extract_coeff_vector ####################################
def extract_coeff_vector(ordered_reactions: Dict[str, Dict[str, int]], reaction_id: str) -> Optional[List[int]]:
    """
    Extract the stoichiometric coefficient vector for a specific reaction index from the ordered reactions.

    Args:
        ordered_reactions (dict): A dictionary containing reaction information.
                                  Keys are reaction IDs and values are dictionaries with metabolite IDs and coefficients.
        reaction_id (str): The name of the reaction for which to extract the coefficient vector.

    Returns:
        list or None: A list containing the coefficient values for the specified reaction index.
                      Returns None if no reaction is found at the given index.
    """
    if reaction_id in ordered_reactions:
        coeff_vector = list(ordered_reactions[reaction_id].values())
        return coeff_vector
    else:
        print(f"No reaction found for '{reaction_id}'.")
        return None


############################ extract_abundance_vector ####################################
def extract_abundance_vector(updated_abundances: pd.Series) -> Optional[List[float]]:
    """
    Extract the abundance vector from an abundance series.

    Args:
        updated_abundances (pandas.Series): A pandas Series containing abundance values for metabolites.

    Returns:
        list or None: A list containing the abundance values extracted from the updated abundance series.
                      Returns None if updated_abundances is None
    """
    if updated_abundances is None:
       print("Error: abundances_series cannot be None.")
       return None

    return updated_abundances.tolist()


############################ rps_calculator ####################################
def rps_calculator(coeff_vector: List[int], abundances_vector: List[float]) -> Optional[List[float]]:
    """
    Calculate the Reaction Propensity scores (RPS) based on the availability of reaction substrates, for (ideally) each input model reaction and for each sample. 
    The score is computed as the product of the concentrations of the reacting substances, with each concentration raised to a power equal to its stoichiometric coefficient 
    for each reaction using the provided coefficient and abundance vectors.

    Args:
        coeff_vector (list): A list containing the coefficient values for each reaction.
        abundances_vector (list): A list containing the abundance values for each metabolite.

    Returns:
        list or None: A list containing the RPS scores calculated for each reaction.
                      Returns None if either the coefficient vector or the abundance vector is None.

    """
    if coeff_vector is None or abundances_vector is None:
        print("Error: input vectors cannot be None.")
        return None
    
    rps_scores = []

    for coeff, abundance in zip(coeff_vector, abundances_vector):
        if abundance == 0:
            rps_scores.append(0)
        else:
            rps_score = np.prod(abundance ** coeff)
            rps_scores.append(rps_score)
    
    return rps_scores



############################ calculate_rps_for_reactions ####################################
def calculate_rps_for_reactions(reactions: Dict[str, Dict[str, int]], abundances_vector: List[float]) -> Optional[List[List[float]]]:
    """
    Calculate RPS scores for each reaction based on the provided reactions and updated abundances.

    Args:
        reactions (dict): A dictionary containing reactions data.
                          Keys are reaction IDs or indices, and values are dictionaries representing stoichiometric coefficients.
        abundances_vector (list): A list containing updated metabolite abundances.

    Returns:
        list: A list of RPS scores for each reaction.
    """
    rps_scores = []

    for reaction_id, stoichiometry in reactions.items():
        coeff_vector = extract_coeff_vector(reactions, reaction_id)

        if coeff_vector is not None:
            rps_score = rps_calculator(coeff_vector, abundances_vector)
            rps_scores.append(rps_score)
        else:
            rps_scores.append(None)

    return rps_scores
############################ MAIN #############################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    