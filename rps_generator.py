import re
import sys
import csv
import math
import argparse

import numpy  as np
import pickle as pk
import pandas as pd

from enum   import Enum
from typing import Optional, List, Dict, Tuple

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
        dataset = pd.read_csv(data, sep = ',', header = 0, engine='python')
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
    updated_abundances = {}
    for name, abundance in abundances_dict.items():
        id = get_metabolite_id(name, syn_dict)
        if id: updated_abundances[id] = abundance

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
        total_contribution = 1
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
    
    for (_, series) in dataframe.iterrows():
        updated_abundances = update_metabolite_names(series.to_dict(), syn_dict) if flag else series.to_dict()
        abundances, missing_list = check_missing_metab(reactions, updated_abundances)
        rps_scores.append(calculate_rps(reactions, abundances, black_list, missing_list))
    
    output_ras = pd.DataFrame.from_dict(resolve_rules)
    
    output_ras.insert(0, 'Reactions', ids)
    output_to_csv = pd.DataFrame.to_csv(output_ras, sep = '\t', index = False)
    
    text_file = open(file, "w")
    
    text_file.write(output_to_csv)
    text_file.close()

    with open(file, "w") as file:
        file.write(rps_scores)



############################ parse_custom_reactions ####################################
# Reaction direction encoding:
class ReactionDir(Enum):
  """
  A reaction can go forwards, backwards or be reversible (able to proceed in both directions).
  Models created / managed with cobrapy encode this information within the reaction's
  formula using the arrows this enum keeps as values.
  """
  FORWARD    = "-->"
  BACKWARD   = "<--"
  REVERSIBLE = "<=>"

  @classmethod
  def fromReaction(cls, reaction :str) -> 'ReactionDir':
    """
    Takes a whole reaction formula string and looks for one of the arrows, returning the
    corresponding reaction direction.

    Args:
      reaction : the reaction's formula.
    
    Raises:
      ValueError : if no valid arrow is found.
    
    Returns:
      ReactionDir : the corresponding reaction direction.
    """
    for member in cls:
      if member.value in reaction: return member

    raise ValueError("No valid arrow found within reaction string.")

ReactionsDict = Dict[str, Dict[str, float]]
def add_custom_reaction(reactionsDict :ReactionsDict, rId :str, reaction :str) -> None:
  """
  Adds an entry to the given reactionsDict. Each entry consists of a given unique reaction id
  (key) and a :dict (value) matching each substrate in the reaction to its stoichiometric coefficient.
  Keys and values are both obtained from the reaction's formula: if a substrate (custom metabolite id)
  appears without an explicit coeff, the value 1.0 will be used instead.

  Args:
    reactionsDict : dictionary encoding custom reactions information.
    rId : unique reaction id.
    reaction : the reaction's formula.
  
  Returns:
    None

  Side effects:
    reactionsDict : mut
  """
  reaction = reaction.strip()
  if not reaction: return

  reactionsDict[rId] = {}
  # We assume the '+' separating consecutive metabs in a reaction is spaced from them,
  # to avoid confusing it for electrical charge:
  for word in reaction.split(" + "):
    metabId, stoichCoeff = word, 1.0
    # Implicit stoichiometric coeff is equal to 1, some coeffs are floats.

    # Accepted coeffs can be integer or floats with a dot (.) decimal separator
    # and must be separated from the metab with a space:
    foundCoeff = re.search(r"\d+(\.\d+)? ", word)
    if foundCoeff:
      wholeMatch  = foundCoeff.group(0)
      metabId     = word[len(wholeMatch) + 1:].strip()
      stoichCoeff = float(wholeMatch.strip())

    reactionsDict[rId][metabId] = stoichCoeff

  if not reactionsDict[rId]: del reactionsDict[rId] # Empty reactions are removed.

def parse_custom_reactions(customReactionsPath :str) -> ReactionsDict:
  """
  Creates a custom dictionary encoding reactions information from a csv file containing
  data about these reactions, the path of which is given as input.

  Args:
    customReactionsPath : path to the reactions information file.
  
  Returns:
    ReactionsDict : dictionary encoding custom reactions information.
  """
  reactionsData :Dict[str, str]
  with open(customReactionsPath, "r") as fd:
    # We expect 2 columns: the first containing reaction ids and another with the actual reactions.
    reactionsData = { row[0] : row[1] for row in csv.reader(fd) }

  reactionsDict :ReactionsDict = {}
  for rId, reaction in reactionsData.items():
    reactionDir = ReactionDir.fromReaction(reaction)
    left, right = reaction.split(f" {reactionDir.value} ")

    # Reversible reactions are split into distinct reactions, one for each direction.
    # In general we only care about substrates, the product information is lost.
    reactionIsReversible = reactionDir is ReactionDir.REVERSIBLE
    if reactionDir is not ReactionDir.BACKWARD:
      add_custom_reaction(reactionsDict, rId + "_F" * reactionIsReversible, left)
    
    if reactionDir is not ReactionDir.FORWARD:
      add_custom_reaction(reactionsDict, rId + "_B" * reactionIsReversible, right)

  return reactionsDict

def parse_custom_reactions_old(reactions_csv_path: str) -> Dict[str, List[str]]:
    """
    Parses custom reactions from a CSV file and generates a dictionary representing the reactions as keys while the values are lists of metabolite
    that take part in those reactions as substrates.

    Parameters:
        reactions_csv_path (str): The file path to the CSV file containing reaction data.

    Returns:
        dict: A dictionary representing the parsed reactions. Keys are reaction IDs and values are lists of metabolites involved in each reaction.
    """
    reaction_dict = {}
    ban_list = ['+', '-->', '<--', '<=>']
    with open(reactions_csv_path, newline='') as csvfile:
      reader = csv.reader(csvfile)
      next(reader, None)

      for row in reader:
          reaction_id = row[1]
          reaction_column = row[2].strip()
          metabolites = row[2].strip()

          if '-->' in reaction_column:
            arrow_index = metabolites.index('-->') if '-->' in metabolites else None
            if arrow_index is not None:
                  left_metabolites = metabolites[:arrow_index]
                  reaction_dict[reaction_id] = [m for m in left_metabolites.split() if m not in ban_list and not m.replace('.', '').isdigit()]

          elif '<---' in reaction_column:
              arrow_index = metabolites.index('<---') if '<---' in metabolites else None
              if arrow_index is not None:
                  right_metabolites = metabolites[arrow_index + 1:]
                  reaction_dict[reaction_id] = [m for m in right_metabolites.split() if m not in ban_list and not m.replace('.', '').isdigit()]
          elif '<=>' in reaction_column:
              reaction_dict[reaction_id] = [m for m in metabolites.split() if m not in ban_list and not m.replace('.', '').isdigit()]
  
    return reaction_dict


############################ count_reaction_number ####################################
def count_reaction_number(reactions_csv_path: str) -> int:
    """
    Counts the total number of reactions from a CSV file. It is assumed that the first row has to be omitted as it does not contain any data

    Parameters:
        reactions_csv_path (str): The file path to the CSV file containing reaction data.

    Returns:
        int: The total number of reactions in the CSV file.
    """
    df = pd.read_csv(reactions_csv_path, skiprows=1) 
    total_reactions_num = len(df)
    return total_reactions_num


############################ get_sorted_dict ####################################
def get_sorted_dict(reaction_dict: Dict[str, List[str]]) -> Dict[str, int]:
    """
    Generates a sorted dictionary based on the count of metabolites in a reaction dictionary.

    Parameters:
        reaction_dict (dict): A dictionary representing reactions, where keys are reaction IDs and values are lists of metabolites involved in each reaction.

    Returns:
        dict: A sorted dictionary where keys are metabolites and values are their counts across all reactions.
    """
    count_dict = {}
    for metabolite_list in reaction_dict.values():
        for metabolite in metabolite_list:
            if metabolite not in count_dict:
                count_dict[metabolite] = 1
            else:
                count_dict[metabolite] += 1

    sorted_dict = dict(sorted(count_dict.items(), key=lambda x: x[1], reverse=True))
    return sorted_dict


############################ get_top_metabolites ####################################
def get_percentages_list(sorted_dict: Dict[str, int], total_reactions_num: int) -> List[float]:
    """
    Calculates the percentages of occurrence for each metabolite based on their counts in a sorted dictionary and the total amount of reactions in the model.

    Parameters:
        sorted_dict (dict): A dictionary representing sorted metabolite counts, where keys are metabolites and values are their counts.
        total_reactions_num (int): The total number of reactions.

    Returns:
        list: A list containing the percentages of occurrence for each metabolite.
    """
    for key in sorted_dict:
        sorted_dict[key] /= total_reactions_num

    percentages = list(sorted_dict.values())
    return percentages 


############################ main ####################################
def get_black_list(percentages: List[float], sorted_dict: Dict[str, int]) -> List[str]:
    """
    Generates a black list of metabolites based on the percentage occurrences and sorted metabolite counts.

    Parameters:
        percentages (list): A list containing the percentages of occurrence for each metabolite.
        sorted_dict (dict): A dictionary representing sorted metabolite counts, where keys are metabolites and values are their counts.

    Returns:
        list: A black list containing metabolites identified as outliers based on log-scaled percentages.
    """
    log_scale_percentages = np.log(percentages)

    data = np.array(log_scale_percentages)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    threshold = 1.5 * iqr
    outliers = np.where((data < q1 - threshold) | (data > q3 + threshold))

    black_list = []
    black_list =[list(sorted_dict.keys())[outlier] for outlier in outliers[0]] 
    return black_list


############################ main ####################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    args = process_args(sys.argv)

    with open(args.tool_dir + '/local/pickle files/black_list.pickle', 'rb') as bl:
        black_list = pk.load(bl)

    with open(args.tool_dir + '/local/pickle files/synonyms.pickle', 'rb') as sd:
        syn_dict = pk.load(sd)

    name = "RPS Dataset"
    dataset = read_dataset(args.input, "dataset")

    flag = True
    if args.reaction_choice == 'default':        
        reactions = pk.load(open(args.tool_dir + '/local/pickle files/reactions.pickle', 'rb'))
    elif args.reaction_choice == 'Custom':
        reactions = parse_custom_reactions(args.custom)   
        count_reaction_number(args.custom)
        sorted_dict = get_sorted_dict(reactions)
        percentages = get_percentages_list(sorted_dict, count_reaction_number(args.custom))
        custom_black_list = get_black_list(percentages, sorted_dict)
        black_list = custom_black_list
        flag = False
      
    rps_for_cell_lines(dataset, reactions, black_list, syn_dict, name, flag)
    
    

    print('Execution succeded')
    return None


##############################################################################
if __name__ == "__main__":
    main()