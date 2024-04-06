from __future__ import division
import sys
import pandas as pd
import collections
import pickle as pk
import math
import argparse
from typing import Union, Optional, List, Dict, Tuple, TypeVar

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
                                     ' genes to create a comparison\'s map.')
    parser.add_argument('-rs', '--rules_selector', 
                        type = str,
                        default = 'HMRcore',
                        choices = ['HMRcore', 'Recon', 'ENGRO2','Custom'], 
                        help = 'chose which type of dataset you want use')
    parser.add_argument('-cr', '--custom',
                        type = str,
                        help='your dataset if you want custom rules')
    parser.add_argument('-n', '--none',
                        type = str,
                        default = 'true',
                        choices = ['true', 'false'], 
                        help = 'compute Nan values')
    parser.add_argument('-td', '--tool_dir',
                        type = str,
                        required = True,
                        help = 'your tool directory')
    parser.add_argument('-ol', '--out_log', 
                        help = "Output log")    
    parser.add_argument('-id', '--input',
                        type = str,
                        help = 'input dataset')
    parser.add_argument('-ra', '--ras_output',
                        type = str,
                        required = True,
                        help = 'ras output')
    
    args = parser.parse_args()
    return args

########################### warning ###########################################
def warning(s :str) -> None:
    """
    Log a warning message to an output log file and print it to the console.

    Args:
        s (str): The warning message to be logged and printed.
    
    Returns:
        None
    """
    args = process_args(sys.argv)
    with open(args.out_log, 'a') as log:
            log.write(s)
            
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
    
############################ load id e rules ##################################
def load_id_rules(reactions :Dict[str, Dict[str, List[str]]]) -> Tuple[str, Dict[str, List[str]]]:
    """
    Load IDs and rules from a dictionary of reactions.

    Args:
        reactions (dict): A dictionary where keys are IDs and values are rules.

    Returns:
        tuple: A tuple containing two lists, the first list containing IDs and the second list containing rules.
    """
    ids, rules = [], []
    for key, value in reactions.items():
            ids.append(key)
            rules.append(value)
    return (ids, rules)

############################ check_methods ####################################
def gene_type(l :str, name :str) -> str:
    """
    Determine the type of gene ID.

    Args:
        l (str): The gene identifier to check.
        name (str): The name of the dataset, used in error messages.

    Returns:
        str: The type of gene ID ('hugo_id', 'ensembl_gene_id', 'symbol', or 'entrez_id').

    Raises:
        sys.exit: If the gene ID type is not supported, the execution is aborted.
    """
    if check_hgnc(l):
        return 'hugo_id'
    elif check_ensembl(l):
        return 'ensembl_gene_id'
    elif check_symbol(l):
        return 'symbol'
    elif check_entrez(l):
        return 'entrez_id'
    else:
        sys.exit('Execution aborted:\n' +
                 'gene ID type in ' + name + ' not supported. Supported ID'+
                 'types are: HUGO ID, Ensemble ID, HUGO symbol, Entrez ID\n')

def check_hgnc(l :str) -> bool:
    """
    Check if a gene identifier follows the HGNC format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the HGNC format, False otherwise.
    """
    if len(l) > 5:
        if (l.upper()).startswith('HGNC:'):
            return l[5:].isdigit()
        else:
            return False
    else:
        return False

def check_ensembl(l :str) -> bool:
    """
    Check if a gene identifier follows the Ensembl format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the Ensembl format, False otherwise.
    """
    if len(l) == 15:
        if (l.upper()).startswith('ENS'):
            return l[4:].isdigit()
        else:  
            return False 
    else: 
        return False 

def check_symbol(l :str) -> bool:
    """
    Check if a gene identifier follows the symbol format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the symbol format, False otherwise.
    """
    if len(l) > 0:
        if l[0].isalpha() and l[1:].isalnum():
            return True
        else:
            return False
    else:
        return False

def check_entrez(l :str) -> bool:
    """
    Check if a gene identifier follows the Entrez ID format.

    Args:
        l (str): The gene identifier to check.

    Returns:
        bool: True if the gene identifier follows the Entrez ID format, False otherwise.
    """ 
    if len(l) > 0:
        return l.isdigit()
    else: 
        return False

def check_bool(b :str) -> Optional[bool]:
    """
    Check if a string represents a boolean value.

    Args:
        b (str): The string to check.

    Returns:
        bool: True if the string represents True, False if the string represents False.
        None: if the input string is anything other than "true" or "false"
    """
    if b == 'true':
        return True
    elif b == 'false':
        return False
    
############################ resolve_methods ##################################
def replace_gene_value(l :str, d :str) -> Tuple[Union[int, float], list]:
    """
    Replace gene identifiers with corresponding values from a dictionary.

    Args:
        l (list): List of gene identifiers.
        d (Any): Dictionary containing gene identifiers as keys and their corresponding values.

    Returns:
        tuple: A tuple containing two lists: the first list contains replaced values, and the second list contains any errors encountered during replacement.
    """
    tmp = []
    err = []
    while l:
        if isinstance(l[0], list):
            tmp_rules, tmp_err = replace_gene_value(l[0], d)
            tmp.append(tmp_rules)
            err.extend(tmp_err)
        else:
            value = replace_gene(l[0], d)
            tmp.append(value)
            if value == None:
                err.append(l[0])
        l = l[1:]
    return (tmp, err)

def replace_gene(l :str, d :str) -> Union[int, float]:
    """
    Replace a single gene identifier with its corresponding value from a dictionary.

    Args:
        l (str): Gene identifier to replace.
        d (Any): Dictionary containing gene identifiers as keys and their corresponding values.

    Returns:
        float/int: Corresponding value from the dictionary if found, None otherwise.

    Raises:
        sys.exit: If the value associated with the gene identifier is not valid.
    """
    if l =='and' or l == 'or':
        return l
    else:
        value = d.get(l, None)
        if not(value == None or isinstance(value, (int, float))):
            sys.exit('Execution aborted: ' + value + ' value not valid\n')
        return value

T = TypeVar("T", bound = Optional[Union[int, float]])
def computes(val1 :T, op :str, val2 :T, cn :bool) -> T:
    """
    Compute the RAS value between two value and an operator ('and' or 'or').

    Args:
        val1(Optional(Union[float, int])): First value.
        op (str): Operator ('and' or 'or').
        val2(Optional(Union[float, int])): Second value.
        cn (bool): Control boolean value.

    Returns:
        Optional(Union[float, int]): Result of the computation.
    """
    if val1 != None and val2 != None:
        if op == 'and':
            return min(val1, val2)
        else:
            return val1 + val2
    elif op == 'and':
        if cn is True:
            if val1 != None:
                return val1
            elif val2 != None:
                return val2
            else:
                return None
        else:
            return None
    else:
        if val1 != None:
            return val1
        elif val2 != None:
            return val2
        else:
            return None

# ris should be Literal[None] but Literal is not supported in Python 3.7
def control(ris, l :List[Union[int, float, list]], cn :bool) -> Union[bool, int, float]: #Union[Literal[False], int, float]:
    """
    Control the format of the expression.

    Args:
        ris: Intermediate result.
        l (list): Expression to control.
        cn (bool): Control boolean value.

    Returns:
        Union[Literal[False], int, float]: Result of the control.
    """
    if len(l) == 1:
        if isinstance(l[0], (float, int)) or l[0] == None:
            return l[0]
        elif isinstance(l[0], list):
            return control(None, l[0], cn)
        else:
            return False
    elif len(l) > 2:
        return control_list(ris, l, cn)
    else:
        return False

def control_list(ris, l :List[Optional[Union[float, int, list]]], cn :bool) -> Optional[bool]: #Optional[Literal[False]]:
    """
    Control the format of a list of expressions.

    Args:
        ris: Intermediate result.
        l (list): List of expressions to control.
        cn (bool): Control boolean value.

    Returns:
        Optional[Literal[False]]: Result of the control.
    """
    while l:
        if len(l) == 1:
            return False
        elif (isinstance(l[0], (float, int)) or
              l[0] == None) and l[1] in ['and', 'or']:
            if isinstance(l[2], (float, int)) or l[2] == None:
                ris = computes(l[0], l[1], l[2], cn)            
            elif isinstance(l[2], list):
                tmp = control(None, l[2], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(l[0], l[1], tmp, cn)
            else:
                return False
            l = l[3:]
        elif l[0] in ['and', 'or']:
            if isinstance(l[1], (float, int)) or l[1] == None:
                ris = computes(ris, l[0], l[1], cn)
            elif isinstance(l[1], list):
                tmp = control(None,l[1], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(ris, l[0], tmp, cn)
            else:
                return False
            l = l[2:]
        elif isinstance(l[0], list) and l[1] in ['and', 'or']:
            if isinstance(l[2], (float, int)) or l[2] == None:
                tmp = control(None, l[0], cn)
                if tmp is False:
                    return False
                else:
                    ris = computes(tmp, l[1], l[2], cn)
            elif isinstance(l[2], list):
                tmp = control(None, l[0], cn)
                tmp2 = control(None, l[2], cn)
                if tmp is False or tmp2 is False:
                    return False
                else:
                    ris = computes(tmp, l[1], tmp2, cn)
            else:
                return False
            l = l[3:]
        else:
            return False
    return ris

############################ make recon #######################################
def check_and_doWord(l :List[str]) -> Union[bool, tuple]: #Union[Literal[False], tuple]:
    """
    Check and parse intems in the input list, removing spaces and checking brackets.

    Args:
        l (list): List of characters representing words.

    Returns:
        tuple: A tuple containing two lists: the first list contains the parsed words and operators, the second list contains only the gene identifiers.
        False: if the brackets are not balanced.
    """
    tmp = []
    tmp_genes = []
    count = 0
    while l:
        if count >= 0:
            if l[0] == '(':
                count += 1
                tmp.append(l[0])
                l.pop(0)
            elif l[0] == ')':
                count -= 1
                tmp.append(l[0])
                l.pop(0)
            elif l[0] == ' ':
                l.pop(0)
            else:
                word = []
                while l:
                    if l[0] in [' ', '(', ')']:
                        break
                    else:
                        word.append(l[0])
                        l.pop(0)
                word = ''.join(word)
                tmp.append(word)
                if not(word in ['or', 'and']):
                    tmp_genes.append(word)
        else:
            return False
    if count == 0:
        return (tmp, tmp_genes)
    else:
        return False

def brackets_to_list(l :List[str]) -> List[str]:
    """
    Convert expression inside brackets to the correct list format.

    Args:
        l (list): List representing an expression containing brackets.

    Returns:
        list: List representing the expression without brackets.
    """
    tmp = []
    while l:
        if l[0] == '(':
            l.pop(0)
            tmp.append(resolve_brackets(l))
        else:
            tmp.append(l[0])
            l.pop(0)
    return tmp

def resolve_brackets(l :List[str]) -> List[str]:
    """
    Resolve expression inside brackets.

    Args:
        l (list): List representing an expression containing brackets.

    Returns:
        list: List representing the resolved expression.
    """
    tmp = []
    while l[0] != ')':
        if l[0] == '(':
            l.pop(0)
            tmp.append(resolve_brackets(l))
        else:
            tmp.append(l[0])
            l.pop(0)
    l.pop(0)
    return tmp

def priorityAND(l :List[str]) -> List[str]:
    """
    Prioritize the 'and' operator over 'or'. It creates a boolean expression assigning priority to 'and' operator
    over the 'or'. It returns a modified list in which 'and' operators are run before the 'or' ones.

    Args:
        l (list): List representing an expression.

    Returns:
        list: List with 'and' operations having higher priority over 'or'.
    """
    tmp = []
    flag = True
    while l:
        if len(l) == 1:
            if isinstance(l[0], list):
                tmp.append(priorityAND(l[0]))
            else:
                tmp.append(l[0])
            l = l[1:]
        elif l[0] == 'or':
            tmp.append(l[0])
            flag = False
            l = l[1:]
        elif l[1] == 'or':
            if isinstance(l[0], list): 
                tmp.append(priorityAND(l[0]))
            else:
                tmp.append(l[0])
            tmp.append(l[1])
            flag = False
            l = l[2:]
        elif l[1] == 'and':
            tmpAnd = []
            if isinstance(l[0], list): 
                tmpAnd.append(priorityAND(l[0]))
            else:
                tmpAnd.append(l[0])
            tmpAnd.append(l[1])
            if isinstance(l[2], list): 
                tmpAnd.append(priorityAND(l[2]))
            else:
                tmpAnd.append(l[2])
            l = l[3:]
            while l:
                if l[0] == 'and':
                    tmpAnd.append(l[0])
                    if isinstance(l[1], list): 
                        tmpAnd.append(priorityAND(l[1]))
                    else:
                        tmpAnd.append(l[1])
                    l = l[2:]
                elif l[0] == 'or':
                    flag = False
                    break
            if flag == True: #when there are only AND in list
                tmp.extend(tmpAnd)
            elif flag == False:
                tmp.append(tmpAnd)
    return tmp

def checkRule(l :List[str]) -> bool:
    """
    Check if the expression follows the rule format.

    Args:
        l (list): List representing an expression.

    Returns:
        bool: True if the expression follows the rule format, False otherwise.
    """
    if len(l) == 1:
        if isinstance(l[0], list):
            if checkRule(l[0]) is False:
                return False
    elif len(l) > 2:
        if checkRule2(l) is False:
            return False
    else:
        return False
    return True

def checkRule2(l :List[str]) -> bool:
    """
    Check if the expression follows the rule format.

    Args:
        l (list): List representing an expression.

    Returns:
        bool: True if the expression follows the rule format, False otherwise.
    """
    while l:
        if len(l) == 1:
            return False
        elif isinstance(l[0], list) and l[1] in ['and', 'or']:
            if checkRule(l[0]) is False:
                return False
            if isinstance(l[2], list):
                if checkRule(l[2]) is False:
                    return False
            l = l[3:]
        elif l[1] in ['and', 'or']:
            if isinstance(l[2], list):
                if checkRule(l[2]) is False:
                    return False
            l = l[3:]
        elif l[0] in ['and', 'or']:
            if isinstance(l[1], list):
                if checkRule(l[1]) is False:
                    return False
            l = l[2:]
        else:
            return False
    return True

def do_rules(rules :List[str]) -> Tuple[List[str], List[str]]:
    """
    Process a list of rules.

    Args:
        rules (list): List of strings representing rules.

    Returns:
        tuple: A tuple containing:
              split_rules (list): List of structured rules.
              gene_in_rule (list): List containing genes found in the rules.
    """
    split_rules = []
    err_rules = []
    tmp_gene_in_rule = []
    for i in range(len(rules)):
        tmp = list(rules[i])
        if tmp:
            tmp, tmp_genes = check_and_doWord(tmp)
            tmp_gene_in_rule.extend(tmp_genes)
            if tmp is False:
                split_rules.append([])
                err_rules.append(rules[i])
            else:
                tmp = brackets_to_list(tmp)
                if checkRule(tmp):
                    split_rules.append(priorityAND(tmp))
                else:
                    split_rules.append([])
                    err_rules.append(rules[i])
        else:
            split_rules.append([])
    if err_rules:
        warning('Warning: wrong format rule in ' + str(err_rules) + '\n')
    return (split_rules, list(set(tmp_gene_in_rule)))

def make_recon(data) -> Tuple[List[str], List[str], Dict[str, str]]:
    """
    Read reaction rules from a given dataset, process them using the `do_rules` function,
    and return the IDs of reactions, structured rules, and a dictionary of genes found in the rules.

    Args:
        data (str): Path to the dataset file.

    Returns:
        tuple: A tuple containing:
            - ids (list): List of reaction IDs.
            - split_rules (list): List of structured rules.
            - gene_in_rule (dict): Dictionary containing genes found in the rules.
    """
    try:
        import cobra as cb
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            recon = cb.io.read_sbml_model(data)
        react = recon.reactions
        rules = [react[i].gene_reaction_rule for i in range(len(react))]
        ids = [react[i].id for i in range(len(react))]
    except:
        try:
            data = (pd.read_csv(data, sep = '\t', dtype = str, engine='python')).fillna('')
            if len(data.columns) < 2:
                sys.exit('Execution aborted: wrong format of '+
                         'custom datarules\n')
            if not len(data.columns) == 2:
                warning('Warning: more than 2 columns in custom datarules.\n' +
                        'Extra columns have been disregarded\n')
            ids = list(data.iloc[:, 0])
            rules = list(data.iloc[:, 1])
        except pd.errors.EmptyDataError:
            sys.exit('Execution aborted: wrong format of custom datarules\n')
        except pd.errors.ParserError:
            sys.exit('Execution aborted: wrong format of custom datarules\n')            
    split_rules, tmp_genes = do_rules(rules)
    gene_in_rule = {}
    for i in tmp_genes:
        gene_in_rule[i] = 'ok'
    return (ids, split_rules, gene_in_rule)

############################ gene #############################################
def data_gene(gene: pd.DataFrame, type_gene: str, name: str, gene_custom: Optional[Dict[str, str]]) -> Dict[str, str]:
    """
    Process gene data to ensure correct formatting and handle duplicates.

    Args:
        gene (DataFrame): DataFrame containing gene data.
        type_gene (str): Type of gene data (e.g., 'hugo_id', 'ensembl_gene_id', 'symbol', 'entrez_id').
        name (str): Name of the dataset.
        gene_custom (dict or None): Custom gene data dictionary if provided.

    Returns:
        dict: A dictionary containing gene data with gene IDs as keys and corresponding values.
    """
    args = process_args(sys.argv)    
    for i in range(len(gene)):
        tmp = gene.iloc[i, 0]
        if tmp.startswith(' ') or tmp.endswith(' '):
            gene.iloc[i, 0] = (tmp.lstrip()).rstrip()
    gene_dup = [item for item, count in 
               collections.Counter(gene[gene.columns[0]]).items() if count > 1]
    pat_dup = [item for item, count in 
               collections.Counter(list(gene.columns)).items() if count > 1]

    if gene_dup:
        if gene_custom == None:
            if args.rules_selector == 'HMRcore':
                gene_in_rule = pk.load(open(args.tool_dir + '/local/pickle files/HMRcore_genes.p', 'rb'))
            
            elif args.rules_selector == 'Recon':
                gene_in_rule = pk.load(open(args.tool_dir + '/local/pickle files/Recon_genes.p', 'rb'))
            
            elif args.rules_selector == 'ENGRO2':
                gene_in_rule = pk.load(open(args.tool_dir + '/local/pickle files/ENGRO2_genes.p', 'rb'))
            
            gene_in_rule = gene_in_rule.get(type_gene)
        
        else:
            gene_in_rule = gene_custom
        tmp = []
        for i in gene_dup:
            if gene_in_rule.get(i) == 'ok':
                tmp.append(i)
        if tmp:
            sys.exit('Execution aborted because gene ID '
                     +str(tmp)+' in '+name+' is duplicated\n')
    if pat_dup:
        warning('Warning: duplicated label\n' + str(pat_dup) + 'in ' + name + 
                '\n')
        
    return (gene.set_index(gene.columns[0])).to_dict()

############################ resolve ##########################################
ResolvedRules = Dict[str, List[Optional[Union[float, int]]]]
def resolve(genes: Dict[str, str], rules: List[str], ids: List[str], resolve_none: bool, name: str) -> Tuple[Optional[ResolvedRules], Optional[list]]:
    """
    Resolve rules using gene data to compute scores for each rule.

    Args:
        genes (dict): Dictionary containing gene data with gene IDs as keys and corresponding values.
        rules (list): List of rules to resolve.
        ids (list): List of IDs corresponding to the rules.
        resolve_none (bool): Flag indicating whether to resolve None values in the rules.
        name (str): Name of the dataset.

    Returns:
        tuple: A tuple containing resolved rules as a dictionary and a list of gene IDs not found in the data.
    """
    resolve_rules = {}
    not_found = []
    flag = False
    for key, value in genes.items():
        tmp_resolve = []
        for i in range(len(rules)):
            tmp = rules[i]
            if tmp:
                tmp, err = replace_gene_value(tmp, value)
                if err:
                    not_found.extend(err)
                ris = control(None, tmp, resolve_none)
                if ris is False or ris == None:
                    tmp_resolve.append(None)
                else:
                    tmp_resolve.append(ris)
                    flag = True
            else:
                tmp_resolve.append(None)    
        resolve_rules[key] = tmp_resolve
    if flag is False:
        warning('Warning: no computable score (due to missing gene values)' +
                'for class ' + name + ', the class has been disregarded\n')
        return (None, None)
    return (resolve_rules, list(set(not_found)))

############################ create_ras #######################################
def create_ras (resolve_rules: Optional[ResolvedRules], dataset_name: str, rules: List[str], ids: List[str], file: str) -> None:
    """
    Create a RAS (Reaction Activity Score) file from resolved rules.

    Args:
        resolve_rules (dict): Dictionary containing resolved rules.
        dataset_name (str): Name of the dataset.
        rules (list): List of rules.
        ids (list): List of IDs corresponding to the rules.
        file (str): Path to the output RAS file.

    Returns:
        None
    """
    if resolve_rules == None:
        warning("Couldn't generate RAS for current dataset: " + dataset_name)

    for geni in resolve_rules.values():
        for i, valori in enumerate(geni):
            if valori == None:
                geni[i] = 'None'
                
    output_ras = pd.DataFrame.from_dict(resolve_rules)
    
    output_ras.insert(0, 'Reactions', ids)
    output_to_csv = pd.DataFrame.to_csv(output_ras, sep = '\t', index = False)
    
    text_file = open(file, "w")
    
    text_file.write(output_to_csv)
    text_file.close()

############################ MAIN #############################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    """
    args = process_args(sys.argv)

    if args.rules_selector == 'HMRcore':        
        recon = pk.load(open(args.tool_dir + '/local/pickle files/HMRcore_rules.p', 'rb'))
    elif args.rules_selector == 'Recon':
        recon = pk.load(open(args.tool_dir + '/local/pickle files/Recon_rules.p', 'rb'))
    elif args.rules_selector == 'ENGRO2':
        recon = pk.load(open(args.tool_dir + '/local/pickle files/ENGRO2_rules.p', 'rb'))
    elif args.rules_selector == 'Custom':
        ids, rules, gene_in_rule = make_recon(args.custom)
        
    resolve_none = check_bool(args.none)
    
    
    name = "RAS Dataset"
    dataset = read_dataset(args.input, "dataset")

    dataset.iloc[:, 0] = (dataset.iloc[:, 0]).astype(str)

    type_gene = gene_type(dataset.iloc[0, 0], name) 
        
    if args.rules_selector != 'Custom':
        genes = data_gene(dataset, type_gene, name, None)
        ids, rules = load_id_rules(recon.get(type_gene))
    elif args.rules_selector == 'Custom':
        genes = data_gene(dataset, type_gene, name, gene_in_rule)
    
    resolve_rules, err = resolve(genes, rules, ids, resolve_none, name)

    create_ras(resolve_rules, name, rules, ids, args.ras_output)
      
    if err != None and err:
        warning('Warning: gene\n' + str(err) + '\nnot found in class '
            + name + ', the expression level for this gene ' +
            'will be considered NaN\n')

    
    print('Execution succeded')

    return None

###############################################################################
if __name__ == "__main__":
    main()
