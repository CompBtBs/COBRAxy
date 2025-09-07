from __future__ import division
# galaxy complains this ^^^ needs to be at the very beginning of the file, for some reason.
import sys
import argparse
import collections
import pandas as pd
import pickle as pk
import utils.general_utils as utils
import utils.rule_parsing as ruleUtils
from typing import Union, Optional, List, Dict, Tuple, TypeVar
import os

ERRORS = []
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
    parser = argparse.ArgumentParser(
        usage = '%(prog)s [options]',
        description = "process some value's genes to create a comparison's map.")
    
    parser.add_argument("-rl", "--model_upload", type = str,
        help = "path to input file with custom rules, if provided")

    parser.add_argument("-rn", "--model_upload_name", type = str, help = "custom rules name")
    # ^ I need this because galaxy converts my files into .dat but I need to know what extension they were in
    
    parser.add_argument(
        '-n', '--none',
        type = utils.Bool("none"), default = True,
        help = 'compute Nan values')
    
    parser.add_argument(
        '-td', '--tool_dir',
        type = str,
        required = True, help = 'your tool directory')
    
    parser.add_argument(
        '-ol', '--out_log',
        type = str,
        help = "Output log")    
    
    parser.add_argument(
        '-in', '--input', #id è diventato in
        type = str,
        help = 'input dataset')
    
    parser.add_argument(
        '-ra', '--ras_output',
        type = str,
        required = True, help = 'ras output')

    
    return parser.parse_args(args)

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

############################ load id e rules ##################################
def load_id_rules(reactions :Dict[str, Dict[str, List[str]]]) -> Tuple[List[str], List[Dict[str, List[str]]]]:
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
    return l.upper().startswith('ENS')
 

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
 
    for i in range(len(gene)):
        tmp = gene.iloc[i, 0]
        gene.iloc[i, 0] = tmp.strip().split('.')[0]

    gene_dup = [item for item, count in 
               collections.Counter(gene[gene.columns[0]]).items() if count > 1]
    pat_dup = [item for item, count in 
               collections.Counter(list(gene.columns)).items() if count > 1]
    
    gene_in_rule = None

    if gene_dup:
        if gene_custom == None:

            if str(ARGS.rules_selector) == 'HMRcore':
                gene_in_rule = pk.load(open(ARGS.tool_dir + '/local/pickle files/HMRcore_genes.p', 'rb'))
            
            elif str(ARGS.rules_selector) == 'Recon':
                gene_in_rule = pk.load(open(ARGS.tool_dir + '/local/pickle files/Recon_genes.p', 'rb'))
            
            elif str(ARGS.rules_selector) == 'ENGRO2':
                gene_in_rule = pk.load(open(ARGS.tool_dir + '/local/pickle files/ENGRO2_genes.p', 'rb'))

            utils.logWarning(f"{ARGS.tool_dir}'/local/pickle files/ENGRO2_genes.p'", ARGS.out_log)

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
    
    if pat_dup: utils.logWarning(f"Warning: duplicated label\n{pat_dup} in {name}", ARGS.out_log)
    return (gene.set_index(gene.columns[0])).to_dict()

############################ resolve ##########################################
def replace_gene_value(l :str, d :str) -> Tuple[Union[int, float], list]:
    """
    Replace gene identifiers with corresponding values from a dictionary.

    Args:
        l (str): String of gene identifier.
        d (str): String corresponding to its value.

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
        d (str): String corresponding to its value.

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
        utils.logWarning(
            f"Warning: no computable score (due to missing gene values) for class {name}, the class has been disregarded",
            ARGS.out_log)
        
        return (None, None)
    
    return (resolve_rules, list(set(not_found)))
############################ create_ras #######################################
def create_ras(resolve_rules: Optional[ResolvedRules], dataset_name: str, rules: List[str], ids: List[str], file: str) -> None:
    """
    Create a RAS (Reaction Activity Score) file from resolved rules.

    Args:
        resolve_rules (dict): Dictionary containing resolved rules.
        dataset_name (str): Name of the dataset.
        rules (list): List of rules.
        file (str): Path to the output RAS file.

    Returns:
        None
    """
    if resolve_rules is None:
        utils.logWarning(f"Couldn't generate RAS for current dataset: {dataset_name}", ARGS.out_log)

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

################################- NEW RAS COMPUTATION -################################
Expr = Optional[Union[int, float]]
Ras  = Expr
def ras_for_cell_lines(dataset: pd.DataFrame, rules: Dict[str, ruleUtils.OpList]) -> Dict[str, Dict[str, Ras]]:
    """
    Generates the RAS scores for each cell line found in the dataset.

    Args:
        dataset (pd.DataFrame): Dataset containing gene values.
        rules (dict): The dict containing reaction ids as keys and rules as values.

    Side effects:
        dataset : mut
    
    Returns:
        dict: A dictionary where each key corresponds to a cell line name and each value is a dictionary
        where each key corresponds to a reaction ID and each value is its computed RAS score.
    """
    ras_values_by_cell_line = {}
    dataset.set_index(dataset.columns[0], inplace=True)
    
    for cell_line_name in dataset.columns: #[1:]:
        cell_line = dataset[cell_line_name].to_dict()
        ras_values_by_cell_line[cell_line_name]= get_ras_values(rules, cell_line)
    return ras_values_by_cell_line

def get_ras_values(value_rules: Dict[str, ruleUtils.OpList], dataset: Dict[str, Expr]) -> Dict[str, Ras]:
    """
    Computes the RAS (Reaction Activity Score) values for each rule in the given dict.

    Args:
        value_rules (dict): A dictionary where keys are reaction ids and values are OpLists.
        dataset : gene expression data of one cell line.

    Returns:
        dict: A dictionary where keys are reaction ids and values are the computed RAS values for each rule.
    """
    return {key: ras_op_list(op_list, dataset) for key, op_list in value_rules.items()}

def get_gene_expr(dataset :Dict[str, Expr], name :str) -> Expr:
    """
    Extracts the gene expression of the given gene from a cell line dataset.

    Args:
        dataset : gene expression data of one cell line.
        name : gene name.
    
    Returns:
        Expr : the gene's expression value.
    """
    expr = dataset.get(name, None)
    if expr is None: ERRORS.append(name)
  
    return expr

def ras_op_list(op_list: ruleUtils.OpList, dataset: Dict[str, Expr]) -> Ras:
    """
    Computes recursively the RAS (Reaction Activity Score) value for the given OpList, considering the specified flag to control None behavior.

    Args:
        op_list (OpList): The OpList representing a rule with gene values.
        dataset : gene expression data of one cell line.

    Returns:
        Ras: The computed RAS value for the given OpList.
    """
    op = op_list.op
    ras_value :Ras = None
    if not op: return get_gene_expr(dataset, op_list[0])
    if op is ruleUtils.RuleOp.AND and not ARGS.none and None in op_list: return None

    for i in range(len(op_list)):
        item = op_list[i]
        if isinstance(item, ruleUtils.OpList):
            item = ras_op_list(item, dataset)

        else:
          item = get_gene_expr(dataset, item)

        if item is None:
          if op is ruleUtils.RuleOp.AND and not ARGS.none: return None
          continue

        if ras_value is None:
          ras_value = item
        else:
          ras_value = ras_value + item if op is ruleUtils.RuleOp.OR else min(ras_value, item)

    return ras_value

def save_as_tsv(rasScores: Dict[str, Dict[str, Ras]], reactions :List[str]) -> None:
    """
    Save computed ras scores to the given path, as a tsv file.

    Args:
        rasScores : the computed ras scores.
        path : the output tsv file's path.
    
    Returns:
        None
    """
    for scores in rasScores.values(): # this is actually a lot faster than using the ootb dataframe metod, sadly
        for reactId, score in scores.items():
            if score is None: scores[reactId] = "None"

    output_ras = pd.DataFrame.from_dict(rasScores)
    output_ras.insert(0, 'Reactions', reactions)
    output_ras.to_csv(ARGS.ras_output, sep = '\t', index = False)

############################ MAIN #############################################
#TODO: not used but keep, it will be when the new translator dicts will be used.
def translateGene(geneName :str, encoding :str, geneTranslator :Dict[str, Dict[str, str]]) -> str:
    """
    Translate gene from any supported encoding to HugoID.

    Args:
        geneName (str): the name of the gene in its current encoding.
        encoding (str): the encoding.
        geneTranslator (Dict[str, Dict[str, str]]): the dict containing all supported gene names
        and encodings in the current model, mapping each to the corresponding HugoID encoding.

    Raises:
        ValueError: When the gene isn't supported in the model.

    Returns:
        str: the gene in HugoID encoding.
    """
    supportedGenesInEncoding = geneTranslator[encoding]
    if geneName in supportedGenesInEncoding: return supportedGenesInEncoding[geneName]
    raise ValueError(f"Gene \"{geneName}\" non trovato, verifica di star utilizzando il modello corretto!")

def load_custom_rules() -> Dict[str, ruleUtils.OpList]:
    """
    Opens custom rules file and extracts the rules. If the file is in .csv format an additional parsing step will be
    performed, significantly impacting the runtime.

    Returns:
        Dict[str, ruleUtils.OpList] : dict mapping reaction IDs to rules.
    """
    datFilePath = utils.FilePath.fromStrPath(ARGS.model_upload) # actual file, stored in galaxy as a .dat

    #try: filenamePath = utils.FilePath.fromStrPath(ARGS.model_upload_name) # file's name in input, to determine its original ext
    #except utils.PathErr as err:      
    #    utils.logWarning(f"Cannot determine file extension from filename '{ARGS.model_upload_name}'. Assuming tabular format.", ARGS.out_log)
    #    filenamePath = None
     
    #if filenamePath.ext is utils.FileFormat.PICKLE: return utils.readPickle(datFilePath)

    dict_rule = {}

    try:
        # Proviamo prima con delimitatore tab
        for line in utils.readCsv(datFilePath, delimiter = "\t"):
            if len(line) < 3:  # Controlliamo che ci siano almeno 3 colonne
                utils.logWarning(f"Skipping malformed line: {line}", ARGS.out_log)
                continue
            
            if line[2] == "":
                dict_rule[line[0]] = ruleUtils.OpList([""])
            else:
                dict_rule[line[0]] = ruleUtils.parseRuleToNestedList(line[2])
                
    except Exception as e:
        # Se fallisce con tab, proviamo con virgola
        try:
            dict_rule = {}
            for line in utils.readCsv(datFilePath, delimiter = ","):
                if len(line) < 3:
                    utils.logWarning(f"Skipping malformed line: {line}", ARGS.out_log)
                    continue
                
                if line[2] == "":
                    dict_rule[line[0]] = ruleUtils.OpList([""])
                else:
                    dict_rule[line[0]] = ruleUtils.parseRuleToNestedList(line[2])
        except Exception as e2:
            raise ValueError(f"Unable to parse rules file. Tried both tab and comma delimiters. Original errors: Tab: {e}, Comma: {e2}")

    if not dict_rule:
            raise ValueError("No valid rules found in the uploaded file. Please check the file format.")
    # csv rules need to be parsed, those in a pickle format are taken to be pre-parsed.
    return dict_rule


def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.
    
    Returns:
        None
    """
    # get args from frontend (related xml)
    global ARGS
    ARGS = process_args(args)

    # read dataset
    dataset = read_dataset(ARGS.input, "dataset")
    dataset.iloc[:, 0] = (dataset.iloc[:, 0]).astype(str)

    # remove versioning from gene names
    dataset.iloc[:, 0] = dataset.iloc[:, 0].str.split('.').str[0]

    rules = load_custom_rules()
    reactions = list(rules.keys())

    save_as_tsv(ras_for_cell_lines(dataset, rules), reactions)
    if ERRORS: utils.logWarning(
        f"The following genes are mentioned in the rules but don't appear in the dataset: {ERRORS}",
        ARGS.out_log)  


    ############

    # handle custom models
    #model :utils.Model = ARGS.rules_selector

    #if model is utils.Model.Custom:
    #    rules = load_custom_rules()
    #    reactions = list(rules.keys())

    #    save_as_tsv(ras_for_cell_lines(dataset, rules), reactions)
    #    if ERRORS: utils.logWarning(
    #        f"The following genes are mentioned in the rules but don't appear in the dataset: {ERRORS}",
    #        ARGS.out_log)
        
    #    return
    
    # This is the standard flow of the ras_generator program, for non-custom models.
    #name = "RAS Dataset"
    #type_gene = gene_type(dataset.iloc[0, 0], name)

    #rules      = model.getRules(ARGS.tool_dir)
    #genes      = data_gene(dataset, type_gene, name, None)
    #ids, rules = load_id_rules(rules.get(type_gene))

    #resolve_rules, err = resolve(genes, rules, ids, ARGS.none, name)
    #create_ras(resolve_rules, name, rules, ids, ARGS.ras_output)
    
    #if err: utils.logWarning(
    #    f"Warning: gene(s) {err} not found in class \"{name}\", " +
    #    "the expression level for this gene will be considered NaN",
    #    ARGS.out_log)
    
    print("Execution succeded")

###############################################################################
if __name__ == "__main__":
    main()
