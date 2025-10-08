"""
Generate Reaction Activity Scores (RAS) from a gene expression dataset and GPR rules.

The script reads a tabular dataset (genes x samples) and a rules file (GPRs),
computes RAS per reaction for each sample/cell line, and writes a tabular output.
"""
from __future__ import division
import sys
import argparse
import pandas as pd
import numpy as np
import utils.general_utils as utils
from typing import  List, Dict
import ast

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
        help = "path to input file containing the rules")

    parser.add_argument("-rn", "--model_upload_name", type = str, help = "custom rules name")
    # Galaxy converts files into .dat, this helps infer the original extension when needed.
    
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
        '-in', '--input',
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
        dataset = pd.read_csv(data, sep = '\t', header = 0, engine='python', index_col=0)
        dataset = dataset.astype(float)
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong file format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong file format of ' + name + '\n')
    return dataset


def load_custom_rules() -> Dict[str,str]:
    """
    Opens custom rules file and extracts the rules. If the file is in .csv format an additional parsing step will be
    performed, significantly impacting the runtime.

    Returns:
        Dict[str, ruleUtils.OpList] : dict mapping reaction IDs to rules.
    """
    datFilePath = utils.FilePath.fromStrPath(ARGS.model_upload)  # actual file, stored in Galaxy as a .dat

    dict_rule = {}

    try:
        rows = utils.readCsv(datFilePath, delimiter = "\t", skipHeader=False)
        if len(rows) <= 1:
            raise ValueError("Model tabular with 1 column is not supported.")

        if not rows:
            raise ValueError("Model tabular is file is empty.")
        
        id_idx, idx_gpr = utils.findIdxByName(rows[0], "GPR")
        
    # First, try using a tab delimiter
        for line in rows[1:]:
            if len(line) <= idx_gpr:
                utils.logWarning(f"Skipping malformed line: {line}", ARGS.out_log)
                continue
            
            dict_rule[line[id_idx]] = line[idx_gpr] 

    except Exception as e:
        # If parsing with tabs fails, try comma delimiter
        try:
            rows = utils.readCsv(datFilePath, delimiter = ",", skipHeader=False)
            
            if len(rows) <= 1:
                raise ValueError("Model tabular with 1 column is not supported.")

            if not rows:
                raise ValueError("Model tabular is file is empty.")
            
            id_idx, idx_gpr = utils.findIdxByName(rows[0], "GPR")
            
            # Try again parsing row content with the GPR column using comma-separated values
            for line in rows[1:]:
                if len(line) <= idx_gpr:
                    utils.logWarning(f"Skipping malformed line: {line}", ARGS.out_log)
                    continue
                
                dict_rule[line[id_idx]] =line[idx_gpr]
                    
        except Exception as e2:
            raise ValueError(f"Unable to parse rules file. Tried both tab and comma delimiters. Original errors: Tab: {e}, Comma: {e2}")

    if not dict_rule:
            raise ValueError("No valid rules found in the uploaded file. Please check the file format.")
    # csv rules need to be parsed, those in a pickle format are taken to be pre-parsed.
    return dict_rule



def computeRAS(
            dataset,gene_rules,rules_total_string,
            or_function=np.sum,    # type of operation to do in case of an or expression (max, sum, mean)
            and_function=np.min,   # type of operation to do in case of an and expression(min, sum)
            ignore_nan = True
            ):


    logic_operators = ['and', 'or', '(', ')']
    reactions=list(gene_rules.keys())

    # Build the dictionary for the GPRs
    df_reactions = pd.DataFrame(index=reactions)
    gene_rules=[gene_rules[reaction_id].replace("OR","or").replace("AND","and").replace("(","( ").replace(")"," )") for reaction_id in reactions]        
    df_reactions['rule'] = gene_rules
    df_reactions = df_reactions.reset_index()
    df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))

    dict_rule_reactions = df_reactions.to_dict()['index']

    # build useful structures for RAS computation
    genes =dataset.index.intersection(rules_total_string)
    
    #check if there is one gene at least 
    if len(genes)==0:
        raise ValueError("ERROR: No genes from the count matrix match the metabolic model. Check that gene annotations are consistent between model and dataset.") 
    
    cell_ids = list(dataset.columns)
    count_df_filtered = dataset.loc[genes]
    count_df_filtered = count_df_filtered.rename(index=lambda x: x.replace("-", "_").replace(":", "_"))

    ras_df=np.full((len(dict_rule_reactions), len(cell_ids)), np.nan)

    # for loop on rules
    genes_not_mapped=[]
    ind = 0       
    for rule, reaction_ids in dict_rule_reactions.items():
        if len(rule) != 0:
            # there is one gene at least in the formula
            warning_rule="_"
            if "-" in rule:
                warning_rule="-"
            if ":" in rule:
                warning_rule=":"
            rule_orig=rule.split().copy()  #original rule in list
            rule = rule.replace(warning_rule,"_")
             #modified rule
            rule_split = rule.split()
            rule_split_elements = list(filter(lambda x: x not in logic_operators, rule_split))  # remove of all logical operators
            rule_split_elements = list(set(rule_split_elements))                                # genes in formula
            
            # which genes are in the count matrix?                
            genes_in_count_matrix = [el for el in rule_split_elements if el in genes]
            genes_notin_count_matrix = []
            for el in rule_split_elements:
                if el not in genes: #not present in original dataset
                    if el.replace("_",warning_rule) in rule_orig: 
                        genes_notin_count_matrix.append(el.replace("_",warning_rule))
                    else:
                        genes_notin_count_matrix.append(el)
            genes_not_mapped.extend(genes_notin_count_matrix)
            
            # add genes not present in the data
            if len(genes_in_count_matrix) > 0: #there is at least one gene in the count matrix                 
                    if len(rule_split) == 1:
                        #one gene --> one reaction
                        ras_df[ind] = count_df_filtered.loc[genes_in_count_matrix]
                    else:    
                        if len(genes_notin_count_matrix) > 0 and ignore_nan == False:
                                ras_df[ind] = np.nan
                        else:                   
                            # more genes in the formula
                            check_only_and=("and" in rule and "or" not in rule) #only and
                            check_only_or=("or" in rule and "and" not in rule)  #only or
                            if check_only_and or check_only_or:
                                #or/and sequence
                                matrix = count_df_filtered.loc[genes_in_count_matrix].values
                                #compute for all cells
                                if check_only_and: 
                                    ras_df[ind] = and_function(matrix, axis=0)
                                else:
                                    ras_df[ind] = or_function(matrix, axis=0)
                            else:
                                # complex expression (e.g. A or (B and C))
                                data = count_df_filtered.loc[genes_in_count_matrix]  # dataframe of genes in the GPRs
                                tree = ast.parse(rule, mode="eval").body
                                values_by_cell = [dict(zip(data.index, data[col].values)) for col in data.columns]
                                for j, values in enumerate(values_by_cell):
                                    ras_df[ind, j] = _evaluate_ast(tree, values, or_function, and_function)
    
        ind +=1
    
    #create the dataframe of ras (rules x samples)
    ras_df= pd.DataFrame(data=ras_df,index=range(len(dict_rule_reactions)), columns=cell_ids)
    ras_df['Reactions'] = [reaction_ids for rule,reaction_ids in dict_rule_reactions.items()]
    
    #create the reaction dataframe for ras (reactions x samples)
    ras_df = ras_df.explode("Reactions").set_index("Reactions")

    #total genes not mapped from the gpr
    genes_not_mapped = sorted(set(genes_not_mapped))

    return ras_df,genes_not_mapped

# function to evalute a complex boolean expression e.g. A or (B and C)
# function to evalute a complex boolean expression e.g. A or (B and C)
def _evaluate_ast( node, values,or_function,and_function):
    if isinstance(node, ast.BoolOp):
        
        vals = [_evaluate_ast(v, values,or_function,and_function) for v in node.values]
       
        vals = [v for v in vals if v is not None]
        if not vals:
            return np.nan
      
        vals = [np.array(v) if isinstance(v, (list, np.ndarray)) else v for v in vals]

        if isinstance(node.op, ast.Or):
            return or_function(vals)
        elif isinstance(node.op, ast.And):
            return and_function(vals)

    elif isinstance(node, ast.Name):
        return values.get(node.id, None)
    elif isinstance(node, ast.Constant):
        key = str(node.value)     #convert in str       
        return values.get(key, None)   
    else:
        raise ValueError(f"Error in boolean expression: {ast.dump(node)}")

def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.
    
    Returns:
        None
    """
    # get args from frontend (related xml)
    global ARGS
    ARGS = process_args(args)

    # read dataset and remove versioning from gene names
    dataset = read_dataset(ARGS.input, "dataset")
    orig_gene_list=dataset.index.copy()
    dataset.index =  [str(el.split(".")[0]) for el in dataset.index]  

    #load GPR rules
    rules = load_custom_rules()
    
    #create a list of all the gpr
    rules_total_string=""
    for id,rule in rules.items():
        rules_total_string+=rule.replace("(","").replace(")","") + " "
    rules_total_string=list(set(rules_total_string.split(" ")))

    if any(dataset.index.duplicated(keep=False)):
        genes_duplicates=orig_gene_list[dataset.index.duplicated(keep=False)]
        genes_duplicates_in_model=[elem for elem in genes_duplicates if elem in rules_total_string]

        if len(genes_duplicates_in_model)>0:#metabolic genes have duplicated entries in the dataset
            list_str=", ".join(genes_duplicates_in_model)
            list_genes=f"ERROR: Duplicate entries in the gene dataset present in one or more GPR. The following metabolic genes are duplicated: "+list_str
            raise ValueError(list_genes)       
        else:
            list_str=", ".join(genes_duplicates)
            list_genes=f"INFO: Duplicate entries in the gene dataset. The following genes are duplicated in the dataset but not mentioned in the GPRs: "+list_str           
            utils.logWarning(list_genes,ARGS.out_log)  

    #check if nan value must be ignored in the GPR 
    if ARGS.none:
    #    #e.g. (A or nan --> A)
        ignore_nan = True
    else:
        #e.g. (A or nan --> nan)
        ignore_nan = False
    
    #compure ras
    ras_df,genes_not_mapped=computeRAS(dataset,rules,
                    rules_total_string,
                    or_function=np.sum,    # type of operation to do in case of an or expression (max, sum, mean)
                    and_function=np.min,
                    ignore_nan=ignore_nan)
    
    #save to csv and replace nan with None
    ras_df.replace([np.nan,None],"None").to_csv(ARGS.ras_output, sep = '\t')

    #report genes not present in the data
    if len(genes_not_mapped)>0: 
        genes_not_mapped_str=", ".join(genes_not_mapped)
        utils.logWarning(
        f"INFO: The following genes are mentioned in the GPR rules but don't appear in the dataset: "+genes_not_mapped_str,
        ARGS.out_log)  
 
    print("Execution succeeded")

###############################################################################
if __name__ == "__main__":
    main()