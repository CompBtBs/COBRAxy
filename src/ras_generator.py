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
try:
    from .utils import general_utils as utils
except:
    import utils.general_utils as utils
from typing import List, Dict
import ast

# Optional imports for AnnData mode (not used in ras_generator.py)
try:
    from progressbar import ProgressBar, Bar, Percentage
    from scanpy import AnnData
    from cobra.flux_analysis.variability import find_essential_reactions, find_essential_genes
except ImportError:
    # These are only needed for AnnData mode, not for ras_generator.py
    pass

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
        default = os.path.dirname(os.path.abspath(__file__)),
        help = 'your tool directory (default: auto-detected package location)')
    
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


"""
Class to compute the RAS values

"""

class RAS_computation:

    def __init__(self, adata=None, model=None, dataset=None, gene_rules=None, rules_total_string=None):
        """
        Initialize RAS computation with two possible input modes:
        
        Mode 1 (Original - for sampling_main.py):
            adata: AnnData object with gene expression (cells × genes)
            model: COBRApy model object with reactions and GPRs
            
        Mode 2 (New - for ras_generator.py):
            dataset: pandas DataFrame with gene expression (genes × samples)
            gene_rules: dict mapping reaction IDs to GPR strings
            rules_total_string: list of all gene names in GPRs (for validation)
        """
        self._logic_operators = ['and', 'or', '(', ')']
        self.val_nan = np.nan
        
        # Determine which mode we're in
        if adata is not None and model is not None:
            # Mode 1: AnnData + COBRApy model (original)
            self._init_from_anndata(adata, model)
        elif dataset is not None and gene_rules is not None:
            # Mode 2: DataFrame + rules dict (ras_generator style)
            self._init_from_dataframe(dataset, gene_rules, rules_total_string)
        else:
            raise ValueError(
                "Invalid initialization. Provide either:\n"
                "  - adata + model (for AnnData input), or\n"
                "  - dataset + gene_rules (for DataFrame input)"
            )
    
    def _normalize_gene_name(self, gene_name):
        """Normalize gene names by replacing special characters."""
        return gene_name.replace("-", "_").replace(":", "_")
    
    def _normalize_rule(self, rule):
        """Normalize GPR rule: lowercase operators, add spaces around parentheses, normalize gene names."""
        rule = rule.replace("OR", "or").replace("AND", "and")
        rule = rule.replace("(", "( ").replace(")", " )")
        # Normalize gene names in the rule
        tokens = rule.split()
        normalized_tokens = [token if token in self._logic_operators else self._normalize_gene_name(token) for token in tokens]
        return " ".join(normalized_tokens)
    
    def _init_from_anndata(self, adata, model):
        """Initialize from AnnData and COBRApy model (original mode)."""
        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=[reaction.id for reaction in model.reactions])
        gene_rules = [self._normalize_rule(reaction.gene_reaction_rule) for reaction in model.reactions]
        df_reactions['rule'] = gene_rules
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))
        
        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = model
        self.count_adata = adata.copy()
        
        # Normalize gene names in both model and dataset
        model_genes = [self._normalize_gene_name(gene.id) for gene in model.genes]
        dataset_genes = [self._normalize_gene_name(gene) for gene in self.count_adata.var.index]
        self.genes = pd.Index(dataset_genes).intersection(model_genes)
        
        if len(self.genes) == 0:
            raise ValueError("ERROR: No genes from the count matrix match the metabolic model. Check that gene annotations are consistent between model and dataset.")
        
        self.cell_ids = list(self.count_adata.obs.index.values)
        # Get expression data with normalized gene names
        self.count_df_filtered = self.count_adata.to_df().T
        self.count_df_filtered.index = [self._normalize_gene_name(g) for g in self.count_df_filtered.index]
        self.count_df_filtered = self.count_df_filtered.loc[self.genes]
    
    def _init_from_dataframe(self, dataset, gene_rules, rules_total_string):
        """Initialize from DataFrame and rules dict (ras_generator mode)."""
        reactions = list(gene_rules.keys())
        
        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=reactions)
        gene_rules_list = [self._normalize_rule(gene_rules[reaction_id]) for reaction_id in reactions]
        df_reactions['rule'] = gene_rules_list
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))

        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = None
        self.count_adata = None
        
        # Normalize gene names in dataset
        dataset_normalized = dataset.copy()
        dataset_normalized.index = [self._normalize_gene_name(g) for g in dataset_normalized.index]
        
        # Determine which genes are in both dataset and GPRs
        if rules_total_string is not None:
            rules_genes = [self._normalize_gene_name(g) for g in rules_total_string]
            self.genes = dataset_normalized.index.intersection(rules_genes)
        else:
            # Extract all genes from rules
            all_genes_in_rules = set()
            for rule in gene_rules_list:
                tokens = rule.split()
                for token in tokens:
                    if token not in self._logic_operators:
                        all_genes_in_rules.add(token)
            self.genes = dataset_normalized.index.intersection(all_genes_in_rules)
        
        if len(self.genes) == 0:
            raise ValueError("ERROR: No genes from the count matrix match the metabolic model. Check that gene annotations are consistent between model and dataset.")
        
        self.cell_ids = list(dataset_normalized.columns)
        self.count_df_filtered = dataset_normalized.loc[self.genes]
 
    def compute(self,
                or_expression=np.sum,       # type of operation to do in case of an or expression (sum, max, mean)
                and_expression=np.min,      # type of operation to do in case of an and expression(min, sum)
                drop_na_rows=False,          # if True remove the nan rows of the ras  matrix
                drop_duplicates=False,      # if true, remove duplicates rows
                ignore_nan=True,            # if True, ignore NaN values in GPR evaluation (e.g., A or NaN -> A)
                print_progressbar=True,     # if True, print the progress bar
                add_count_metadata=True,    # if True add metadata of cells in the ras adata
                add_met_metadata=True,      # if True add metadata from the metabolic model (gpr and compartments of reactions)
                add_essential_reactions=False,
                add_essential_genes=False
                ):

        self.or_function = or_expression
        self.and_function = and_expression
        
        ras_df = np.full((len(self.dict_rule_reactions), len(self.cell_ids)), np.nan)
        genes_not_mapped = set()  # Track genes not in dataset
        
        if print_progressbar:
            pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(self.dict_rule_reactions)).start()
        
        # Process each unique GPR rule
        for ind, (rule, reaction_ids) in enumerate(self.dict_rule_reactions.items()):
            if len(rule) == 0:
                # Empty rule - keep as NaN
                pass
            else:
                # Extract genes from rule
                rule_genes = [token for token in rule.split() if token not in self._logic_operators]
                rule_genes_unique = list(set(rule_genes))
                
                # Which genes are in the dataset?
                genes_present = [g for g in rule_genes_unique if g in self.genes]
                genes_missing = [g for g in rule_genes_unique if g not in self.genes]
                
                if genes_missing:
                    genes_not_mapped.update(genes_missing)
                
                if len(genes_present) == 0:
                    # No genes in dataset - keep as NaN
                    pass
                elif len(genes_missing) > 0 and not ignore_nan:
                    # Some genes missing and we don't ignore NaN - set to NaN
                    pass
                else:
                    # Evaluate the GPR expression using AST
                    # For single gene, AST handles it fine: ast.parse("GENE_A") works
                    # more genes in the formula
                    check_only_and=("and" in rule and "or" not in rule) #only and
                    check_only_or=("or" in rule and "and" not in rule)  #only or
                    if check_only_and or check_only_or:
                        #or/and sequence
                        matrix = self.count_df_filtered.loc[genes_present].values
                        #compute for all cells
                        if check_only_and: 
                            ras_df[ind] = self.and_function(matrix, axis=0)
                        else:
                            ras_df[ind] = self.or_function(matrix, axis=0)
                    else:
                        # complex expression (e.g. A or (B and C))
                        data = self.count_df_filtered.loc[genes_present]  # dataframe of genes in the GPRs
                        tree = ast.parse(rule, mode="eval").body
                        values_by_cell = [dict(zip(data.index, data[col].values)) for col in data.columns]
                        for j, values in enumerate(values_by_cell):
                            ras_df[ind, j] =self._evaluate_ast(tree, values, self.or_function, self.and_function, ignore_nan)

            if print_progressbar:
                pbar.update(ind + 1)
        
        if print_progressbar:
            pbar.finish()
        
        # Store genes not mapped for later use
        self.genes_not_mapped = sorted(genes_not_mapped)
        
        # create the dataframe of ras (rules x samples)
        ras_df = pd.DataFrame(data=ras_df, index=range(len(self.dict_rule_reactions)), columns=self.cell_ids)
        ras_df['Reactions'] = [reaction_ids for rule, reaction_ids in self.dict_rule_reactions.items()]
        
        reactions_common = pd.DataFrame()
        reactions_common["Reactions"] = ras_df['Reactions']
        reactions_common["proof2"] = ras_df['Reactions']
        reactions_common = reactions_common.explode('Reactions')
        reactions_common = reactions_common.set_index("Reactions")

        ras_df = ras_df.explode("Reactions")
        ras_df = ras_df.set_index("Reactions")

        if drop_na_rows:
            ras_df = ras_df.dropna(how="all")
            
        if drop_duplicates:
            ras_df = ras_df.drop_duplicates()
        
        # If initialized from DataFrame (ras_generator mode), return DataFrame instead of AnnData
        if self.count_adata is None:
            return ras_df, self.genes_not_mapped
        
        # Original AnnData mode: create AnnData structure for RAS
        ras_adata = AnnData(ras_df.T)

        #add metadata
        if add_count_metadata:
            ras_adata.var["common_gprs"] = reactions_common.loc[ras_df.index]
            ras_adata.var["common_gprs"] = ras_adata.var["common_gprs"].apply(lambda x: ",".join(x))
            for el in self.count_adata.obs.columns:
                ras_adata.obs["countmatrix_"+el]=self.count_adata.obs[el]

        if add_met_metadata:
            if self.model is not None and len(self.model.compartments)>0:
                  ras_adata.var['compartments']=[list(self.model.reactions.get_by_id(reaction).compartments) for reaction in ras_adata.var.index]  
                  ras_adata.var['compartments']=ras_adata.var["compartments"].apply(lambda x: ",".join(x))
            
            if self.model is not None:
                ras_adata.var['GPR rule'] = [self.model.reactions.get_by_id(reaction).gene_reaction_rule for reaction in ras_adata.var.index]

        if add_essential_reactions:
            if self.model is not None:
                essential_reactions=find_essential_reactions(self.model)
                essential_reactions=[el.id for el in essential_reactions]            
                ras_adata.var['essential reactions']=["yes" if el in essential_reactions else "no" for el in ras_adata.var.index]
        
        if add_essential_genes:
            if self.model is not None:
                essential_genes=find_essential_genes(self.model)
                essential_genes=[el.id for el in essential_genes]
                ras_adata.var['essential genes']=[" ".join([gene for gene in genes.split()  if gene in essential_genes]) for genes in ras_adata.var["GPR rule"]]
        
        return ras_adata

    def _evaluate_ast(self, node, values, or_function, and_function, ignore_nan):
        """
        Evaluate a boolean expression using AST (Abstract Syntax Tree).
        Handles all GPR types: single gene, simple (A and B), nested (A or (B and C)).
        
        Args:
            node: AST node to evaluate
            values: Dictionary mapping gene names to their expression values
            or_function: Function to apply for OR operations
            and_function: Function to apply for AND operations
            ignore_nan: If True, ignore None/NaN values (e.g., A or None -> A)
            
        Returns:
            Evaluated expression result (float or np.nan)
        """
        if isinstance(node, ast.BoolOp):
            # Boolean operation (and/or)
            vals = [self._evaluate_ast(v, values, or_function, and_function, ignore_nan) for v in node.values]
            
            if ignore_nan:
                # Filter out None/NaN values
                vals = [v for v in vals if v is not None and not (isinstance(v, float) and np.isnan(v))]
            
            if not vals:
                return np.nan
            
            if isinstance(node.op, ast.Or):
                return or_function(vals)
            elif isinstance(node.op, ast.And):
                return and_function(vals)

        elif isinstance(node, ast.Name):
            # Variable (gene name)
            return values.get(node.id, None)
        elif isinstance(node, ast.Constant):
            # Constant (shouldn't happen in GPRs, but handle it)
            return values.get(str(node.value), None)
        else:
            raise ValueError(f"Unexpected node type in GPR: {ast.dump(node)}")


# ============================================================================
# STANDALONE FUNCTION FOR RAS_GENERATOR COMPATIBILITY
# ============================================================================

def computeRAS(
    dataset, 
    gene_rules, 
    rules_total_string,
    or_function=np.sum,
    and_function=np.min,
    ignore_nan=True
):
    """
    Compute RAS from tabular data and GPR rules (ras_generator.py compatible).
    
    This is a standalone function that wraps the RAS_computation class
    to provide the same interface as ras_generator.py.
    
    Args:
        dataset: pandas DataFrame with gene expression (genes × samples)
        gene_rules: dict mapping reaction IDs to GPR strings
        rules_total_string: list of all gene names in GPRs
        or_function: function for OR operations (default: np.sum)
        and_function: function for AND operations (default: np.min)
        ignore_nan: if True, ignore NaN in GPR evaluation (default: True)
    
    Returns:
        tuple: (ras_df, genes_not_mapped)
            - ras_df: DataFrame with RAS values (reactions × samples)
            - genes_not_mapped: list of genes in GPRs not found in dataset
    """
    # Create RAS computation object in DataFrame mode
    ras_obj = RAS_computation(
        dataset=dataset,
        gene_rules=gene_rules,
        rules_total_string=rules_total_string
    )
    
    # Compute RAS
    result = ras_obj.compute(
        or_expression=or_function,
        and_expression=and_function,
        ignore_nan=ignore_nan,
        print_progressbar=False,  # No progress bar for ras_generator
        add_count_metadata=False,  # No metadata in DataFrame mode
        add_met_metadata=False,
        add_essential_reactions=False,
        add_essential_genes=False
    )
    
    # Result is a tuple (ras_df, genes_not_mapped) in DataFrame mode
    return result

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