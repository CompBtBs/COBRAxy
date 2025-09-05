import os
import csv
import cobra
import pickle
import argparse
import pandas as pd
import utils.general_utils as utils
import utils.rule_parsing  as rulesUtils
from typing import Optional, Tuple, Union, List, Dict
import utils.reaction_parsing as reactionUtils
import openpyxl

ARGS : argparse.Namespace
def process_args(args: List[str] = None) -> argparse.Namespace:
    """
    Parse command-line arguments for CustomDataGenerator.
    """

    parser = argparse.ArgumentParser(
        usage="%(prog)s [options]",
        description="Generate custom data from a given model"
    )

    parser.add_argument("--out_log", type=str, required=True,
                        help="Output log file")

    parser.add_argument("--model", type=str,
                        help="Built-in model identifier (e.g., ENGRO2, Recon, HMRcore)")
    parser.add_argument("--input", type=str,
                        help="Custom model file (JSON or XML)")
    parser.add_argument("--name", type=str, required=True,
                        help="Model name (default or custom)")
    
    parser.add_argument("--medium_selector", type=str, required=True,
                        help="Medium selection option (default/custom)")
    parser.add_argument("--medium", type=str,
                        help="Custom medium file if medium_selector=Custom")
    
    parser.add_argument("--out_tabular", type=str,
                        help="Output file for the merged dataset (CSV or XLSX)")
    
    parser.add_argument("--tool_dir", type=str, default=os.path.dirname(__file__),
                        help="Tool directory (passed from Galaxy as $__tool_directory__)")


    return parser.parse_args(args)

################################- INPUT DATA LOADING -################################
def load_custom_model(file_path :utils.FilePath, ext :Optional[utils.FileFormat] = None) -> cobra.Model:
    """
    Loads a custom model from a file, either in JSON or XML format.

    Args:
        file_path : The path to the file containing the custom model.
        ext : explicit file extension. Necessary for standard use in galaxy because of its weird behaviour.

    Raises:
        DataErr : if the file is in an invalid format or cannot be opened for whatever reason.    
    
    Returns:
        cobra.Model : the model, if successfully opened.
    """
    ext = ext if ext else file_path.ext
    try:
        if ext is utils.FileFormat.XML:
            return cobra.io.read_sbml_model(file_path.show())
        
        if ext is utils.FileFormat.JSON:
            return cobra.io.load_json_model(file_path.show())

    except Exception as e: raise utils.DataErr(file_path, e.__str__())
    raise utils.DataErr(file_path,
        f"Formato \"{file_path.ext}\" non riconosciuto, sono supportati solo file JSON e XML")

################################- DATA GENERATION -################################
ReactionId = str
def generate_rules(model: cobra.Model, *, asParsed = True) -> Union[Dict[ReactionId, rulesUtils.OpList], Dict[ReactionId, str]]:
    """
    Generates a dictionary mapping reaction ids to rules from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the rules to an optimized runtime format, otherwise leaves them as strings.

    Returns:
        Dict[ReactionId, rulesUtils.OpList] : the generated dictionary of parsed rules.
        Dict[ReactionId, str] : the generated dictionary of raw rules.
    """
    # Is the below approach convoluted? yes
    # Ok but is it inefficient? probably
    # Ok but at least I don't have to repeat the check at every rule (I'm clinically insane)
    _ruleGetter   =  lambda reaction : reaction.gene_reaction_rule
    ruleExtractor = (lambda reaction :
        rulesUtils.parseRuleToNestedList(_ruleGetter(reaction))) if asParsed else _ruleGetter

    return {
        reaction.id : ruleExtractor(reaction)
        for reaction in model.reactions
        if reaction.gene_reaction_rule }

def generate_reactions(model :cobra.Model, *, asParsed = True) -> Dict[ReactionId, str]:
    """
    Generates a dictionary mapping reaction ids to reaction formulas from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the reactions to an optimized runtime format, otherwise leaves them as they are.

    Returns:
        Dict[ReactionId, str] : the generated dictionary.
    """

    unparsedReactions = {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction 
    }

    if not asParsed: return unparsedReactions
    
    return reactionUtils.create_reaction_dict(unparsedReactions)

def get_medium(model:cobra.Model) -> pd.DataFrame:
    trueMedium=[]
    for r in model.reactions:
        positiveCoeff=0
        for m in r.metabolites:
            if r.get_coefficient(m.id)>0:
                positiveCoeff=1;
        if (positiveCoeff==0 and r.lower_bound<0):
            trueMedium.append(r.id)

    df_medium = pd.DataFrame()
    df_medium["reaction"] = trueMedium
    return df_medium

def generate_bounds(model:cobra.Model) -> pd.DataFrame:

    rxns = []
    for reaction in model.reactions:
        rxns.append(reaction.id)

    bounds = pd.DataFrame(columns = ["lower_bound", "upper_bound"], index=rxns)

    for reaction in model.reactions:
        bounds.loc[reaction.id] = [reaction.lower_bound, reaction.upper_bound]
    return bounds


###############################- FILE SAVING -################################
def save_as_csv_filePath(data :dict, file_path :utils.FilePath, fieldNames :Tuple[str, str]) -> None:
    """
    Saves any dictionary-shaped data in a .csv file created at the given file_path as FilePath.

    Args:
        data : the data to be written to the file.
        file_path : the path to the .csv file.
        fieldNames : the names of the fields (columns) in the .csv file.
    
    Returns:
        None
    """
    with open(file_path.show(), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fieldNames, dialect="excel-tab")
        writer.writeheader()

        for key, value in data.items():
            writer.writerow({ fieldNames[0] : key, fieldNames[1] : value })

def save_as_csv(data :dict, file_path :str, fieldNames :Tuple[str, str]) -> None:
    """
    Saves any dictionary-shaped data in a .csv file created at the given file_path as string.

    Args:
        data : the data to be written to the file.
        file_path : the path to the .csv file.
        fieldNames : the names of the fields (columns) in the .csv file.
    
    Returns:
        None
    """
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fieldNames, dialect="excel-tab")
        writer.writeheader()

        for key, value in data.items():
            writer.writerow({ fieldNames[0] : key, fieldNames[1] : value })

def save_as_tabular_df(df: pd.DataFrame, path: str) -> None:
    try:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        df.to_csv(path, sep="\t", index=False)
    except Exception as e:
        raise utils.DataErr(path, f"failed writing tabular output: {e}")


###############################- ENTRY POINT -################################
def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.
    
    Returns:
        None
    """
    # get args from frontend (related xml)
    global ARGS
    ARGS = process_args(args)


    if ARGS.input:
        # load custom model
        model = load_custom_model(
            utils.FilePath.fromStrPath(ARGS.input), utils.FilePath.fromStrPath(ARGS.name).ext)
    else:
        # load built-in model

        try:
            model_enum = utils.Model[ARGS.model]  # e.g., Model['ENGRO2']
        except KeyError:
            raise utils.ArgsErr("model", "one of Recon/ENGRO2/HMRcore/Custom_model", ARGS.model)

        # Load built-in model (Model.getCOBRAmodel uses tool_dir to locate local models)
        try:
            model = model_enum.getCOBRAmodel(toolDir=ARGS.tool_dir)
        except Exception as e:
            # Wrap/normalize load errors as DataErr for consistency
            raise utils.DataErr(ARGS.model, f"failed loading built-in model: {e}")

    # Determine final model name: explicit --name overrides, otherwise use the model id
    model_name = ARGS.name if ARGS.name else ARGS.model

    # generate data
    rules = generate_rules(model, asParsed = False)
    reactions = generate_reactions(model, asParsed = False)
    bounds = generate_bounds(model)
    medium = get_medium(model)

    df_rules = pd.DataFrame(list(rules.items()), columns = ["ReactionID", "Rule"])
    df_reactions = pd.DataFrame(list(reactions.items()), columns = ["ReactionID", "Reaction"])

    df_bounds = bounds.reset_index().rename(columns = {"index": "ReactionID"})
    df_medium = medium.rename(columns = {"reaction": "ReactionID"})
    df_medium["InMedium"] = True # flag per indicare la presenza nel medium

    merged = df_reactions.merge(df_rules, on = "ReactionID", how = "outer")
    merged = merged.merge(df_bounds, on = "ReactionID", how = "outer")

    merged = merged.merge(df_medium, on = "ReactionID", how = "left")

    merged["InMedium"] = merged["InMedium"].fillna(False)

    merged = merged.sort_values(by = "InMedium", ascending = False)

    #out_file = os.path.join(ARGS.output_path, f"{os.path.basename(ARGS.name).split('.')[0]}_custom_data")

    #merged.to_csv(out_file, sep = '\t', index = False)


    ####


    if not ARGS.out_tabular:
        raise utils.ArgsErr("out_tabular", "output path (--out_tabular) is required when output_format == tabular", ARGS.out_tabular)
    save_as_tabular_df(merged, ARGS.out_tabular)
    expected = ARGS.out_tabular

    # verify output exists and non-empty
    if not expected or not os.path.exists(expected) or os.path.getsize(expected) == 0:
        raise utils.DataErr(expected, "Output non creato o vuoto")

print("CustomDataGenerator: completed successfully")

if __name__ == '__main__':
    main()