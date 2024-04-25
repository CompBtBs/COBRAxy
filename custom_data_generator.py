import os
import csv
import cobra
import pickle
import argparse
import utils.general_utils as utils
import utils.rule_parsing  as rulesUtils
from typing import Optional, Tuple, Union, Dict

ARGS : argparse.Namespace
def process_args() -> argparse.Namespace:
    """
    Interfaces the script of a module with its frontend, making the user's choices for
    various parameters available as values in code.

    Args:
        args : Always obtained (in file) from sys.argv

    Returns:
        Namespace : An object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        usage = "%(prog)s [options]",
        description = "generate custom data from a given model")
    
    parser.add_argument("-ol", "--out_log",    type = str, required = True, help = "Output log")
    parser.add_argument("-id", "--input",      type = str, required = True, help = "Input model")
    parser.add_argument("-mn", "--name",       type = str, required = True, help = "Input model name")
    # ^ I need this because galaxy converts my files into .dat but I need to know what extension they were in

    parser.add_argument(
        "-of", "--output_format",
        type = utils.FileFormat,
        default = utils.FileFormat.PICKLE,
        choices = [utils.FileFormat.CSV, utils.FileFormat.PICKLE],
        # ^^^ Not all variants are valid here, otherwise list(utils.FileFormat) would be best.
        required = True,
        help = "Extension of all output files")
    
    argsNamespace = parser.parse_args()
    argsNamespace.out_dir = "result"
    # ^ can't get this one to work from xml, there doesn't seem to be a way to get the directory attribute from the collection

    return argsNamespace

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

def generate_reactions(model :cobra.Model) -> Dict[ReactionId, str]:
    """
    Generates a dictionary mapping reaction ids to reaction formulas from the model.

    Args:
        model : the model to derive data from.

    Returns:
        Dict[ReactionId, str] : the generated dictionary.
    """
    return {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction }

###############################- FILE SAVING -################################
def save_as_csv(data :dict, file_path :utils.FilePath, fieldNames :Tuple[str, str]) -> None:
    """
    Saves any dictionary-shaped data in a .csv file created at the given file_path.

    Args:
        data : the data to be written to the file.
        file_path : the path to the .csv file.
        fieldNames : the names of the fields (columns) in the .csv file.
    
    Returns:
        None
    """
    with open(file_path.show(), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fieldNames)
        writer.writeheader()

        for key, value in data.items():
            writer.writerow({ fieldNames[0] : key, fieldNames[1] : value })



###############################- ENTRY POINT -################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.
    
    Returns:
        None
    """
    # get args from frontend (related xml)
    global ARGS
    ARGS = process_args()

    # this is the worst thing I've seen so far, congrats to the former MaREA devs for suggesting this!
    if os.path.isdir(ARGS.out_dir) == False: os.makedirs(ARGS.out_dir)

    # load custom model
    model = load_custom_model(
        utils.FilePath.fromStrPath(ARGS.input), utils.FilePath.fromStrPath(ARGS.name).ext)
    
    # generate data and save it in the desired format and in a location galaxy understands
    # (it should show up as a collection in the history)
    reactions     = generate_reactions(model)
    rulesPath     = utils.FilePath("rules",     ARGS.output_format, prefix = ARGS.out_dir)
    reactionsPath = utils.FilePath("reactions", ARGS.output_format, prefix = ARGS.out_dir)

    if ARGS.output_format is utils.FileFormat.PICKLE:
        rules = generate_rules(model, asParsed = True)
        utils.writePickle(rulesPath,     rules)
        utils.writePickle(reactionsPath, reactions)
    
    elif ARGS.output_format is utils.FileFormat.CSV:
        rules = generate_rules(model, asParsed = False)
        save_as_csv(rules,     rulesPath,     ("ReactionID", "Rule"))
        save_as_csv(reactions, reactionsPath, ("ReactionID", "Reaction"))

    # ^ Please if anyone works on this after updating python to 3.12 change the if/elif into a match statement!!

if __name__ == '__main__':
    main()