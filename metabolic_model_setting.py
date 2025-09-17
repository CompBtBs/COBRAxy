"""
Scripts to generate a tabular file of a metabolic model (built-in or custom).

This script loads a COBRA model (built-in or custom), optionally applies
medium and gene nomenclature settings, derives reaction-related metadata
(GPR rules, formulas, bounds, objective coefficients, medium membership,
and compartments for ENGRO2), and writes a tabular summary.
"""

import os
import csv
import cobra
import argparse
import pandas as pd
import utils.general_utils as utils
from typing import Optional, Tuple, List
import utils.model_utils as modelUtils
import logging

ARGS : argparse.Namespace
def process_args(args: List[str] = None) -> argparse.Namespace:
    """
    Parse command-line arguments for metabolic_model_setting.
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
                        help="Medium selection option")

    parser.add_argument("--gene_format", type=str, default="Default",
                        help="Gene nomenclature format: Default (original), ENSNG, HGNC_SYMBOL, HGNC_ID, ENTREZ")
    
    parser.add_argument("--out_tabular", type=str,
                        help="Output file for the merged dataset (CSV or XLSX)")
    
    parser.add_argument("--tool_dir", type=str, default=os.path.dirname(__file__),
                        help="Tool directory (passed from Galaxy as $__tool_directory__)")


    return parser.parse_args(args)

################################- INPUT DATA LOADING -################################
def load_custom_model(file_path :utils.FilePath, ext :Optional[utils.FileFormat] = None) -> cobra.Model:
    """
    Loads a custom model from a file, either in JSON, XML, MAT, or YML format.

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

        if ext is utils.FileFormat.MAT:
            return cobra.io.load_matlab_model(file_path.show())

        if ext is utils.FileFormat.YML:
            return cobra.io.load_yaml_model(file_path.show())

    except Exception as e: raise utils.DataErr(file_path, e.__str__())
    raise utils.DataErr(
        file_path,
        f"Unrecognized format '{file_path.ext}'. Only JSON, XML, MAT, YML are supported."
    )


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
    """
    Save a pandas DataFrame as a tab-separated file, creating directories as needed.

    Args:
        df: The DataFrame to write.
        path: Destination file path (will be written as TSV).

    Raises:
        DataErr: If writing the output fails for any reason.

    Returns:
        None
    """
    try:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        df.to_csv(path, sep="\t", index=False)
    except Exception as e:
        raise utils.DataErr(path, f"failed writing tabular output: {e}")


###############################- ENTRY POINT -################################
def main(args:List[str] = None) -> None:
    """
    Initialize and generate custom data based on the frontend input arguments.
    
    Returns:
        None
    """
    # Parse args from frontend (Galaxy XML)
    global ARGS
    ARGS = process_args(args)


    if ARGS.input:
        # Load a custom model from file
        model = load_custom_model(
            utils.FilePath.fromStrPath(ARGS.input), utils.FilePath.fromStrPath(ARGS.name).ext)
    else:
        # Load a built-in model

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
    
    if ARGS.name == "ENGRO2" and ARGS.medium_selector != "Default":
        df_mediums = pd.read_csv(ARGS.tool_dir + "/local/medium/medium.csv", index_col = 0)
        ARGS.medium_selector = ARGS.medium_selector.replace("_", " ")
        medium = df_mediums[[ARGS.medium_selector]]
        medium = medium[ARGS.medium_selector].to_dict()

        # Reset all medium reactions lower bound to zero
        for rxn_id, _ in model.medium.items():
            model.reactions.get_by_id(rxn_id).lower_bound = float(0.0)
        
        # Apply selected medium uptake bounds (negative for uptake)
        for reaction, value in medium.items():
            if value is not None:
                model.reactions.get_by_id(reaction).lower_bound = -float(value)

    if (ARGS.name == "Recon" or ARGS.name == "ENGRO2") and ARGS.gene_format != "Default":
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        model = modelUtils.translate_model_genes(
            model=model,
            mapping_df= pd.read_csv(ARGS.tool_dir + "/local/mappings/genes_human.csv", dtype={'entrez_id': str}),
            target_nomenclature=ARGS.gene_format,
            source_nomenclature='HGNC_symbol',
            logger=logger
        )

    # generate data
    rules = modelUtils.generate_rules(model, asParsed = False)
    reactions = modelUtils.generate_reactions(model, asParsed = False)
    bounds = modelUtils.generate_bounds(model)
    medium = modelUtils.get_medium(model)
    objective_function = modelUtils.extract_objective_coefficients(model)
    
    if ARGS.name == "ENGRO2":
        compartments = modelUtils.generate_compartments(model)

    df_rules = pd.DataFrame(list(rules.items()), columns = ["ReactionID", "GPR"])
    df_reactions = pd.DataFrame(list(reactions.items()), columns = ["ReactionID", "Formula"])

    df_bounds = bounds.reset_index().rename(columns = {"index": "ReactionID"})
    df_medium = medium.rename(columns = {"reaction": "ReactionID"})
    df_medium["InMedium"] = True

    merged = df_reactions.merge(df_rules, on = "ReactionID", how = "outer")
    merged = merged.merge(df_bounds, on = "ReactionID", how = "outer")
    merged = merged.merge(objective_function, on = "ReactionID", how = "outer")
    if ARGS.name == "ENGRO2": 
        merged = merged.merge(compartments, on = "ReactionID", how = "outer")
    merged = merged.merge(df_medium, on = "ReactionID", how = "left")

    merged["InMedium"] = merged["InMedium"].fillna(False)

    merged = merged.sort_values(by = "InMedium", ascending = False)

    if not ARGS.out_tabular:
        raise utils.ArgsErr("out_tabular", "output path (--out_tabular) is required when output_format == tabular", ARGS.out_tabular)
    save_as_tabular_df(merged, ARGS.out_tabular)
    expected = ARGS.out_tabular

    # verify output exists and non-empty
    if not expected or not os.path.exists(expected) or os.path.getsize(expected) == 0:
        raise utils.DataErr(expected, "Output not created or empty")

    print("Metabolic_model_setting: completed successfully")

if __name__ == '__main__':

    main()
