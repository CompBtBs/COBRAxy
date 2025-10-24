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
from pathlib import Path


ARGS : argparse.Namespace
def process_args(args: List[str] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
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
                        help="Custom model file (JSON, XML, MAT, YAML)")
    parser.add_argument("--name", nargs='*', required=True,
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
def detect_file_format(file_path: str) -> utils.FileFormat:
    """
    Detect file format by examining file content and extension.
    Handles Galaxy .dat files by looking at content.
    """
    try:
        with open(file_path, 'r') as f:
            first_lines = ''.join([f.readline() for _ in range(5)])
        
        # Check for XML (SBML)
        if '<?xml' in first_lines or '<sbml' in first_lines:
            return utils.FileFormat.XML
        
        # Check for JSON
        if first_lines.strip().startswith('{'):
            return utils.FileFormat.JSON
            
        # Check for YAML
        if any(line.strip().endswith(':') for line in first_lines.split('\n')[:3]):
            return utils.FileFormat.YML
            
    except:
        pass
    
    # Fall back to extension-based detection
    if file_path.endswith('.xml') or file_path.endswith('.sbml'):
        return utils.FileFormat.XML
    elif file_path.endswith('.json'):
        return utils.FileFormat.JSON
    elif file_path.endswith('.mat'):
        return utils.FileFormat.MAT
    elif file_path.endswith('.yml') or file_path.endswith('.yaml'):
        return utils.FileFormat.YML
    
    # Default to XML for unknown extensions
    return utils.FileFormat.XML

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
    
def is_placeholder(gid) -> bool:
    """Return True if the gene id looks like a placeholder (e.g., 0/NA/NAN/empty)."""
    if gid is None:
        return True
    s = str(gid).strip().lower()
    return s in {"0", "", "na", "nan"}  # lowercase for simple matching

def sample_valid_gene_ids(genes, limit=10):
    """Yield up to `limit` valid gene IDs, skipping placeholders (e.g., the first 0 in RECON)."""
    out = []
    for g in genes:
        gid = getattr(g, "id", getattr(g, "gene_id", g))
        if not is_placeholder(gid):
            out.append(str(gid))
            if len(out) >= limit:
                break
    return out


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

    # Convert name from list to string (handles names with spaces)
    if isinstance(ARGS.name, list):
        ARGS.name = ' '.join(ARGS.name)

    if ARGS.input:
        # Load a custom model from file with auto-detected format
        detected_format = detect_file_format(ARGS.input)
        model = load_custom_model(utils.FilePath.fromStrPath(ARGS.input), detected_format)
    else:
        # Load a built-in model
        if not ARGS.model:
            raise utils.ArgsErr("model", "either --model or --input must be provided", "None")

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
    
    if ARGS.name == "ENGRO2" and ARGS.medium_selector != "Default":
        df_mediums = pd.read_csv(ARGS.tool_dir + "/local/medium/medium.csv", index_col = 0)
        #ARGS.medium_selector = ARGS.medium_selector.replace("_", " ") medium.csv uses underscores now
        medium = df_mediums[[ARGS.medium_selector]]
        medium = medium[ARGS.medium_selector].to_dict()

        # Reset all medium reactions lower bound to zero
        for rxn_id, _ in model.medium.items():
            model.reactions.get_by_id(rxn_id).lower_bound = float(0.0)
        
        # Apply selected medium uptake bounds (negative for uptake)
        for reaction, value in medium.items():
            if value is not None:
                model.reactions.get_by_id(reaction).lower_bound = -float(value)

    # Initialize translation_issues dictionary
    translation_issues = {}
    
    if (ARGS.name == "Recon" or ARGS.name == "ENGRO2") and ARGS.gene_format != "Default":
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        model, translation_issues = modelUtils.translate_model_genes(
            model=model,
            mapping_df= pd.read_csv(ARGS.tool_dir + "/local/mappings/genes_human.csv", dtype={'entrez_id': str}),
            target_nomenclature=ARGS.gene_format,
            source_nomenclature='HGNC_symbol',
            logger=logger
        )

    if ARGS.input and ARGS.gene_format != "Default":
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        # Take a small, clean sample of gene IDs (skipping placeholders like 0)
        ids_sample = sample_valid_gene_ids(model.genes, limit=10)
        if not ids_sample:
            raise utils.DataErr(
                "Custom_model",
                "No valid gene IDs found (many may be placeholders like 0)."
            )

        # Detect source nomenclature on the sample
        types = []
        for gid in ids_sample:
            try:
                t = modelUtils.gene_type(gid, "Custom_model")
            except Exception as e:
                # Keep it simple: skip problematic IDs
                logger.debug(f"gene_type failed for {gid}: {e}")
                t = None
            if t:
                types.append(t)

        if not types:
            raise utils.DataErr(
                "Custom_model",
                "Could not detect a known gene nomenclature from the sample."
            )

        unique_types = set(types)
        if len(unique_types) > 1:
            raise utils.DataErr(
                "Custom_model",
                "Mixed or inconsistent gene nomenclatures detected. "
                "Please unify them before converting."
            )

        source_nomenclature = types[0]

        # Convert only if needed
        if source_nomenclature != ARGS.gene_format:
            model, translation_issues = modelUtils.translate_model_genes(
                model=model,
                mapping_df= pd.read_csv(ARGS.tool_dir + "/local/mappings/genes_human.csv", dtype={'entrez_id': str}),
                target_nomenclature=ARGS.gene_format,
                source_nomenclature=source_nomenclature,
                logger=logger
            )

    # generate data using unified function
    if not ARGS.out_tabular:
        raise utils.ArgsErr("out_tabular", "output path (--out_tabular) is required when output_format == tabular", ARGS.out_tabular)
    
    merged = modelUtils.export_model_to_tabular(
        model=model,
        output_path=ARGS.out_tabular,
        translation_issues=translation_issues,
        include_objective=True,
        save_function=save_as_tabular_df
    )
    expected = ARGS.out_tabular

    # verify output exists and non-empty
    if not expected or not os.path.exists(expected) or os.path.getsize(expected) == 0:
        raise utils.DataErr(expected, "Output not created or empty")

    print("Completed successfully")

if __name__ == '__main__':

    main()
