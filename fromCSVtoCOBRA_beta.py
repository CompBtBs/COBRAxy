import os
import csv
import cobra
import pickle
import argparse
import pandas as pd
import utils.general_utils as utils
from typing import Optional, Tuple, Union, List, Dict
import logging
import utils.rule_parsing  as rulesUtils
import utils.reaction_parsing as reactionUtils
import utils.model_utils as modelUtils

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

    parser.add_argument("--input", type=str,
                        help="Input tabular file")
    
    parser.add_argument("--format", type=str, required=True, choices=["sbml", "json", "mat", "yaml"],
                        help="Model format (SBML, JSON, MATLAB, YAML)")

    parser.add_argument("--output", type=str,
                        help="Output model file")
    
    parser.add_argument("--tool_dir", type=str, default=os.path.dirname(__file__),
                        help="Tool directory (passed from Galaxy as $__tool_directory__)")
    


    return parser.parse_args(args)

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

    model = modelUtils.build_cobra_model_from_csv(ARGS.model_upload)


    if ARGS.format == "sbml":
        cobra.io.write_sbml_model(model, ARGS.output)
    elif ARGS.format == "json":
        cobra.io.save_json_model(model, ARGS.output)
    elif ARGS.format == "mat":
        cobra.io.save_matlab_model(model, ARGS.output)
    elif ARGS.format == "yaml":
        cobra.io.save_yaml_model(model, ARGS.output)

if __name__ == '__main__':
    main()