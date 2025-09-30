"""
Convert a tabular (CSV/TSV/Tabular) description of a COBRA model into a COBRA file.

Supported output formats: SBML, JSON, MATLAB (.mat), YAML.
The script logs to a user-provided file for easier debugging in Galaxy.
"""

import os
import cobra
import argparse
from typing import List
import logging
import utils.model_utils as modelUtils

ARGS : argparse.Namespace
def process_args(args: List[str] = None) -> argparse.Namespace:
    """
    Parse command-line arguments for the CSV-to-COBRA conversion tool.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
    usage="%(prog)s [options]",
    description="Convert a tabular/CSV file to a COBRA model"
    )


    parser.add_argument("--out_log", type=str, required=True,
    help="Output log file")


    parser.add_argument("--input", type=str, required=True,
    help="Input tabular file (CSV/TSV)")


    parser.add_argument("--format", type=str, required=True, choices=["sbml", "json", "mat", "yaml"],
    help="Model format (SBML, JSON, MATLAB, YAML)")


    parser.add_argument("--output", type=str, required=True,
    help="Output model file path")

    parser.add_argument("--tool_dir", type=str, default=os.path.dirname(__file__),
    help="Tool directory (passed from Galaxy as $__tool_directory__)")


    return parser.parse_args(args)


###############################- ENTRY POINT -################################

def main(args: List[str] = None) -> None:
    """
    Entry point: parse arguments, build the COBRA model from a CSV/TSV file,
    and save it in the requested format.

    Returns:
        None
    """
    global ARGS
    ARGS = process_args(args)

    # configure logging to the requested log file (overwrite each run)
    logging.basicConfig(filename=ARGS.out_log,
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        filemode='w')

    logging.info('Starting fromCSVtoCOBRA tool')
    logging.debug('Args: input=%s format=%s output=%s tool_dir=%s', ARGS.input, ARGS.format, ARGS.output, ARGS.tool_dir)

    try:
        # Basic sanity checks
        if not os.path.exists(ARGS.input):
            logging.error('Input file not found: %s', ARGS.input)

        out_dir = os.path.dirname(os.path.abspath(ARGS.output))
        
        if out_dir and not os.path.isdir(out_dir):
            try:
                os.makedirs(out_dir, exist_ok=True)
                logging.info('Created missing output directory: %s', out_dir)
            except Exception as e:
                logging.exception('Cannot create output directory: %s', out_dir)

        model = modelUtils.build_cobra_model_from_csv(ARGS.input)
        
        
        logging.info('Created model with name: %s (ID: %s)', model.name, model.id)

        # Save model in requested format - Galaxy handles the filename
        if ARGS.format == "sbml":
            cobra.io.write_sbml_model(model, ARGS.output)
        elif ARGS.format == "json":
            cobra.io.save_json_model(model, ARGS.output)
        elif ARGS.format == "mat":
            cobra.io.save_matlab_model(model, ARGS.output)
        elif ARGS.format == "yaml":
            cobra.io.save_yaml_model(model, ARGS.output)
        else:
            logging.error('Unknown format requested: %s', ARGS.format)
            raise ValueError(f"Unknown format: {ARGS.format}")


        logging.info('Model successfully written to %s (format=%s)', ARGS.output, ARGS.format)
        print(f"Model created successfully in {ARGS.format.upper()} format")

    except Exception as e:
        # Log full traceback to the out_log so Galaxy users/admins can see what happened
        logging.exception('Unhandled exception in fromCSVtoCOBRA')
        print(f"ERROR: {str(e)}")
        raise


if __name__ == '__main__':
    main()
