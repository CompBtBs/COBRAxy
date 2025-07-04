from __future__ import division
import csv
from enum import Enum
import re
import sys
import numpy as np
import pandas as pd
import itertools as it
import scipy.stats as st
import lxml.etree as ET
import math
import utils.general_utils as utils
from PIL import Image
import os
import argparse
import pyvips
from typing import Tuple, Union, Optional, List, Dict
import copy

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

ERRORS = []
########################## argparse ##########################################
ARGS :argparse.Namespace
def process_args(args:List[str] = None) -> argparse.Namespace:
    """
    Interfaces the script of a module with its frontend, making the user's choices for various parameters available as values in code.

    Args:
        args : Always obtained (in file) from sys.argv

    Returns:
        Namespace : An object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        usage = "%(prog)s [options]",
        description = "process some value's genes to create a comparison's map.")
    
    #General:
    parser.add_argument(
        '-td', '--tool_dir',
        type = str,
        required = True,
        help = 'your tool directory')
    
    parser.add_argument('-on', '--control', type = str)
    parser.add_argument('-ol', '--out_log', help = "Output log")

    #Computation details:
    parser.add_argument(
        '-co', '--comparison',
        type = str, 
        default = 'manyvsmany',
        choices = ['manyvsmany', 'onevsrest', 'onevsmany'])

    parser.add_argument(
        '-te' ,'--test',
        type = str, 
        default = 'ks', 
        choices = ['ks', 'ttest_p', 'ttest_ind', 'wilcoxon', 'mw', 'DESeq'],
        help = 'Statistical test to use (default: %(default)s)')
    
    parser.add_argument(
        '-pv' ,'--pValue',
        type = float, 
        default = 0.1, 
        help = 'P-Value threshold (default: %(default)s)')

    parser.add_argument(
        '-adj' ,'--adjusted',
        type = utils.Bool("adjusted"), default = False, 
        help = 'Apply the FDR (Benjamini-Hochberg) correction (default: %(default)s)')
    
    parser.add_argument(
        '-fc', '--fChange',
        type = float, 
        default = 1.5, 
        help = 'Fold-Change threshold (default: %(default)s)')
    
    parser.add_argument(
        "-ne", "--net",
        type = utils.Bool("net"), default = False,
        help = "choose if you want net enrichment for RPS")

    parser.add_argument(
        '-op', '--option',
        type = str, 
        choices = ['datasets', 'dataset_class'],
        help='dataset or dataset and class')
    
    #RAS:
    parser.add_argument(
        "-ra", "--using_RAS",
        type = utils.Bool("using_RAS"), default = True,
        help = "choose whether to use RAS datasets.")

    parser.add_argument(
        '-id', '--input_data',
        type = str,
        help = 'input dataset')
    
    parser.add_argument(
        '-ic', '--input_class',
        type = str, 
        help = 'sample group specification')
    
    parser.add_argument(
        '-ids', '--input_datas',
        type = str,
        nargs = '+', 
        help = 'input datasets')
    
    parser.add_argument(
        '-na', '--names',
        type = str,
        nargs = '+', 
        help = 'input names')
    
    #RPS:
    parser.add_argument(
        "-rp", "--using_RPS",
        type = utils.Bool("using_RPS"), default = False,
        help = "choose whether to use RPS datasets.")
    
    parser.add_argument(
        '-idr', '--input_data_rps',
        type = str,
        help = 'input dataset rps')
    
    parser.add_argument(
        '-icr', '--input_class_rps', 
        type = str,
        help = 'sample group specification rps')
    
    parser.add_argument(
        '-idsr', '--input_datas_rps', 
        type = str,
        nargs = '+', 
        help = 'input datasets rps')
    
    parser.add_argument(
        '-nar', '--names_rps', 
        type = str,
        nargs = '+', 
        help = 'input names rps')
    
    #Output:
    parser.add_argument(
        "-gs", "--generate_svg",
        type = utils.Bool("generate_svg"), default = True,
        help = "choose whether to use RAS datasets.")
    
    parser.add_argument(
        "-gp", "--generate_pdf",
        type = utils.Bool("generate_pdf"), default = True,
        help = "choose whether to use RAS datasets.")
    
    parser.add_argument(
        '-cm', '--custom_map',
        type = str,
        help='custom map to use')
    
    parser.add_argument(
        '-idop', '--output_path', 
        type = str,
        default='result',
        help = 'output path for maps')
    
    parser.add_argument(
        '-mc',  '--choice_map',
        type = utils.Model, default = utils.Model.HMRcore,
        choices = [utils.Model.HMRcore, utils.Model.ENGRO2, utils.Model.Custom])

    args :argparse.Namespace = parser.parse_args(args)
    if args.using_RAS and not args.using_RPS: args.net = False

    return args
          
############################ dataset input ####################################
def read_dataset(data :str, name :str) -> pd.DataFrame:
    """
    Tries to read the dataset from its path (data) as a tsv and turns it into a DataFrame.

    Args:
        data : filepath of a dataset (from frontend input params or literals upon calling)
        name : name associated with the dataset (from frontend input params or literals upon calling)

    Returns:
        pd.DataFrame : dataset in a runtime operable shape
    
    Raises:
        sys.exit : if there's no data (pd.errors.EmptyDataError) or if the dataset has less than 2 columns
    """
    try:
        dataset = pd.read_csv(data, sep = '\t', header = 0, engine='python')
    except pd.errors.EmptyDataError:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    if len(dataset.columns) < 2:
        sys.exit('Execution aborted: wrong format of ' + name + '\n')
    return dataset

############################ map_methods ######################################
FoldChange = Union[float, int, str] # Union[float, Literal[0, "-INF", "INF"]]
def fold_change(avg1 :float, avg2 :float) -> FoldChange:
    """
    Calculates the fold change between two gene expression values.

    Args:
        avg1 : average expression value from one dataset avg2 : average expression value from the other dataset

    Returns:
        FoldChange :
            0 : when both input values are 0
            "-INF" : when avg1 is 0
            "INF" : when avg2 is 0
            float : for any other combination of values
    """
    if avg1 == 0 and avg2 == 0:
        return 0
    
    if avg1 == 0:
        return '-INF' # TODO: maybe fix
    
    if avg2 == 0:
        return 'INF'
    
    # (threshold_F_C - 1) / (abs(threshold_F_C) + 1) con threshold_F_C > 1
    return (avg1 - avg2) / (abs(avg1) + abs(avg2))

# TODO: I would really like for this one to get the Thanos treatment
def fix_style(l :str, col :Optional[str], width :str, dash :str) -> str:
    """
    Produces a "fixed" style string to assign to a reaction arrow in the SVG map, assigning style properties to the corresponding values passed as input params.

    Args:
        l : current style string of an SVG element
        col : new value for the "stroke" style property
        width : new value for the "stroke-width" style property
        dash : new value for the "stroke-dasharray" style property

    Returns:
        str : the fixed style string
    """
    tmp = l.split(';')
    flag_col = False
    flag_width = False
    flag_dash = False
    for i in range(len(tmp)):
        if tmp[i].startswith('stroke:'):
            tmp[i] = 'stroke:' + col
            flag_col = True
        if tmp[i].startswith('stroke-width:'):
            tmp[i] = 'stroke-width:' + width
            flag_width = True
        if tmp[i].startswith('stroke-dasharray:'):
            tmp[i] = 'stroke-dasharray:' + dash
            flag_dash = True
    if not flag_col:
        tmp.append('stroke:' + col)
    if not flag_width:
        tmp.append('stroke-width:' + width)
    if not flag_dash:
        tmp.append('stroke-dasharray:' + dash)
    return ';'.join(tmp)

# TODO: remove, there's applyRPS whatever
# The type of d values is collapsed, losing precision, because the dict containst lists instead of tuples, please fix!
def fix_map(d :Dict[str, List[Union[float, FoldChange]]], core_map :ET.ElementTree, threshold_P_V :float, threshold_F_C :float, max_z_score :float) -> ET.ElementTree:
    """
    Edits the selected SVG map based on the p-value and fold change data (d) and some significance thresholds also passed as inputs.

    Args:
        d : dictionary mapping a p-value and a fold-change value (values) to each reaction ID as encoded in the SVG map (keys)
        core_map : SVG map to modify
        threshold_P_V : threshold for a p-value to be considered significant
        threshold_F_C : threshold for a fold change value to be considered significant
        max_z_score : highest z-score (absolute value)
    
    Returns:
        ET.ElementTree : the modified core_map

    Side effects:
        core_map : mut
    """
    maxT = 12
    minT = 2
    grey = '#BEBEBE'
    blue = '#6495ed'
    red = '#ecac68'
    for el in core_map.iter():
        el_id = str(el.get('id'))
        if el_id.startswith('R_'):
            tmp = d.get(el_id[2:])
            if tmp != None:
                p_val, f_c, z_score, avg1, avg2 = tmp
                
                if math.isnan(p_val) or (isinstance(f_c, float) and math.isnan(f_c)): continue

                if p_val <= threshold_P_V: # p-value is OK
                    if not isinstance(f_c, str): # FC is finite
                        if abs(f_c) < ((threshold_F_C - 1) / (abs(threshold_F_C) + 1)): # FC is not OK
                            col = grey
                            width = str(minT)
                        else: # FC is OK
                            if f_c < 0:
                                col = blue
                            elif f_c > 0:
                                col = red
                            width = str(
                                min(
                                    max(abs(z_score * maxT) / max_z_score, minT),
                                    maxT))
                    
                    else: # FC is infinite
                        if f_c == '-INF':
                            col = blue
                        elif f_c == 'INF':
                            col = red
                        width = str(maxT)
                    dash = 'none'
                else: # p-value is not OK
                    dash = '5,5'
                    col = grey
                    width = str(minT)
                el.set('style', fix_style(el.get('style', ""), col, width, dash))
    return core_map

def getElementById(reactionId :str, metabMap :ET.ElementTree) -> utils.Result[ET.Element, utils.Result.ResultErr]:
    """
    Finds any element in the given map with the given ID. ID uniqueness in an svg file is recommended but
    not enforced, if more than one element with the exact ID is found only the first will be returned.

    Args:
        reactionId (str): exact ID of the requested element.
        metabMap (ET.ElementTree): metabolic map containing the element.

    Returns:
        utils.Result[ET.Element, ResultErr]: result of the search, either the first match found or a ResultErr.
    """
    return utils.Result.Ok(
        f"//*[@id=\"{reactionId}\"]").map(
        lambda xPath : metabMap.xpath(xPath)[0]).mapErr(
        lambda _ : utils.Result.ResultErr(f"No elements with ID \"{reactionId}\" found in map"))
        # ^^^ we shamelessly ignore the contents of the IndexError, it offers nothing to the user.

def styleMapElement(element :ET.Element, styleStr :str) -> None:
    currentStyles :str = element.get("style", "")
    if re.search(r";stroke:[^;]+;stroke-width:[^;]+;stroke-dasharray:[^;]+$", currentStyles):
        currentStyles = ';'.join(currentStyles.split(';')[:-3]) # TODO: why the last 3? Are we sure?

    #TODO: this is attempting to solve the styling override problem, not sure it does tho

    element.set("style", currentStyles + styleStr)

# TODO: maybe remove vvv
class ReactionDirection(Enum):
    Unknown = ""
    Direct  = "_F"
    Inverse = "_B"

    @classmethod
    def fromDir(cls, s :str) -> "ReactionDirection":
        # vvv as long as there's so few variants I actually condone the if spam:
        if s == ReactionDirection.Direct.value:  return ReactionDirection.Direct
        if s == ReactionDirection.Inverse.value: return ReactionDirection.Inverse
        return ReactionDirection.Unknown

    @classmethod
    def fromReactionId(cls, reactionId :str) -> "ReactionDirection":
        return ReactionDirection.fromDir(reactionId[-2:])

def getArrowBodyElementId(reactionId :str) -> str:
    if reactionId.endswith("_RV"): reactionId = reactionId[:-3] #TODO: standardize _RV
    elif ReactionDirection.fromReactionId(reactionId) is not ReactionDirection.Unknown: reactionId = reactionId[:-2]
    return f"R_{reactionId}"

def getArrowHeadElementId(reactionId :str) -> Tuple[str, str]:
    """
    We attempt extracting the direction information from the provided reaction ID, if unsuccessful we provide the IDs of both directions.

    Args:
        reactionId : the provided reaction ID.

    Returns:
        Tuple[str, str]: either a single str ID for the correct arrow head followed by an empty string or both options to try.
    """
    if reactionId.endswith("_RV"): reactionId = reactionId[:-3] #TODO: standardize _RV
    elif ReactionDirection.fromReactionId(reactionId) is not ReactionDirection.Unknown:
        return reactionId[:-3:-1] + reactionId[:-2], "" # ^^^ Invert _F to F_

    return f"F_{reactionId}", f"B_{reactionId}"

class ArrowColor(Enum):
    """
    Encodes possible arrow colors based on their meaning in the enrichment process.
    """
    Invalid       = "#BEBEBE" # gray, fold-change under treshold or not significant p-value
    Transparent   = "#ffffff00" # transparent, to make some arrow segments disappear
    UpRegulated   = "#ecac68" # orange, up-regulated reaction
    DownRegulated = "#6495ed" # lightblue, down-regulated reaction

    UpRegulatedInv = "#FF0000"
    # ^^^ bright red, up-regulated net value for a reversible reaction with
    # conflicting enrichment in the two directions.

    DownRegulatedInv = "#0000FF"
    # ^^^ bright blue, down-regulated net value for a reversible reaction with
    # conflicting enrichment in the two directions.

    @classmethod
    def fromFoldChangeSign(cls, foldChange :float, *, useAltColor = False) -> "ArrowColor":
        colors = (cls.DownRegulated, cls.DownRegulatedInv) if foldChange < 0 else (cls.UpRegulated, cls.UpRegulatedInv)
        return colors[useAltColor]

    def __str__(self) -> str: return self.value

class Arrow:
    """
    Models the properties of a reaction arrow that change based on enrichment.
    """
    MIN_W = 2
    MAX_W = 12

    def __init__(self, width :int, col: ArrowColor, *, isDashed = False) -> None:
        """
        (Private) Initializes an instance of Arrow.

        Args:
            width : width of the arrow, ideally to be kept within Arrow.MIN_W and Arrow.MAX_W (not enforced).
            col : color of the arrow.
            isDashed : whether the arrow should be dashed, meaning the associated pValue resulted not significant.
        
        Returns:
            None : practically, a Arrow instance.
        """
        self.w    = width
        self.col  = col
        self.dash = isDashed
    
    def applyTo(self, reactionId :str, metabMap :ET.ElementTree, styleStr :str) -> None:
        if getElementById(reactionId, metabMap).map(lambda el : styleMapElement(el, styleStr)).isErr:
            ERRORS.append(reactionId)

    def styleReactionElements(self, metabMap :ET.ElementTree, reactionId :str, *, mindReactionDir = True) -> None:
        # If We're dealing with RAS data or in general don't care about the direction of the reaction we only style the arrow body
        if not mindReactionDir:
            return self.applyTo(getArrowBodyElementId(reactionId), metabMap, self.toStyleStr())
        
        # Now we style the arrow head(s):
        idOpt1, idOpt2 = getArrowHeadElementId(reactionId)
        self.applyTo(idOpt1, metabMap, self.toStyleStr(downSizedForTips = True))
        if idOpt2: self.applyTo(idOpt2, metabMap, self.toStyleStr(downSizedForTips = True))
    
    # TODO: this seems to be unused, remove
    def getMapReactionId(self, reactionId :str, mindReactionDir :bool) -> str:
        """
        Computes the reaction ID as encoded in the map for a given reaction ID from the dataset.

        Args:
            reactionId: the reaction ID, as encoded in the dataset.
            mindReactionDir: if True forward (F_) and backward (B_) directions will be encoded in the result.
    
        Returns:
            str : the ID of an arrow's body or tips in the map.
        """
        # we assume the reactionIds also don't encode reaction dir if they don't mind it when styling the map.
        if not mindReactionDir: return "R_" + reactionId

        #TODO: this is clearly something we need to make consistent in RPS
        return (reactionId[:-3:-1] + reactionId[:-2]) if reactionId[:-2] in ["_F", "_B"] else f"F_{reactionId}" # "Pyr_F" --> "F_Pyr"

    def toStyleStr(self, *, downSizedForTips = False) -> str:
        """
        Collapses the styles of this Arrow into a str, ready to be applied as part of the "style" property on an svg element.

        Returns:
            str : the styles string.
        """
        width = self.w
        if downSizedForTips: width *= 0.8
        return f";stroke:{self.col};stroke-width:{width};stroke-dasharray:{'5,5' if self.dash else 'none'}"

# vvv These constants could be inside the class itself a static properties, but python
# was built by brainless organisms so here we are!
INVALID_ARROW       = Arrow(Arrow.MIN_W, ArrowColor.Invalid)
INSIGNIFICANT_ARROW = Arrow(Arrow.MIN_W, ArrowColor.Invalid, isDashed = True)
TRANSPARENT_ARROW   = Arrow(Arrow.MIN_W, ArrowColor.Transparent) # Who cares how big it is if it's transparent

# TODO: A more general version of this can be used for RAS as well, we don't need "fix map" or whatever
def applyRpsEnrichmentToMap(rpsEnrichmentRes :Dict[str, Union[Tuple[float, FoldChange], Tuple[float, FoldChange, float, float]]], metabMap :ET.ElementTree, maxNumericZScore :float) -> None:
    """
    Applies RPS enrichment results to the provided metabolic map.

    Args:
        rpsEnrichmentRes : RPS enrichment results.
        metabMap : the metabolic map to edit.
        maxNumericZScore : biggest finite z-score value found.
    
    Side effects:
        metabMap : mut
    
    Returns:
        None
    """
    for reactionId, values in rpsEnrichmentRes.items():
        pValue = values[0]
        foldChange = values[1]
        z_score = values[2]

        if math.isnan(pValue) or (isinstance(foldChange, float) and math.isnan(foldChange)): continue

        if isinstance(foldChange, str): foldChange = float(foldChange)
        if pValue >= ARGS.pValue: # pValue above tresh: dashed arrow
            INSIGNIFICANT_ARROW.styleReactionElements(metabMap, reactionId)
            continue

        if abs(foldChange) < (ARGS.fChange - 1) / (abs(ARGS.fChange) + 1):
            INVALID_ARROW.styleReactionElements(metabMap, reactionId)
            continue
        
        width = Arrow.MAX_W
        if not math.isinf(z_score):
            try: width = min(
                max(abs(z_score * Arrow.MAX_W) / maxNumericZScore, Arrow.MIN_W),
                Arrow.MAX_W)
            
            except ZeroDivisionError: pass
        
        if not reactionId.endswith("_RV"): # RV stands for reversible reactions
            Arrow(width, ArrowColor.fromFoldChangeSign(foldChange)).styleReactionElements(metabMap, reactionId)
            continue
        
        reactionId = reactionId[:-3] # Remove "_RV"
        
        inversionScore = (values[3] < 0) + (values[4] < 0) # Compacts the signs of averages into 1 easy to check score
        if inversionScore == 2: foldChange *= -1
        
        # If the score is 1 (opposite signs) we use alternative colors vvv
        arrow = Arrow(width, ArrowColor.fromFoldChangeSign(foldChange, useAltColor = inversionScore == 1))
        
        # vvv These 2 if statements can both be true and can both happen
        if ARGS.net: # style arrow head(s):
            arrow.styleReactionElements(metabMap, reactionId + ("_B" if inversionScore == 2 else "_F"))
        
        if not ARGS.using_RAS: # style arrow body
            arrow.styleReactionElements(metabMap, reactionId, mindReactionDir = False)

############################ split class ######################################
def split_class(classes :pd.DataFrame, dataset_values :Dict[str, List[float]]) -> Dict[str, List[List[float]]]:
    """
    Generates a :dict that groups together data from a :DataFrame based on classes the data is related to.

    Args:
        classes : a :DataFrame of only string values, containing class information (rows) and keys to query the resolve_rules :dict
        dataset_values : a :dict containing :float data

    Returns:
        dict : the dict with data grouped by class

    Side effects:
        classes : mut
    """
    class_pat :Dict[str, List[List[float]]] = {}
    for i in range(len(classes)):
        classe :str = classes.iloc[i, 1]
        if pd.isnull(classe): continue

        l :List[List[float]] = []
        sample_ids: List[str] = []

        for j in range(i, len(classes)):
            if classes.iloc[j, 1] == classe:
                pat_id :str = classes.iloc[j, 0] # sample name
                values = dataset_values.get(pat_id, None) # the column of values for that sample
                if values != None:
                    l.append(values)
                    sample_ids.append(pat_id)
                classes.iloc[j, 1] = None # TODO: problems?
        
        if l:
            class_pat[classe] = {
                "values": list(map(list, zip(*l))),  # trasposta
                "samples": sample_ids
            }
            continue
        
        utils.logWarning(
            f"Warning: no sample found in class \"{classe}\", the class has been disregarded", ARGS.out_log)
    
    return class_pat

############################ conversion ##############################################
#conversion from svg to png 
def svg_to_png_with_background(svg_path :utils.FilePath, png_path :utils.FilePath, dpi :int = 72, scale :int = 1, size :Optional[float] = None) -> None:
    """
    Internal utility to convert an SVG to PNG (forced opaque) to aid in PDF conversion.

    Args:
        svg_path : path to SVG file
        png_path : path for new PNG file
        dpi : dots per inch of the generated PNG
        scale : scaling factor for the generated PNG, computed internally when a size is provided
        size : final effective width of the generated PNG

    Returns:
        None
    """
    if size:
        image = pyvips.Image.new_from_file(svg_path.show(), dpi=dpi, scale=1)
        scale = size / image.width
        image = image.resize(scale)
    else:
        image = pyvips.Image.new_from_file(svg_path.show(), dpi=dpi, scale=scale)

    white_background = pyvips.Image.black(image.width, image.height).new_from_image([255, 255, 255])
    white_background = white_background.affine([scale, 0, 0, scale])

    if white_background.bands != image.bands:
        white_background = white_background.extract_band(0)

    composite_image = white_background.composite2(image, 'over')
    composite_image.write_to_file(png_path.show())

def convert_to_pdf(file_svg :utils.FilePath, file_png :utils.FilePath, file_pdf :utils.FilePath) -> None:
    """
    Converts the SVG map at the provided path to PDF.

    Args:
        file_svg : path to SVG file
        file_png : path to PNG file
        file_pdf : path to new PDF file

    Returns:
        None
    """
    svg_to_png_with_background(file_svg, file_png)
    try:
        image = Image.open(file_png.show())
        image = image.convert("RGB")
        image.save(file_pdf.show(), "PDF", resolution=100.0)
        print(f'PDF file {file_pdf.filePath} successfully generated.')
    
    except Exception as e:
        raise utils.DataErr(file_pdf.show(), f'Error generating PDF file: {e}')

############################ map ##############################################
def buildOutputPath(dataset1Name :str, dataset2Name = "rest", *, details = "", ext :utils.FileFormat) -> utils.FilePath:
    """
    Builds a FilePath instance from the names of confronted datasets ready to point to a location in the
    "result/" folder, used by this tool for output files in collections.

    Args:
        dataset1Name : _description_
        dataset2Name : _description_. Defaults to "rest".
        details : _description_
        ext : _description_

    Returns:
        utils.FilePath : _description_
    """
    # This function returns a util data structure but is extremely specific to this module.
    # RAS also uses collections as output and as such might benefit from a method like this, but I'd wait
    # TODO: until a third tool with multiple outputs appears before porting this to utils.
    return utils.FilePath(
        f"{dataset1Name}_vs_{dataset2Name}" + (f" ({details})" if details else ""),
        # ^^^ yes this string is built every time even if the form is the same for the same 2 datasets in
        # all output files: I don't care, this was never the performance bottleneck of the tool and
        # there is no other net gain in saving and re-using the built string.
        ext,
        prefix = ARGS.output_path)

FIELD_NOT_AVAILABLE = '/'
def writeToCsv(rows: List[list], fieldNames :List[str], outPath :utils.FilePath) -> None:
    fieldsAmt = len(fieldNames)
    with open(outPath.show(), "w", newline = "") as fd:
        writer = csv.DictWriter(fd, fieldnames = fieldNames, delimiter = '\t')
        writer.writeheader()
        
        for row in rows:
            sizeMismatch = fieldsAmt - len(row)
            if sizeMismatch > 0: row.extend([FIELD_NOT_AVAILABLE] * sizeMismatch)
            writer.writerow({ field : data for field, data in zip(fieldNames, row) })

OldEnrichedScores = Dict[str, List[Union[float, FoldChange]]] #TODO: try to use Tuple whenever possible
def temp_thingsInCommon(tmp :OldEnrichedScores, core_map :ET.ElementTree, max_z_score :float, dataset1Name :str, dataset2Name = "rest", ras_enrichment = True) -> None:
    # this function compiles the things always in common between comparison modes after enrichment.
    # TODO: organize, name better.
    suffix = "RAS" if ras_enrichment else "RPS"
    writeToCsv(
        [ [reactId] + values for reactId, values in tmp.items() ],
        ["ids", "P_Value", "fold change", "z-score", "average_1", "average_2"],
        buildOutputPath(dataset1Name, dataset2Name, details = f"Tabular Result ({suffix})", ext = utils.FileFormat.TSV))
    
    if ras_enrichment:
        fix_map(tmp, core_map, ARGS.pValue, ARGS.fChange, max_z_score)
        return

    for reactId, enrichData in tmp.items(): tmp[reactId] = tuple(enrichData)
    applyRpsEnrichmentToMap(tmp, core_map, max_z_score)

def computePValue(dataset1Data: List[float], dataset2Data: List[float]) -> Tuple[float, float]:
    """
    Computes the statistical significance score (P-value) of the comparison between coherent data
    from two datasets. The data is supposed to, in both datasets:
    - be related to the same reaction ID;
    - be ordered by sample, such that the item at position i in both lists is related to the
      same sample or cell line.

    Args:
        dataset1Data : data from the 1st dataset.
        dataset2Data : data from the 2nd dataset.

    Returns:
        tuple: (P-value, Z-score)
            - P-value from the selected test on the provided data.
            - Z-score of the difference between means of the two datasets.
    """
    match ARGS.test:
        case "ks":
            # Perform Kolmogorov-Smirnov test
            _, p_value = st.ks_2samp(dataset1Data, dataset2Data)
        case "ttest_p":
            # Datasets should have same size
            if len(dataset1Data) != len(dataset2Data):
                raise ValueError("Datasets must have the same size for paired t-test.")
            # Perform t-test for paired samples
            _, p_value = st.ttest_rel(dataset1Data, dataset2Data)
        case "ttest_ind":
            # Perform t-test for independent samples
            _, p_value = st.ttest_ind(dataset1Data, dataset2Data)
        case "wilcoxon":
            # Datasets should have same size
            if len(dataset1Data) != len(dataset2Data):
                raise ValueError("Datasets must have the same size for Wilcoxon signed-rank test.")
            # Perform Wilcoxon signed-rank test
            np.random.seed(42) # Ensure reproducibility since zsplit method is used
            _, p_value = st.wilcoxon(dataset1Data, dataset2Data, zero_method='zsplit')
        case "mw":
            # Perform Mann-Whitney U test
            _, p_value = st.mannwhitneyu(dataset1Data, dataset2Data)
        case _:
            p_value = np.nan # Default value if no valid test is selected
    
    # Calculate means and standard deviations
    mean1 = np.mean(dataset1Data)
    mean2 = np.mean(dataset2Data)
    std1 = np.std(dataset1Data, ddof=1)
    std2 = np.std(dataset2Data, ddof=1)
    
    n1 = len(dataset1Data)
    n2 = len(dataset2Data)
    
    # Calculate Z-score
    z_score = (mean1 - mean2) / np.sqrt((std1**2 / n1) + (std2**2 / n2))
    
    return p_value, z_score


def DESeqPValue(comparisonResult :Dict[str, List[Union[float, FoldChange]]], dataset1Data :List[List[float]], dataset2Data :List[List[float]], ids :List[str]) -> None:
    """
    Computes the p-value for each reaction in the comparisonResult dictionary using DESeq2.

    Args:
        comparisonResult : dictionary mapping a p-value and a fold-change value (values) to each reaction ID as encoded in the SVG map (keys)
        dataset1Data : data from the 1st dataset.
        dataset2Data : data from the 2nd dataset.
        ids : list of reaction IDs.

    Returns:
        None : mutates the comparisonResult dictionary in place with the p-values.
    """

    # pyDESeq2 needs at least 2 replicates per sample so I check this
    if len(dataset1Data[0]) < 2 or len(dataset2Data[0]) < 2:
        raise ValueError("Datasets must have at least 2 replicates each")

    # pyDESeq2 is based on pandas, so we need to convert the data into a DataFrame and clean it from NaN values
    dataframe1 = pd.DataFrame(dataset1Data, index=ids)
    dataframe2 = pd.DataFrame(dataset2Data, index=ids)
    
    # pyDESeq2 requires datasets to be samples x reactions and integer values
    dataframe1_clean = dataframe1.dropna(axis=0, how="any").T.astype(int)
    dataframe2_clean = dataframe2.dropna(axis=0, how="any").T.astype(int)
    dataframe1_clean.index = [f"ds1_rep{i+1}" for i in range(dataframe1_clean.shape[0])]
    dataframe2_clean.index = [f"ds2_rep{j+1}" for j in range(dataframe2_clean.shape[0])]

    # pyDESeq2 works on a DataFrame with values and another with infos about how samples are split (like dataset class)
    dataframe = pd.concat([dataframe1_clean, dataframe2_clean], axis=0)
    metadata = pd.DataFrame({"dataset": (["dataset1"]*dataframe1_clean.shape[0] + ["dataset2"]*dataframe2_clean.shape[0])}, index=dataframe.index)

    # Ensure the index of the metadata matches the index of the dataframe
    if not dataframe.index.equals(metadata.index):
        raise ValueError("The index of the metadata DataFrame must match the index of the counts DataFrame.")

    # Prepare and run pyDESeq2
    inference = DefaultInference()
    dds = DeseqDataSet(counts=dataframe, metadata=metadata, design="~dataset", inference=inference, quiet=True, low_memory=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["dataset", "dataset1", "dataset2"], inference=inference, quiet=True)
    ds.summary()

    # Retrieve the p-values from the DESeq2 results
    for reactId in ds.results_df.index:
        comparisonResult[reactId][0] = ds.results_df["pvalue"][reactId]


# TODO: the net RPS computation should be done in the RPS module
def compareDatasetPair(dataset1Data :List[List[float]], dataset2Data :List[List[float]], ids :List[str]) -> Tuple[Dict[str, List[Union[float, FoldChange]]], float, Dict[str, Tuple[np.ndarray, np.ndarray]]]:

    #TODO: the following code still suffers from "dumbvarnames-osis"
    netRPS :Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    comparisonResult :Dict[str, List[Union[float, FoldChange]]] = {}
    count   = 0
    max_z_score = 0

    for l1, l2 in zip(dataset1Data, dataset2Data):
        reactId = ids[count]
        count += 1
        if not reactId: continue # we skip ids that have already been processed

        try: #TODO: identify the source of these errors and minimize code in the try block
            reactDir = ReactionDirection.fromReactionId(reactId)
            # Net score is computed only for reversible reactions when user wants it on arrow tips or when RAS datasets aren't used
            if (ARGS.net or not ARGS.using_RAS) and reactDir is not ReactionDirection.Unknown:
                try: position = ids.index(reactId[:-1] + ('B' if reactDir is ReactionDirection.Direct else 'F'))
                except ValueError: continue # we look for the complementary id, if not found we skip

                nets1 = np.subtract(l1, dataset1Data[position])
                nets2 = np.subtract(l2, dataset2Data[position])
                netRPS[reactId] = (nets1, nets2)

                # Compute p-value and z-score for the RPS scores, if the pyDESeq option is set, p-values will be computed after and this function will return p_value = 0
                p_value, z_score = computePValue(nets1, nets2)
                avg1 = sum(nets1)   / len(nets1)
                avg2 = sum(nets2)   / len(nets2)
                net = fold_change(avg1, avg2)
                
                if math.isnan(net): continue
                comparisonResult[reactId[:-1] + "RV"] = [p_value, net, z_score, avg1, avg2]
                
                # vvv complementary directional ids are set to None once processed if net is to be applied to tips
                if ARGS.net: # If only using RPS, we cannot delete the inverse, as it's needed to color the arrows
                    ids[position] = None
                    continue

            # fallthrough is intended, regular scores need to be computed when tips aren't net but RAS datasets aren't used
            # Compute p-value and z-score for the RAS scores, if the pyDESeq option is set, p-values will be computed after and this function will return p_value = 0
            p_value, z_score = computePValue(l1, l2)
            avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
            # vvv TODO: Check numpy version compatibility
            if np.isfinite(z_score) and max_z_score < abs(z_score): max_z_score = abs(z_score)
            comparisonResult[reactId] = [float(p_value), avg, z_score, sum(l1) / len(l1), sum(l2) / len(l2)]
        
        except (TypeError, ZeroDivisionError): continue
    
    if ARGS.test == "DESeq":
        # Compute p-values using DESeq2
        DESeqPValue(comparisonResult, dataset1Data, dataset2Data, ids)

    # Apply multiple testing correction if set by the user
    if ARGS.adjusted:

        # Retrieve the p-values from the comparisonResult dictionary, they have to be different from NaN
        validPValues = [(reactId, result[0]) for reactId, result in comparisonResult.items() if not np.isnan(result[0])]
        # Unpack the valid p-values
        reactIds, pValues = zip(*validPValues)
        # Adjust the p-values using the Benjamini-Hochberg method
        adjustedPValues = st.false_discovery_control(pValues)
        # Update the comparisonResult dictionary with the adjusted p-values
        for reactId , adjustedPValue in zip(reactIds, adjustedPValues):
            comparisonResult[reactId][0] = adjustedPValue

    return comparisonResult, max_z_score, netRPS

def computeEnrichment(class_pat: Dict[str, List[List[float]]], ids: List[str], *, fromRAS=True) -> Tuple[List[Tuple[str, str, dict, float]], dict]:
    """
    Compares clustered data based on a given comparison mode and applies enrichment-based styling on the
    provided metabolic map.

    Args:
        class_pat : the clustered data.
        ids : ids for data association.
        fromRAS : whether the data to enrich consists of RAS scores.

    Returns:
        tuple: A tuple containing:
        - List[Tuple[str, str, dict, float]]: List of tuples with pairs of dataset names, comparison dictionary and max z-score.
        - dict : net RPS values for each dataset's reactions
    
    Raises:
        sys.exit : if there are less than 2 classes for comparison
    """
    class_pat = {k.strip(): v for k, v in class_pat.items()}
    if (not class_pat) or (len(class_pat.keys()) < 2):
        sys.exit('Execution aborted: classes provided for comparisons are less than two\n')
    
    # { datasetName : { reactId : netRPS, ... }, ... }
    netRPSResults :Dict[str, Dict[str, np.ndarray]] = {}
    enrichment_results = []

    if ARGS.comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            comparisonDict, max_z_score, netRPS = compareDatasetPair(class_pat.get(i), class_pat.get(j), ids)
            enrichment_results.append((i, j, comparisonDict, max_z_score))
            netRPSResults[i] = { reactId : net[0] for reactId, net in netRPS.items() }
            netRPSResults[j] = { reactId : net[1] for reactId, net in netRPS.items() }
    
    elif ARGS.comparison == "onevsrest":
        for single_cluster in class_pat.keys():
            rest = [item for k, v in class_pat.items() if k != single_cluster for item in v]
            comparisonDict, max_z_score, netRPS = compareDatasetPair(class_pat.get(single_cluster), rest, ids)
            enrichment_results.append((single_cluster, "rest", comparisonDict, max_z_score))
            netRPSResults[single_cluster] = { reactId : net[0] for reactId, net in netRPS.items() }
            netRPSResults["rest"]         = { reactId : net[1] for reactId, net in netRPS.items() }
    
    elif ARGS.comparison == "onevsmany":
        controlItems = class_pat.get(ARGS.control)
        for otherDataset in class_pat.keys():
            if otherDataset == ARGS.control:
                continue
            
            #comparisonDict, max_z_score, netRPS = compareDatasetPair(controlItems, class_pat.get(otherDataset), ids)
            comparisonDict, max_z_score, netRPS = compareDatasetPair(class_pat.get(otherDataset),controlItems, ids)
            #enrichment_results.append((ARGS.control, otherDataset, comparisonDict, max_z_score))
            enrichment_results.append(( otherDataset,ARGS.control, comparisonDict, max_z_score))
            netRPSResults[otherDataset] = { reactId : net[0] for reactId, net in netRPS.items() }
            netRPSResults[ARGS.control] = { reactId : net[1] for reactId, net in netRPS.items() }
    
    return enrichment_results, netRPSResults

def createOutputMaps(dataset1Name: str, dataset2Name: str, core_map: ET.ElementTree) -> None:
    svgFilePath = buildOutputPath(dataset1Name, dataset2Name, details="SVG Map", ext=utils.FileFormat.SVG)
    utils.writeSvg(svgFilePath, core_map)

    if ARGS.generate_pdf:
        pngPath = buildOutputPath(dataset1Name, dataset2Name, details="PNG Map", ext=utils.FileFormat.PNG)
        pdfPath = buildOutputPath(dataset1Name, dataset2Name, details="PDF Map", ext=utils.FileFormat.PDF)
        svg_to_png_with_background(svgFilePath, pngPath)
        try:
            image = Image.open(pngPath.show())
            image = image.convert("RGB")
            image.save(pdfPath.show(), "PDF", resolution=100.0)
            print(f'PDF file {pdfPath.filePath} successfully generated.')
        
        except Exception as e:
            raise utils.DataErr(pdfPath.show(), f'Error generating PDF file: {e}')

    if not ARGS.generate_svg: # This argument is useless, who cares if the user wants the svg or not
        os.remove(svgFilePath.show())

ClassPat = Dict[str, List[List[float]]]
def getClassesAndIdsFromDatasets(datasetsPaths :List[str], datasetPath :str, classPath :str, names :List[str]) -> Tuple[List[str], ClassPat, Dict[str, List[str]]]:
    # TODO: I suggest creating dicts with ids as keys instead of keeping class_pat and ids separate,
    # for the sake of everyone's sanity.
    columnNames :Dict[str, List[str]] = {} # { datasetName : [ columnName, ... ], ... }
    class_pat :ClassPat = {}
    if ARGS.option == 'datasets':
        num = 1
        for path, name in zip(datasetsPaths, names):
            name = str(name)
            if name == 'Dataset':
                name += '_' + str(num)
            
            values, ids = getDatasetValues(path, name)
            if values != None:
                class_pat[name]   = list(map(list, zip(*values.values()))) # TODO: ???
                columnNames[name] = ["Reactions", *values.keys()]
            
            num += 1
    
    elif ARGS.option == "dataset_class":
        classes = read_dataset(classPath, "class")
        classes = classes.astype(str)

        values, ids = getDatasetValues(datasetPath, "Dataset Class (not actual name)")
        if values != None:
            class_pat_with_samples_id = split_class(classes, values)

            for clas, values_and_samples_id in class_pat_with_samples_id.items():
                class_pat[clas] = values_and_samples_id["values"]
                columnNames[clas] = ["Reactions", *values_and_samples_id["samples"]]
    
    return ids, class_pat, columnNames
    #^^^ TODO: this could be a match statement over an enum, make it happen future marea dev with python 3.12! (it's why I kept the ifs)

#TODO: create these damn args as FilePath objects
def getDatasetValues(datasetPath :str, datasetName :str) -> Tuple[ClassPat, List[str]]:
    """
    Opens the dataset at the given path and extracts the values (expected nullable numerics) and the IDs.

    Args:
        datasetPath : path to the dataset
        datasetName (str): dataset name, used in error reporting

    Returns:
        Tuple[ClassPat, List[str]]: values and IDs extracted from the dataset
    """
    dataset = read_dataset(datasetPath, datasetName)
    IDs = pd.Series.tolist(dataset.iloc[:, 0].astype(str))

    dataset = dataset.drop(dataset.columns[0], axis = "columns").to_dict("list")
    return { id : list(map(utils.Float("Dataset values, not an argument"), values)) for id, values in dataset.items() }, IDs

############################ MAIN #############################################
def main(args:List[str] = None) -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    
    Raises:
        sys.exit : if a user-provided custom map is in the wrong format (ET.XMLSyntaxError, ET.XMLSchemaParseError)
    """
    global ARGS
    ARGS = process_args(args)

    # Create output folder
    if not os.path.isdir(ARGS.output_path):
        os.makedirs(ARGS.output_path, exist_ok=True)
    
    core_map: ET.ElementTree = ARGS.choice_map.getMap(
        ARGS.tool_dir,
        utils.FilePath.fromStrPath(ARGS.custom_map) if ARGS.custom_map else None)
    
    # TODO: in the future keep the indices WITH the data and fix the code below.

    # Prepare enrichment results containers
    ras_results = []
    rps_results = []

    # Compute RAS enrichment if requested
    if ARGS.using_RAS: #       vvv columnNames only matter with RPS data
        ids_ras, class_pat_ras, _ = getClassesAndIdsFromDatasets(
            ARGS.input_datas, ARGS.input_data, ARGS.input_class, ARGS.names)
        ras_results, _ = computeEnrichment(class_pat_ras, ids_ras, fromRAS=True)
        #           ^^^ netRPS only matter with RPS data

    # Compute RPS enrichment if requested
    if ARGS.using_RPS:
        ids_rps, class_pat_rps, columnNames = getClassesAndIdsFromDatasets(
            ARGS.input_datas_rps, ARGS.input_data_rps, ARGS.input_class_rps, ARGS.names_rps)
        
        rps_results, netRPS = computeEnrichment(class_pat_rps, ids_rps, fromRAS=False)

    # Organize by comparison pairs
    comparisons: Dict[Tuple[str, str], Dict[str, Tuple]] = {}
    for i, j, comparison_data, max_z_score in ras_results:
        comparisons[(i, j)] = {'ras': (comparison_data, max_z_score), 'rps': None}
    
    for i, j, comparison_data, max_z_score,  in rps_results:
        comparisons.setdefault((i, j), {}).update({'rps': (comparison_data, max_z_score)})

    # For each comparison, create a styled map with RAS bodies and RPS heads
    for (i, j), res in comparisons.items():
        map_copy = copy.deepcopy(core_map)

        # Apply RAS styling to arrow bodies
        if res.get('ras'):
            tmp_ras, max_z_ras = res['ras']
            temp_thingsInCommon(tmp_ras, map_copy, max_z_ras, i, j, ras_enrichment=True)

        # Apply RPS styling to arrow heads
        if res.get('rps'):
            tmp_rps, max_z_rps = res['rps']
            # applyRpsEnrichmentToMap styles only heads unless only RPS are used
            temp_thingsInCommon(tmp_rps, map_copy, max_z_rps, i, j, ras_enrichment=False)

        # Output both SVG and PDF/PNG as configured
        createOutputMaps(i, j, map_copy)
    
    # Add net RPS output file
    if ARGS.net or not ARGS.using_RAS:
        for datasetName, rows in netRPS.items():
            writeToCsv(
                [[reactId, *netValues] for reactId, netValues in rows.items()],
                # vvv In weird comparison modes the dataset names are not recorded properly..
                columnNames.get(datasetName, ["Reactions"]),
                utils.FilePath(
                    "Net_RPS_" + datasetName,
                    ext = utils.FileFormat.CSV,
                    prefix = ARGS.output_path))

    print('Execution succeeded')
###############################################################################
if __name__ == "__main__":
    main()
