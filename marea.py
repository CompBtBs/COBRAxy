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
import os
import argparse
import pyvips
import utils.general_utils as utils
from PIL import Image
from typing import Tuple, Union, Optional, List, Dict

ERRORS = []
########################## argparse ##########################################
ARGS :argparse.Namespace
def process_args() -> argparse.Namespace:
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
        default = '1vs1',
        choices = ['manyvsmany', 'onevsrest', 'onevsmany'])
    
    parser.add_argument(
        '-pv' ,'--pValue',
        type = float, 
        default = 0.1, 
        help = 'P-Value threshold (default: %(default)s)')
    
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
        '-mc',  '--choice_map',
        type = utils.Model, default = utils.Model.HMRcore,
        choices = [utils.Model.HMRcore, utils.Model.ENGRO2, utils.Model.Custom])

    args :argparse.Namespace = parser.parse_args()
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

############################ dataset name #####################################
def name_dataset(name_data :str, count :int) -> str:
    """
    Produces a unique name for a dataset based on what was provided by the user. The default name for any dataset is "Dataset", thus if the user didn't change it this function appends f"_{count}" to make it unique.

    Args:
        name_data : name associated with the dataset (from frontend input params)
        count : counter from 1 to make these names unique (external)

    Returns:
        str : the name made unique
    """
    if str(name_data) == 'Dataset':
        return str(name_data) + '_' + str(count)
    else:
        return str(name_data)

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
    elif avg1 == 0:
        return '-INF'
    elif avg2 == 0:
        return 'INF'
    else:
        return math.log(avg1 / avg2, 2)
    
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

# The type of d values is collapsed, losing precision, because the dict containst lists instead of tuples, please fix!
def fix_map(d :Dict[str, List[Union[float, FoldChange]]], core_map :ET.ElementTree, threshold_P_V :float, threshold_F_C :float, max_F_C :float) -> ET.ElementTree:
    """
    Edits the selected SVG map based on the p-value and fold change data (d) and some significance thresholds also passed as inputs.

    Args:
        d : dictionary mapping a p-value and a fold-change value (values) to each reaction ID as encoded in the SVG map (keys)
        core_map : SVG map to modify
        threshold_P_V : threshold for a p-value to be considered significant
        threshold_F_C : threshold for a fold change value to be considered significant
        max_F_C : highest fold change (absolute value)
    
    Returns:
        ET.ElementTree : the modified core_map

    Side effects:
        core_map : mut
    """
    maxT = 12
    minT = 2
    grey = '#BEBEBE'
    blue = '#0000FF'
    red = '#E41A1C'
    for el in core_map.iter():
        el_id = str(el.get('id'))
        if el_id.startswith('R_'):
            tmp = d.get(el_id[2:])
            if tmp != None:
                p_val :float = tmp[0]
                f_c = tmp[1]
                if p_val < threshold_P_V:
                    if not isinstance(f_c, str):
                        if abs(f_c) < math.log(threshold_F_C, 2):
                            col = grey
                            width = str(minT)
                        else:
                            if f_c < 0:
                                col = blue
                            elif f_c > 0:
                                col = red
                            width = str(max((abs(f_c) * maxT) / max_F_C, minT))
                    else:
                        if f_c == '-INF':
                            col = blue
                        elif f_c == 'INF':
                            col = red
                        width = str(maxT)
                    dash = 'none'
                else:
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
        currentStyles = ';'.join(currentStyles.split(';')[:-3])

    element.set("style", currentStyles + styleStr)

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
    elif ReactionDirection.fromReactionId(reactionId) is not ReactionDirection.Unknown: return reactionId[:-3:-1] + reactionId[:-2], ""
    return f"F_{reactionId}", f"B_{reactionId}"

class ArrowColor(Enum):
    """
    Encodes possible arrow colors based on their meaning in the enrichment process.
    """
    Invalid       = "#BEBEBE" # gray, fold-change under treshold
    UpRegulated   = "#E41A1C" # red, up-regulated reaction
    DownRegulated = "#0000FF" # blue, down-regulated reaction

    UpRegulatedInv = "#FF7A00"
    # ^^^ different shade of red (actually orange), up-regulated net value for a reversible reaction with
    # conflicting enrichment in the two directions.

    DownRegulatedInv = "#B22CF1"
    # ^^^ different shade of blue (actually purple), down-regulated net value for a reversible reaction with
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
        if downSizedForTips: width *= 0.15
        return f";stroke:{self.col};stroke-width:{width};stroke-dasharray:{'5,5' if self.dash else 'none'}"

# vvv These constants could be inside the class itself a static properties, but python
# was built by brainless organisms so here we are!
INVALID_ARROW = Arrow(Arrow.MIN_W, ArrowColor.Invalid)
INSIGNIFICANT_ARROW = Arrow(Arrow.MIN_W, ArrowColor.Invalid, isDashed = True)

def applyRpsEnrichmentToMap(rpsEnrichmentRes :Dict[str, Union[Tuple[float, FoldChange], Tuple[float, FoldChange, float, float]]], metabMap :ET.ElementTree, maxNumericFoldChange :float) -> None:
    """
    Applies RPS enrichment results to the provided metabolic map.

    Args:
        rpsEnrichmentRes : RPS enrichment results.
        metabMap : the metabolic map to edit.
        maxNumericFoldChange : biggest finite fold-change value found.
    
    Side effects:
        metabMap : mut
    
    Returns:
        None
    """
    for reactionId, values in rpsEnrichmentRes.items():
        pValue = values[0]
        foldChange = values[1]

        if isinstance(foldChange, str): foldChange = float(foldChange)
        if pValue >= ARGS.pValue: # pValue above tresh: dashed arrow
            INSIGNIFICANT_ARROW.styleReactionElements(metabMap, reactionId)
            continue

        if abs(foldChange) < math.log(ARGS.fChange, 2):
            INVALID_ARROW.styleReactionElements(metabMap, reactionId)
            continue
        
        width = Arrow.MAX_W
        if not math.isinf(foldChange):
            try: width = max(abs(foldChange * Arrow.MAX_W) / maxNumericFoldChange, Arrow.MIN_W)
            except ZeroDivisionError: pass
        
        if not reactionId.endswith("_RV"): # RV stands for reversible reactions
            Arrow(width, ArrowColor.fromFoldChangeSign(foldChange)).styleReactionElements(metabMap, reactionId)
            continue
        
        reactionId = reactionId[:-3] # Remove "_RV"
        
        inversionScore = (values[2] < 0) + (values[3] < 0) # Compacts the signs of averages into 1 easy to check score
        if inversionScore == 2: foldChange *= -1
        # ^^^ Style the inverse direction with the opposite sign netValue
        
        # If the score is 1 (opposite signs) we use alternative colors vvv
        arrow = Arrow(width, ArrowColor.fromFoldChangeSign(foldChange, useAltColor = inversionScore == 1))
        
        # vvv These 2 if statements can both be true and can both happen
        if ARGS.net: # style arrow head(s):
            arrow.styleReactionElements(metabMap, reactionId + ("_B" if inversionScore == 2 else "_F"))
        
        if not ARGS.using_RAS: # style arrow body
            arrow.styleReactionElements(metabMap, reactionId, mindReactionDir = False)

############################ split class ######################################
def split_class(classes :pd.DataFrame, resolve_rules :Dict[str, List[float]]) -> Dict[str, List[List[float]]]:
    """
    Generates a :dict that groups together data from a :DataFrame based on classes the data is related to.

    Args:
        classes : a :DataFrame of only string values, containing class information (rows) and keys to query the resolve_rules :dict
        resolve_rules : a :dict containing :float data

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
        for j in range(i, len(classes)):
            if classes.iloc[j, 1] == classe:
                pat_id :str = classes.iloc[j, 0]
                tmp = resolve_rules.get(pat_id, None)
                if tmp != None:
                    l.append(tmp)
                classes.iloc[j, 1] = None
        
        if l:
            class_pat[classe] = list(map(list, zip(*l)))
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

#funzione unica, lascio fuori i file e li passo in input
#conversion from png to pdf
def convert_png_to_pdf(png_file :utils.FilePath, pdf_file :utils.FilePath) -> None:
    """
    Internal utility to convert a PNG to PDF to aid from SVG conversion.

    Args:
        png_file : path to PNG file
        pdf_file : path to new PDF file

    Returns:
        None
    """
    image = Image.open(png_file.show())
    image = image.convert("RGB")
    image.save(pdf_file.show(), "PDF", resolution=100.0)

#function called to reduce redundancy in the code
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
        convert_png_to_pdf(file_png, file_pdf)
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
        prefix = "result")

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
def writeTabularResult(enrichedScores : OldEnrichedScores, ras_enrichment: bool, outPath :utils.FilePath) -> None:
    fieldNames = ["ids", "P_Value", "Log2(fold change)"]
    if not ras_enrichment: fieldNames.extend(["average_1", "average_2"])

    writeToCsv([ [reactId] + values for reactId, values in enrichedScores.items() ], fieldNames, outPath)

def temp_thingsInCommon(tmp :Dict[str, List[Union[float, FoldChange]]], core_map :ET.ElementTree, max_F_C :float, dataset1Name :str, dataset2Name = "rest", ras_enrichment = True) -> None:
    # this function compiles the things always in common between comparison modes after enrichment.
    # TODO: organize, name better.
    writeTabularResult(tmp, ras_enrichment, buildOutputPath(dataset1Name, dataset2Name, details = "Tabular Result", ext = utils.FileFormat.TSV))
    
    if ras_enrichment:
        fix_map(tmp, core_map, ARGS.pValue, ARGS.fChange, max_F_C)
        return

    for reactId, enrichData in tmp.items(): tmp[reactId] = tuple(enrichData)
    applyRpsEnrichmentToMap(tmp, core_map, max_F_C)

def computePValue(dataset1Data :List[float], dataset2Data :List[float]) -> float:
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
        float: P-value from a Kolmogorov-Smirnov test on the provided data.
    """
    return st.ks_2samp(dataset1Data, dataset2Data)[1]

def compareDatasetPair(dataset1Data :List[List[float]], dataset2Data :List[List[float]], ids :List[str]) -> Tuple[Dict[str, List[Union[float, FoldChange]]], float]:
    #TODO: the following code still suffers from "dumbvarnames-osis"
    tmp :Dict[str, List[Union[float, FoldChange]]] = {}
    count   = 0
    max_F_C = 0

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

                p_value = computePValue(nets1, nets2)
                avg1 = sum(nets1)   / len(nets1)
                avg2 = sum(nets2)   / len(nets2)
                net = (avg1 - avg2) / abs(avg2)
                
                if math.isnan(net): continue
                tmp[reactId[:-1] + "RV"] = [p_value, net, avg1, avg2]
                
                # vvv complementary directional ids are set to None once processed if net is to be applied to tips
                if ARGS.net:
                    ids[position] = None
                    continue

            # fallthrough is intended, regular scores need to be computed when tips aren't net but RAS datasets aren't used
            p_value = computePValue(l1, l2)
            avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
            if not isinstance(avg, str) and max_F_C < abs(avg): max_F_C = abs(avg)
            tmp[reactId] = [float(p_value), avg]
        
        except (TypeError, ZeroDivisionError): continue
    
    return tmp, max_F_C

def computeEnrichment(metabMap :ET.ElementTree, class_pat :Dict[str, List[List[float]]], ids :List[str], *, fromRAS = True) -> None:
    """
    Compares clustered data based on a given comparison mode and applies enrichment-based styling on the
    provided metabolic map.

    Args:
        metabMap : SVG map to modify.
        class_pat : the clustered data.
        ids : ids for data association.
        fromRAS : whether the data to enrich consists of RAS scores.

    Returns:
        None

    Raises:
        sys.exit : if there are less than 2 classes for comparison
    
    Side effects:
        metabMap : mut
        ids : mut
    """
    class_pat = { k.strip() : v for k, v in class_pat.items() }
    #TODO: simplfy this stuff vvv and stop using sys.exit (raise the correct utils error)
    if (not class_pat) or (len(class_pat.keys()) < 2): sys.exit('Execution aborted: classes provided for comparisons are less than two\n')

    if ARGS.comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            #TODO: these 2 functions are always called in pair and in this order and need common data,
            # some clever refactoring would be appreciated.
            comparisonDict, max_F_C = compareDatasetPair(class_pat.get(i), class_pat.get(j), ids)
            temp_thingsInCommon(comparisonDict, metabMap, max_F_C, i, j, fromRAS)
    
    elif ARGS.comparison == "onevsrest":
        for single_cluster in class_pat.keys():
            t :List[List[List[float]]] = []
            for k in class_pat.keys():
                if k != single_cluster:
                   t.append(class_pat.get(k))
            
            rest :List[List[float]] = []
            for i in t:
                rest = rest + i
            
            comparisonDict, max_F_C = compareDatasetPair(class_pat.get(single_cluster), rest, ids)
            temp_thingsInCommon(comparisonDict, metabMap, max_F_C, single_cluster, fromRAS)
    
    elif ARGS.comparison == "onevsmany":
        controlItems = class_pat.get(ARGS.control)
        for otherDataset in class_pat.keys():
            if otherDataset == ARGS.control: continue
            
            comparisonDict, max_F_C = compareDatasetPair(controlItems, class_pat.get(otherDataset), ids)
            temp_thingsInCommon(comparisonDict, metabMap, max_F_C, ARGS.control, otherDataset, fromRAS)

def createOutputMaps(dataset1Name :str, dataset2Name :str, core_map :ET.ElementTree) -> None:
    svgFilePath = buildOutputPath(dataset1Name, dataset2Name, details = "SVG Map", ext = utils.FileFormat.SVG)
    utils.writeSvg(svgFilePath, core_map)

    if ARGS.generate_pdf:
        pngPath = buildOutputPath(dataset1Name, dataset2Name, details = "PNG Map", ext = utils.FileFormat.PNG)
        pdfPath = buildOutputPath(dataset1Name, dataset2Name, details = "PDF Map", ext = utils.FileFormat.PDF)
        convert_to_pdf(svgFilePath, pngPath, pdfPath)                     

    if not ARGS.generate_svg: os.remove(svgFilePath.show())

ClassPat = Dict[str, List[List[float]]]
def getClassesAndIdsFromDatasets(datasetsPaths :List[str], datasetPath :str, classPath :str, names :List[str]) -> Tuple[List[str], ClassPat]:
    # TODO: I suggest creating dicts with ids as keys instead of keeping class_pat and ids separate,
    # for the sake of everyone's sanity.
    class_pat :ClassPat = {}
    if ARGS.option == 'datasets':
        num = 1 #TODO: the dataset naming function could be a generator
        for path, name in zip(datasetsPaths, names):
            name = name_dataset(name, num)
            resolve_rules_float, ids = getDatasetValues(path, name)
            if resolve_rules_float != None:
                class_pat[name] = list(map(list, zip(*resolve_rules_float.values())))
        
            num += 1
    
    elif ARGS.option == "dataset_class":
        classes = read_dataset(classPath, "class")
        classes = classes.astype(str)

        resolve_rules_float, ids = getDatasetValues(datasetPath, "Dataset Class (not actual name)")
        if resolve_rules_float != None: class_pat = split_class(classes, resolve_rules_float)
    
    return ids, class_pat
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
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.

    Returns:
        None
    
    Raises:
        sys.exit : if a user-provided custom map is in the wrong format (ET.XMLSyntaxError, ET.XMLSchemaParseError)
    """
    global ARGS
    ARGS = process_args()

    if os.path.isdir('result') == False: os.makedirs('result')
    
    core_map :ET.ElementTree = ARGS.choice_map.getMap(
        ARGS.tool_dir,
        utils.FilePath.fromStrPath(ARGS.custom_map) if ARGS.custom_map else None)
    # TODO: ^^^ ugly but fine for now, the argument is None if the model isn't custom because no file was given.
    # getMap will None-check the customPath and panic when the model IS custom but there's no file (good). A cleaner
    # solution can be derived from my comment in FilePath.fromStrPath

    if ARGS.using_RAS:
        ids, class_pat = getClassesAndIdsFromDatasets(ARGS.input_datas, ARGS.input_data, ARGS.input_class, ARGS.names)
        computeEnrichment(core_map, class_pat, ids)
    
    if ARGS.using_RPS:
        ids, class_pat = getClassesAndIdsFromDatasets(ARGS.input_datas_rps, ARGS.input_data_rps, ARGS.input_class_rps, ARGS.names_rps)
        computeEnrichment(core_map, class_pat, ids, fromRAS = False)
    
    # create output files: TODO: this is the same comparison happening in "maps", find a better way to organize this
    if ARGS.comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2): createOutputMaps(i, j, core_map)
        return
    
    if ARGS.comparison == "onevsrest":
        for single_cluster in class_pat.keys(): createOutputMaps(single_cluster, "rest", core_map)
        return
    
    for otherDataset in class_pat.keys():
        if otherDataset != ARGS.control: createOutputMaps(i, j, core_map)

    if not ERRORS: return
    utils.logWarning(
        f"The following reaction IDs were mentioned in the dataset but weren't found in the map: {ERRORS}",
        ARGS.out_log)
    
    print('Execution succeded')

###############################################################################
if __name__ == "__main__":
    main()