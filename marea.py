from __future__ import division
import csv
from enum import Enum
import re
import sys
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
def process_args(args :List[str]) -> argparse.Namespace:
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
    return args

########################### warning ###########################################
def warning(s :str) -> None:
    """
    Gets a logger (.out_log) in the :Namespace object obtained from calling process_args, opens it in append mode and writes the input param (s) to it.

    Args:
        s : Always starting (in file) with "Warning: "

    Returns:
        None

    Side Effects:
        Edits logger file in the :Namespace object
    """
    args = process_args(sys.argv)
    with open(args.out_log, 'a') as log:
            log.write(s)
            
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

class ArrowColor(Enum):
    """
    Encodes possible arrow colors based on their meaning in the enrichment process.
    """
    Invalid       = "#BEBEBE" # gray, fold-change under treshold
    UpRegulated   = "#E41A1C" # red, up-regulated reaction
    DownRegulated = "#0000FF" # blue, down-regulated reaction

    def __str__(self) -> str: return self.value

class Arrow:
    """
    Models the properties of a reaction arrow that change based on enrichment.
    """
    MIN_W = 2
    MAX_W = 12
    
    def __init__(self, width :int, col: ArrowColor, isDashed = False) -> None:
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
    
    def applyTo(self, metabMap :ET.ElementTree, reactionId :str, mindReactionDir = True) -> None:
        """
        Applies the Arrow's styles to the actual reaction arrow in the given metabolic map.

        Args:
            metabMap : the metabolic map to edit.
            reactionId : the reaction ID associated with the arrow to style, as encoded in the dataset.
            mindReactionDir: if True the arrow's tips corresponding to a specific direction in the reaction will be styled,
            otherwise the body will.
        
        Side effects:
            metabMap : mut
        
        Returns:
            None
        """
        try: arrowEl :ET.Element = metabMap.xpath(
            f"//*[@id=\"{self.getMapReactionId(reactionId, mindReactionDir)}\"]")[0]
        
        except IndexError as err:
            ERRORS.append(reactionId)
            return
        
        currentStyles :str = arrowEl.get("style", "")
        if not re.search(r";stroke:[^;]+;stroke-width:[^;]+;stroke-dasharray:[^;]+$", currentStyles):
            arrowEl.set("style", currentStyles + self.toStyleStr())
            return # I don't bother to check much, as long as I don't extend the style tag at every edit
        
        arrowEl.set("style", ';'.join(currentStyles.split(';')[:-3]) + self.toStyleStr())
    
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

    def toStyleStr(self) -> str:
        """
        Collapses the styles of this Arrow into a str, ready to be applied as part of the "style" property on an svg element.

        Returns:
            str : the styles string.
        """
        return f";stroke:{self.col};stroke-width:{self.w};stroke-dasharray:{'5,5' if self.dash else 'none'}"

def applyRpsEnrichmentToMap(rpsEnrichmentRes :Dict[str, Tuple[float, FoldChange]], metabMap :ET.ElementTree, maxNumericFoldChange :float) -> None:
    """
    (Temporary) Applies RPS enrichment results to the provided metabolic map.

    Args:
        rpsEnrichmentRes : RPS enrichment results.
        metabMap : the metabolic map to edit.
        maxNumericFoldChange : biggest finite fold-change value found.
    
    Side effects:
        metabMap : mut
    
    Returns:
        None
    """
    for reactionId, (pValue, foldChange) in rpsEnrichmentRes.items():
        if isinstance(foldChange,str): foldChange = float(foldChange)
        if pValue >= ARGS.pValue: # pValue above tresh: dashed arrow
            Arrow(Arrow.MIN_W, ArrowColor.Invalid, isDashed = True).applyTo(metabMap, reactionId)
            continue

        if abs(foldChange) < math.log(ARGS.fChange, 2):
            Arrow(Arrow.MIN_W, ArrowColor.Invalid).applyTo(metabMap, reactionId)
            continue
        
        width = Arrow.MAX_W
        if not math.isinf(foldChange):
            try: width = max(abs(foldChange * Arrow.MAX_W) / maxNumericFoldChange, Arrow.MIN_W)
            except ZeroDivisionError: pass
        
        color = ArrowColor.DownRegulated if foldChange < 0 else ArrowColor.UpRegulated
        Arrow(width, color).applyTo(metabMap, reactionId)

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
        if not pd.isnull(classe):
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
            else:
                warning('Warning: no sample found in class ' + classe +
                        ', the class has been disregarded\n')
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

OldEnrichedScores = Dict[str, List[Union[float, FoldChange]]] #TODO: try to use Tuple whenever possible
def old_writeTabularResult(enrichedScores :OldEnrichedScores, outPath :utils.FilePath) -> None:
    tmp_csv = pd.DataFrame.from_dict(enrichedScores, orient = "index")
    tmp_csv = tmp_csv.reset_index()
    header = ['ids', 'P_Value', 'Log2(fold change)']
    tmp_csv.to_csv(outPath.show(), sep = '\t', index = False, header = header)

def new_writeTabularResult(enrichedScores : OldEnrichedScores, outPath :utils.FilePath) -> None:
    fieldNames = ("ids", "P_Value", "Log2(fold change)")
    with open(outPath.show(), "w", newline = "") as fd:
        writer = csv.DictWriter(fd, fieldnames = fieldNames)
        writer.writeheader()

        for reactId, [pValue, foldChange] in enrichedScores.items():
            writer.writerow({
                fieldNames[0] : reactId,
                fieldNames[1] : pValue,
                fieldNames[2] : foldChange
            }) #TODO: if you know a prettier way go for it! :D

def temp_enrichmentUpdate(tmp :Dict[str, List[Union[float, FoldChange]]], core_map :ET.ElementTree, max_F_C :float, dataset1Name :str, dataset2Name = "rest") -> None:
    new_writeTabularResult(
        tmp, buildOutputPath(dataset1Name, dataset2Name, details = "Tabular Result", ext = utils.FileFormat.TSV))
    
    if ARGS.using_RAS:
        fix_map(tmp, core_map, ARGS.pValue, ARGS.fChange, max_F_C)

    if not ARGS.using_RPS: return
    for reactId, enrichData in tmp.items(): tmp[reactId] = tuple(enrichData)
    applyRpsEnrichmentToMap(tmp, core_map, max_F_C)

def maps(core_map :ET.ElementTree, class_pat :Dict[str, List[List[float]]], ids :List[str]) -> None:
    """
    Compares clustered data based on a given comparison mode and generates metabolic maps visualizing the results.

    Args:
        core_map : SVG map to modify
        class_pat : the clustered data
        ids : ids for data association
        threshold_P_V : threshold for a p-value to be considered significant
        threshold_F_C : threshold for a fold change value to be considered significant
        comparison : comparison mode between clusters ("manyvsmany", "onevsrest", "onevsmany")
        control : another frontend-derived input parameter, identifying (I assume) the control sample

    Returns:
        None

    Raises:
        sys.exit : if there are less than 2 classes for comparison
    
    Side effects:
        core_map : mut (passed to fix_map)
    """
    class_pat = { k.strip() : v for k, v in class_pat.items() }
    if (not class_pat) or (len(class_pat.keys()) < 2): sys.exit('Execution aborted: classes provided for comparisons are less than two\n')
    #TODO: below is some repeated code that should be abstracted out
    if ARGS.comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            tmp :Dict[str, List[Union[float, FoldChange]]] = {}
            count = 0
            max_F_C = 0
            for l1, l2 in zip(class_pat.get(i), class_pat.get(j)):
                try:
                    stat_D, p_value = st.ks_2samp(l1, l2)
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            
            temp_enrichmentUpdate(tmp, core_map, max_F_C, i, j)
    
    elif ARGS.comparison == "onevsrest":
        for single_cluster in class_pat.keys():
            t :List[List[List[float]]] = []
            for k in class_pat.keys():
                if k != single_cluster:
                   t.append(class_pat.get(k))
            
            rest :List[List[float]] = []
            for i in t:
                rest = rest + i
            
            tmp = {}
            count = 0
            max_F_C = 0
            
            primo = -1
            for l1, l2 in zip(class_pat.get(single_cluster), rest):
                try:
                    stat_D, p_value = st.ks_2samp(l1, l2)
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if primo == -1:
                        primo = 0
                        print(avg)
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            
            temp_enrichmentUpdate(tmp, core_map, max_F_C, single_cluster)
                        
    elif ARGS.comparison == "onevsmany":
        controlItems = class_pat.get(ARGS.control)
        for otherDataset in class_pat.keys():
            if otherDataset == ARGS.control: continue
            
            tmp = {}
            count = 0
            max_F_C = 0
            for l1, l2 in zip(controlItems, class_pat.get(otherDataset)):
                try:
                    _, p_value = st.ks_2samp(l1, l2)
                    #sum(l1) da errore secondo me perchè ha null
                    avg = fold_change(sum(l1) / len(l1), sum(l2) / len(l2))
                    if not isinstance(avg, str):
                        if max_F_C < abs(avg):
                            max_F_C = abs(avg)
                    tmp[ids[count]] = [float(p_value), avg]
                    count += 1
                except (TypeError, ZeroDivisionError):
                    count += 1
            
            temp_enrichmentUpdate(tmp, core_map, max_F_C, ARGS.control, otherDataset)

def createOutputMaps(dataset1Name :str, dataset2Name :str, core_map :ET.ElementTree) -> None:
    svgFilePath = buildOutputPath(dataset1Name, dataset2Name, details = "SVG Map", ext = utils.FileFormat.SVG)
    utils.writeSvg(svgFilePath, core_map)

    if ARGS.generate_pdf:
        pngPath = buildOutputPath(dataset1Name, dataset2Name, details = "PNG Map", ext = utils.FileFormat.PNG)
        pdfPath = buildOutputPath(dataset1Name, dataset2Name, details = "PDF Map", ext = utils.FileFormat.PDF)
        convert_to_pdf(svgFilePath, pngPath, pdfPath)                     

    if not ARGS.generate_svg: os.remove(svgFilePath.show())

ClassPat = Dict[str, List[List[float]]]
def temp_RASorRPS(datasets, dataset, names) -> Tuple[List[str], ClassPat]:
    class_pat :ClassPat = {}
    if ARGS.option == 'datasets': return temp_optionDatasets(class_pat, datasets, names), class_pat
    return temp_optionDatasetClasses(class_pat, dataset, "RAS")

def temp_optionDatasets(class_pat, datasets, names) -> List[str]:
    num = 1
    for i, j in zip(datasets, names):
        name = name_dataset(j, num)
        resolve_rules_float, ids = temp_doCommons(i, name)
        
        if resolve_rules_float != None:
            class_pat[name] = list(map(list, zip(*resolve_rules_float.values())))
    
        num += 1
    
    return ids
    
def temp_optionDatasetClasses(class_pat, dataset, name) -> Tuple[List[str], ClassPat]:
    resolve_rules_float, ids = temp_doCommons(dataset, name)

    classes = read_dataset(ARGS.input_class, 'class')
    classes = classes.astype(str)
    if resolve_rules_float != None: class_pat = split_class(classes, resolve_rules_float)
    return ids, class_pat

def temp_doCommons(datasetPath :str, datasetName :str) -> Tuple[ClassPat, List[str]]:
    resolve_rules = read_dataset(datasetPath, datasetName)   
    resolve_rules.iloc[:, 0] = (resolve_rules.iloc[:, 0]).astype(str)
    ids = pd.Series.tolist(resolve_rules.iloc[:, 0])

    resolve_rules = resolve_rules.drop(resolve_rules.columns[[0]], axis=1)
    try: resolve_rules = resolve_rules.replace({'None': math.nan})
    except: pass #TODO: dataframes are acoustic but this is still bad, solution: opt out of dataframes before converting

    resolve_rules = resolve_rules.to_dict('list')
    return { k : list(map(float, v)) for k, v in resolve_rules.items() }, ids

def temp_writeAllFiles(core_map :ET.ElementTree, class_pat :ClassPat) -> None:
    if ARGS.comparison == "manyvsmany":
        for i, j in it.combinations(class_pat.keys(), 2):
            createOutputMaps(i, j, core_map)
        return
    
    if ARGS.comparison == "onevsrest":
        for single_cluster in class_pat.keys():
            createOutputMaps(single_cluster, "rest", core_map)
        
        return
    
    for otherDataset in class_pat.keys():
        if otherDataset == ARGS.control: continue
        createOutputMaps(i, j, core_map)

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
    ARGS = process_args(sys.argv)

    if os.path.isdir('result') == False: os.makedirs('result')
    
    core_map :ET.ElementTree = ARGS.choice_map.getMap(
        ARGS.tool_dir,
        utils.FilePath.fromStrPath(ARGS.custom_map) if ARGS.custom_map else None)
    # TODO: ^^^ ugly but fine for now, the argument is None if the model isn't custom because no file was given.
    # getMap will None-check the customPath and panic when the model IS custom but there's no file (good). A cleaner
    # solution can be derived from my comment in FilePath.fromStrPath

    if ARGS.using_RAS:
        ids, class_pat = temp_RASorRPS(ARGS.input_datas, ARGS.input_data, ARGS.names)
        maps(core_map, class_pat, ids)
    
    if ARGS.using_RPS:
        ids, class_pat = temp_RASorRPS(ARGS.input_datas_rps, ARGS.input_data_rps, ARGS.names_rps)
        maps(core_map, class_pat, ids)
    
    temp_writeAllFiles(core_map, class_pat)
    print('Execution succeded')

    if not ERRORS: return
    utils.logWarning(
        f"The following reaction IDs were mentioned in the dataset but weren't found in the map: {ERRORS}",
        ARGS.out_log)

###############################################################################

if __name__ == "__main__":
    main()
