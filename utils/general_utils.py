import math
import re
import sys
import csv
import pickle
import lxml.etree as ET

from enum import Enum
from itertools import count
from typing import Any, Callable, Dict, Generic, List, Literal, Optional, TypeVar, Union, Set, Tuple

import pandas as pd
import cobra
from cobra import Model as cobraModel, Reaction, Metabolite

import zipfile
import gzip
import bz2
from io import StringIO
import utils.rule_parsing  as rulesUtils
import utils.reaction_parsing as reactionUtils



class ValueErr(Exception):
    def __init__(self, param_name, expected, actual):
        super().__init__(f"Invalid value for {param_name}: expected {expected}, got {actual}")

class PathErr(Exception):
    def __init__(self, path, message):
        super().__init__(f"Path error for '{path}': {message}")

class FileFormat(Enum):
    """
    Encodes possible file extensions to conditionally save data in a different format.
    """
    DAT    = ("dat",) # this is how galaxy treats all your files!
    CSV    = ("csv",) # this is how most editable input data is written
    TSV    = ("tsv",) # this is how most editable input data is ACTUALLY written TODO:more support pls!!
    SVG    = ("svg",) # this is how most metabolic maps are written
    PNG    = ("png",) # this is a common output format for images (such as metabolic maps)
    PDF    = ("pdf",) # this is also a common output format for images, as it's required in publications.
    
    # Updated to include compressed variants
    XML    = ("xml", "xml.gz", "xml.zip", "xml.bz2") # SBML files are XML files, sometimes compressed
    JSON   = ("json", "json.gz", "json.zip", "json.bz2") # COBRA models can be stored as JSON files, sometimes compressed
    MAT    = ("mat", "mat.gz", "mat.zip", "mat.bz2") # COBRA models can be stored as MAT files, sometimes compressed
    YML    = ("yml", "yml.gz", "yml.zip", "yml.bz2") # COBRA models can be stored as YML files, sometimes compressed

    TXT    = ("txt",) # this is how most output data is written
    PICKLE = ("pickle", "pk", "p") # this is how all runtime data structures are saved

    def __init__(self, *extensions):
        self.extensions = extensions
        # Store original extension when set via fromExt
        self._original_extension = None

    @classmethod
    def fromExt(cls, ext: str) -> "FileFormat":
        """
        Converts a file extension string to a FileFormat instance.
        Args:
            ext : The file extension as a string.
        Returns:
            FileFormat: The FileFormat instance corresponding to the file extension.
        """
        variantName = ext.upper()
        if variantName in FileFormat.__members__: 
            instance = FileFormat[variantName]
            instance._original_extension = ext
            return instance
        
        variantName = ext.lower()
        for member in cls:
            if variantName in member.value: 
                # Create a copy-like behavior by storing the original extension
                member._original_extension = ext
                return member
        
        raise ValueErr("ext", "a valid FileFormat file extension", ext)

    def __str__(self) -> str:
        """
        (Private) converts to str representation. Good practice for usage with argparse.
        Returns:
            str : the string representation of the file extension.
        """
        # If we have an original extension stored (for compressed files only), use it
        if hasattr(self, '_original_extension') and self._original_extension:
            return self._original_extension
        
        # For XML, JSON, MAT and YML without original extension, use the base extension
        if self == FileFormat.XML:
            return "xml"
        elif self == FileFormat.JSON:
            return "json"
        elif self == FileFormat.MAT:
            return "mat"
        elif self == FileFormat.YML:
            return "yml"
        
        return self.value[-1]

class FilePath():
    """
    Represents a file path. View this as an attempt to standardize file-related operations by expecting
    values of this type in any process requesting a file path.
    """
    def __init__(self, filePath: str, ext: FileFormat, *, prefix="") -> None:
        """
        (Private) Initializes an instance of FilePath.
        Args:
            path : the end of the path, containing the file name.
            ext : the file's extension.
            prefix : anything before path, if the last '/' isn't there it's added by the code.
        Returns:
            None : practically, a FilePath instance.
        """
        self.ext = ext
        self.filePath = filePath

        if prefix and prefix[-1] != '/': 
            prefix += '/'
        self.prefix = prefix
    
    @classmethod
    def fromStrPath(cls, path: str) -> "FilePath":
        """
        Factory method to parse a string from which to obtain, if possible, a valid FilePath instance.
        It detects double extensions such as .json.gz and .xml.bz2, which are common in COBRA models.
        These double extensions are not supported for other file types such as .csv.
        Args:
            path : the string containing the path
        Raises:
            PathErr : if the provided string doesn't represent a valid path.
        Returns:
            FilePath : the constructed instance.
        """
        result = re.search(r"^(?P<prefix>.*\/)?(?P<name>.*)\.(?P<ext>[^.]*)$", path)
        if not result or not result["name"] or not result["ext"]:
            raise PathErr(path, "cannot recognize folder structure or extension in path")

        prefix = result["prefix"] if result["prefix"] else ""
        name, ext = result["name"], result["ext"]

        # Check for double extensions (json.gz, xml.zip, etc.)
        parts = path.split(".")
        if len(parts) >= 3:  
            penultimate = parts[-2]
            last = parts[-1]
            double_ext = f"{penultimate}.{last}"
            
            # Try the double extension first
            try:
                ext_format = FileFormat.fromExt(double_ext)
                name = ".".join(parts[:-2])
                # Extract prefix if it exists
                if '/' in name:
                    prefix = name[:name.rfind('/') + 1]
                    name = name[name.rfind('/') + 1:]
                return cls(name, ext_format, prefix=prefix)
            except ValueErr:
                # If double extension doesn't work, fall back to single extension
                pass

        # Single extension fallback (original logic)
        try:
            ext_format = FileFormat.fromExt(ext)
            return cls(name, ext_format, prefix=prefix)
        except ValueErr:
            raise PathErr(path, f"unsupported file extension: {ext}")

    def show(self) -> str:
        """
        Shows the path as a string.
        Returns:
            str : the path shown as a string.
        """
        return f"{self.prefix}{self.filePath}.{self.ext}"
    
    def __str__(self) -> str: 
        return self.show()

# ERRORS
def terminate(msg :str) -> None:
    """
    Terminate the execution of the script with an error message.
    
    Args:
        msg (str): The error message to be displayed.
    
    Returns:
        None
    """
    sys.exit(f"Execution aborted: {msg}\n")

def logWarning(msg :str, loggerPath :str) -> None:
    """
    Log a warning message to an output log file and print it to the console. The final period and a
    newline is added by the function.

    Args:
        s (str): The warning message to be logged and printed.
        loggerPath : The file path of the output log file. Given as a string, parsed to a FilePath and
        immediately read back (beware relative expensive operation, log with caution).

    Returns:
        None
    """
    # building the path and then reading it immediately seems useless, but it's actually a way of
    # validating that reduces repetition on the caller's side. Besides, logging a message by writing
    # to a file is supposed to be computationally expensive anyway, so this is also a good deterrent from
    # mindlessly logging whenever something comes up, log at the very end and tell the user everything
    # that went wrong. If you don't like it: implement a persistent runtime buffer that gets dumped to
    # the file only at the end of the program's execution.
    with open(FilePath.fromStrPath(loggerPath).show(), 'a') as log: log.write(f"{msg}.\n")

class CustomErr(Exception):
    """
    Custom error class to handle exceptions in a structured way, with a unique identifier and a message.
    """
    __idGenerator = count()
    errName = "Custom Error"
    def __init__(self, msg :str, details = "", explicitErrCode = -1) -> None:
        """
        (Private) Initializes an instance of CustomErr.

        Args:
            msg (str): Error message to be displayed.
            details (str): Informs the user more about the error encountered. Defaults to "".
            explicitErrCode (int): Explicit error code to be used. Defaults to -1.
        
        Returns:
            None : practically, a CustomErr instance.
        """
        self.msg     = msg
        self.details = details

        self.id = max(explicitErrCode, next(CustomErr.__idGenerator))

    def throw(self, loggerPath = "") -> None:
        """
        Raises the current CustomErr instance, logging a warning message before doing so.

        Raises:
            self: The current CustomErr instance.
        
        Returns:
            None
        """
        if loggerPath: logWarning(str(self), loggerPath)
        raise self

    def abort(self) -> None:
        """
        Aborts the execution of the script.
        
        Returns:
            None
        """
        terminate(str(self))

    def __str__(self) -> str:
        """
        (Private) Returns a string representing the current CustomErr instance.

        Returns:
            str: A string representing the current CustomErr instance.
        """
        return f"{CustomErr.errName} #{self.id}: {self.msg}, {self.details}."

class ArgsErr(CustomErr):
    """
    CustomErr subclass for UI arguments errors.
    """
    errName = "Args Error"
    def __init__(self, argName :str, expected :Any, actual :Any, msg = "no further details provided") -> None:
        super().__init__(f"argument \"{argName}\" expected {expected} but got {actual}", msg)

class DataErr(CustomErr):
    """
    CustomErr subclass for data formatting errors.
    """
    errName = "Data Format Error"
    def __init__(self, fileName :str, msg = "no further details provided") -> None:
        super().__init__(f"file \"{fileName}\" contains malformed data", msg)

class PathErr(CustomErr):
    """
    CustomErr subclass for filepath formatting errors.
    """
    errName = "Path Error"
    def __init__(self, path :FilePath, msg = "no further details provided") -> None:
        super().__init__(f"path \"{path}\" is invalid", msg)

class ValueErr(CustomErr):
    """
    CustomErr subclass for any value error.
    """
    errName = "Value Error"
    def __init__(self, valueName: str, expected :Any, actual :Any, msg = "no further details provided") -> None:
        super().__init__("value " + f"\"{valueName}\" " * bool(valueName) + f"was supposed to be {expected}, but got {actual} instead", msg)

# RESULT
T = TypeVar('T')
E = TypeVar('E', bound = CustomErr) # should bind to Result.ResultErr but python happened!
class Result(Generic[T, E]):
    class ResultErr(CustomErr):
        """
        CustomErr subclass for all Result errors.
        """
        errName = "Result Error"
        def __init__(self, msg = "no further details provided") -> None:
            super().__init__(msg)
    """
    Class to handle the result of an operation, with a value and a boolean flag to indicate
    whether the operation was successful or not.
    """
    def __init__(self, value :Union[T, E], isOk :bool) -> None:
        """
        (Private) Initializes an instance of Result.

        Args:
            value (Union[T, E]): The value to be stored in the Result instance.
            isOk (bool): A boolean flag to indicate whether the operation was successful or not.
        
            Returns:
                None : practically, a Result instance.
        """
        self.isOk  = isOk
        self.isErr = not isOk
        self.value = value

    @classmethod
    def Ok(cls,  value :T) -> "Result":
        """
        Constructs a new Result instance with a successful operation.

        Args:
            value (T): The value to be stored in the Result instance, set as successful.

        Returns:
            Result: A new Result instance with a successful operation.
        """
        return Result(value, isOk = True)
    
    @classmethod
    def Err(cls, value :E) -> "Result": 
        """
        Constructs a new Result instance with a failed operation.

        Args:
            value (E): The value to be stored in the Result instance, set as failed.

        Returns:
            Result: A new Result instance with a failed operation.
        """
        return Result(value, isOk = False)

    def unwrap(self) -> T:
        """
        Unwraps the value of the Result instance, if the operation was successful.

        Raises:
            ResultErr: If the operation was not successful.

        Returns:
            T: The value of the Result instance, if the operation was successful.
        """
        if self.isOk: return self.value
        raise Result.ResultErr(f"Unwrapped Result.Err : {self.value}")

    def unwrapOr(self, default :T) -> T:
        """
        Unwraps the value of the Result instance, if the operation was successful, otherwise
        it returns a default value.

        Args:
            default (T): The default value to be returned if the operation was not successful.

        Returns:
            T: The value of the Result instance, if the operation was successful,
            otherwise the default value.
        """
        return self.value if self.isOk else default
    
    def expect(self, err :"Result.ResultErr") -> T:
        """
        Expects that the value of the Result instance is successful, otherwise it raises an error.

        Args:
            err (Exception): The error to be raised if the operation was not successful.

        Raises:
            err: The error raised if the operation was not successful.

        Returns:
            T: The value of the Result instance, if the operation was successful.
        """
        if self.isOk: return self.value
        raise err

    U = TypeVar("U")
    def map(self, mapper: Callable[[T], U]) -> "Result[U, E]":
        """
        Maps the value of the current Result to whatever is returned by the mapper function.
        If the Result contained an unsuccessful operation to begin with it remains unchanged
        (a reference to the current instance is returned).
        If the mapper function panics the returned result instance will be of the error kind.

        Args:
            mapper (Callable[[T], U]): The mapper operation to be applied to the Result value.

        Returns:
            Result[U, E]: The result of the mapper operation applied to the Result value.
        """
        if self.isErr: return self
        try: return Result.Ok(mapper(self.value))
        except Exception as e: return Result.Err(e)
    
    D = TypeVar("D", bound = "Result.ResultErr")
    def mapErr(self, mapper :Callable[[E], D]) -> "Result[T, D]":
        """
        Maps the error of the current Result to whatever is returned by the mapper function.
        If the Result contained a successful operation it remains unchanged
        (a reference to the current instance is returned).
        If the mapper function panics this method does as well.

        Args:
            mapper (Callable[[E], D]): The mapper operation to be applied to the Result error.

        Returns:
            Result[U, E]: The result of the mapper operation applied to the Result error.
        """
        if self.isOk: return self
        return Result.Err(mapper(self.value))

    def __str__(self):
        return f"Result::{'Ok' if self.isOk else 'Err'}({self.value})"

# FILES
def read_dataset(path :FilePath, datasetName = "Dataset (not actual file name!)") -> pd.DataFrame:
    """
    Reads a .csv or .tsv file and returns it as a Pandas DataFrame.

    Args:
        path : the path to the dataset file.
        datasetName : the name of the dataset.

    Raises:
        DataErr: If anything goes wrong when trying to open the file, if pandas thinks the dataset is empty or if
        it has less than 2 columns.
    
    Returns:
        pandas.DataFrame: The dataset loaded as a Pandas DataFrame.
    """
    # I advise against the use of this function. This is an attempt at standardizing bad legacy code rather than
    # removing / replacing it to avoid introducing as many bugs as possible in the tools still relying on this code.
    # First off, this is not the best way to distinguish between .csv and .tsv files and Galaxy itself makes it really
    # hard to implement anything better. Also, this function's name advertizes it as a dataset-specific operation and
    # contains dubious responsibility (how many columns..) while being a file-opening function instead. My suggestion is
    # TODO: stop using dataframes ever at all in anything and find a way to have tight control over file extensions.
    try: dataset = pd.read_csv(path.show(), sep = '\t', header = None, engine = "python")
    except:
        try: dataset = pd.read_csv(path.show(), sep = ',', header = 0, engine = "python")
        except Exception as err: raise DataErr(datasetName, f"encountered empty or wrongly formatted data: {err}")
    
    if len(dataset.columns) < 2: raise DataErr(datasetName, "a dataset is always meant to have at least 2 columns")
    return dataset

def readPickle(path :FilePath) -> Any:
    """
    Reads the contents of a .pickle file, which needs to exist at the given path.

    Args:
        path : the path to the .pickle file.
    
    Returns:
        Any : the data inside a pickle file, could be anything.
    """
    with open(path.show(), "rb") as fd: return pickle.load(fd)

def writePickle(path :FilePath, data :Any) -> None:
    """
    Saves any data in a .pickle file, created at the given path.

    Args:
        path : the path to the .pickle file.
        data : the data to be written to the file.
    
    Returns:
        None
    """
    with open(path.show(), "wb") as fd: pickle.dump(data, fd)

def readCsv(path :FilePath, delimiter = ',', *, skipHeader = True) -> List[List[str]]:
    """
    Reads the contents of a .csv file, which needs to exist at the given path.

    Args:
        path : the path to the .csv file.
        delimiter : allows other subformats such as .tsv to be opened by the same method (\\t delimiter).
        skipHeader : whether the first row of the file is a header and should be skipped.
    
    Returns:
        List[List[str]] : list of rows from the file, each parsed as a list of strings originally separated by commas.
    """
    with open(path.show(), "r", newline = "") as fd: return list(csv.reader(fd, delimiter = delimiter))[skipHeader:]

def readSvg(path :FilePath, customErr :Optional[Exception] = None) -> ET.ElementTree:
    """
    Reads the contents of a .svg file, which needs to exist at the given path.

    Args:
        path : the path to the .svg file.
    
    Raises:
        DataErr : if the map is malformed.
    
    Returns:
        Any : the data inside a svg file, could be anything.
    """
    try: return ET.parse(path.show())
    except (ET.XMLSyntaxError, ET.XMLSchemaParseError) as err:
        raise customErr if customErr else err

def writeSvg(path :FilePath, data:ET.ElementTree) -> None:
    """
    Saves svg data opened with lxml.etree in a .svg file, created at the given path.

    Args:
        path : the path to the .svg file.
        data : the data to be written to the file.
    
    Returns:
        None
    """
    with open(path.show(), "wb") as fd: fd.write(ET.tostring(data))

# UI ARGUMENTS
class Bool:
    def __init__(self, argName :str) -> None:
        self.argName = argName

    def __call__(self, s :str) -> bool: return self.check(s)

    def check(self, s :str) -> bool:
        s = s.lower()
        if s == "true" : return True
        if s == "false": return False
        raise ArgsErr(self.argName, "boolean string (true or false, not case sensitive)", f"\"{s}\"")

class Float:
    def __init__(self, argName = "Dataset values, not an argument") -> None:
        self.argName = argName
    
    def __call__(self, s :str) -> float: return self.check(s)

    def check(self, s :str) -> float:
        try: return float(s)
        except ValueError:
            s = s.lower()
            if s == "nan" or s == "none": return math.nan
            raise ArgsErr(self.argName, "numeric string or \"None\" or \"NaN\" (not case sensitive)", f"\"{s}\"")

# MODELS
OldRule = List[Union[str, "OldRule"]]
class Model(Enum):
    """
    Represents a metabolic model, either custom or locally supported. Custom models don't point
    to valid file paths.
    """

    Recon   = "Recon"
    ENGRO2  = "ENGRO2"
    ENGRO2_no_legend = "ENGRO2_no_legend"
    HMRcore = "HMRcore"
    HMRcore_no_legend = "HMRcore_no_legend"
    Custom  = "Custom" # Exists as a valid variant in the UI, but doesn't point to valid file paths.

    def __raiseMissingPathErr(self, path :Optional[FilePath]) -> None:
        if not path: raise PathErr("<<MISSING>>", "it's necessary to provide a custom path when retrieving files from a custom model")

    def getRules(self, toolDir :str, customPath :Optional[FilePath] = None) -> Dict[str, Dict[str, OldRule]]:
        """
        Open "rules" file for this model.

        Returns:
            Dict[str, Dict[str, OldRule]] : the rules for this model.
        """
        path = customPath if self is Model.Custom else FilePath(f"{self.name}_rules", FileFormat.PICKLE, prefix = f"{toolDir}/local/pickle files/")
        self.__raiseMissingPathErr(path)
        return readPickle(path)
    
    def getTranslator(self, toolDir :str, customPath :Optional[FilePath] = None) -> Dict[str, Dict[str, str]]:
        """
        Open "gene translator (old: gene_in_rule)" file for this model.

        Returns:
            Dict[str, Dict[str, str]] : the translator dict for this model.
        """
        path = customPath if self is Model.Custom else FilePath(f"{self.name}_genes", FileFormat.PICKLE, prefix = f"{toolDir}/local/pickle files/")
        self.__raiseMissingPathErr(path)
        return readPickle(path)
    
    def getMap(self, toolDir = ".", customPath :Optional[FilePath] = None) -> ET.ElementTree:
        path = customPath if self is Model.Custom else FilePath(f"{self.name}_map", FileFormat.SVG, prefix = f"{toolDir}/local/svg metabolic maps/")
        self.__raiseMissingPathErr(path)
        return readSvg(path, customErr = DataErr(path, f"custom map in wrong format"))
    
    def getCOBRAmodel(self, toolDir = ".", customPath :Optional[FilePath] = None, customExtension :Optional[FilePath]=None)->cobra.Model:
        if(self is Model.Custom):
            return self.load_custom_model(customPath, customExtension)
        else:
            return cobra.io.read_sbml_model(FilePath(f"{self.name}", FileFormat.XML, prefix = f"{toolDir}/local/models/").show())
        
    def load_custom_model(self, file_path :FilePath, ext :Optional[FileFormat] = None) -> cobra.Model:
        ext = ext if ext else file_path.ext
        try:
            if str(ext) in FileFormat.XML.value:
                return cobra.io.read_sbml_model(file_path.show())
            
            if str(ext) in FileFormat.JSON.value:
                # Compressed files are not automatically handled by cobra
                if(ext == "json"):
                    return cobra.io.load_json_model(file_path.show())
                else: 
                    return self.extract_model(file_path, ext, "json")

            if str(ext) in FileFormat.MAT.value:
                # Compressed files are not automatically handled by cobra
                if(ext == "mat"):
                    return cobra.io.load_matlab_model(file_path.show())
                else: 
                    return self.extract_model(file_path, ext, "mat")

            if str(ext) in FileFormat.YML.value:
                # Compressed files are not automatically handled by cobra
                if(ext == "yml"):
                    return cobra.io.load_yaml_model(file_path.show())
                else: 
                    return self.extract_model(file_path, ext, "yml")

        except Exception as e: raise DataErr(file_path, e.__str__())
        raise DataErr(file_path,
            f"Fomat \"{file_path.ext}\" is not recognized, only JSON, XML, MAT and YAML (.yml) files are supported.")
    

    def extract_model(self, file_path:FilePath, ext :FileFormat, model_encoding:Literal["json", "mat", "yml"]) -> cobra.Model:
        """
        Extract JSON, MAT and YAML COBRA model from a compressed file (zip, gz, bz2).
        
        Args:
            file_path: File path of the model
            ext: File extensions of class FileFormat (should be .zip, .gz or .bz2)
            
        Returns:
            cobra.Model: COBRApy model 
            
        Raises:
            Exception: Extraction errors
        """
        ext_str = str(ext)

        try:
            if '.zip' in ext_str:
                with zipfile.ZipFile(file_path.show(), 'r') as zip_ref:
                    with zip_ref.open(zip_ref.namelist()[0]) as json_file:
                        content = json_file.read().decode('utf-8')
                        if model_encoding == "json":
                            return cobra.io.load_json_model(StringIO(content))
                        elif model_encoding == "mat":
                            return cobra.io.load_matlab_model(StringIO(content))
                        elif model_encoding == "yml":
                            return cobra.io.load_yaml_model(StringIO(content))
                        else:
                            raise ValueError(f"Unsupported model encoding: {model_encoding}. Supported: json, mat, yml")
            elif '.gz' in ext_str:
                with gzip.open(file_path.show(), 'rt', encoding='utf-8') as gz_ref:
                    if model_encoding == "json":
                        return cobra.io.load_json_model(gz_ref)
                    elif model_encoding == "mat":
                        return cobra.io.load_matlab_model(gz_ref)
                    elif model_encoding == "yml":
                        return cobra.io.load_yaml_model(gz_ref)
                    else:
                        raise ValueError(f"Unsupported model encoding: {model_encoding}. Supported: json, mat, yml")
            elif '.bz2' in ext_str:
                with bz2.open(file_path.show(), 'rt', encoding='utf-8') as bz2_ref:
                    if model_encoding == "json":
                        return cobra.io.load_json_model(bz2_ref)
                    elif model_encoding == "mat":
                        return cobra.io.load_matlab_model(bz2_ref)
                    elif model_encoding == "yml":
                        return cobra.io.load_yaml_model(bz2_ref)
                    else:
                        raise ValueError(f"Unsupported model encoding: {model_encoding}. Supported: json, mat, yml")
            else:
                raise ValueError(f"Compression format not supported: {ext_str}. Supported: .zip, .gz and .bz2")
            
        except Exception as e:
            raise Exception(f"Error during model extraction: {str(e)}")
        


    def __str__(self) -> str: return self.value


def convert_genes(model,annotation):
    from cobra.manipulation import rename_genes
    model2=model.copy()
    try:
        dict_genes={gene.id:gene.notes[annotation]  for gene in model2.genes}
    except:
        print("No annotation in gene dict!")
        return -1
    rename_genes(model2,dict_genes)

    return model2


def build_cobra_model_from_csv(csv_path: str, model_id: str = "new_model") -> cobra.Model:
    """
    Costruisce un modello COBRApy a partire da un file CSV con i dati delle reazioni.
    
    Args:
        csv_path: Path al file CSV (separato da tab)
        model_id: ID del modello da creare
        
    Returns:
        cobra.Model: Il modello COBRApy costruito
    """
    
    # Leggi i dati dal CSV
    df = pd.read_csv(csv_path, sep='\t')
    
    # Crea il modello vuoto
    model = cobraModel(model_id)
    
    # Dict per tenere traccia di metaboliti e compartimenti
    metabolites_dict = {}
    compartments_dict = {}
    
    print(f"Costruendo modello da {len(df)} reazioni...")
    
    # Prima passata: estrai metaboliti e compartimenti dalle formule delle reazioni
    for idx, row in df.iterrows():
        reaction_formula = str(row['Reaction']).strip()
        if not reaction_formula or reaction_formula == 'nan':
            continue
            
        # Estrai metaboliti dalla formula della reazione
        metabolites = extract_metabolites_from_reaction(reaction_formula)
        
        for met_id in metabolites:
            compartment = extract_compartment_from_metabolite(met_id)
            
            # Aggiungi compartimento se non esiste
            if compartment not in compartments_dict:
                compartments_dict[compartment] = compartment
            
            # Aggiungi metabolita se non esiste
            if met_id not in metabolites_dict:
                metabolites_dict[met_id] = Metabolite(
                    id=met_id,
                    compartment=compartment,
                    name=met_id.replace(f"_{compartment}", "").replace("__", "_")
                )
    
    # Aggiungi compartimenti al modello
    model.compartments = compartments_dict
    
    # Aggiungi metaboliti al modello  
    model.add_metabolites(list(metabolites_dict.values()))
    
    print(f"Aggiunti {len(metabolites_dict)} metaboliti e {len(compartments_dict)} compartimenti")
    
    # Seconda passata: aggiungi le reazioni
    reactions_added = 0
    reactions_skipped = 0
    
    for idx, row in df.iterrows():
        try:
            reaction_id = str(row['ReactionID']).strip()
            reaction_formula = str(row['Reaction']).strip()
            
            # Salta reazioni senza formula
            if not reaction_formula or reaction_formula == 'nan':
                reactions_skipped += 1
                continue
            
            # Crea la reazione
            reaction = Reaction(reaction_id)
            reaction.name = reaction_id
            
            # Imposta bounds
            reaction.lower_bound = float(row['lower_bound']) if pd.notna(row['lower_bound']) else -1000.0
            reaction.upper_bound = float(row['upper_bound']) if pd.notna(row['upper_bound']) else 1000.0
            
            # Aggiungi gene rule se presente
            if pd.notna(row['Rule']) and str(row['Rule']).strip():
                reaction.gene_reaction_rule = str(row['Rule']).strip()
            
            # Parse della formula della reazione
            try:
                parse_reaction_formula(reaction, reaction_formula, metabolites_dict)
            except Exception as e:
                print(f"Errore nel parsing della reazione {reaction_id}: {e}")
                reactions_skipped += 1
                continue
            
            # Aggiungi la reazione al modello
            model.add_reactions([reaction])
            reactions_added += 1
            
        except Exception as e:
            print(f"Errore nell'aggiungere la reazione {reaction_id}: {e}")
            reactions_skipped += 1
            continue
    
    print(f"Aggiunte {reactions_added} reazioni, saltate {reactions_skipped} reazioni")
    
    # Imposta l'obiettivo di biomassa
    set_biomass_objective(model)
    
    # Imposta il medium
    set_medium_from_data(model, df)
    
    print(f"Modello completato: {len(model.reactions)} reazioni, {len(model.metabolites)} metaboliti")
    
    return model


# Estrae tutti gli ID metaboliti nella formula (gestisce prefissi numerici + underscore)
def extract_metabolites_from_reaction(reaction_formula: str) -> Set[str]:
    """
    Estrae gli ID dei metaboliti da una formula di reazione.
    Pattern robusto: cattura token che terminano con _<compartimento> (es. _c, _m, _e)
    e permette che comincino con cifre o underscore.
    """
    metabolites = set()
    # coefficiente opzionale seguito da un token che termina con _<letters>
    pattern = r'(?:\d+(?:\.\d+)?\s+)?([A-Za-z0-9_]+_[a-z]+)'
    matches = re.findall(pattern, reaction_formula)
    metabolites.update(matches)
    return metabolites


def extract_compartment_from_metabolite(metabolite_id: str) -> str:
    """
    Estrae il compartimento dall'ID del metabolita.
    """
    # Il compartimento è solitamente l'ultima lettera dopo l'underscore
    if '_' in metabolite_id:
        return metabolite_id.split('_')[-1]
    return 'c'  # default cytoplasm


def parse_reaction_formula(reaction: Reaction, formula: str, metabolites_dict: Dict[str, Metabolite]):
    """
    Parsa una formula di reazione e imposta i metaboliti con i loro coefficienti.
    """

    if reaction.id == 'EX_thbpt_e':
        print(reaction.id)
        print(formula)
    # Dividi in parte sinistra e destra
    if '<=>' in formula:
        left, right = formula.split('<=>')
        reversible = True
    elif '<--' in formula:
        left, right = formula.split('<--')
        reversible = False
        left, right = left, right
    elif '-->' in formula:
        left, right = formula.split('-->')
        reversible = False
    elif '<-' in formula:
        left, right = formula.split('<-')
        reversible = False
        left, right = left, right
    else:
        raise ValueError(f"Formato reazione non riconosciuto: {formula}")
    
    # Parse dei metaboliti e coefficienti
    reactants = parse_metabolites_side(left.strip())
    products = parse_metabolites_side(right.strip())
    
    # Aggiungi metaboliti alla reazione
    metabolites_to_add = {}
    
    # Reagenti (coefficienti negativi)
    for met_id, coeff in reactants.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = -coeff
    
    # Prodotti (coefficienti positivi)
    for met_id, coeff in products.items():
        if met_id in metabolites_dict:
            metabolites_to_add[metabolites_dict[met_id]] = coeff
    
    reaction.add_metabolites(metabolites_to_add)


def parse_metabolites_side(side_str: str) -> Dict[str, float]:
    """
    Parsa un lato della reazione per estrarre metaboliti e coefficienti.
    """
    metabolites = {}
    if not side_str or side_str.strip() == '':
        return metabolites

    terms = side_str.split('+')
    for term in terms:
        term = term.strip()
        if not term:
            continue

        # pattern allineato: coefficiente opzionale + id che termina con _<compartimento>
        match = re.match(r'(?:(\d+\.?\d*)\s+)?([A-Za-z0-9_]+_[a-z]+)', term)
        if match:
            coeff_str, met_id = match.groups()
            coeff = float(coeff_str) if coeff_str else 1.0
            metabolites[met_id] = coeff

    return metabolites



def set_biomass_objective(model: Model):
    """
    Imposta la reazione di biomassa come obiettivo.
    """
    biomass_reactions = [r for r in model.reactions if 'biomass' in r.id.lower()]
    
    if biomass_reactions:
        model.objective = biomass_reactions[0].id
        print(f"Obiettivo impostato su: {biomass_reactions[0].id}")
    else:
        print("Nessuna reazione di biomassa trovata")


def set_medium_from_data(model: Model, df: pd.DataFrame):
    """
    Imposta il medium basato sulla colonna InMedium.
    """
    medium_reactions = df[df['InMedium'] == True]['ReactionID'].tolist()
    
    medium_dict = {}
    for rxn_id in medium_reactions:
        if rxn_id in [r.id for r in model.reactions]:
            reaction = model.reactions.get_by_id(rxn_id)
            if reaction.lower_bound < 0:  # Solo reazioni di uptake
                medium_dict[rxn_id] = abs(reaction.lower_bound)
    
    if medium_dict:
        model.medium = medium_dict
        print(f"Medium impostato con {len(medium_dict)} componenti")


def validate_model(model: Model) -> Dict[str, any]:
    """
    Valida il modello e fornisce statistiche di base.
    """
    validation = {
        'num_reactions': len(model.reactions),
        'num_metabolites': len(model.metabolites),
        'num_genes': len(model.genes),
        'num_compartments': len(model.compartments),
        'objective': str(model.objective),
        'medium_size': len(model.medium),
        'reversible_reactions': len([r for r in model.reactions if r.reversibility]),
        'exchange_reactions': len([r for r in model.reactions if r.id.startswith('EX_')]),
    }
    
    try:
        # Test di crescita
        solution = model.optimize()
        validation['growth_rate'] = solution.objective_value
        validation['status'] = solution.status
    except Exception as e:
        validation['growth_rate'] = None
        validation['status'] = f"Error: {e}"
    
    return validation


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



def generate_compartments(model: cobra.Model) -> pd.DataFrame:
    """
    Generates a DataFrame containing compartment information for each reaction.
    Creates columns for each compartment position (Compartment_1, Compartment_2, etc.)
    
    Args:
        model: the COBRA model to extract compartment data from.
        
    Returns:
        pd.DataFrame: DataFrame with ReactionID and compartment columns
    """
    pathway_data = []

    # First pass: determine the maximum number of pathways any reaction has
    max_pathways = 0
    reaction_pathways = {}

    for reaction in model.reactions:
        # Get unique pathways from all metabolites in the reaction
        if type(reaction.annotation['pathways']) == list:
            reaction_pathways[reaction.id] = reaction.annotation['pathways']
            max_pathways = max(max_pathways, len(reaction.annotation['pathways']))
        else:
            reaction_pathways[reaction.id] = [reaction.annotation['pathways']]

    # Create column names for pathways
    pathway_columns = [f"Pathway_{i+1}" for i in range(max_pathways)]

    # Second pass: create the data
    for reaction_id, pathways in reaction_pathways.items():
        row = {"ReactionID": reaction_id}
        
        # Fill pathway columns
        for i in range(max_pathways):
            col_name = pathway_columns[i]
            if i < len(pathways):
                row[col_name] = pathways[i]
            else:
                row[col_name] = None  # or "" if you prefer empty strings

        pathway_data.append(row)

    return pd.DataFrame(pathway_data)