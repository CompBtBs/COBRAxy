"""
General utilities for COBRAxy.

This module provides:
- File and path helpers (FileFormat, FilePath)
- Error and result handling utilities (CustomErr, Result)
- Basic I/O helpers (CSV/TSV, pickle, SVG)
- Lightweight CLI argument parsers (Bool, Float)
- Model loader utilities for COBRA models, including compressed formats
"""
import math
import re
import sys
import csv
import pickle
import lxml.etree as ET

from enum import Enum
from itertools import count
from typing import Any, Callable, Dict, Generic, List, Literal, Optional, TypeVar, Union, Tuple

import pandas as pd
import cobra

import zipfile
import gzip
import bz2
from io import StringIO


from typing import Any, Callable, Dict, Generic, List, Literal, Optional, TypeVar, Union, Tuple
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
    DAT    = ("dat",)
    CSV    = ("csv",)
    TSV    = ("tsv",)
    SVG    = ("svg",)
    PNG    = ("png",)
    PDF    = ("pdf",)

    # Compressed variants for common model formats
    XML    = ("xml", "xml.gz", "xml.zip", "xml.bz2")
    JSON   = ("json", "json.gz", "json.zip", "json.bz2")
    MAT    = ("mat", "mat.gz", "mat.zip", "mat.bz2")
    YML    = ("yml", "yml.gz", "yml.zip", "yml.bz2")

    TXT    = ("txt",)
    PICKLE = ("pickle", "pk", "p")

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
        if hasattr(self, '_original_extension') and self._original_extension:
            return self._original_extension
        
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
    Represents a file path with format-aware helpers.
    """
    def __init__(self, filePath: str, ext: FileFormat, *, prefix="") -> None:
        """
        Initialize FilePath.
        Args:
            path: File name stem.
            ext: File extension (FileFormat).
            prefix: Optional directory path (trailing '/' auto-added).
        """
        self.ext = ext
        self.filePath = filePath

        if prefix and prefix[-1] != '/': 
            prefix += '/'
        self.prefix = prefix
    
    @classmethod
    def fromStrPath(cls, path: str) -> "FilePath":
        """
    Parse a string path into a FilePath, supporting double extensions for models (e.g., .json.gz).
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

        parts = path.split(".")
        if len(parts) >= 3:  
            penultimate = parts[-2]
            last = parts[-1]
            double_ext = f"{penultimate}.{last}"
            
            try:
                ext_format = FileFormat.fromExt(double_ext)
                name = ".".join(parts[:-2])
                if '/' in name:
                    prefix = name[:name.rfind('/') + 1]
                    name = name[name.rfind('/') + 1:]
                return cls(name, ext_format, prefix=prefix)
            except ValueErr:
                pass

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
        msg (str): The warning message to be logged and printed.
        loggerPath : The file path of the output log file. Given as a string, parsed to a FilePath and
        immediately read back (beware relative expensive operation, log with caution).

    Returns:
        None
    """
    # Note: validates path via FilePath; keep logging minimal to avoid overhead.
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
        Raises the current CustomErr instance, optionally logging it first.

        Args:
            loggerPath (str): Optional path to a log file to append this error before raising.

        Raises:
            self: The current CustomErr instance.

        Returns:
            None
        """
        if loggerPath:
            logWarning(str(self), loggerPath)
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
        Initialize an instance of Result.

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
        Construct a successful Result.

        Args:
            value (T): The value to be stored in the Result instance, set as successful.

        Returns:
            Result: A new Result instance with a successful operation.
        """
        return Result(value, isOk = True)
    
    @classmethod
    def Err(cls, value :E) -> "Result": 
        """
        Construct a failed Result.

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

def findIdxByName(header: List[str], name: str, colName="name") -> Optional[int]:
    """
    Find the indices of the 'ReactionID' column and a user-specified column name
    within the header row of a tabular file.

    Args:
        header (List[str]): The header row, as a list of column names.
        name (str): The name of the column to look for (e.g. 'GPR').
        colName (str, optional): Label used in error messages for clarity. Defaults to "name".

    Returns:
        Tuple[int, int]: A tuple containing:
            - The index of the 'ReactionID' column.
            - The index of the requested column `name`.

    Raises:
        ValueError: If 'ReactionID' or the requested column `name` is not found in the header.

    Notes:
        Both 'ReactionID' and the requested column are mandatory for downstream processing.
    """

    col_index = {col_name: idx for idx, col_name in enumerate(header)}

    if name not in col_index or "ReactionID" not in col_index:
        raise ValueError(f"Tabular file must contain 'ReactionID' and {name} columns.")

    id_idx = col_index["ReactionID"]
    idx_gpr = col_index[name]

    return id_idx, idx_gpr


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
    """Simple boolean CLI argument parser accepting 'true' or 'false' (case-insensitive)."""
    def __init__(self, argName :str) -> None:
        self.argName = argName

    def __call__(self, s :str) -> bool: return self.check(s)

    def check(self, s :str) -> bool:
        s = s.lower()
        if s == "true" : return True
        if s == "false": return False
        raise ArgsErr(self.argName, "boolean string (true or false, not case sensitive)", f"\"{s}\"")

class Float:
    """Float CLI argument parser supporting NaN and None keywords (case-insensitive)."""
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
    Custom  = "Custom" 

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
        """Open the SVG metabolic map for this model."""
        path = customPath if self is Model.Custom else FilePath(f"{self.name}_map", FileFormat.SVG, prefix = f"{toolDir}/local/svg metabolic maps/")
        self.__raiseMissingPathErr(path)
        return readSvg(path, customErr = DataErr(path, f"custom map in wrong format"))
    
    def getCOBRAmodel(self, toolDir = ".", customPath :Optional[FilePath] = None, customExtension :Optional[FilePath]=None)->cobra.Model:
        """Load the COBRA model for this enum variant (supports Custom with explicit path/extension)."""
        if(self is Model.Custom):
            return self.load_custom_model(customPath, customExtension)
        else:
            return cobra.io.read_sbml_model(FilePath(f"{self.name}", FileFormat.XML, prefix = f"{toolDir}/local/models/").show())
        
    def load_custom_model(self, file_path :FilePath, ext :Optional[FileFormat] = None) -> cobra.Model:
        """Load a COBRA model from a custom path, supporting XML, JSON, MAT, and YML (compressed or not)."""
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


