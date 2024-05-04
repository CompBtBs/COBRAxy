import math
import re
import sys
import csv
import pickle
import lxml.etree as ET

from enum import Enum
from itertools import count
from typing import Any, Callable, Dict, Generic, List, Optional, TypeVar, Union

import pandas as pd

# FILES
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

    XML    = ("xml",) # this is one main way cobra models appear in
    JSON   = ("json",) # this is the other
    
    PICKLE = ("pickle", "pk", "p") # this is how all runtime data structures are saved
    #TODO: we're in a pickle (ba dum tss), there's no point in supporting many extensions internally. The
    # issue will never be solved for user-uploaded files and those are saved as .dat by galaxy anyway so it
    # doesn't matter as long as we CAN recognize these 3 names as valid pickle extensions. We must however
    # agree on an internal standard and use only that one, otherwise constructing usable paths becomes a nightmare.
    @classmethod
    def fromExt(cls, ext :str) -> "FileFormat":
        """
        Converts a file extension string to a FileFormat instance.

        Args:
            ext : The file extension as a string.

        Returns:
            FileFormat: The FileFormat instance corresponding to the file extension.
        """
        variantName = ext.upper()
        if variantName in FileFormat.__members__: return FileFormat[variantName]
        
        variantName = variantName.lower()
        for member in cls:
            if variantName in member.value: return member
        
        raise ValueErr("ext", "a valid FileFormat file extension", ext)

    def __str__(self) -> str:
        """
        (Private) converts to str representation. Good practice for usage with argparse.

        Returns:
            str : the string representation of the file extension.
        """
        return self.value[-1] #TODO: fix, it's the dumb pickle thing

class FilePath():
    """
    Represents a file path. View this as an attempt to standardize file-related operations by expecting
    values of this type in any process requesting a file path.
    """
    def __init__(self, filePath :str, ext :FileFormat, *, prefix = "") -> None:
        """
        (Private) Initializes an instance of FilePath.

        Args:
            path : the end of the path, containing the file name.
            ext : the file's extension.
            prefix : anything before path, if the last '/' isn't there it's added by the code.
        
        Returns:
            None : practically, a FilePath instance.
        """
        self.ext      = ext
        self.filePath = filePath

        if prefix and prefix[-1] != '/': prefix += '/'
        self.prefix = prefix
    
    @classmethod
    def fromStrPath(cls, path :str) -> "FilePath":
        """
        Factory method to parse a string from which to obtain, if possible, a valid FilePath instance.

        Args:
            path : the string containing the path
        
        Raises:
            PathErr : if the provided string doesn't represent a valid path.
        
        Returns:
            FilePath : the constructed instance.
        """
        # This method is often used to construct FilePath instances from ARGS UI arguments. These arguments *should*
        # always be correct paths and could be used as raw strings, however most if not all functions that work with
        # file paths request the FilePath objects specifically, which is a very good thing in any case other than this.
        # What ends up happening is we spend time parsing a string into a FilePath so that the function accepts it, only
        # to call show() immediately to bring back the string and open the file it points to.
        # TODO: this is an indication that the arguments SHOULD BE OF TYPE FilePath if they are filepaths, this ENSURES
        # their correctness when modifying the UI and avoids the pointless back-and-forth.
        result = re.search(r"^(?P<prefix>.*\/)?(?P<name>.*)\.(?P<ext>[^.]*)$", path)
        if not result or not result["name"] or not result["ext"]:
            raise PathErr(path, "cannot recognize folder structure or extension in path")

        prefix = result["prefix"] if result["prefix"] else ""
        return cls(result["name"], FileFormat.fromExt(result["ext"]), prefix = prefix)

    def show(self) -> str:
        """
        Shows the path as a string.

        Returns:
            str : the path shown as a string.
        """
        return f"{self.prefix}{self.filePath}.{self.ext}"
    
    def __str__(self) -> str: return self.show()

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
    def __init__(self, argName :str) -> None:
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
    HMRcore = "HMRcore"
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

    def __str__(self) -> str: return self.value