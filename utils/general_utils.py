import pickle
import re
import sys
from enum import Enum
from itertools import count
from argparse import Namespace
from typing import Any, Callable, Generic, Optional, Tuple, TypeVar, Union

# RESULT
T = TypeVar('T')
E = TypeVar('E', bound = Exception)
class Result(Generic[T, E]):
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
            ValueError: If the operation was not successful.

        Returns:
            T: The value of the Result instance, if the operation was successful.
        """
        if self.isOk: return self.value
        raise ValueError(f"Unwrapped Result.Err : {self.value}")

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
    
    def expect(self, err :Exception) -> T:
        """
        Expects that the value of the Result instance is successful, otherwise it raises an error.

        Args:
            err (Exception): The error to be raised if the operation was not successful.

        Raises:
            err: The error raised if the operation was not successful

        Returns:
            T: The value of the Result instance, if the operation was successful.
        """
        if self.isOk: return self.value
        raise err

    U = TypeVar("U")
    def map(self, mapper: Callable[[T], U]) -> "Result[U, E]":
        """
        Pipes the Result value into the mapper operation if successful and it returns
        the Result of that operation.

        Args:
            mapper (Callable[[T], U]): The mapper operation to be applied to the Result value.

        Returns:
            Result[U, E]: The result of the mapper operation applied to the Result value.
        """
        if self.isErr: return self
        try: return Result.Ok(mapper(self.value))
        except Exception as e: return Result.Err(e)

# FILES
class FileFormat(Enum):
    """
    Encodes possible file extensions to conditionally save data in a different format.
    """
    DAT    = "dat" # this is how galaxy treats all your files!
    CSV    = "csv"
    XML    = "xml"
    JSON   = "json"
    PICKLE = "p"#"pickle"

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
        if variantName not in FileFormat.__members__:
            raise PathErr(ext, f"\"{ext}\" is not a valid file extension")
        
        return FileFormat[variantName]

    def __str__(self) -> str:
        """
        (Private) converts to str representation. Good practice for usage with argparse.

        Returns:
            str : the string representation of the file extension.
        """
        return self.value

class FilePath():
    """
    Represents a file path. View this as an attempt to standardize file-related operations by expecting values of this type
    in any process requesting a file path.
    """
    def __init__(self, filePath :str, ext :FileFormat, *, prefix = "") -> None:
        """
        (Private) Initializes an instance of FilePath.

        Args:
            path : the end of the path, containing the file name.
            ext : the file's extension.
            prefix : anything before path, already parsed (ending in '/')
        
        Returns:
            None : practically, a FilePath instance.
        """
        self.ext      = ext
        self.prefix   = prefix
        self.filePath = filePath
    
    @classmethod
    def fromStrPath(cls, path :str) -> "FilePath":
        result = re.search(r"^(?:(?P<prefix>.*)\/)?(?P<name>.*)\.(?P<ext>[^.]*)$", path)
        if not result or not result["name"] or not result["ext"]:
            raise PathErr(path, f"cannot recognize folder structure or extension in path \"{path}\"")

        prefix = result["prefix"] if result["prefix"] else ""
        return cls(result["name"], FileFormat.fromExt(result["ext"]), prefix = prefix)
    
    def show(self) -> str:
        """
        Shows the path as a string.

        Returns:
            str : the path.
        """
        return str(self)
    
    def __str__(self) -> str: return f"{self.prefix}/{self.filePath}.{self.ext}"

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

def logWarning(msg :str, loggerPath :FilePath) -> None:
    """
    Log a warning message to an output log file and print it to the console. The final period and a
    newline is added by the function.

    Args:
        s (str): The warning message to be logged and printed.
        loggerPath : The file path of the output log file.

    Returns:
        None
    """
    with open(loggerPath.show(), 'a') as log: log.write(f"{msg}.\n")

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

    def throw(self, loggerPath :Optional[FilePath] = None) -> None:
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

class DataErr(CustomErr):
    """
    utils.CustomErr subclass for data formatting errors.
    """
    errName = "Data Format Error"
    def __init__(self, fileName :str, msg = "no further details provided") -> None:
        super().__init__(f"file \"{fileName}\" contains malformed data", msg)

class PathErr(CustomErr):
    """
    utils.CustomErr subclass for path, filename and extension formatting errors.
    """
    errName = "Path Error"
    def __init__(self, path :str, msg = "no further details provided") -> None:
        super().__init__(f"path \"{path}\" is invalid", msg)

# PICKLE FILES
def readPickle(path :FilePath) -> Any:
    with open(path.show(), "rb") as fd: return pickle.load(fd)

def writePickle(path :FilePath, data :Any) -> None:
    with open(path.show(), "wb") as fd: pickle.dump(data, fd)