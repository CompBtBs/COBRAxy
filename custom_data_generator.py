import os
import csv
import sys
import cobra
import pickle
import argparse
from enum import Enum, auto
from itertools import count
from typing import Tuple, Union, Optional, List, Dict, TypeVar, Generic, Callable

class FileFormat(Enum):
    """
    Encodes possible file extensions to conditionally save data in a different format.
    """
    CSV    = "csv"
    PICKLE = "pickle"

    def __str__(self) -> str:
        """
        (Private) converts to str representation. Good practice for usage with argparse.
        """
        return self.value

ARGS : argparse.Namespace
def process_args() -> argparse.Namespace:
    """
    Interfaces the script of a module with its frontend, making the user's choices for
    various parameters available as values in code.

    Args:
        args : Always obtained (in file) from sys.argv

    Returns:
        Namespace : An object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        usage = "%(prog)s [options]",
        description = "generate custom data from a given model")
    
    parser.add_argument("-ol", "--out_log",    type = str, required = True, help = "Output log")
    parser.add_argument("-id", "--input",      type = str, required = True, help = "Input model")
    parser.add_argument("-mn", "--name",       type = str, required = True, help = "Input model name")
    # ^ I need this because galaxy converts my files into .dat but I need to know what extension they were in

    parser.add_argument(
        "-of", "--output_format",
        type = FileFormat, default = FileFormat.PICKLE, choices = list(FileFormat),
        required = True,
        help = "Extension of all output files")
    
    argsNamespace = parser.parse_args()
    argsNamespace.out_dir = "result"
    # ^ can't get this one to work from xml, there doesn't seem to be a way to get the directory attribute from the collection

    return argsNamespace

################################- ERROR HANDLING -################################
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

    def throw(self) -> None:
        """
        Raises the current CustomErr instance, logging a warning message before doing so.

        Raises:
            self: The current CustomErr instance.
        
        Returns:
            None
        """
        logWarning(str(self))
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
        Returns a string representing the current CustomErr instance.

        Returns:
            str: A string representing the current CustomErr instance.
        """
        return f"{CustomErr.errName} #{self.id}: {self.msg}, {self.details}."

class RuleErr(CustomErr):
    """
    CustomErr subclass for rule syntax errors.
    """
    errName = "Rule Syntax Error"
    def __init__(self, rule :str, msg = "no further details provided") -> None:
        super().__init__(
            f"rule \"{rule}\" is malformed, {msg}",
            "please verify your input follows the validity guidelines")

class DataErr(CustomErr):
    """
    CustomErr subclass for data formatting errors.
    """
    errName = "Data Format Error"
    def __init__(self, fileName :str, msg = "no further details provided") -> None:
        super().__init__(f"file \"{fileName}\" contains malformed data", msg)

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

def terminate(msg :str) -> None:
    """
    Terminate the execution of the script with an error message.
    
    Args:
        msg (str): The error message to be displayed.
    
    Returns:
        None
    """
    sys.exit(f"Execution aborted: {msg}\n")

def logWarning(msg :str) -> None:
    """
    Log a warning message to an output log file and print it to the console.

    Args:
        s (str): The warning message to be logged and printed.

    Returns:
        None
    """
    with open(ARGS.out_log, 'a') as log: log.write(f"{msg}\n")

################################- INPUT DATA LOADING -################################
def load_custom_model(file_path :str, ext = "") -> cobra.Model:
    """
    Loads a custom model from a file, either in JSON or XML format.

    Args:
        file_path : The path to the file containing the custom model.
        ext : explicit file extension. Necessary for standard use in galaxy because of its weird behaviour.

    Raises:
        DataErr : if the file is in an invalid format or cannot be opened for whatever reason.    
    
    Returns:
        Result[cobra.Model, DataErr]: A Result instance containing the custom model if the
        operation was successful, otherwise a DataErr instance.
    """
    try:
        if ext.lower() == "xml"  or file_path.lower().endswith(".xml"):  return cobra.io.read_sbml_model(file_path)
        if ext.lower() == "json" or file_path.lower().endswith(".json"): return cobra.io.load_json_model(file_path)

    except Exception as e: raise DataErr(file_path, e.__str__())
    raise DataErr(file_path,
        f"Formato \"{file_path.split('.')[-1]}\" non riconosciuto, sono supportati solo file JSON e XML")

################################- RULE PARSING -################################
class OpList(List[Union[str, "OpList"]]):
    """
    Represents a parsed rule and each of its nesting levels, including the operator that level uses.
    """
    def __init__(self, op = "") -> None:
        """
        (Private) Initializes an instance of OpList.

        Args:
            op (str): Operator to be assigned to the OpList. Defaults to "".
        
        Returns:
            None : practically, an OpList instance.
        """
        self.op = op

    def setOpIfMissing(self, op :str) -> None:
        """
        Sets the operator of the OpList if it's missing.

        Args:
            op (str): Operator to be assigned to the OpList.
        
        Returns:
            None
        """
        if not self.op: self.op = op

    def __repr__(self, indent = "") -> str:
        """
        (Private) Returns a string representation of the current OpList instance.

        Args:
            indent (str): Indentation level . Defaults to "".

        Returns:
            str: A string representation of the current OpList instance.
        """
        nextIndent = indent + "  "
        return f"<{self.op}>[\n" + ",\n".join([
            f"{nextIndent}{item.__repr__(nextIndent) if isinstance(item, OpList) else item}"
            for item in self ]) + f"\n{indent}]"

class RuleStack:
    """
    FILO stack structure to save the intermediate representation of a Rule during parsing, with the
    current nesting level at the top of the stack.
    """
    def __init__(self) -> None:
        """
        (Private) initializes an instance of RuleStack.

        Returns:
            None : practically, a RuleStack instance.
        """
        self.__stack = [OpList()] # the stack starts out with the result list already allocated
        self.__updateCurrent()

    def pop(self) -> None:
        """
        Removes the OpList on top of the stack, also flattening it once when possible.

        Side Effects:
            self : mut

        Returns:
            None
        """
        oldTop = self.__stack.pop()
        if len(oldTop) == 1 and isinstance(oldTop[0], OpList): self.__stack[-1][-1] = oldTop[0]
        self.__updateCurrent()

    def push(self, operator = "") -> None:
        """
        Adds a new nesting level, in the form of a new OpList on top of the stack.

        Args:
            operator : the operator assigned to the new OpList.

        Side Effects:
            self : mut
        
        Returns:
            None
        """
        newLevel = OpList(operator)
        self.current.append(newLevel)
        self.__stack.append(newLevel)
        self.__updateCurrent()

    def popForward(self) -> None:
        """
        Moves the last "actual" item from the 2nd list to the beginning of the top list, as per
        the example below:
        stack  : [list_a, list_b]
        list_a : [item1, item2, list_b] --> [item1, list_b]
        list_b : [item3, item4]         --> [item2, item3, item4]

        This is essentially a "give back as needed" operation.

        Side Effects:
            self : mut
        
        Returns:
            None
        """
        self.current.insert(0, self.__stack[-2].pop(-2))

    def currentIsAnd(self) -> bool:
        """
        Checks if the current OpList's assigned operator is "and".

        Returns:
            bool : True if the current OpList's assigned operator is "and", False otherwise.
        """
        return self.current.op == "and"

    def obtain(self, err :Optional[CustomErr] = None) -> Optional[OpList]:
        """
        Obtains the first OpList on the stack, only if it's the only element.

        Args:
            err : The error to raise if obtaining the result is not possible.

        Side Effects:
            self : mut    
        
        Raises:
            err: If given, otherwise None is returned.

        Returns:
            Optional[OpList]: The first OpList on the stack, only if it's the only element.
        """

        if len(self.__stack) == 1: return self.__stack.pop()
        if err: raise err
        return None

    def __updateCurrent(self) -> None:
        """
        Updates the current OpList to the one on top of the stack.

        Side Effects:
            self : mut
        
        Returns:
            None
        """
        self.current = self.__stack[-1]

def parseRuleToNestedList(rule :str) -> OpList:
    """
    Parse a single rule from its string representation to an OpList, making all priority explicit
    through nesting levels.

    Args:
        rule : the string representation of a rule to be parsed.
    
    Raises:
        RuleErr : whenever something goes wrong during parsing.
    
    Returns:
        OpList : the parsed rule.
    """
    source = iter(rule
        .replace("(", "( ").replace(")", " )") # Single out parens as words
        .strip()  # remove whitespace at extremities
        .split()) # split by spaces

    stack = RuleStack()
    nestingErr = RuleErr(rule, "mismatch between open and closed parentheses")
    try:
        while True: # keep reading until source ends
            while True:
                operand = next(source, None) # expected name or rule opening
                if operand is None: raise RuleErr(rule, "found trailing open parentheses")
                if operand == "and" or operand == "or" or operand == ")": # found operator instead, panic
                    raise RuleErr(rule, f"found \"{operand}\" in unexpected position")

                if operand != "(": break # found name

                # found rule opening, we add new nesting level but don't know the operator
                stack.push()

            stack.current.append(operand)

            while True: # keep reading until operator is found or source ends
                operator = next(source, None) # expected operator or rule closing
                if operator and operator != ")": break # found operator

                if stack.currentIsAnd(): stack.pop() # we close the "and" chain

                if not operator: break
                stack.pop() # we close the parentheses

            # we proceed with operator:
            if not operator: break # there is no such thing as a double loop break.. yet
            if operator == "or" and stack.currentIsAnd():
                stack.pop()

            elif operator == "and" and not stack.currentIsAnd():
                stack.push(operator)
                stack.popForward()

            stack.current.setOpIfMissing(operator) # buffer now knows what operator its data had

    except: raise nestingErr

    parsedRule = stack.obtain(nestingErr)
    return parsedRule[0] if len(parsedRule) == 1 and isinstance(parsedRule[0], list) else parsedRule

################################- DATA GENERATION -################################
ReactionId = str
def generate_rules(model: cobra.Model, asParsed = True) -> Union[Dict[ReactionId, OpList], Dict[ReactionId, str]]:
    """
    Generates a dictionary mapping reaction ids to rules from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the rules to an optimized runtime format, otherwise leaves them as strings.

    Returns:
        Dict[ReactionId, OpList] : the generated dictionary of parsed rules.
        Dict[ReactionId, str] : the generated dictionary of raw rules.
    """
    # Is the below approach convoluted? yes
    # Ok but is it inefficient? probably
    # Ok but at least I don't have to repeat the check at every rule (I'm clinically insane)
    _ruleGetter   = lambda reaction : reaction.gene_reaction_rule
    ruleExtractor = (lambda reaction : parseRuleToNestedList(_ruleGetter(reaction))) if asParsed else _ruleGetter

    return {
        reaction.id : ruleExtractor(reaction)
        for reaction in model.reactions
        if reaction.gene_reaction_rule }

def generate_reactions(model :cobra.Model) -> Dict[ReactionId, str]:
    """
    Generates a dictionary mapping reaction ids to reaction formulas from the model.

    Args:
        model : the model to derive data from.

    Returns:
        Dict[ReactionId, str] : the generated dictionary.
    """
    return {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction }

###############################- FILE SAVING -################################
class FilePath():
    """
    Represents a file path. View this as an attempt to standardize file-related operations by expecting values of this type
    in any process requesting a file path.
    """
    def __init__(self, path :str, ext :FileFormat, prefix = "") -> None:
        """
        (Private) Initializes an instance of FilePath.

        Args:
            path : the end of the path, containing the file name.
            ext : the file's extension.
            prefix : anything before path, already parsed (ending in '/')
        
        Returns:
            None : practically, a FilePath instance.
        """
        self.path = f"{prefix}{path}.{ext}"
    
    @classmethod
    def Absolute(cls, path :str, ext :FileFormat) -> "FilePath":
        """
        Factory method for an absolute FilePath instance, where the given path is already complete.

        Args:
            path : the file path.
            ext : the file's extension.
        
        Returns:
            FilePath : the created instance.
        """
        return cls(path, ext)

    @classmethod
    def Output(cls, path :str, ext :FileFormat) -> "FilePath":
        """
        Factory method for an output file's FilePath instance, where the given path is prefixed with the designed output
        folder path in this tool contained in ARGS.out_dir.

        Args:
            path : the file path.
            ext : the file's extension.
        
        Returns:
            FilePath : the created instance.
        """
        return cls(path, ext, prefix = f"{ARGS.out_dir}/")
    
    def show(self) -> str:
        """
        Shows the path as a string.

        Returns:
            str : the path.
        """
        return str(self)
    
    def __str__(self) -> str: return self.path

def save_as_csv(data :dict, file_path :FilePath, fieldNames :Tuple[str, str]) -> None:
    """
    Saves any dictionary-shaped data in a .csv file created at the given file_path.

    Args:
        data : the data to be written to the file.
        file_path : the path to the .csv file.
        fieldNames : the names of the fields (columns) in the .csv file.
    
    Returns:
        None
    """
    with open(file_path.show(), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = fieldNames)
        writer.writeheader()

        for key, value in data.items():
            writer.writerow({ fieldNames[0] : key, fieldNames[1] : value })

def save_as_pickle(data :dict, file_path :FilePath) -> None:
    """
    Saves any dictionary-shaped data in a .pickle file created at the given file_path.

    Args:
        data : the data to be written to the file.
        file_path : the path to the .pickle file
    
    Returns:
        None
    """
    with open(file_path.show(), "wb") as fd: pickle.dump(data, fd)

###############################- ENTRY POINT -################################
def main() -> None:
    """
    Initializes everything and sets the program in motion based on the fronted input arguments.
    
    Returns:
        None
    """
    # get args from frontend (related xml)
    global ARGS
    ARGS = process_args()

    # this is the worst thing I've seen so far, congrats to the former MaREA devs for suggesting this!
    if os.path.isdir(ARGS.out_dir) == False: os.makedirs(ARGS.out_dir)

    # load custom model
    model = load_custom_model(ARGS.input, ARGS.name.split('.')[-1])
    
    # generate data and save it in the desired format and in a location galaxy understands
    # (it should show up as a collection in the history)
    reactions     = generate_reactions(model)
    rulesPath     = FilePath.Output("rules",     ARGS.output_format)
    reactionsPath = FilePath.Output("reactions", ARGS.output_format)

    if ARGS.output_format is FileFormat.PICKLE:
        rules = generate_rules(model, asParsed = True)
        save_as_pickle(rules,     rulesPath)
        save_as_pickle(reactions, reactionsPath)
    
    elif ARGS.output_format is FileFormat.CSV:
        rules = generate_rules(model, asParsed = False)
        save_as_csv(rules,     rulesPath,     ("ReactionID", "Rule"))
        save_as_csv(reactions, reactionsPath, ("ReactionID", "Reaction"))

    # ^ Please if anyone works on this after updating python to 3.12 change the if/elif into a match statement!!

if __name__ == '__main__':
    main()