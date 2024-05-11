import utils.general_utils as utils
from typing import List, Union, Optional

class RuleErr(utils.CustomErr):
    """
    CustomErr subclass for rule syntax errors.
    """
    errName = "Rule Syntax Error"
    def __init__(self, rule :str, msg = "no further details provided") -> None:
        super().__init__(
            f"rule \"{rule}\" is malformed, {msg}",
            "please verify your input follows the validity guidelines")

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
        Moves the last "actual" item from the 2nd to last list to the beginning of the top list, as per
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

    def obtain(self, err :Optional[utils.CustomErr] = None) -> Optional[OpList]:
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
        (Private) Updates the current OpList to the one on top of the stack.

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