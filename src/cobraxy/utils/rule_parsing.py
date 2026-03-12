"""
Parsing utilities for gene rules (GPRs).

This module provides:
- RuleErr: structured errors for malformed rules
- RuleOp: valid logical operators (AND/OR)
- OpList: nested list structure representing parsed rules with explicit operator
- RuleStack: helper stack to build nested OpLists during parsing
- parseRuleToNestedList: main entry to parse a rule string into an OpList
"""
from enum import Enum
from typing import List, Union, Optional

try:
    from . import general_utils as utils
except:
    import general_utils as utils

class RuleErr(utils.CustomErr):
    """
    Error type for rule syntax errors.
    """
    errName = "Rule Syntax Error"
    def __init__(self, rule :str, msg = "no further details provided") -> None:
        super().__init__(
            f"rule \"{rule}\" is malformed, {msg}",
            "please verify your input follows the validity guidelines")

class RuleOp(Enum):
    """
    Valid logical operators for gene rules.
    """
    OR  = "or"
    AND = "and"

    @classmethod
    def isOperator(cls, op :str) -> bool:
        return op.upper() in cls.__members__

    def __str__(self) -> str: return self.value

class OpList(List[Union[str, "OpList"]]):
    """
    Parsed rule structure: a list with an associated operator for that level.
    """
    def __init__(self, op :Optional[RuleOp] = None) -> None:
        """
        (Private) Initializes an instance of OpList.

        Args:
            op (str): Operator to be assigned to the OpList. Defaults to "".
        
        Returns:
            None : practically, an OpList instance.
        """
        self.op = op

    def setOpIfMissing(self, op :RuleOp) -> None:
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
    FILO stack used during parsing to build nested OpLists; the top is the current level.
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
        return self.current.op is RuleOp.AND

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
    Parse a rule string into an OpList, making operator precedence explicit via nesting.

    Args:
        rule: Rule string to parse (supports parentheses, 'and', 'or').
    
    Raises:
        RuleErr: If the rule is malformed (e.g., mismatched parentheses or misplaced operators).
    
    Returns:
        OpList: Parsed rule as an OpList structure.
    """
    source = iter(rule
        .replace("(", "( ").replace(")", " )") # single out parentheses as words
        .strip()  # trim edges
        .split()) # split by spaces

    stack = RuleStack()
    nestingErr = RuleErr(rule, "mismatch between open and closed parentheses")
    try:
        while True: # read until source ends
            while True:
                operand = next(source, None) # expect operand or '('
                if operand is None: raise RuleErr(rule, "found trailing open parentheses")
                if operand in ("and", "or", ")"): # unexpected operator position
                    raise RuleErr(rule, f"found \"{operand}\" in unexpected position")

                if operand != "(": break # got a name

                # found rule opening: add a new nesting level
                stack.push()

            stack.current.append(operand)

            while True: # read until operator found or source ends
                operator = next(source, None) # expect operator or ')'
                if operator and operator != ")": break # got operator

                if stack.currentIsAnd(): stack.pop() # close current AND chain

                if not operator: break
                stack.pop() # close parentheses

            if not operator: break
            
            if not RuleOp.isOperator(operator): raise RuleErr(
                rule, f"found \"{operator}\" in unexpected position, expected operator")
            
            operator = RuleOp(operator)
            if operator is RuleOp.OR and stack.currentIsAnd():
                stack.pop()

            elif operator is RuleOp.AND and not stack.currentIsAnd():
                stack.push(operator)
                stack.popForward()

            stack.current.setOpIfMissing(operator)

    except RuleErr as err: raise err # bubble up proper errors
    except: raise nestingErr # everything else is interpreted as a nesting error.

    parsedRule = stack.obtain(nestingErr)
    return parsedRule[0] if len(parsedRule) == 1 and isinstance(parsedRule[0], list) else parsedRule