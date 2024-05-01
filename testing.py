# This is a general-purpose "testing utilities" module for the MaREA tool.

## Imports:
from typing import Dict, Callable, Type, List
from enum import Enum, auto
from collections.abc import Iterable

# Imports needed by the modules I'm testing:
import pickle as pk
import math

## Generic utilities:
class TestResult:
    """
    Represents the result of a test and contains all the relevant information about it. Loosely models two variants:
    - Ok: The test passed, no further information is saved besides the target's name.
    - Err: The test failed, an error message and further contextual details are also saved.

    This class does not ensure a static proof of the two states' behaviour, their meaning or mutual exclusivity outside
    of the :bool property "isPass", meant for outside reads.
    """
    def __init__(self, isPass :bool, targetName :str, errMsg = "", details = "") -> None:
        """
        (Private) Initializes an instance of TestResult.

        Args:
            isPass : distinction between TestResult.Ok (True) and TestResult.Err (False).
            targetName : the name of the target object / property / function / module being tested, not always set
            to a meaningful value at this stage.
            
            errMsg : concise error message explaining the test's failure.
            details : contextual details about the error.
        
        Returns:
            None : practically, a TestResult instance.
        """
        self.isPass = isPass
        self.isFail = not isPass # Convenience above all

        self.targetName = targetName
        if isPass: return

        self.errMsg   = errMsg
        self.details  = details

    @classmethod
    def Ok(cls, targetName = "") -> "TestResult":
        """
        Factory method for TestResult.Ok, where all we need to know is that our test passed.

        Args:
            targetName : the name of the target object / property / function / module being tested, not always set
            to a meaningful value at this stage.
        
        Returns:
            TestResult : a new Ok instance.
        """
        return cls(True, targetName)

    @classmethod
    def Err(cls, errMsg :str, details :str, targetName = "") -> "TestResult":
        """
        Factory method for TestResult.Err, where we store relevant error information.

        Args:
            errMsg : concise error message explaining the test's failure.
            details : contextual details about the error.
            targetName : the name of the target object / property / function / module being tested, not always set
            to a meaningful value at this stage.
        
        Returns:
            TestResult : a new Err instance.
        """
        return cls(False, targetName, errMsg, details)

    def log(self, isCompact = True) -> str:
        """
        Dumps all the available information in a :str, ready for logging.

        Args:
            isCompact : if True limits the amount of information displayed to the targetName.
        
        Returns:
            str : information about this test result.

        """
        if isCompact:
            return f"{TestResult.__name__}::{'Ok' if self.isPass else 'Err'}(Unit test on {self.targetName})"
        
        logMsg = f"Unit test on {self.targetName} {'passed' if self.isPass else f'failed because {self.errMsg}'}"
        if self.details: logMsg += f", {self.details}"
        return logMsg

    def throw(self) -> None:
        #TODO: finer Exception typing would be desirable
        """
        Logs the result information and panics.

        Raises:
            Exception : an error containing log information about the test result.

        Returns:
            None

        """
        raise Exception(self.log())

class CheckingMode:
    """
    (Private) Represents a way to check a value for correctness, in the context of "testing" it.
    """

    def __init__(self) -> None:
        """
        (Private) Implemented on child classes, initializes an instance of CheckingMode.

        Returns:
            None : practically, a CheckingMode instance.
        """
        self.logMsg = "CheckingMode base class should not be used directly"

    def __checkPasses__(self, _) -> bool:
        """
        (Private) Implemented on child classes, performs the actual correctness check on a received value.

        Returns:
            bool : True if the check passed, False if it failed.
        """
        return True
    
    def check(self, value) -> TestResult:
        """
        Converts the :bool evaluation of the value's correctness to a TestResult.

        Args:
            value : the value to check.

        Returns:
            TestResult : the result of the check.
        """
        return TestResult.Ok() if self.__checkPasses__(value) else TestResult.Err(self.logMsg, f"got {value} instead")
    
    def __repr__(self) -> str:
        """
        (Private) Implemented on child classes, formats :object as :str.
        """
        return self.__class__.__name__

class ExactValue(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value needs to match another exactly.
    """
    
    #I suggest solving the more complex equality checking edge cases with the "Satisfies" and "MatchingShape" variants.
    def __init__(self, value) -> None:
        self.value  = value
        self.logMsg = f"value needed to match {value} exactly"
    
    def __checkPasses__(self, value) -> bool:
        return self.value == value

    def __repr__(self) -> str:
        return f"{super().__repr__()}({self.value})"

class AcceptedValues(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value needs to appear in a list of accepted values.
    """
    def __init__(self, *values) -> None:
        self.values = values
        self.logMsg = f"value needed to be one of these: {values}"
    
    def __checkPasses__(self, value) -> bool:
        return value in self.values

    def __repr__(self) -> str:
        return f"{super().__repr__()}{self.values}"

class SatisfiesPredicate(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value needs to verify a given predicate, as in
    the predicate accepts it as input and returns True.
    """
    def __init__(self, pred :Callable[..., bool], predName = "") -> None:
        self.pred = pred
        self.logMsg = f"value needed to verify a predicate{bool(predName) * f' called {predName}'}"
    
    def __checkPasses__(self, *params) -> bool:
        return self.pred(*params)

    def __repr__(self) -> str:
        return f"{super().__repr__()}(T) -> bool"

class IsOfType(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value needs to be of a certain type.
    """
    def __init__(self, type :Type) -> None:
        self.type = type
        self.logMsg = f"value needed to be of type {type.__name__}"
    
    def __checkPasses__(self, value :Type) -> bool:
        return isinstance(value, self.type)

    def __repr__(self) -> str:
        return f"{super().__repr__()}:{self.type.__name__}"

class Exists(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value needs to exist (or not!). Mainly employed as a quick default
    check that always passes, it still upholds its contract when it comes to checking for existing properties in objects
    without much concern on what value they contain.
    """
    def __init__(self, exists = True) -> None:
        self.exists = exists
        self.logMsg = f"value needed to {(not exists) * 'not '}exist"
    
    def __checkPasses__(self, _) -> bool: return self.exists

    def __repr__(self) -> str:
        return f"{super().__repr__() if self.exists else 'IsMissing'}"

class MatchingShape(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value is an object that needs to have a certain shape,
    as in to posess properties with a given name and value. Each property is checked for existance and correctness with
    its own given CheckingMode.
    """
    def __init__(self, props :Dict[str, CheckingMode], objName = "") -> None:
        """
        (Private) Initializes an instance of MatchingShape.

        Args:
            props : :dict using property names as keys and checking modes for the property's value as values.
            objName : label for the object we're testing the shape of.

        Returns:
            None : practically, a MatchingShape instance.
        """
        self.props   = props
        self.objName = objName

        self.shapeRepr = " {\n" + "\n".join([f"  {propName} : {prop}" for propName, prop in props.items()]) + "\n}"
    
    def check(self, obj :object) -> TestResult:
        objIsDict = isinstance(obj, dict) # Python forces us to distinguish between object properties and dict keys
        for propName, checkingMode in self.props.items():
            # Checking if the property exists:
            if (not objIsDict and not hasattr(obj, propName)) or (objIsDict and propName not in obj):
                if not isinstance(checkingMode, Exists): return TestResult.Err(
                    f"property \"{propName}\" doesn't exist on object {self.objName}", "", self.objName)
                
                if not checkingMode.exists: return TestResult.Ok(self.objName)
                # Either the property value is meant to be checked (checkingMode is anything but Exists)
                # or we want the property to not exist, all other cases are handled correctly ahead
        
            checkRes = checkingMode.check(obj[propName] if objIsDict else getattr(obj, propName))
            if checkRes.isPass: continue

            checkRes.targetName = self.objName
            return TestResult.Err(
                f"property \"{propName}\" failed check {checkingMode} on shape {obj}",
                checkRes.log(isCompact = False),
                self.objName)

        return TestResult.Ok(self.objName)

    def __repr__(self) -> str:
        return super().__repr__() + self.shapeRepr

class Many(CheckingMode):
    """
    CheckingMode subclass variant to be used when the checked value is an Iterable we want to check item by item.
    """
    def __init__(self, *values :CheckingMode) -> None:
        self.values = values
        self.shapeRepr = " [\n" + "\n".join([f"  {value}" for value in values]) + "\n]"

    def check(self, coll :Iterable) -> TestResult:
        amt         = len(coll)
        expectedAmt = len(self.values)
        # Length equality is forced:
        if amt != expectedAmt: return TestResult.Err(
            "items' quantities don't match", f"expected {expectedAmt} items, but got {amt}")

        # Items in the given collection value are paired in order with the corresponding checkingMode meant for each of them
        for item, checkingMode in zip(coll, self.values):
            checkRes = checkingMode.check(item)
            if checkRes.isFail: return TestResult.Err(
                f"item in list failed check {checkingMode}",
                checkRes.log(isCompact = False))

        return TestResult.Ok()

    def __repr__(self) -> str:
        return super().__repr__() + self.shapeRepr

class LogMode(Enum):
    """
    Represents the level of detail of a logged message. Models 4 variants, in order of increasing detail:
    - Minimal  : Logs the overall test result for the entire module.
    - Default  : Also logs all single test fails, in compact mode.
    - Detailed : Logs all function test results, in compact mode.
    - Pedantic : Also logs all single test results in detailed mode.
    """
    Minimal  = auto()
    Default  = auto()
    Detailed = auto()
    Pedantic = auto()

    def isMoreVerbose(self, requiredMode :"LogMode") -> bool:
        """
        Compares the instance's level of detail with that of another.

        Args:
            requiredMode : the other instance.
        
        Returns:
            bool : True if the caller instance is a more detailed variant than the other.
        """
        return self.value >= requiredMode.value

## Specific Unit Testing utilities:
class UnitTest:
    """
    Represents a unit test, the test of a single function's isolated correctness.
    """
    def __init__(self, func :Callable, inputParams :list, expectedRes :CheckingMode) -> None:
        """
        (Private) Initializes an instance of UnitTest.

        Args:
            func : the function to test.
            inputParams : list of parameters to pass as inputs to the function, in order.
            expectedRes : checkingMode to test the function's return value for correctness.

        Returns:
            None : practically, a UnitTest instance.
        """
        self.func        = func
        self.inputParams = inputParams
        self.expectedRes = expectedRes

        self.funcName = func.__name__
    
    def test(self) -> TestResult:
        """
        Tests the function.

        Returns:
            TestResult : the test's result.
        """
        result = None
        try: result = self.func(*self.inputParams)
        except Exception as e: return TestResult.Err("the function panicked at runtime", e, self.funcName)

        checkRes = self.expectedRes.check(result)
        checkRes.targetName = self.funcName
        return checkRes

class UnitTester:
    """
    Manager class for unit testing an entire module, groups single UnitTests together and executes them in order on a
    per-function basis (tests about the same function are executed consecutively) giving back as much information as
    possible depending on the selected logMode. More customization options are available.
    """
    def __init__(self, moduleName :str, logMode = LogMode.Default, stopOnFail = True, *funcTests :'UnitTest') -> None:
        """
        (Private) initializes an instance of UnitTester.

        Args:
            moduleName : name of the tested module.
            logMode : level of detail applied to all messages logged during the test.
            stopOnFail : if True, the test stops entirely after one unit test fails.
            funcTests : the unit tests to perform on the module.

        Returns:
            None : practically, a UnitTester instance.
        """
        self.logMode    = logMode
        self.moduleName = moduleName
        self.stopOnFail = stopOnFail

        # This ensures the per-function order:
        self.funcTests :Dict[str, List[UnitTest]]= {}
        for test in funcTests:
            if test.funcName in self.funcTests: self.funcTests[test.funcName].append(test)
            else: self.funcTests[test.funcName] = [test]

    def logTestResult(self, testRes :TestResult) -> None:
        """
        Prints the formatted result information of a unit test.

        Args:
            testRes : the result of the test.
        
        Returns:
            None
        """
        if testRes.isPass: return self.log("Passed!", LogMode.Detailed, indent = 2)
        
        failMsg = "Failed! "
        # Doing it this way prevents .log computations when not needed
        if self.logMode.isMoreVerbose(LogMode.Detailed):
            # Given that Pedantic is the most verbose variant, there's no point in comparing with LogMode.isMoreVerbose
            failMsg += testRes.log(self.logMode is not LogMode.Pedantic)
        
        self.log(failMsg, indent = 2)

    def log(self, msg :str, minRequiredMode = LogMode.Default, indent = 0) -> None:
        """
        Prints and formats a message only when the UnitTester instance is set to a level of detail at least equal
        to a minimum requirement, given as input.

        Args:
            msg : the message to print.
            minRequiredMode : minimum detail requirement.
            indent : formatting information, counter from 0 that adds 2 spaces each number up
        
        Returns:
            None
        """
        if self.logMode.isMoreVerbose(minRequiredMode): print("  " * indent + msg)

    def testFunction(self, name :str) -> TestResult:
        """
        Perform all unit tests relative to the same function, plus the surrounding logs and checks.

        Args:
            name : the name of the tested function.
        
        Returns :
            TestResult : the overall Ok result of all the tests passing or the first Err. This behaviour is unrelated
            to that of the overall testing procedure (stopOnFail), it always works like this for tests about the
            same function.
        """
        self.log(f"Unit testing {name}...", indent = 1)

        allPassed = True
        for unitTest in self.funcTests[name]:
            testRes = unitTest.test()
            self.logTestResult(testRes)
            if testRes.isPass: continue
            
            allPassed = False
            if self.stopOnFail: break
        
        self.log("", LogMode.Detailed) # Provides one extra newline of space when needed, to better format the output
        if allPassed: return TestResult.Ok(name)

        if self.logMode is LogMode.Default: self.log("")        
        return TestResult.Err(f"Unlogged err", "unit test failed", name)
    
    def testModule(self) -> None:
        """
        Runs all the provided unit tests in order but on a per-function basis.

        Returns:
            None
        """
        self.log(f"Unit testing module {self.moduleName}...", LogMode.Minimal)
        
        fails = 0
        testStatusMsg = "complete"
        for funcName in self.funcTests.keys():
            if self.testFunction(funcName).isPass: continue
            fails += 1

            if self.stopOnFail:
                testStatusMsg = "interrupted"
                break

        self.log(f"Testing {testStatusMsg}: {fails} problem{'s' * (fails != 1)} found.\n", LogMode.Minimal)
        # ^^^ Manually applied an extra newline of space.

## Unit testing all the modules:
def unit_marea() -> None:
    import marea as m
    import pandas as pd

    unitTester = UnitTester("marea", LogMode.Pedantic, False,
        UnitTest(m.read_dataset, ["", ""], IsOfType(pd.DataFrame)),

        UnitTest(m.name_dataset, ["customName", 12], ExactValue("customName")),
        UnitTest(m.name_dataset, ["Dataset", 12], ExactValue("Dataset_12")),

        UnitTest(m.fold_change, [0.5, 0.5], ExactValue(0.0)),
        UnitTest(m.fold_change, [0, 0.35], ExactValue("-INF")),
        UnitTest(m.fold_change, [0.5, 0], ExactValue("INF")),
        UnitTest(m.fold_change, [0, 0], ExactValue(0)),

        UnitTest(
            m.Arrow(m.Arrow.MAX_W, m.ArrowColor.DownRegulated, True).toStyleStr, [],
            ExactValue(";stroke:#0000FF;stroke-width:12;stroke-dasharray:5,5")),
    ).testModule()

def unit_rps_generator() -> None:
    import rps_generator as rps
    import pandas as pd
    dataset = pd.DataFrame({
        "cell lines" : ["normal", "cancer"],
        "pyru_vate"  : [5.3, 7.01],
        "glu,cose"   : [8.2, 4.0],
        "unknown"    : [3.0, 3.97],
        "()atp"      : [7.05, 8.83],
    })

    abundancesNormalRaw = {
        "pyru_vate" : 5.3,
        "glu,cose"  : 8.2,
        "unknown"   : 3.0,
        "()atp"     : 7.05,
    }

    abundancesNormal = {
        "pyr"    : 5.3,
        "glc__D" : 8.2,
        "atp"    : 7.05,
    }

    with open('/home/vboxuser/galaxy/tools/marea/local/pickle files/synonyms.pickle', 'rb') as sd:
        synsDict = pk.load(sd)

    reactionsDict = {
        "r1" : {
            "glc__D" : 1
        },
    
        "r2" : {
            "co2" : 2,
            "pyr" : 3,
        },

        "r3" : {
            "atp"    : 2,
            "glc__D" : 4,
        },
        
        "r4" : {
            "atp" : 3,
        }
    }

    abundancesNormalEdited = {
        "pyr"    : 5.3,
        "glc__D" : 8.2,
        "atp"    : 7.05,
        "co2"    : 1,
    }

    blackList = ["atp"] # No jokes allowed!
    missingInDataset = ["co2"]

    normalRpsShape = MatchingShape({
        "r1" : ExactValue(8.2 ** 1),
        "r2" : ExactValue((1 ** 2) * (5.3 ** 3)),
        "r3" : ExactValue((8.2 ** 4) * (7.05 ** 2)),
        "r4" : SatisfiesPredicate(lambda n : math.isnan(n))
    }, "rps dict")

    UnitTester("rps_generator", LogMode.Pedantic, False,
        UnitTest(rps.get_abund_data, [dataset, 0], MatchingShape({
            "pyru_vate" : ExactValue(5.3),
            "glu,cose"  : ExactValue(8.2),
            "unknown"   : ExactValue(3.0),
            "()atp"     : ExactValue(7.05),
            "name"      : ExactValue("normal")
        }, "abundance series")),
        
        UnitTest(rps.get_abund_data, [dataset, 1], MatchingShape({
            "pyru_vate" : ExactValue(7.01),
            "glu,cose"  : ExactValue(4.0),
            "unknown"   : ExactValue(3.97),
            "()atp"     : ExactValue(8.83),
            "name"      : ExactValue("cancer")
        }, "abundance series")),

        UnitTest(rps.get_abund_data, [dataset, -1], ExactValue(None)),

        UnitTest(rps.update_metabolite_names, [abundancesNormalRaw, synsDict], MatchingShape({
            "pyr"     : ExactValue(5.3),
            "glc__D"  : ExactValue(8.2),
            "atp"     : ExactValue(7.05),
            "unknown" : Exists(False)
        }, "abundance dict")),

        UnitTest(rps.check_missing_metab, [reactionsDict, abundancesNormal.copy()], Many(MatchingShape({
            "pyr"    : ExactValue(5.3),
            "glc__D" : ExactValue(8.2),
            "atp"    : ExactValue(7.05),
            "co2"    : ExactValue(1)
        }, "updated abundances"), Many(ExactValue("co2")))),

        UnitTest(rps.clean_metabolite_name, ["4,4'-diphenylmethane diisocyanate"], ExactValue("44diphenylmethanediisocyanate")),

        UnitTest(rps.get_metabolite_id, ["tryptophan", synsDict], ExactValue("trp__L")),

        UnitTest(rps.ReactionDir.fromReaction, ["atp <=> adp + pi"], ExactValue(rps.ReactionDir.REVERSIBLE)),
        UnitTest(rps.ReactionDir.fromReaction, ["atp --> adp + pi"], ExactValue(rps.ReactionDir.FORWARD)),
        UnitTest(rps.ReactionDir.fromReaction, ["atp <-- adp + pi"], ExactValue(rps.ReactionDir.BACKWARD)),
        UnitTest(rps.ReactionDir.fromReaction, ["atp ??? adp + pi"], Exists(False)), # should panic

        UnitTest(rps.calculate_rps, [reactionsDict, abundancesNormalEdited, blackList, missingInDataset], normalRpsShape),

        UnitTest(rps.rps_for_cell_lines, [dataset, reactionsDict, blackList, synsDict, "", True], Many(normalRpsShape, MatchingShape({
            "r1" : ExactValue(4.0 ** 1),
            "r2" : ExactValue((1 ** 2) * (7.01 ** 3)),
            "r3" : ExactValue((4.0 ** 4) * (8.83 ** 2)),
            "r4" : SatisfiesPredicate(lambda n : math.isnan(n))
        }, "rps dict"))),

        #UnitTest(rps.main, [], ExactValue(None)) # Complains about sys argvs
    ).testModule()

def unit_custom_data_generator() -> None:
    import custom_data_generator as cdg

    UnitTester("custom data generator", LogMode.Pedantic, False,
        UnitTest(lambda :True, [], ExactValue(True)), # No tests can be done without a model at hand!
    ).testModule()

def unit_utils() -> None:
    import utils.general_utils as utils
    import utils.rule_parsing as ruleUtils
    import utils.reaction_parsing as reactionUtils

    UnitTester("utils", LogMode.Pedantic, False,
        UnitTest(utils.CustomErr, ["myMsg", "more details"], MatchingShape({
            "details" : ExactValue("more details"),
            "msg"     : ExactValue("myMsg"),
            "id"      : ExactValue(0) # this will fail if any custom errors happen anywhere else before!
        })),

        UnitTest(utils.CustomErr, ["myMsg", "more details", 42], MatchingShape({
            "details" : ExactValue("more details"),
            "msg"     : ExactValue("myMsg"),
            "id"      : ExactValue(42)
        })),

        UnitTest(utils.Bool("someArg").check, ["TrUe"],  ExactValue(True)),
        UnitTest(utils.Bool("someArg").check, ["FALse"], ExactValue(False)),
        UnitTest(utils.Bool("someArg").check, ["foo"],   Exists(False)), # should panic!

        UnitTest(utils.Model.ENGRO2.getRules, ["."], IsOfType(dict)),
        UnitTest(utils.Model.Custom.getRules, [".", ""], Exists(False)), # expected panic

        # rule utilities tests:
        UnitTest(
            ruleUtils.parseRuleToNestedList, ["A or B and C or (D and E)"],
            Many(
                ExactValue("A"), 
                Many(ExactValue("B"), ExactValue("C")),
                Many(ExactValue("D"), ExactValue("E"))
        )),

        UnitTest(
            reactionUtils.create_reaction_dict,
            [{'shdgd': '2 pyruvate + 1 h2o <=> 1 h2o + 2 acetate', 'sgwrw': '2 co2 + 6 h2o --> 3 atp'}], 
            MatchingShape({
                "shdgd_B" : MatchingShape({
                    "acetate" : ExactValue(2),
                    "h2o" : ExactValue(1),
                }),

                "shdgd_F" : MatchingShape({
                    "pyruvate" : ExactValue(2),
                    "h2o" : ExactValue(1)
                }),

                "sgwrw" : MatchingShape({
                    "co2" : ExactValue(2),
                    "h2o" : ExactValue(6),
                })
            }, "reaction dict")),   
    ).testModule()

def unit_ras_generator() -> None:
    import ras_generator as ras
    import utils.rule_parsing as ruleUtils

    # Making an alias to mask the name of the inner function and separate the 2 tests:
    def opListAlias(op_list, dataset):
        ras.ARGS.none = False
        return ras.ras_op_list(op_list, dataset)
    
    ras.ARGS = ras.process_args()
    rule = ruleUtils.OpList("and")
    rule.extend(["foo", "bar", "baz"])

    dataset = { "foo" : 5, "bar" : 2, "baz" : None }
    
    UnitTester("ras generator", LogMode.Pedantic, False,
        UnitTest(ras.ras_op_list, [rule, dataset], ExactValue(2)),
        UnitTest(opListAlias, [rule, dataset], ExactValue(None)),
    ).testModule()

if __name__ == "__main__":
    unit_marea()
    unit_custom_data_generator()
    unit_utils()
    unit_ras_generator()