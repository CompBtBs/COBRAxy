"""
Unit tests for the utils modules (general_utils, rule_parsing, reaction_parsing, model_utils, CBS_backend).

Run with: python -m pytest test_utils.py -v
Or: python test_utils.py
"""

import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

# Try to import pytest, but don't fail if not available
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    # Create a dummy pytest.raises for compatibility
    class _DummyPytest:
        class raises:
            def __init__(self, *args, **kwargs):
                self.expected_exceptions = args
            def __enter__(self):
                return self
            def __exit__(self, exc_type, exc_val, exc_tb):
                if exc_type is None:
                    raise AssertionError("Expected an exception but none was raised")
                if not any(issubclass(exc_type, e) for e in self.expected_exceptions):
                    return False
                return True
    pytest = _DummyPytest()

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import utils.general_utils as utils
import utils.rule_parsing as ruleUtils
import utils.reaction_parsing as reactionUtils
import utils.model_utils as modelUtils

# Get the tool directory (one level up from test directory)
TOOL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

class TestGeneralUtils:
    """Tests for utils.general_utils module"""
    
    def test_bool_check_true(self):
        """Test Bool validator with true values"""
        bool_checker = utils.Bool("testArg")
        assert bool_checker.check("true") == True
        assert bool_checker.check("True") == True
        assert bool_checker.check("TRUE") == True
    
    def test_bool_check_false(self):
        """Test Bool validator with false values"""
        bool_checker = utils.Bool("testArg")
        assert bool_checker.check("false") == False
        assert bool_checker.check("False") == False
        assert bool_checker.check("FALSE") == False
    
    def test_bool_check_invalid(self):
        """Test Bool validator with invalid values"""
        bool_checker = utils.Bool("testArg")
        with pytest.raises((ValueError, utils.ArgsErr)):
            bool_checker.check("foo")
        with pytest.raises((ValueError, utils.ArgsErr)):
            bool_checker.check("1")
    
    def test_custom_error(self):
        """Test CustomErr class"""
        err = utils.CustomErr("Test message", "Test details")
        assert err.msg == "Test message"
        assert err.details == "Test details"
        assert isinstance(err.id, int)
    
    def test_custom_error_with_id(self):
        """Test CustomErr class with custom ID"""
        err = utils.CustomErr("Test message", "Test details", explicitErrCode=42)
        assert err.msg == "Test message"
        assert err.details == "Test details"
        assert err.id == 42
    
    def test_filepath_creation(self):
        """Test FilePath creation"""
        fp = utils.FilePath("test", utils.FileFormat.CSV)
        assert "test" in fp.show()
        assert ".csv" in fp.show() or ".tsv" in fp.show()
    
    def test_filepath_with_prefix(self):
        """Test FilePath with prefix"""
        fp = utils.FilePath("test", utils.FileFormat.CSV, prefix="/tmp")
        path = fp.show()
        assert "/tmp" in path
        assert "test" in path
    
    def test_model_enum(self):
        """Test Model enum"""
        assert utils.Model.ENGRO2 is not None
        assert utils.Model.Recon is not None
        assert utils.Model.Custom is not None


class TestRuleParsing:
    """Tests for utils.rule_parsing module"""
    
    def test_parse_single_gene(self):
        """Test parsing a single gene"""
        result = ruleUtils.parseRuleToNestedList("GENE1")
        assert "GENE1" in str(result)
    
    def test_parse_or_rule(self):
        """Test parsing OR rule"""
        result = ruleUtils.parseRuleToNestedList("A or B")
        assert result.op == ruleUtils.RuleOp.OR
        assert len(result) == 2  # OpList is a list itself
    
    def test_parse_and_rule(self):
        """Test parsing AND rule"""
        result = ruleUtils.parseRuleToNestedList("A and B")
        assert result.op == ruleUtils.RuleOp.AND
        assert len(result) == 2  # OpList is a list itself
    
    def test_parse_complex_rule(self):
        """Test parsing complex nested rule"""
        result = ruleUtils.parseRuleToNestedList("A or (B and C)")
        assert result.op == ruleUtils.RuleOp.OR
    
    def test_parse_invalid_rule(self):
        """Test parsing invalid rule"""
        with pytest.raises(Exception):
            ruleUtils.parseRuleToNestedList("A foo B")
    
    def test_parse_mismatched_parentheses(self):
        """Test parsing rule with mismatched parentheses"""
        with pytest.raises(Exception):
            ruleUtils.parseRuleToNestedList("A)")
    
    def test_rule_op_enum(self):
        """Test RuleOp enum"""
        assert ruleUtils.RuleOp("or") == ruleUtils.RuleOp.OR
        assert ruleUtils.RuleOp("and") == ruleUtils.RuleOp.AND
    
    def test_rule_op_is_operator(self):
        """Test RuleOp.isOperator"""
        assert ruleUtils.RuleOp.isOperator("or") == True
        assert ruleUtils.RuleOp.isOperator("and") == True
        assert ruleUtils.RuleOp.isOperator("foo") == False


class TestReactionParsing:
    """Tests for utils.reaction_parsing module"""
    
    def test_reaction_dir_reversible(self):
        """Test ReactionDir detection for reversible reactions"""
        result = reactionUtils.ReactionDir.fromReaction("atp <=> adp + pi")
        assert result == reactionUtils.ReactionDir.REVERSIBLE
    
    def test_reaction_dir_forward(self):
        """Test ReactionDir detection for forward reactions"""
        result = reactionUtils.ReactionDir.fromReaction("atp --> adp + pi")
        assert result == reactionUtils.ReactionDir.FORWARD
    
    def test_reaction_dir_backward(self):
        """Test ReactionDir detection for backward reactions"""
        result = reactionUtils.ReactionDir.fromReaction("atp <-- adp + pi")
        assert result == reactionUtils.ReactionDir.BACKWARD
    
    def test_reaction_dir_invalid(self):
        """Test ReactionDir with invalid arrow"""
        with pytest.raises(Exception):
            reactionUtils.ReactionDir.fromReaction("atp ??? adp + pi")
    
    def test_create_reaction_dict(self):
        """Test creating reaction dictionary"""
        reactions = {
            'r1': '2 pyruvate + 1 h2o <=> 1 h2o + 2 acetate',
            'r2': '2 co2 + 6 h2o --> 3 atp'
        }
        result = reactionUtils.create_reaction_dict(reactions)
        
        # Check that we have the expected reactions
        assert 'r1_B' in result or 'r1_F' in result
        assert 'r2' in result


class TestModelUtils:
    """Tests for utils.model_utils module"""
    
    def test_gene_type_detection(self):
        """Test gene type detection"""
        # Test with entrez ID (numeric)
        assert modelUtils.gene_type("123456", "test") == "entrez_id"
        
        # Test with Ensembl ID
        assert modelUtils.gene_type("ENSG00000123456", "test") == "ENSG"
        
        # Test with symbol
        assert modelUtils.gene_type("TP53", "test") == "HGNC_symbol"


class TestModelLoading:
    """Tests for model loading functionality"""
    
    def test_engro2_model_exists(self):
        """Test that ENGRO2 model files exist"""
        model_path = os.path.join(TOOL_DIR, "local", "models", "ENGRO2.xml")
        assert os.path.exists(model_path), f"ENGRO2 model not found at {model_path}"
    
    def test_recon_model_exists(self):
        """Test that Recon model files exist"""
        model_path = os.path.join(TOOL_DIR, "local", "models", "Recon.xml")
        assert os.path.exists(model_path), f"Recon model not found at {model_path}"
    
    def test_pickle_files_exist(self):
        """Test that pickle files exist"""
        pickle_dir = os.path.join(TOOL_DIR, "local", "pickle files")
        assert os.path.exists(pickle_dir), f"Pickle directory not found at {pickle_dir}"
        
        # Check for some expected pickle files
        expected_files = ["synonyms.pickle", "black_list.pickle"]
        for fname in expected_files:
            fpath = os.path.join(pickle_dir, fname)
            assert os.path.exists(fpath), f"Expected pickle file not found: {fpath}"
    
    def test_map_files_exist(self):
        """Test that SVG map files exist"""
        map_dir = os.path.join(TOOL_DIR, "local", "svg metabolic maps")
        assert os.path.exists(map_dir), f"Map directory not found at {map_dir}"
    
    def test_medium_file_exists(self):
        """Test that medium file exists"""
        medium_path = os.path.join(TOOL_DIR, "local", "medium", "medium.csv")
        assert os.path.exists(medium_path), f"Medium file not found at {medium_path}"
    
    def test_mapping_file_exists(self):
        """Test that mapping file exists"""
        mapping_path = os.path.join(TOOL_DIR, "local", "mappings", "genes_human.csv")
        assert os.path.exists(mapping_path), f"Mapping file not found at {mapping_path}"


if __name__ == "__main__":
    # Run tests with pytest if available, otherwise run basic checks
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        print("pytest not available, running basic tests...")
        
        # Run basic tests manually
        test_classes = [
            TestGeneralUtils(),
            TestRuleParsing(),
            TestReactionParsing(),
            TestModelUtils(),
            TestModelLoading()
        ]
        
        failed = 0
        passed = 0
        
        for test_class in test_classes:
            class_name = test_class.__class__.__name__
            print(f"\n{class_name}:")
            
            for method_name in dir(test_class):
                if method_name.startswith("test_"):
                    try:
                        method = getattr(test_class, method_name)
                        method()
                        print(f"  ✓ {method_name}")
                        passed += 1
                    except Exception as e:
                        print(f"  ✗ {method_name}: {str(e)}")
                        failed += 1
        
        print(f"\n{'='*60}")
        print(f"Results: {passed} passed, {failed} failed")
        if failed > 0:
            sys.exit(1)
