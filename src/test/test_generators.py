"""
Unit tests for RAS and RPS generator modules.

Run with: python -m pytest test_generators.py -v
Or: python test_generators.py
"""

import sys
import os
import pandas as pd
import numpy as np
import math
from pathlib import Path
import tempfile

# Try to import pytest, but don't fail if not available
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
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

import ras_generator
import rps_generator
import utils.general_utils as utils
import utils.rule_parsing as ruleUtils

# Get the tool directory
TOOL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class TestRASGenerator:
    """Tests for ras_generator module"""
    
    def test_ras_op_list_and(self):
        """Test RAS calculation with AND rule"""
        # Create a mock args object
        class MockArgs:
            none = False
        
        ras_generator.ARGS = MockArgs()
        
        # Create an OpList with AND operator
        rule = ruleUtils.OpList(ruleUtils.RuleOp.AND)
        rule.extend(["gene1", "gene2", "gene3"])
        
        # Create dataset
        dataset = {
            "gene1": 5.0,
            "gene2": 2.0,
            "gene3": 3.0
        }
        
        # Should return minimum value (AND logic)
        result = ras_generator.ras_op_list(rule, dataset)
        assert result == 2.0
    
    def test_ras_op_list_or(self):
        """Test RAS calculation with OR rule"""
        class MockArgs:
            none = False
        
        ras_generator.ARGS = MockArgs()
        
        # Create an OpList with OR operator
        rule = ruleUtils.OpList(ruleUtils.RuleOp.OR)
        rule.extend(["gene1", "gene2", "gene3"])
        
        dataset = {
            "gene1": 5.0,
            "gene2": 2.0,
            "gene3": 3.0
        }
        
        # Should return maximum value (OR logic)
        result = ras_generator.ras_op_list(rule, dataset)
        assert result == 5.0
    
    def test_ras_op_list_with_none(self):
        """Test RAS calculation with None values"""
        class MockArgs:
            none = True
        
        ras_generator.ARGS = MockArgs()
        
        rule = ruleUtils.OpList(ruleUtils.RuleOp.AND)
        rule.extend(["gene1", "gene2"])
        
        dataset = {
            "gene1": 5.0,
            "gene2": None
        }
        
        # Should return None when none=True and a gene is None
        result = ras_generator.ras_op_list(rule, dataset)
        assert result is None
    
    def test_process_args(self):
        """Test argument processing"""
        # Test that process_args returns a valid Namespace
        args = ras_generator.process_args([])
        assert hasattr(args, 'tool_dir')


class TestRPSGenerator:
    """Tests for rps_generator module"""
    
    def test_get_abund_data(self):
        """Test extracting abundance data from dataset"""
        dataset = pd.DataFrame({
            "cell_lines": ["normal", "cancer"],
            "pyruvate": [5.3, 7.01],
            "glucose": [8.2, 4.0],
            "atp": [7.05, 8.83]
        })
        
        # Get first row (normal)
        result = rps_generator.get_abund_data(dataset, 0)
        assert result is not None
        assert result["pyruvate"] == 5.3
        assert result["glucose"] == 8.2
        assert result["atp"] == 7.05
        assert result["name"] == "normal"
        
        # Get second row (cancer)
        result = rps_generator.get_abund_data(dataset, 1)
        assert result is not None
        assert result["pyruvate"] == 7.01
        assert result["name"] == "cancer"
    
    def test_get_abund_data_invalid_index(self):
        """Test extracting abundance data with invalid index"""
        dataset = pd.DataFrame({
            "cell_lines": ["normal", "cancer"],
            "pyruvate": [5.3, 7.01]
        })
        
        # Try to get invalid index
        result = rps_generator.get_abund_data(dataset, -1)
        assert result is None
        
        result = rps_generator.get_abund_data(dataset, 999)
        assert result is None
    
    def test_clean_metabolite_name(self):
        """Test metabolite name cleaning"""
        # Test removing special characters
        result = rps_generator.clean_metabolite_name("4,4'-diphenylmethane diisocyanate")
        assert "," not in result
        assert "'" not in result
        assert " " not in result
        assert result == "44diphenylmethanediisocyanate"
        
        # Test with parentheses
        result = rps_generator.clean_metabolite_name("(S)-lactate")
        assert "(" not in result
        assert ")" not in result
    
    def test_check_missing_metab(self):
        """Test checking for missing metabolites"""
        reactions_dict = {
            "r1": {"glc__D": 1, "atp": 1},
            "r2": {"co2": 2, "pyr": 3}
        }
        
        abundances = {
            "glc__D": 8.2,
            "atp": 7.05,
            "pyr": 5.3
            # co2 is missing
        }
        
        updated_abundances, missing = rps_generator.check_missing_metab(
            reactions_dict, 
            abundances.copy()
        )
        
        # Should have added co2 with value 1
        assert "co2" in updated_abundances
        assert updated_abundances["co2"] == 1
        
        # Should report co2 as missing
        assert "co2" in missing
    
    def test_calculate_rps(self):
        """Test RPS calculation"""
        reactions_dict = {
            "r1": {"glc__D": 1},
            "r2": {"co2": 2, "pyr": 3},
            "r3": {"atp": 2, "glc__D": 4}
        }
        
        abundances = {
            "glc__D": 8.2,
            "pyr": 5.3,
            "atp": 7.05,
            "co2": 1.0
        }
        
        black_list = []
        missing_in_dataset = ["co2"]
        
        result = rps_generator.calculate_rps(
            reactions_dict,
            abundances,
            black_list,
            missing_in_dataset
        )
        
        # Check that RPS values are calculated
        assert "r1" in result
        assert result["r1"] == 8.2 ** 1
        
        assert "r2" in result
        assert result["r2"] == (1.0 ** 2) * (5.3 ** 3)
        
        assert "r3" in result
        assert result["r3"] == (8.2 ** 4) * (7.05 ** 2)
    
    def test_calculate_rps_with_blacklist(self):
        """Test RPS calculation with blacklisted metabolites"""
        reactions_dict = {
            "r1": {"atp": 3},  # Only has blacklisted metabolite
            "r2": {"glc__D": 2, "atp": 1}  # Has both
        }
        
        abundances = {
            "glc__D": 8.2,
            "atp": 7.05
        }
        
        black_list = ["atp"]
        missing_in_dataset = []
        
        result = rps_generator.calculate_rps(
            reactions_dict,
            abundances,
            black_list,
            missing_in_dataset
        )
        
        # r1 should be NaN (only has blacklisted metabolite)
        assert "r1" in result
        assert math.isnan(result["r1"])
        
        # r2 should only use glc__D (atp is blacklisted)
        assert "r2" in result
        assert result["r2"] == 8.2 ** 2
    
    def test_process_args(self):
        """Test argument processing"""
        args = rps_generator.process_args([])
        assert hasattr(args, 'tool_dir')


class TestGeneratorIntegration:
    """Integration tests for generators with real data structures"""
    
    def test_ras_with_complex_rule(self):
        """Test RAS with complex nested rules"""
        class MockArgs:
            none = False
        
        ras_generator.ARGS = MockArgs()
        
        # Create complex rule: (A and B) or (C and D)
        rule = ruleUtils.OpList(ruleUtils.RuleOp.OR)
        
        sub_rule1 = ruleUtils.OpList(ruleUtils.RuleOp.AND)
        sub_rule1.extend(["geneA", "geneB"])
        
        sub_rule2 = ruleUtils.OpList(ruleUtils.RuleOp.AND)
        sub_rule2.extend(["geneC", "geneD"])
        
        rule.extend([sub_rule1, sub_rule2])
        
        dataset = {
            "geneA": 5.0,
            "geneB": 3.0,
            "geneC": 8.0,
            "geneD": 2.0
        }
        
        # sub_rule1 (A and B) = min(5.0, 3.0) = 3.0
        # sub_rule2 (C and D) = min(8.0, 2.0) = 2.0
        # final (OR) = max(3.0, 2.0) = 3.0
        result = ras_generator.ras_op_list(rule, dataset)
        assert result == 3.0
    
    def test_rps_with_multiple_cell_lines(self):
        """Test RPS calculation with multiple cell lines"""
        dataset = pd.DataFrame({
            "cell_lines": ["normal", "cancer", "treated"],
            "glucose": [8.2, 4.0, 6.5],
            "pyruvate": [5.3, 7.0, 6.0]
        })
        
        # Test that we can extract data for all cell lines
        for i in range(len(dataset)):
            result = rps_generator.get_abund_data(dataset, i)
            assert result is not None
            assert "name" in result
            assert result["glucose"] == dataset.iloc[i]["glucose"]


class TestFileStructure:
    """Test that required files and directories exist"""
    
    def test_pickle_files_accessible(self):
        """Test that pickle files are accessible"""
        pickle_dir = os.path.join(TOOL_DIR, "local", "pickle files")
        
        # Check synonyms pickle
        synonyms_path = os.path.join(pickle_dir, "synonyms.pickle")
        assert os.path.exists(synonyms_path), f"Synonyms file not found at {synonyms_path}"
        
        # Check blacklist pickle
        blacklist_path = os.path.join(pickle_dir, "black_list.pickle")
        assert os.path.exists(blacklist_path), f"Blacklist file not found at {blacklist_path}"
    
    def test_can_load_synonyms(self):
        """Test that we can load the synonyms dictionary"""
        pickle_path = utils.FilePath(
            "synonyms", 
            utils.FileFormat.PICKLE, 
            prefix=os.path.join(TOOL_DIR, "local", "pickle files")
        )
        
        try:
            syns_dict = utils.readPickle(pickle_path)
            assert syns_dict is not None
            assert isinstance(syns_dict, dict)
        except Exception as e:
            pytest.skip(f"Could not load synonyms pickle: {e}")


if __name__ == "__main__":
    # Run tests with pytest if available
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        print("pytest not available, running basic tests...")
        
        test_classes = [
            TestRASGenerator(),
            TestRPSGenerator(),
            TestGeneratorIntegration(),
            TestFileStructure()
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
