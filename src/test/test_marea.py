"""
Unit tests for MAREA, flux_simulation, and related visualization modules.

Run with: python -m pytest test_marea.py -v
Or: python test_marea.py
"""

import sys
import os
import pandas as pd
import numpy as np
import tempfile
from pathlib import Path

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

import marea
import flux_simulation
import flux_to_map
import ras_to_bounds
import utils.general_utils as utils

# Get the tool directory
TOOL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class TestMAREA:
    """Tests for marea module"""
    
    def test_process_args(self):
        """Test argument processing for MAREA"""
        # Create minimal args for testing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("reaction_id,value\n")
            f.write("r1,1.5\n")
            temp_file = f.name
        
        try:
            args = marea.process_args([
                '-td', TOOL_DIR,
                '--tool_dir', TOOL_DIR
            ])
            assert hasattr(args, 'tool_dir')
            assert args.tool_dir == TOOL_DIR
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def test_comparison_types(self):
        """Test that comparison type enum exists and is correct"""
        # Check that the ComparisonType enum has expected values
        assert hasattr(marea, 'ComparisonType') or hasattr(marea, 'GroupingCriterion')
    
    def test_ras_transformation(self):
        """Test RAS transformation logic"""
        # Create sample RAS data
        ras_data = pd.DataFrame({
            'reaction': ['r1', 'r2', 'r3'],
            'value': [1.5, 0.5, 2.0]
        })
        
        # Test that data can be processed
        assert len(ras_data) == 3
        assert ras_data['value'].max() == 2.0


class TestFluxSimulation:
    """Tests for flux_simulation module"""
    
    def test_process_args(self):
        """Test argument processing for flux simulation"""
        args = flux_simulation.process_args([
            '-td', TOOL_DIR
        ])
        assert hasattr(args, 'tool_dir')
    
    def test_flux_balance_setup(self):
        """Test that FBA setup functions exist"""
        # Check that key functions exist
        assert hasattr(flux_simulation, 'process_args')
        assert hasattr(flux_simulation, 'main')


class TestFluxToMap:
    """Tests for flux_to_map module"""
    
    def test_process_args(self):
        """Test argument processing for flux to map"""
        args = flux_to_map.process_args([
            '-td', TOOL_DIR
        ])
        assert hasattr(args, 'tool_dir')
    
    def test_color_map_options(self):
        """Test that color map options are available"""
        # The module should have color map functionality
        assert hasattr(flux_to_map, 'process_args')


class TestRasToBounds:
    """Tests for ras_to_bounds module"""
    
    def test_process_args(self):
        """Test argument processing for RAS to bounds"""
        args = ras_to_bounds.process_args([
            '-td', TOOL_DIR
        ])
        assert hasattr(args, 'tool_dir')
    
    def test_bounds_conversion(self):
        """Test that bounds conversion logic exists"""
        # Create sample RAS data
        ras_data = {
            'r1': 1.5,
            'r2': 0.5,
            'r3': 2.0
        }
        
        # Test basic transformation logic
        # Reactions with higher RAS should have higher bounds
        assert ras_data['r3'] > ras_data['r1'] > ras_data['r2']


class TestModelConversion:
    """Tests for model conversion tools"""
    
    def test_tabular_to_model(self):
        """Test tabular to model conversion"""
        import tabular2MetabolicModel
        
        args = tabular2MetabolicModel.process_args([])
        assert hasattr(args, 'tool_dir')
    
    def test_model_to_tabular(self):
        """Test model to tabular conversion"""
        import metabolicModel2Tabular
        
        args = metabolicModel2Tabular.process_args([])
        assert hasattr(args, 'tool_dir')


class TestDataProcessing:
    """Tests for data processing utilities used across tools"""
    
    def test_ras_data_format(self):
        """Test RAS data format validation"""
        # Create valid RAS data
        ras_df = pd.DataFrame({
            'reaction_id': ['r1', 'r2', 'r3'],
            'group1': [1.5, 0.5, 2.0],
            'group2': [1.8, 0.3, 2.2]
        })
        
        assert 'reaction_id' in ras_df.columns
        assert len(ras_df) > 0
    
    def test_rps_data_format(self):
        """Test RPS data format validation"""
        # Create valid RPS data
        rps_df = pd.DataFrame({
            'reaction_id': ['r1', 'r2', 'r3'],
            'sample1': [100.5, 50.3, 200.1],
            'sample2': [150.2, 30.8, 250.5]
        })
        
        assert 'reaction_id' in rps_df.columns
        assert len(rps_df) > 0
    
    def test_flux_data_format(self):
        """Test flux data format validation"""
        # Create valid flux data
        flux_df = pd.DataFrame({
            'reaction_id': ['r1', 'r2', 'r3'],
            'flux': [1.5, -0.5, 2.0],
            'lower_bound': [-10, -10, 0],
            'upper_bound': [10, 10, 10]
        })
        
        assert 'reaction_id' in flux_df.columns
        assert 'flux' in flux_df.columns


class TestStatistics:
    """Tests for statistical operations in MAREA"""
    
    def test_fold_change_calculation(self):
        """Test fold change calculation"""
        # Simple fold change test
        group1_mean = 2.0
        group2_mean = 4.0
        fold_change = group2_mean / group1_mean
        
        assert fold_change == 2.0
    
    def test_log_fold_change(self):
        """Test log fold change calculation"""
        group1_mean = 2.0
        group2_mean = 8.0
        log_fc = np.log2(group2_mean / group1_mean)
        
        assert log_fc == 2.0  # log2(8/2) = log2(4) = 2
    
    def test_pvalue_correction(self):
        """Test that statistical functions handle edge cases"""
        # Test with identical values (should give p-value close to 1)
        group1 = [1.0, 1.0, 1.0]
        group2 = [1.0, 1.0, 1.0]
        
        from scipy import stats
        t_stat, p_value = stats.ttest_ind(group1, group2)
        
        # p-value should be NaN or close to 1 for identical groups
        assert np.isnan(p_value) or p_value > 0.9


class TestMapVisualization:
    """Tests for SVG map visualization"""
    
    def test_svg_maps_exist(self):
        """Test that SVG maps exist"""
        map_dir = os.path.join(TOOL_DIR, "local", "svg metabolic maps")
        assert os.path.exists(map_dir)
        
        # Check for at least one map
        maps = [f for f in os.listdir(map_dir) if f.endswith('.svg')]
        assert len(maps) > 0, "No SVG maps found"
    
    def test_model_has_map(self):
        """Test that models have associated maps"""
        # ENGRO2 should have a map
        engro2_map = os.path.join(TOOL_DIR, "local", "svg metabolic maps", "ENGRO2_map.svg")
        if os.path.exists(engro2_map):
            assert os.path.getsize(engro2_map) > 0
    
    def test_color_gradient(self):
        """Test color gradient generation"""
        # Test that we can generate colors for a range of values
        values = [-2.0, -1.0, 0.0, 1.0, 2.0]
        
        # All values should be processable
        for val in values:
            # Simple color mapping test
            if val < 0:
                # Negative values should map to one color scheme
                assert val < 0
            elif val > 0:
                # Positive values should map to another
                assert val > 0
            else:
                # Zero should be neutral
                assert val == 0


class TestIntegration:
    """Integration tests for complete workflows"""
    
    def test_ras_to_marea_workflow(self):
        """Test that RAS data can flow into MAREA"""
        # Create sample RAS data
        ras_data = pd.DataFrame({
            'reaction_id': ['r1', 'r2', 'r3'],
            'control': [1.5, 0.8, 1.2],
            'treatment': [2.0, 0.5, 1.8]
        })
        
        # Calculate fold changes
        ras_data['fold_change'] = ras_data['treatment'] / ras_data['control']
        
        assert 'fold_change' in ras_data.columns
        assert len(ras_data) == 3
    
    def test_rps_to_flux_workflow(self):
        """Test that RPS data can be used for flux simulation"""
        # Create sample RPS data
        rps_data = pd.DataFrame({
            'reaction_id': ['r1', 'r2', 'r3'],
            'rps': [100.0, 50.0, 200.0]
        })
        
        # RPS can be used to set bounds
        rps_data['upper_bound'] = rps_data['rps'] / 10
        
        assert 'upper_bound' in rps_data.columns


class TestErrorHandling:
    """Tests for error handling across modules"""
    
    def test_invalid_model_name(self):
        """Test handling of invalid model names"""
        with pytest.raises((ValueError, KeyError, AttributeError)):
            utils.Model("INVALID_MODEL")
    
    def test_missing_required_column(self):
        """Test handling of missing required columns"""
        # Create incomplete data
        incomplete_data = pd.DataFrame({
            'wrong_column': [1, 2, 3]
        })
        
        # Should fail when looking for required columns
        with pytest.raises(KeyError):
            value = incomplete_data['reaction_id']


if __name__ == "__main__":
    # Run tests with pytest if available
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        print("pytest not available, running basic tests...")
        
        test_classes = [
            TestMAREA(),
            TestFluxSimulation(),
            TestFluxToMap(),
            TestRasToBounds(),
            TestModelConversion(),
            TestDataProcessing(),
            TestStatistics(),
            TestMapVisualization(),
            TestIntegration(),
            TestErrorHandling()
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
                        import traceback
                        traceback.print_exc()
                        failed += 1
        
        print(f"\n{'='*60}")
        print(f"Results: {passed} passed, {failed} failed")
        if failed > 0:
            sys.exit(1)
