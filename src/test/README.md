# COBRAxy Test Suite

This directory contains comprehensive unit tests for all COBRAxy modules.

## Test Files

- **test_utils.py** - Tests for utility modules (general_utils, rule_parsing, reaction_parsing, model_utils, CBS_backend)
- **test_generators.py** - Tests for RAS and RPS generator modules
- **test_marea.py** - Tests for MAREA, flux_simulation, flux_to_map, and visualization modules
- **test_clustering.py** - Tests for marea_cluster and clustering algorithms
- **testing.py** - Original testing framework (legacy)

## Running Tests

### Option 1: Using pytest (Recommended)

Install pytest if not already installed:
```bash
pip install pytest pytest-cov
```

Run all tests:
```bash
cd /hdd/home/flapi/COBRAxy/src/test
pytest -v
```

Run specific test file:
```bash
pytest test_utils.py -v
pytest test_generators.py -v
pytest test_marea.py -v
pytest test_clustering.py -v
```

Run tests with coverage:
```bash
pytest --cov=../ --cov-report=html
```

### Option 2: Run individual test files

Each test file can be run standalone:
```bash
python test_utils.py
python test_generators.py
python test_marea.py
python test_clustering.py
```

### Option 3: Run all tests with the run script

```bash
python run_all_tests.py
```

## Test Structure

Each test file is organized into classes that group related tests:

```python
class TestModuleName:
    """Tests for module_name"""
    
    def test_specific_feature(self):
        """Test description"""
        # Test code
        assert result == expected
```

## Adding New Tests

To add new tests:

1. Choose the appropriate test file or create a new one
2. Create a test class if needed
3. Add test methods (must start with `test_`)
4. Use assertions to verify behavior

Example:
```python
class TestMyFeature:
    def test_my_new_function(self):
        """Test my new function"""
        result = my_function(input_data)
        assert result == expected_output
```

## Test Coverage

Current test coverage includes:

### Utils Module
- ✓ Bool validator
- ✓ CustomErr class
- ✓ FilePath creation and validation
- ✓ Model enum
- ✓ Rule parsing
- ✓ Reaction parsing
- ✓ Model utilities

### Generators
- ✓ RAS calculation with AND/OR rules
- ✓ RPS calculation with metabolite abundances
- ✓ Missing metabolite handling
- ✓ Blacklist functionality
- ✓ Complex nested rules

### MAREA and Visualization
- ✓ Argument processing
- ✓ Data format validation
- ✓ Statistical operations (fold change, p-values)
- ✓ SVG map visualization
- ✓ Model conversion tools

### Clustering
- ✓ K-means clustering
- ✓ DBSCAN clustering
- ✓ Hierarchical clustering
- ✓ Data scaling/normalization
- ✓ Cluster evaluation metrics
- ✓ Visualization preparation

## Dependencies

Required packages for running tests:
- pytest (optional but recommended)
- pandas
- numpy
- scipy
- scikit-learn
- cobra
- lxml

All dependencies are listed in the main setup.py file.

## Continuous Integration

Tests can be integrated into CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest pytest-cov
      - name: Run tests
        run: |
          cd src/test
          pytest -v --cov
```

## Troubleshooting

### Import Errors
If you get import errors, make sure you're running from the test directory:
```bash
cd /hdd/home/flapi/COBRAxy/src/test
python test_utils.py
```

### Missing Files
If tests fail due to missing pickle files or models, verify that the `local/` directory structure is intact:
```
src/
├── local/
│   ├── pickle files/
│   ├── svg metabolic maps/
│   ├── models/
│   ├── mappings/
│   └── medium/
```

### Path Issues
Tests automatically set up paths relative to the test directory. If you encounter path issues, check that `TOOL_DIR` is set correctly in the test file.

## Contributing

When adding new features to COBRAxy:
1. Write tests for the new functionality
2. Ensure all existing tests still pass
3. Aim for >80% code coverage
4. Document any new test files in this README

## License

Same license as COBRAxy main project.
