# Troubleshooting

Common issues and solutions when using COBRAxy.

## Installation Issues

### Missing Build Tools

**Problem**: `gcc: command not found` or compilation errors (Linux/macOS)
```bash
# Ubuntu/Debian
sudo apt-get install build-essential cmake pkg-config

# macOS
xcode-select --install
brew install cmake pkg-config
```

**Problem**: `CMake not found`
```bash
# Ubuntu/Debian
sudo apt-get install cmake

# macOS
brew install cmake

# Or via conda
conda install -c conda-forge cmake
```

### Python Import Errors

**Problem**: `ModuleNotFoundError: No module named 'cobra'`
```bash
# Solution: Reinstall COBRAxy with dependencies
cd COBRAxy/src
pip install .

# Or install missing dependency directly
pip install cobra
```

**Problem**: `ImportError: No module named 'cobraxy'`  
```python
# Solution: Ensure COBRAxy is installed
pip install /path/to/COBRAxy/src/

# Or add to Python path temporarily
import sys
sys.path.insert(0, '/path/to/COBRAxy/src')
```

### System Dependencies

**Problem**: GLPK solver not found
```bash
# Ubuntu/Debian
sudo apt-get install libglpk40 glpk-utils
pip install swiglpk

# macOS  
brew install glpk
pip install swiglpk

# Windows (using conda)
conda install -c conda-forge glpk swiglpk
```

**Problem**: SVG processing errors
```bash
# Install libvips for image processing
# Ubuntu/Debian: sudo apt-get install libvips
# macOS: brew install vips
```

## Data Format Issues

### Gene Expression Problems

**Problem**: "No computable scores" error
```
Cause: Gene IDs don't match between data and model
Solution: 
1. Check gene ID format (HGNC vs symbols vs Ensembl)
2. Verify first column contains gene identifiers
3. Ensure tab-separated format
4. Try different built-in model
```

**Problem**: Many "gene not found" warnings
```python
# Check gene overlap with model
import pickle
genes_dict = pickle.load(open('src/local/pickle files/ENGRO2_genes.p', 'rb'))
model_genes = set(genes_dict['hugo_id'].keys())

import pandas as pd
data_genes = set(pd.read_csv('expression.tsv', sep='\t').iloc[:, 0])

overlap = len(model_genes.intersection(data_genes))
print(f"Gene overlap: {overlap}/{len(data_genes)} ({overlap/len(data_genes)*100:.1f}%)")
```

**Problem**: File format not recognized
```tsv
# Correct format - tab-separated:
Gene_ID	Sample_1	Sample_2
HGNC:5	10.5	11.2
HGNC:10	3.2	4.1

# Wrong - comma-separated or spaces will fail
```

### Model Issues

**Problem**: Custom model not loading
```
Solution:
1. Check TSV format with "GPR" column header
2. Verify reaction IDs are unique
3. Test GPR syntax (use 'and'/'or', proper parentheses)
4. Check file permissions and encoding (UTF-8)
```

## Tool Execution Errors



### File Path Problems

**Problem**: "File not found" errors
```python
# Use absolute paths
from pathlib import Path

input_file = str(Path('expression.tsv').absolute())

args = ['-in', input_file, ...]
```

**Problem**: Permission denied
```bash
# Check write permissions
ls -la output_directory/

# Fix permissions
chmod 755 output_directory/
chmod 644 input_files/*
```

### Galaxy Integration Issues

**Problem**: COBRAxy tools not appearing in Galaxy
```xml
<!-- Check tool_conf.xml syntax -->
<section id="cobraxy" name="COBRAxy">
  <tool file="cobraxy/ras_generator.xml" />
</section>

<!-- Verify file paths are correct -->
ls tools/cobraxy/ras_generator.xml
```

**Problem**: Tool execution fails in Galaxy
```
Check Galaxy logs:
- main.log: General Galaxy issues
- handler.log: Job execution problems  
- uwsgi.log: Web server issues

Common fixes:
1. Restart Galaxy after adding tools
2. Check Python environment has COBRApy installed
3. Verify file permissions on tool files
```



**Problem**: Flux sampling hangs
```bash
# Check solver availability
python -c "import cobra; print(cobra.Configuration().solver)"

# Should show: glpk, cplex, or gurobi
# Install GLPK if missing:
pip install swiglpk
```

### Large Dataset Handling

**Problem**: Cannot process large expression matrices
```python
# Process in chunks
def process_large_dataset(expression_file, chunk_size=1000):
    df = pd.read_csv(expression_file, sep='\t')
    
    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        chunk_file = f'chunk_{i}.tsv'
        chunk.to_csv(chunk_file, sep='\t', index=False)
        
        # Process chunk
        ras_generator.main(['-in', chunk_file, ...])
```

## Output Validation

### Unexpected Results

**Problem**: All RAS values are zero or null
```python
# Debug gene mapping
import pandas as pd
ras_df = pd.read_csv('ras_output.tsv', sep='\t', index_col=0)

# Check data quality
print(f"Null percentage: {ras_df.isnull().sum().sum() / ras_df.size * 100:.1f}%")
print(f"Zero percentage: {(ras_df == 0).sum().sum() / ras_df.size * 100:.1f}%")

# Check expression data preprocessing
expr_df = pd.read_csv('expression.tsv', sep='\t', index_col=0)
print(f"Expression range: {expr_df.min().min():.2f} to {expr_df.max().max():.2f}")
```

**Problem**: RAS values seem too high/low
```
Possible causes:
1. Expression data not log-transformed
2. Wrong normalization method
3. Incorrect gene ID mapping
4. GPR rule interpretation issues

Solutions:
1. Check expression data preprocessing
2. Validate against known control genes
3. Compare with published metabolic activity patterns
```

### Missing Pathway Maps

**Problem**: MAREA generates no output maps
```
Debug steps:
1. Check RAS input has non-null values
2. Verify model choice matches RAS generation
3. Check statistical significance thresholds
4. Look at log files for specific errors
```

## Environment Issues

### Conda/Virtual Environment Problems

**Problem**: Tool import fails in virtual environment
```bash
# Activate environment properly
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate  # Windows

# Verify COBRAxy installation
pip list | grep cobra
python -c "import cobra; print('COBRApy version:', cobra.__version__)"
```

**Problem**: Version conflicts
```bash
# Create clean environment
conda create -n cobraxy python=3.9
conda activate cobraxy

# Install COBRAxy fresh
cd COBRAxy/src
pip install -e .
```

### Cross-Platform Issues

**Problem**: Windows path separator issues
```python
# Use pathlib for cross-platform paths
from pathlib import Path

# Instead of: '/path/to/file'  
# Use: str(Path('path') / 'to' / 'file')
```

**Problem**: Line ending issues (Windows/Unix)
```bash
# Convert line endings if needed
dos2unix input_file.tsv  # Unix
unix2dos input_file.tsv  # Windows
```

## Debugging Strategies

### Enable Detailed Logging

```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Many tools accept log file parameter
args = [..., '--out_log', 'detailed.log']
```

### Test with Small Datasets

```python
# Create minimal test case
test_data = """Gene_ID	Sample1	Sample2
HGNC:5	10.0	15.0
HGNC:10	5.0	8.0"""

with open('test_input.tsv', 'w') as f:
    f.write(test_data)

# Test basic functionality
ras_generator.main(['-in', 'test_input.tsv', 
                   '-ra', 'test_output.tsv', '-rs', 'ENGRO2'])
```

### Check Dependencies

```python
# Verify all required packages
required_packages = ['cobra', 'pandas', 'numpy', 'scipy']

for package in required_packages:
    try:
        __import__(package)
        print(f"✓ {package}")
    except ImportError:
        print(f"✗ {package} - MISSING")
```

## Getting Help

### Information to Include in Bug Reports

When reporting issues, include:

1. **System information**:
   ```bash
   python --version
   pip list | grep cobra
   uname -a  # Linux/macOS
   ```

2. **Complete error messages**: Copy full traceback
3. **Input file format**: First few lines of input data
4. **Command/parameters used**: Exact command or Python code
5. **Expected vs actual behavior**: What should happen vs what happens

### Community Resources

- **GitHub Issues**: [Report bugs and ask questions](https://github.com/CompBtBs/COBRAxy/issues)
- **COBRApy Community**: [General metabolic modeling help](https://github.com/opencobra/cobrapy)

### Self-Help Checklist

Before reporting issues:

- Checked this troubleshooting guide
- Verified installation completeness
- Tested with built-in example data
- Searched existing GitHub issues
- Tried alternative models/parameters
- Checked file formats and permissions

## Prevention Tips

### Best Practices

1. **Use virtual environments** to avoid conflicts
2. **Validate input data** before processing
3. **Start with small datasets** for testing
4. **Keep backups** of working configurations
5. **Document successful workflows** for reuse
6. **Test after updates** to catch regressions

### Data Quality Checks

```python
def validate_expression_data(filename):
    """Validate gene expression file format."""
    df = pd.read_csv(filename, sep='\t')
    
    # Check basic format
    assert df.shape[0] > 0, "Empty file"
    assert df.shape[1] > 1, "Need at least 2 columns"
    
    # Check numeric data  
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    assert len(numeric_cols) > 0, "No numeric expression data"
    
    # Check for missing values
    null_pct = df.isnull().sum().sum() / df.size * 100
    if null_pct > 50:
        print(f"Warning: {null_pct:.1f}% missing values")
    
    print(f"✓ File valid: {df.shape[0]} genes × {df.shape[1]-1} samples")
```

This troubleshooting guide covers the most common issues. For tool-specific problems, check the individual tool documentation pages.