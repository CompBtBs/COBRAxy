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


## Galaxy Tool Issues

### 1. Flux simulation 

**Error message**: Execution aborted: wrong format of bounds dataset

**Meaning:**  
Flux simulation cannot read the bounds of the metabolic model for the constrained simulation problem (optimization or sampling).  
This usually happens if the input “Bound file(s): *” is incorrect. For example, it occurs when the **RasToBounds - Cell Class** file is passed instead of the collection of bound files named **"RAS to bounds"**.

**Suggested Action:**  
Check the input files and ensure the correct bounds collection is used.

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


This troubleshooting guide covers the most common issues. For tool-specific problems, check the individual tool documentation pages.
