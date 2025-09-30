# Installation

This guide walks you through installing COBRAxy on your system.

## System Requirements

- **Python**: 3.8-3.11
- **Operating System**: Linux (recommended), macOS, Windows
- **Storage**: 2GB free space for installation and temporary files

## Quick Install

The fastest way to install COBRAxy:

```bash
# Clone the repository
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy

# Install COBRAxy
pip install .
```

## Development Install

For development or if you want to modify COBRAxy:

```bash
# Clone and install in development mode
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy
pip install -e .
```

## Dependencies

COBRAxy automatically installs its Python dependencies:

- **COBRApy** - Core metabolic modeling
- **pandas** - Data manipulation
- **numpy** - Numerical computations
- **scipy** - Scientific computing

## Optional System Libraries

Install additional libraries for enhanced features:

### Ubuntu/Debian

```bash
# Install GLPK solver
sudo apt-get update
sudo apt-get install libglpk40 glpk-utils

# Install libvips for SVG processing
sudo apt-get install libvips

# Install Python GLPK bindings
pip install swiglpk
```

### macOS

```bash
# Using Homebrew
brew install glpk vips

# Install Python bindings
pip install swiglpk
```

### Windows

```bash
# Using conda (recommended for Windows)
conda install -c conda-forge glpk

# Or using pip
pip install swiglpk
```

## Verify Installation

Test your installation:

```bash
# Check if COBRAxy tools are available
ras_generator --help
flux_simulation --help

# Test with example data (if available)
cd COBRAxy
python testing.py
```

## Troubleshooting Installation

### Common Issues

**Import Error: No module named 'cobra'**
```bash
# Install COBRApy manually
pip install cobra
```

**GLPK solver not found**
```bash
# Install GLPK solver
# Ubuntu/Debian: sudo apt-get install glpk-utils
# macOS: brew install glpk
# Then: pip install swiglpk
```

**Permission denied errors**
```bash
# Use user installation
pip install --user .
# Or use virtual environment (recommended)
python -m venv cobraxy-env
source cobraxy-env/bin/activate  # Linux/macOS
# cobraxy-env\Scripts\activate  # Windows
pip install .
```

## Virtual Environment (Recommended)

Using a virtual environment prevents conflicts with other Python packages:

```bash
# Create virtual environment
python -m venv cobraxy-env

# Activate environment
source cobraxy-env/bin/activate  # Linux/macOS
# cobraxy-env\Scripts\activate  # Windows

# Install COBRAxy
pip install .

# When done, deactivate
deactivate
```

## Next Steps

After successful installation:

1. **[Quick Start Guide](quickstart.md)** - Run your first analysis
2. **[Tutorial: Python API](tutorials/python-api.md)** - Learn programmatic usage
3. **[Tutorial: Galaxy Setup](tutorials/galaxy-setup.md)** - Set up web interface

## Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Search [existing issues](https://github.com/CompBtBs/COBRAxy/issues)
3. Create a [new issue](https://github.com/CompBtBs/COBRAxy/issues/new) with:
   - Your operating system
   - Python version (`python --version`)
   - Complete error message
   - Installation method used