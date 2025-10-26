# Installation

This guide walks you through installing COBRAxy on your system.

## System Requirements

- **Python**: 3.8-3.13
- **Operating System**: Linux (recommended), macOS, Windows
- **Build tools**: C/C++ compiler (gcc, clang, or MSVC), CMake, pkg-config

## System Dependencies

Install required build tools before installing COBRAxy:

```bash
# Ubuntu/Debian
sudo apt-get install build-essential cmake pkg-config libvips libglpk40 glpk-utils

# macOS
xcode-select --install
brew install cmake pkg-config vips glpk

# Windows (with Chocolatey)
choco install cmake visualstudio2022buildtools pkgconfiglite
```

## Installation Methods

### Recommended: Using Conda

Create an isolated environment with all dependencies:

```bash
# Create a new conda environment
conda create -n cobraxy python=3.13 -y
conda activate cobraxy

# Install build tools via conda (optional, if not using system packages)
conda install -c conda-forge cmake pkg-config -y

# Clone and install COBRAxy
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy/src
pip install .

# Optional: Install GLPK bindings for enhanced performance
pip install swiglpk
```

### Alternative: Direct Installation

If you have system dependencies already installed:

```bash
# Clone the repository
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy

# Install COBRAxy
pip install src/
```

### Development Install

For development or if you want to modify COBRAxy:

```bash
# Clone and install in editable mode
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy
pip install -e src/
```

## What Gets Installed

COBRAxy automatically installs Python dependencies including:
- **COBRApy**: Core metabolic modeling
- **pandas, numpy**: Data handling
- **matplotlib, seaborn**: Visualization
- **scikit-learn, scipy**: Statistical analysis
- **pyvips**: Image processing
- **anndata, pydeseq2**: Advanced analysis

## Verify Installation

Test your installation:

```bash
# Check if COBRAxy tools are available
ras_generator --help
flux_simulation --help
marea --help

# Check Python can import COBRAxy modules
python -c "import ras_generator; print('COBRAxy installed successfully!')"
```

## Troubleshooting Installation

### Missing Compiler Errors

If you see errors about missing compilers during installation:

```bash
# Ubuntu/Debian
sudo apt-get install build-essential

# macOS
xcode-select --install
```

### CMake Not Found

```bash
# Ubuntu/Debian
sudo apt-get install cmake

# macOS
brew install cmake

# Or via conda
conda install -c conda-forge cmake
```

### pkg-config Issues

```bash
# Ubuntu/Debian
sudo apt-get install pkg-config

# macOS
brew install pkg-config

# Or via conda
conda install -c conda-forge pkg-config
```

## Alternative: Virtual Environment (without Conda)

Using a virtual environment prevents conflicts with other Python packages:

```bash
# Create virtual environment
python -m venv cobraxy-env

# Activate environment
source cobraxy-env/bin/activate  # Linux/macOS
# cobraxy-env\Scripts\activate  # Windows

# Install COBRAxy
cd COBRAxy
pip install src/

# When done, deactivate
deactivate
```

## Next Steps

After successful installation:

1. **[Quick Start Guide](quickstart.md)** - Run your first analysis
2. **[Tutorial: Galaxy Setup](tutorials/galaxy-setup.md)** - Set up web interface

## Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Search [existing issues](https://github.com/CompBtBs/COBRAxy/issues)
3. Create a [new issue](https://github.com/CompBtBs/COBRAxy/issues/new) with:
   - Your operating system
   - Python version (`python --version`)
   - Complete error message
   - Installation method used