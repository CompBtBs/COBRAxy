# Installation

This guide walks you through installing COBRAxy on your system.

## System Requirements

- **Python**: 3.8-3.13
- **Operating System**: Linux (recommended), macOS, Windows

## Quick Install

The fastest way to install COBRAxy:

```bash
# Clone the repository
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy

# Install COBRAxy
pip install src/
```

## Development Install

For development or if you want to modify COBRAxy:

```bash
# Clone and install in development mode
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy
pip install -e src/
```

## Dependencies

COBRAxy automatically installs its Python dependencies (COBRApy, pandas, numpy, etc.)

## Verify Installation

Test your installation:

```bash
# Check if COBRAxy tools are available
ras_generator --help
flux_simulation --help
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