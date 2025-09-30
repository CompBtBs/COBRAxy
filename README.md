# COBRAxy

A Python toolkit for metabolic flux analysis and visualization, with Galaxy integration.

COBRAxy transforms gene expression and metabolite data into meaningful metabolic insights through flux sampling and interactive pathway maps.
DOC: https://compbtbs.github.io/COBRAxy
## Features

- **Reaction Activity Scores (RAS)** from gene expression data
- **Reaction Propensity Scores (RPS)** from metabolite abundance
- **Flux sampling** with CBS or OptGP algorithms  
- **Statistical analysis** with pFBA, FVA, and sensitivity analysis
- **Interactive maps** with SVG/PDF export and custom styling
- **Galaxy tools** for web-based analysis
- **Built-in models** including ENGRO2 and Recon

## Quick Start

### Installation

```bash
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy
pip install .
```

### Basic Workflow

```bash
# 1. Generate RAS from expression data
ras_generator -td $(pwd) -in expression.tsv -ra ras_output.tsv -rs ENGRO2

# 2. Generate RPS from metabolite data (optional)
rps_generator -td $(pwd) -id metabolites.tsv -rp rps_output.tsv

# 3. Create enriched pathway maps with statistical analysis
marea -td $(pwd) -using_RAS true -input_data ras_output.tsv -choice_map ENGRO2 -gs true -idop base_maps

# 4. Apply RAS constraints to model for flux simulation
ras_to_bounds -td $(pwd) -ms ENGRO2 -ir ras_output.tsv -rs true -idop bounds_output

# 5. Sample metabolic fluxes with constrained model
flux_simulation -td $(pwd) -ms ENGRO2 -in bounds_output/*.tsv -a CBS -ns 1000 -idop flux_results

# 6. Add flux data to enriched maps
flux_to_map -td $(pwd) -if flux_results/*.tsv -mp base_maps/*.svg -idop final_maps
```

## Tools

| Tool | Purpose | Input | Output |
|------|---------|--------|---------|
| `metabolic_model_setting` | Extract model components | SBML model | Rules, reactions, bounds, medium |
| `ras_generator` | Compute reaction activity scores | Gene expression data | RAS values |
| `rps_generator` | Compute reaction propensity scores | Metabolite abundance | RPS values |
| `marea` | Statistical pathway analysis | RAS + RPS data | Enrichment + base maps |
| `ras_to_bounds` | Apply RAS constraints to model | RAS + SBML model | Constrained bounds |
| `flux_simulation` | Sample metabolic fluxes | Constrained model | Flux distributions |
| `flux_to_map` | Add fluxes to enriched maps | Flux samples + base maps | Final styled maps |
| `marea_cluster` | Cluster analysis | Expression/flux data | Sample clusters |

## Requirements

- **Python**: 3.8-3.11
- **OS**: Linux, macOS, Windows (Linux recommended)
- **Dependencies**: Automatically installed via pip (COBRApy, pandas, numpy, etc.)

**Optional system libraries** (for enhanced features):
```bash
# Ubuntu/Debian
sudo apt-get install libvips libglpk40 glpk-utils

# For Python GLPK bindings
pip install swiglpk
```

## Data Flow

```
Gene Expression    Metabolite Data    SBML Model
      ↓                   ↓               ↓
  RAS Generator      RPS Generator   Model Tables
      ↓                   ↓               
    RAS Values       RPS Values           
    | ↓                   ↓               
    | └─────────┬─────────┘               
    |           ↓                         
    |        MAREA                        
    |    (Enrichment +                    
    |     Base Maps)                      
    ↓                
    RAS Values  →  RAS to Bounds  ←── Model Tables
                        ↓
                  Constrained Model
                        ↓
                  Flux Simulation
                        ↓
                   Flux Samples
                        ↓
                   Flux to Map  ←── Maps (ENGRO2)
                        ↓
               Final Enriched Maps
```

## Built-in Models & Data

COBRAxy includes ready-to-use resources:

- **Models**: ENGRO2, Recon (human metabolism)
- **Gene mappings**: HGNC, Ensembl, Entrez ID conversions
- **Pathway maps**: Pre-styled SVG templates
- **Medium compositions**: Standard growth conditions

Located in `local/` directory for immediate use.

## Command Line Usage

All tools support `--help` for detailed options. Key commands:

### Generate RAS/RPS scores
```bash
# From gene expression
ras_generator -td $(pwd) -in expression.tsv -ra ras_output.tsv -rs ENGRO2

# From metabolite data  
rps_generator -td $(pwd) -id metabolites.tsv -rp rps_output.tsv
```

### Flux sampling
```bash
flux_simulation -td $(pwd) -ms ENGRO2 -in bounds/*.tsv -a CBS -ns 1000 -idop results/
```

### Statistical analysis & visualization
```bash
marea -td $(pwd) -using_RAS true -input_data ras.tsv -choice_map ENGRO2 -gs true -idop maps/
```

## Galaxy Integration

COBRAxy provides Galaxy tool wrappers (`.xml` files) for web-based analysis:

- Upload data through Galaxy interface
- Chain tools in visual workflows  
- Share and reproduce analyses
- Access via Galaxy ToolShed

## Tutorials

### Local Galaxy Installation

To set up a local Galaxy instance with COBRAxy tools:

1. **Install Galaxy**:
   ```bash
   # Clone Galaxy repository
   git clone -b release_23.1 https://github.com/galaxyproject/galaxy.git
   cd galaxy
   
   # Install dependencies and start Galaxy
   sh run.sh
   ```

2. **Install COBRAxy tools**:
   ```bash
   # Add COBRAxy tools to Galaxy
   mkdir -p tools/cobraxy
   cp path/to/COBRAxy/Galaxy_tools/*.xml tools/cobraxy/
   
   # Update tool_conf.xml to include COBRAxy tools
   # Add section in config/tool_conf.xml:
   # <section id="cobraxy" name="COBRAxy">
   #   <tool file="cobraxy/ras_generator.xml" />
   #   <tool file="cobraxy/rps_generator.xml" />
   #   <tool file="cobraxy/marea.xml" />
   #   <!-- Add other tools -->
   # </section>
   ```

3. **Galaxy Tutorial Resources**:
   - [Galaxy Installation Guide](https://docs.galaxyproject.org/en/master/admin/)
   - [Tool Development Tutorial](https://training.galaxyproject.org/training-material/topics/dev/)
   - [Galaxy Admin Training](https://training.galaxyproject.org/training-material/topics/admin/)

### Python Direct Usage

For programmatic use of COBRAxy tools in Python scripts:

1. **Installation for Development**:
   ```bash
   # Clone and install in development mode
   git clone https://github.com/CompBtBs/COBRAxy.git
   cd COBRAxy
   pip install -e .
   ```

2. **Python API Usage**:
   ```python
   import sys
   import os
   
   # Add COBRAxy to Python path
   sys.path.append('/path/to/COBRAxy')
   
   # Import tool modules
   import ras_generator
   import rps_generator
   import flux_simulation
   import marea
   import ras_to_bounds
   
   # Set working directory
   tool_dir = "/path/to/COBRAxy"
   os.chdir(tool_dir)
   
   # Generate RAS scores
   ras_args = [
       '-td', tool_dir,
       '-in', 'data/expression.tsv',
       '-ra', 'output/ras_values.tsv',
       '-rs', 'ENGRO2'
   ]
   ras_generator.main(ras_args)
   
   # Generate RPS scores (optional)
   rps_args = [
       '-td', tool_dir,
       '-id', 'data/metabolites.tsv',
       '-rp', 'output/rps_values.tsv'
   ]
   rps_generator.main(rps_args)
   
   # Create enriched pathway maps
   marea_args = [
       '-td', tool_dir,
       '-using_RAS', 'true',
       '-input_data', 'output/ras_values.tsv',
       '-choice_map', 'ENGRO2',
       '-gs', 'true',
       '-idop', 'maps'
   ]
   marea.main(marea_args)
   
   # Apply RAS constraints to model
   bounds_args = [
       '-td', tool_dir,
       '-ms', 'ENGRO2',
       '-ir', 'output/ras_values.tsv',
       '-rs', 'true',
       '-idop', 'bounds'
   ]
   ras_to_bounds.main(bounds_args)
   
   # Sample metabolic fluxes
   flux_args = [
       '-td', tool_dir,
       '-ms', 'ENGRO2',
       '-in', 'bounds/bounds_output.tsv',
       '-a', 'CBS',
       '-ns', '1000',
       '-idop', 'flux_results'
   ]
   flux_simulation.main(flux_args)
   ```

3. **Python Tutorial Resources**:
   - [COBRApy Documentation](https://cobrapy.readthedocs.io/)
   - [Metabolic Modeling with Python](https://opencobra.github.io/cobrapy/building_model.html)
   - [Flux Sampling Tutorial](https://cobrapy.readthedocs.io/en/stable/sampling.html)
   - [Jupyter Notebooks Examples](examples/) (included in repository)

## Input/Output Formats

| Data Type | Format | Description |
|-----------|---------|-------------|
| Gene expression | TSV | Genes (rows) × Samples (columns) |
| Metabolites | TSV | Metabolites (rows) × Samples (columns) |  
| Models | SBML | Standard metabolic model format |
| Results | TSV/CSV | Tabular flux/score data |
| Maps | SVG/PDF | Styled pathway visualizations |

## Troubleshooting

**Common issues:**

- **Missing GLPK**: Install `glpk-utils` and `swiglpk` for optimal CBS performance
- **SVG errors**: Install `libvips` system library
- **Memory issues**: Reduce sampling count (`-ns`) or use fewer batches (`-nb`)

## Contributing

Contributions welcome! Please:
- Follow existing code style
- Add documentation for new features
- Test with provided example data
- Submit focused pull requests

## Citation

If you use COBRAxy in research, please cite:
- [COBRApy](https://opencobra.github.io/cobrapy/) for core metabolic modeling
- [MaREA](https://galaxyproject.org/use/marea4galaxy/) for enrichment methods
- This repository for integrated workflow

## Links

- [COBRApy Documentation](https://opencobra.github.io/cobrapy/)
- [Galaxy Project](https://usegalaxy.org/)
- [GSoC 2024 Project](https://summerofcode.withgoogle.com/programs/2024/projects/LSrCKfq7)
