<div align="center">
  <img src="docs/_media/logo.png" alt="COBRAxy Logo" width="200"/>
</div>

# COBRAxy

A Python-based command-line suite for metabolic flux analysis and visualization, with [Galaxy](http://marea4galaxy.cloud.ba.infn.it/galaxy) integration.

COBRAxy enables the integration of transcriptomics data with COBRA-based metabolic models, offering a comprehensive framework for studying metabolism in both health and disease. With COBRAxy, users can load and enrich metabolic models by incorporating transcriptomic data and adjusting the model's medium conditions.
DOC: https://compbtbs.github.io/COBRAxy
## Features

- **Galaxy Tools** - Web-based analysis with intuitive interface
- **Import/Export** of metabolic models in multiple formats (SBML, JSON, MAT, YAML)
- **Reaction Activity Scores (RAS)** - Compute metabolic activity from gene expression data
- **Reaction Propensity Scores (RPS)** - Infer metabolic preferences from metabolite abundance
- **Flux computation** - Compute metabolic flux distributions using different optimization or sampling algorithms
- **Statistical Analysis** - Perform statistically significant flux differences between groups of samples and report on an enriched metabolic map
- **Built-in Models** - Ready-to-use models including ENGRO2 and Recon3D


## Requirements

- **Python**: 3.8-3.13
- **OS**: Linux, macOS, Windows (Linux/macOS recommended)
- **Dependencies**: Automatically installed via pip (COBRApy, pandas, numpy, etc.)
- **Build tools**: C/C++ compiler (gcc, clang, or MSVC), CMake for compiling Python extensions, pkg-config

### Installation

**Recommended: Using Conda**

```bash
# Create a new conda environment
conda create -n cobraxy python=3.13 -y
conda activate cobraxy

# Install system dependencies via conda (optional, if not using system packages)
conda install -c conda-forge gcc cmake pkg-config swiglpk -y

# Clone and install COBRAxy
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy/src
pip install .
```

### Basic Workflow

```bash
# 1. Generate RAS from expression data
ras_generator -in expression.tsv -ra ras_output.tsv -rs ENGRO2

# 2. Generate RPS from metabolite data (optional)
rps_generator -id metabolites.tsv -rp rps_output.tsv

# 3. Create enriched pathway maps with statistical analysis
marea -using_RAS true -input_data ras_output.tsv -choice_map ENGRO2 -gs true -idop base_maps

# 4. Apply RAS constraints to model for flux simulation
ras_to_bounds -ms ENGRO2 -ir ras_output.tsv -rs true -idop bounds_output

# 5. Sample metabolic fluxes with constrained model
flux_simulation -ms ENGRO2 -in bounds_output/*.tsv -a CBS -ns 1000 -idop flux_results

# 6. Add flux data to enriched maps
flux_to_map -if flux_results/*.tsv -mp base_maps/*.svg -idop final_maps
```

## Tools

| Tool | Purpose | Input | Output |
|------|---------|--------|---------|
| [`importMetabolicModel`](https://compbtbs.github.io/COBRAxy/#/tools/import-metabolic-model) | Import and extract model components | SBML/JSON/MAT/YAML model | Tabular model data |
| [`exportMetabolicModel`](https://compbtbs.github.io/COBRAxy/#/tools/export-metabolic-model) | Export tabular data to model format | Tabular model data | SBML/JSON/MAT/YAML model |
| [`ras_generator`](https://compbtbs.github.io/COBRAxy/#/tools/ras-generator) | Compute reaction activity scores | Gene expression data | RAS values |
| [`rps_generator`](https://compbtbs.github.io/COBRAxy/#/tools/rps-generator) | Compute reaction propensity scores | Metabolite abundance | RPS values |
| [`marea`](https://compbtbs.github.io/COBRAxy/#/tools/marea) | Statistical pathway analysis | RAS + RPS data | Enrichment + base maps |
| [`ras_to_bounds`](https://compbtbs.github.io/COBRAxy/#/tools/ras-to-bounds) | Apply RAS constraints to model | RAS + SBML model | Constrained bounds |
| [`flux_simulation`](https://compbtbs.github.io/COBRAxy/#/tools/flux-simulation) | Sample metabolic fluxes | Constrained model | Flux distributions |
| [`flux_to_map`](https://compbtbs.github.io/COBRAxy/#/tools/flux-to-map) | Add fluxes to enriched maps | Flux samples + base maps | Final styled maps |
| [`marea_cluster`](https://compbtbs.github.io/COBRAxy/#/tools/marea-cluster) | Cluster analysis | Expression/flux data | Sample clusters |


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

- **Models**: ENGRO2 (recommended), Recon (comprehensive) - see [Built-in Models Reference](https://compbtbs.github.io/COBRAxy/#/reference/built-in-models)
- **Gene mappings**: HGNC, Ensembl, Entrez ID conversions
- **Pathway maps**: Pre-styled SVG templates for ENGRO2 
- **Medium compositions**: Standard growth conditions

Located in `src/local/` directory for immediate use.

## Command Line Usage

All tools support `--help` for detailed options. Key commands:

### Generate RAS/RPS scores
```bash
# From gene expression
ras_generator -in expression.tsv -ra ras_output.tsv -rs ENGRO2

# From metabolite data  
rps_generator -id metabolites.tsv -rp rps_output.tsv
```

### Flux sampling
```bash
flux_simulation -ms ENGRO2 -in bounds/*.tsv -a CBS -ns 1000 -idop results/
```

### Statistical analysis & visualization
```bash
marea -using_RAS true -input_data ras.tsv -choice_map ENGRO2 -gs true -idop maps/
```

## Galaxy Integration

COBRAxy provides Galaxy tool wrappers (`.xml` files) for web-based analysis:

- Upload data through Galaxy interface
- Chain tools in visual workflows  
- Share and reproduce analyses

For Galaxy installation and setup, refer to the [official Galaxy documentation](https://docs.galaxyproject.org/).

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
