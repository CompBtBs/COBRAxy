# Tools Reference

Complete documentation for all COBRAxy tools with parameters, examples, and usage guidelines.

## Tool Overview

| Tool | Category | Purpose | Input | Output |
|------|----------|---------|--------|--------|
| [RAS Generator](ras-generator.md) | Core | Compute reaction activity scores | Gene expression + GPR rules | RAS values |
| [RPS Generator](rps-generator.md) | Core | Compute reaction propensity scores | Metabolite abundance | RPS values |
| [MAREA](marea.md) | Analysis | Statistical pathway enrichment | RAS/RPS data | Enriched maps + statistics |
| [RAS to Bounds](ras-to-bounds.md) | Modeling | Apply RAS constraints to model | RAS + SBML model | Constrained bounds |
| [Flux Simulation](flux-simulation.md) | Modeling | Sample metabolic fluxes | Constrained model | Flux distributions |
| [Flux to Map](flux-to-map.md) | Visualization | Add flux data to maps | Flux samples + maps | Final visualizations |
| [Model Setting](metabolic-model-setting.md) | Utility | Extract model components | SBML model | Rules, reactions, bounds |
| [MAREA Cluster](marea-cluster.md) | Analysis | Cluster analysis | Expression/flux data | Sample clusters |

## Tool Categories

### Core Tools
Generate activity/propensity scores from omics data:
- **[RAS Generator](ras-generator.md)** - Gene expression → Reaction activity
- **[RPS Generator](rps-generator.md)** - Metabolites → Reaction propensity

### Analysis Tools  
Perform statistical analysis and enrichment:
- **[MAREA](marea.md)** - Pathway enrichment and visualization
- **[MAREA Cluster](marea-cluster.md)** - Sample clustering and classification

### Modeling Tools
Constraint-based flux analysis:
- **[RAS to Bounds](ras-to-bounds.md)** - Apply activity constraints
- **[Flux Simulation](flux-simulation.md)** - Sample flux distributions

### Visualization Tools
Create publication-ready figures:
- **[Flux to Map](flux-to-map.md)** - Combine flux data with pathway maps

### Utility Tools
Model manipulation and setup:
- **[Model Setting](metabolic-model-setting.md)** - Extract model information

## Common Parameters

Many tools share these parameters:

### Input/Output Parameters
- `-td, --tool_dir`: COBRAxy installation directory (required for all tools)
- `-in, --input`: Input dataset file
- `-idop, --output_dir`: Output directory for results

### Model Selection  
- `-rs, --rules_selector`: Choose built-in model (ENGRO2, Recon, HMRcore)
- `-ms, --model_selector`: Model for flux analysis

### Analysis Options
- `-gs, --gene_set_analysis`: Enable statistical enrichment testing
- `-a, --algorithm`: Sampling algorithm (CBS, OptGP)
- `-ns, --n_samples`: Number of flux samples to generate

## Usage Patterns

### Command Line Usage
```bash
# Basic pattern for all tools
tool_name -td $(pwd) [tool-specific options]

# Example: Generate RAS scores
ras_generator -td $(pwd) -in expression.tsv -ra ras_output.tsv -rs ENGRO2
```

### Python API Usage
```python
import tool_module

# All tools accept argument lists
args = ['-td', '/path/to/cobraxy', '-in', 'input.tsv', '-out', 'output.tsv']
tool_module.main(args)
```

### Galaxy Integration
All tools include Galaxy XML wrappers for web-based usage through the Galaxy interface.

## Parameter Reference

### File Format Requirements

**Gene Expression Files**
- Format: TSV (tab-separated values)
- Structure: Genes (rows) × Samples (columns)
- First column: Gene IDs (HGNC, Ensembl, etc.)
- Subsequent columns: Expression values

**Metabolite Files**
- Format: TSV (tab-separated values)  
- Structure: Metabolites (rows) × Samples (columns)
- First column: Metabolite names
- Subsequent columns: Abundance values

**Model Files**
- Format: SBML (.xml) or tabular rules (.tsv)
- Content: Metabolic network with GPR rules

### Built-in Models

| Model | Organism | Reactions | Genes | Best For |
|-------|----------|-----------|-------|----------|
| **ENGRO2** | Human | ~2,000 | ~500 | Focused analysis, faster computation |
| **Recon** | Human | ~10,000 | ~2,000 | Comprehensive metabolism |
| **HMRcore** | Human | ~5,000 | ~1,000 | Balanced coverage |

## Tool Selection Guide

### Choose Your Analysis Path

**For Pathway Enrichment**
1. [RAS Generator](ras-generator.md) → Generate activity scores
2. [RPS Generator](rps-generator.md) → Generate propensity scores (optional)
3. [MAREA](marea.md) → Statistical analysis and visualization

**For Flux Analysis**  
1. [RAS Generator](ras-generator.md) → Generate activity scores
2. [RAS to Bounds](ras-to-bounds.md) → Apply constraints
3. [Flux Simulation](flux-simulation.md) → Sample fluxes
4. [Flux to Map](flux-to-map.md) → Create visualizations

**For Model Exploration**
1. [Model Setting](metabolic-model-setting.md) → Extract model info
2. Analyze model structure and gene coverage

**For Sample Classification**
1. Generate RAS/RPS scores
2. [MAREA Cluster](marea-cluster.md) → Cluster samples

## Performance Considerations

### Memory Requirements
- **Small datasets** (< 1000 genes): 2-4 GB RAM
- **Medium datasets** (1000-5000 genes): 4-8 GB RAM  
- **Large datasets** (> 5000 genes): 8+ GB RAM
- **Flux sampling**: Additional 2-4 GB per 1000 samples

### Speed Optimization
- Use **ENGRO2** model for faster analysis
- Install **GLPK** solver for CBS sampling performance
- Reduce **sampling count** for initial testing
- Use **SSD storage** for temporary files

### Batch Processing
- Process multiple datasets in parallel
- Use appropriate **chunk sizes** for large sample sets
- Monitor **memory usage** during batch runs

## Troubleshooting

### Common Issues Across Tools

**File Path Problems**
- Use absolute paths when possible
- Ensure all input files exist before starting
- Check write permissions for output directories

**Memory Errors**
- Reduce dataset size or sampling parameters
- Close other applications to free memory
- Use more powerful hardware for large analyses

**Model Issues**
- Verify model file format and gene ID consistency
- Check gene ID mapping between data and model
- Use built-in models to avoid compatibility issues

### Getting Help

For tool-specific issues:
1. Check individual tool documentation
2. Review parameter requirements and formats
3. Test with smaller datasets first
4. Consult [troubleshooting guide](../troubleshooting.md)

## Contributing

Help improve tool documentation:
- Report unclear instructions
- Suggest additional examples
- Contribute usage patterns
- Fix documentation errors

Each tool page includes detailed parameter descriptions, examples, and troubleshooting tips. Select a tool from the sidebar to get started!