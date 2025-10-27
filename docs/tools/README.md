Complete reference for all COBRAxy tools with parameters, examples, and usage guidelines.

## Available Tools

| Tool | Purpose | Input | Output |
|------|---------|--------|--------|
| [Import Metabolic Model](tools/import-metabolic-model) | Import and extract model components | SBML/JSON/MAT/YAML model | Tabular model data |
| [Export Metabolic Model](tools/export-metabolic-model) | Export tabular data to model format | Tabular model data | SBML/JSON/MAT/YAML model |
| [RAS Generator](tools/ras-generator) | Compute reaction activity scores | Gene expression + GPR rules | RAS values |
| [RPS Generator](tools/rps-generator) | Compute reaction propensity scores | Metabolite abundance | RPS values |
| [MAREA](tools/marea) | Statistical pathway enrichment | RAS/RPS data | Enriched maps + statistics |
| [RAS to Bounds](tools/ras-to-bounds) | Apply RAS constraints to model | RAS + SBML model | Constrained bounds |
| [Flux Simulation](tools/flux-simulation) | Sample metabolic fluxes | Constrained model | Flux distributions |
| [Flux to Map](tools/flux-to-map) | Visualize flux data on maps | Flux samples + statistical comparison | Color-coded pathway maps |
| [MAREA Cluster](tools/marea-cluster) | Cluster analysis | Expression/RAS/RPS/flux data | Sample clusters + validation plots |


## Analysis Workflows

**Enrichment Analysis**: Gene Expression → RAS Generator → MAREA → Pathway Maps

**Flux Simulation**: Gene Expression → RAS Generator → RAS to Bounds → Flux Simulation → Flux to Map

## Usage Patterns

### Galaxy Integration
All tools include Galaxy XML wrappers for web-based usage through the Galaxy interface.

### Command Line Usage
```bash
# Basic pattern for all tools
tool_name [tool-specific options]

# Example: Generate RAS scores
ras_generator -in expression.tsv -ra ras_output.tsv -rs ENGRO2
```

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

COBRAxy includes two pre-installed models. See [Built-in Models Reference](reference/built-in-models) for details.

| Model | Size | Best For |
|-------|------|----------|
| **ENGRO2** | ~2,000 reactions | General analysis (recommended) |
| **Recon** | ~10,000 reactions | Comprehensive studies |

## Tool Selection Guide

### Choose Your Analysis Path

**For Pathway Enrichment**
1. [RAS Generator](tools/ras-generator) → Generate activity scores
2. [RPS Generator](tools/rps-generator) → Generate propensity scores (optional)
3. [MAREA](tools/marea) → Statistical analysis and visualization

**For Flux Analysis**  
1. [RAS Generator](tools/ras-generator) → Generate activity scores
2. [RAS to Bounds](tools/ras-to-bounds) → Apply constraints
3. [Flux Simulation](tools/flux-simulation) → Sample fluxes
4. [Flux to Map](tools/flux-to-map) → Create visualizations

**For Model Exploration**
1. [Import Metabolic Model](tools/import-metabolic-model) → Extract model info
2. Analyze model structure and gene coverage

**For Model Creation**
1. Create/edit tabular model data
2. [Export Metabolic Model](tools/export-metabolic-model) → Create SBML/JSON/MAT/YAML model

**For Sample Classification**
1. Generate RAS/RPS scores
2. [MAREA Cluster](tools/marea-cluster) → Cluster samples



## Troubleshooting

### Common Issues Across Tools

**Model Issues**
- Verify model file format and gene ID consistency
- Check gene ID mapping between data and model
- Use built-in models to avoid compatibility issues

### Getting Help

For tool-specific issues:
1. Check individual tool documentation
2. Review parameter requirements and formats
3. Test with smaller datasets first
4. Consult [troubleshooting guide](/troubleshooting.md)

## Contributing

Help improve tool documentation:
- Report unclear instructions
- Suggest additional examples
- Contribute usage patterns
- Fix documentation errors

Each tool page includes detailed parameter descriptions, examples, and troubleshooting tips. Select a tool from the sidebar to get started!