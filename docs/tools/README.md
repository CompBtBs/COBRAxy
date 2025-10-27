Complete reference for all COBRAxy tools with parameters, examples, and usage guidelines.

## Available Tools

| Galaxy Tool | Python script | Purpose | Input | Output |
|------|---------|---------|--------|--------|
| [Import Metabolic Model](tools/import-metabolic-model) | importMetabolicModel | Import and extract model components | SBML/JSON/MAT/YAML model | Tabular model data |
| [Export Metabolic Model](tools/export-metabolic-model) | exportMetabolicModel |Export tabular data to model format | Tabular model data | SBML/JSON/MAT/YAML model |
| [Expression2RAS](tools/ras-generator) | ras_generator | Compute reaction activity scores | Gene expression + GPR rules | RAS values |
| [Expression2RPS](tools/rps-generator) | rps_generator | Compute reaction propensity scores | Metabolite abundance | RPS values |
| [Metabolic Reaction Enrichment Analysis](tools/marea) | marea | Statistical pathway enrichment | RAS/RPS data | Enriched maps + statistics |
| [RAS2Bounds](tools/ras-to-bounds) | ras_to_bounds | Apply RAS constraints to model | RAS + SBML model | Constrained bounds |
| [Flux Simulation](tools/flux-simulation) | flux_simulation | Sample metabolic fluxes | Constrained model | Flux distributions |
| [Metabolic Flux Enrichment Analysis](tools/flux-to-map) | flux_to_map | Visualize flux data on maps | Flux samples + statistical comparison | Color-coded pathway maps |
| [Cluster Analysis](tools/marea-cluster) | marea_cluster | Cluster analysis | Expression/RAS/RPS/flux data | Sample clusters + validation plots |

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

**For Model Creation**
1. Create/edit tabular model data
2. [Export Metabolic Model](tools/export-metabolic-model) → Create SBML/JSON/MAT/YAML model

## Troubleshooting

### Common Issues Across Tools

**Model Issues**
- Verify model file format and gene ID consistency
- Check gene ID mapping between data and model

### Getting Help

For tool-specific issues:
1. Check individual tool documentation
2. Review parameter requirements and formats
3. Test with smaller datasets first
4. Consult [troubleshooting guide](troubleshooting)

## Contributing

Help improve tool documentation:
- Report unclear instructions
- Suggest additional examples
- Contribute usage patterns
- Fix documentation errors

Each tool page includes detailed parameter descriptions, examples, and troubleshooting tips.
