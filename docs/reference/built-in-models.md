# Built-in Models

COBRAxy includes two pre-installed metabolic models for human metabolism analysis.

## ENGRO2 (Recommended)

**Best for**: General metabolic analysis

- ~2,000 reactions, ~1,500 metabolites, ~500 genes
- Balanced coverage with fast computation
- Core metabolic pathways well-represented
- **Use for**: Tissue profiling, disease comparisons, time-series analysis

## Recon (Comprehensive)

**Best for**: Genome-wide studies

- ~10,000 reactions, ~5,000 metabolites, ~2,000 genes
- Most complete human metabolic network
- Includes rare and specialized pathways
- **Use for**: Comprehensive studies, rare diseases (requires >16 GB RAM)

## Usage

```bash
# Specify model name in any COBRAxy tool
tool_name --model ENGRO2 [options]
```

## Technical Notes

**Gene IDs**: All models support HGNC ID (e.g., `HGNC:5`), HGNC Symbol (e.g., `ALDOA`), Ensembl, and Entrez formats.

**GPR Rules**: Gene-Protein-Reaction associations use Boolean logic (`and`, `or`).

**Custom Models**: Use [Import Metabolic Model](../tools/import-metabolic-model.md) to extract and [Export Metabolic Model](../tools/export-metabolic-model.md) to create custom models.

## See Also

- [RAS Generator](../tools/ras-generator.md) - Uses model GPR rules
- [RAS to Bounds](../tools/ras-to-bounds.md) - Applies flux constraints to models
- [Flux Simulation](../tools/flux-simulation.md) - Samples model flux distributions
- [Import Metabolic Model](../tools/import-metabolic-model.md) - Extract model data
- [Export Metabolic Model](../tools/export-metabolic-model.md) - Create custom models
