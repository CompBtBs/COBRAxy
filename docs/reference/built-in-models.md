# Built-in Models

COBRAxy includes two pre-installed metabolic models for human metabolism analysis.

## ENGRO2 (Recommended)

**Best for**: Core metabolic analysis

- ~500 reactions, ~400 metabolites, ~500 genes
- Core metabolic model

## Recon (Comprehensive)

**Best for**: Genome-wide studies

- ~10,000 reactions, ~5,000 metabolites, ~2,000 genes
- Most complete human metabolic network, including all metabolic pathways

## Usage

```bash
# Specify model name in any COBRAxy tool
tool_name --model ENGRO2 [options]
```

## Technical Notes

**Gene IDs**: All models support HGNC ID (e.g., `HGNC:5`), HGNC Symbol (e.g., `ALDOA`), Ensembl, and Entrez formats.

**GPR Rules**: Gene-Protein-Reaction associations use Boolean logic (`and`, `or`).

**Custom Models**: Use [Import Metabolic Model](../tools/import-metabolic-model) to extract and [Export Metabolic Model](../tools/export-metabolic-model) to create custom models.

## See Also

- [RAS Generator](../tools/ras-generator) - Uses model GPR rules
- [RAS to Bounds](../tools/ras-to-bounds) - Applies flux constraints to models
- [Flux Simulation](../tools/flux-simulation) - Samples model flux distributions
- [Import Metabolic Model](../tools/import-metabolic-model) - Extract model data
- [Export Metabolic Model](../tools/export-metabolic-model) - Create custom models
