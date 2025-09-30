# RPS Generator

Generate Reaction Propensity Scores (RPS) from metabolite abundance data.

## Overview

The RPS Generator computes reaction propensity scores based on metabolite abundance measurements. RPS values indicate how likely metabolic reactions are to be active based on the availability of their substrate and product metabolites.

## Usage

### Command Line

```bash
rps_generator -td /path/to/COBRAxy \
              -id metabolite_abundance.tsv \
              -rp output_rps.tsv \
              -ol log.txt
```

### Galaxy Interface

Select "RPS Generator" from the COBRAxy tool suite and upload your metabolite abundance file.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |
| Input Dataset | `-id, --input` | Metabolite abundance TSV file (rows=metabolites, cols=samples) |
| RPS Output | `-rp, --rps_output` | Output file path for RPS scores |

### Optional Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Custom Reactions | `-rl, --model_upload` | Path to custom reactions file | Built-in reactions |
| Output Log | `-ol, --out_log` | Log file for warnings/errors | Standard output |

## Input Format

### Metabolite Abundance File

Tab-separated values (TSV) format:

```
Metabolite	Sample1	Sample2	Sample3
glucose	100.5	85.2	92.7
pyruvate	45.3	38.9	41.2
lactate	15.8	22.1	18.5
```

**Requirements:**
- First column: metabolite names (case-insensitive)
- Subsequent columns: abundance values for each sample
- Missing values: use 0 or leave empty
- File encoding: UTF-8

### Custom Reactions File (Optional)

If using custom reactions instead of built-in ones:

```
ReactionID	Reaction
R00001	glucose + ATP -> glucose-6-phosphate + ADP
R00002	glucose-6-phosphate <-> fructose-6-phosphate
```

## Output Format

### RPS Scores File

```
Reaction	Sample1	Sample2	Sample3
R00001	0.85	0.72	0.79
R00002	0.45	0.38	0.52
R00003	0.12	0.28	0.21
```

- Values range from 0 (low propensity) to 1 (high propensity)
- NaN values indicate insufficient metabolite data for that reaction

## Algorithm

1. **Metabolite Matching**: Input metabolite names are matched against internal synonyms
2. **Abundance Normalization**: Raw abundances are normalized per sample
3. **Reaction Scoring**: For each reaction, RPS is computed based on:
   - Substrate availability (geometric mean of substrate abundances)
   - Product formation potential
   - Stoichiometric coefficients

## Examples

### Basic Usage

```bash
# Generate RPS from metabolite data
rps_generator -td /opt/COBRAxy \
              -id /data/metabolomics.tsv \
              -rp /results/rps_scores.tsv
```

### With Custom Reactions

```bash
# Use custom reaction set
rps_generator -td /opt/COBRAxy \
              -id /data/metabolomics.tsv \
              -rl /custom/reactions.tsv \
              -rp /results/custom_rps.tsv \
              -ol /logs/rps.log
```

## Tips and Best Practices

### Data Preparation

- **Metabolite Names**: Use standard nomenclature (KEGG, ChEBI, or common names)
- **Missing Data**: Remove samples with >50% missing metabolites
- **Outliers**: Consider log-transformation for highly variable metabolites
- **Replicates**: Average technical replicates before analysis

### Quality Control

- Check log file for unmatched metabolite names
- Verify RPS score distributions (should span 0-1 range)
- Compare results with expected pathway activities

### Integration with Other Tools

RPS scores are typically used with:
- [MAREA](marea.md) for pathway enrichment analysis
- [Flux to Map](flux-to-map.md) for metabolic map visualization

## Troubleshooting

### Common Issues

**No RPS scores generated**
- Check metabolite name format and spelling
- Verify input file has correct TSV format
- Ensure tool directory contains reaction databases

**Many NaN values in output**
- Insufficient metabolite coverage for reactions
- Consider using a smaller, more focused reaction set

**Memory errors**
- Reduce dataset size or split into batches
- Increase available system memory

### Error Messages

| Error | Cause | Solution |
|-------|--------|----------|
| "File not found" | Missing input file | Check file path and permissions |
| "Invalid format" | Malformed TSV | Verify column headers and data types |
| "No metabolites matched" | Name mismatch | Check metabolite nomenclature |

## See Also

- [RAS Generator](ras-generator.md) - Generate reaction activity scores from gene expression
- [MAREA](marea.md) - Statistical analysis and visualization
- [Flux Simulation](flux-simulation.md) - Constraint-based modeling