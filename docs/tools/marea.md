# MAREA

Metabolic Enrichment Analysis and Visualization.

## Overview

MAREA performs statistical comparison of metabolic scores (RAS/RPS) and visualizes results on pathway maps.

## Galaxy Interface

In Galaxy: **COBRAxy → MAREA**

1. Upload RAS/RPS scores and sample class file
2. Select map and configure statistical parameters
3. Click **Execute**

## Usage

```bash
marea -input_data scores.tsv \
      -input_class classes.csv \
      -choice_map ENGRO2 \
      -comparison manyvsmany \
      -pvalue 0.05 \
      -idop output/
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Input Data | `-input_data` | RAS/RPS scores file | - |
| Input Class | `-input_class` | Sample class definitions | - |
| Map Choice | `-choice_map` | ENGRO2, Recon, or Custom | ENGRO2 |
| Custom Map | `-custom_map` | Path to custom SVG map | - |
| Comparison | `-comparison` | manyvsmany, onevsrest, onevsmany | manyvsmany |
| P-value | `-pvalue` | Significance threshold | 0.05 |
| FDR Correction | `-fdr` | Apply FDR correction | true |
| Test Type | `-test_type` | t, wilcoxon, ks | t |
| Output Path | `-idop` | Output directory | marea/ |

## Input Formats

### Metabolic Scores

```
Reaction	Sample1	Sample2	Sample3
R00001	1.25	0.85	1.42
R00002	0.65	1.35	0.72
```

### Sample Classes

```
SampleID	Class
Sample1	Control
Sample2	Treatment
Sample3	Treatment
```

## Statistical Tests

- **t**: Student's t-test (parametric)
- **wilcoxon**: Wilcoxon/Mann-Whitney (non-parametric)
- **ks**: Kolmogorov-Smirnov (distribution-free)

## Comparison Types

- **manyvsmany**: All pairwise comparisons
- **onevsrest**: Each class vs all others
- **onevsmany**: One reference vs multiple classes

## Output

- `*_map.svg`: Annotated pathway maps
- `comparison_results.tsv`: Statistical results
- `*.log`: Processing log

## Examples

### Basic Analysis

```bash
marea -input_data ras_scores.tsv \
      -input_class classes.csv \
      -choice_map ENGRO2 \
      -comparison manyvsmany \
      -pvalue 0.05 \
      -idop results/
```

### Custom Map

```bash
marea -input_data rps_scores.tsv \
      -input_class classes.csv \
      -choice_map Custom \
      -custom_map pathway.svg \
      -comparison onevsrest \
      -idop results/
```

### Non-parametric Test

```bash
marea -input_data scores.tsv \
      -input_class classes.csv \
      -choice_map ENGRO2 \
      -test_type wilcoxon \
      -pvalue 0.01 \
      -fdr true \
      -idop results/
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "No matching reactions" | Verify reaction IDs |
| "Insufficient samples" | Increase sample sizes per group |

## See Also

- [RAS Generator](ras-generator.md)
- [RPS Generator](rps-generator.md)
- [Flux to Map](flux-to-map.md)
- [Built-in Models](../reference/built-in-models.md)
