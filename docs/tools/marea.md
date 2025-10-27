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
| Test Type | `-test_type` | t, wilcoxon, ks, DESeq | t |
| Net RPS | `--net` | Use net contribution for reversible reactions (RPS only) | false |
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

**File Format Notes:**
- Use **tab-separated** values (TSV) or **comma-separated** (CSV)
- First row must contain column headers
- Sample names must match between scores and class file
- Class names should not contain spaces

## Statistical Tests

- **t**: Student's t-test (parametric, assumes normality)
- **wilcoxon**: Wilcoxon/Mann-Whitney (non-parametric)
- **ks**: Kolmogorov-Smirnov (distribution-free)
- **DESeq**: DESeq2-like test (**RAS only**, requires ≥2 replicates per sample)

**Note on DESeq**: The DESeq2-like statistical test is specifically designed for RAS data and implements variance modeling similar to RNA-seq differential expression analysis. It requires at least 2 biological replicates per sample group and should NOT be used with RPS data.

## Comparison Types

- **manyvsmany**: All pairwise comparisons
- **onevsrest**: Each class vs all others
- **onevsmany**: One reference vs multiple classes

## Advanced Options

### Net RPS Values

When analyzing RPS data with reversible reactions, the `--net` parameter controls arrow coloring:

**When `--net false` (default):**
- Each direction of a reversible reaction colored independently
- Forward and backward contributions shown separately

**When `--net true` (RPS only):**
- Arrow tips colored with net contribution of both directions
- Useful for visualizing overall metabolite flow direction

**Note**: This option only applies to RPS analysis and affects visualization of reversible reactions on metabolic maps.

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

- [RAS Generator](tools/ras-generator)
- [RPS Generator](tools/rps-generator)
- [Flux to Map](tools/flux-to-map)
- [Built-in Models](reference/built-in-models)
