# Flux to Map

Visualize flux distributions on metabolic pathway maps.

## Overview

This tool analyzes and visualizes statistical differences in reaction fluxes of groups of samples, returned by the Flux Simulation tool. The results can be visualized on s SVG metabolic map.

## Galaxy Interface

In Galaxy: **COBRAxy → Metabolic Flux Enrichment Analysis**

1. Upload flux data and sample class file
2. Select the map and configure the comparison type
3. Click **Run tool**

## Command-line console

```bash
flux_to_map -input_data fluxes.csv \
            -input_class classes.csv \
            -choice_map ENGRO2 \
            -comparison manyvsmany \
            -pvalue 0.05 \
            -idop output/
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Input Data | `-input_data` | Flux data file | - |
| Input Class | `-input_class` | Sample class definitions | - |
| Map Choice | `-choice_map` | ENGRO2, Recon, or Custom | ENGRO2 |
| Custom Map | `-custom_map` | Path to custom SVG map | - |
| Comparison | `-comparison` | manyvsmany, onevsrest, onevsmany | manyvsmany |
| P-value | `-pvalue` | Significance threshold | 0.05 |
| FDR Correction | `-fdr` | Apply FDR correction | true |
| Test Type | `-test_type` | t, wilcoxon, ks | t |
| Color Map | `--color_map` | Color scheme: viridis or jet | viridis |
| Output Path | `-idop` | Output directory | flux_to_map/ |

## Input Formats

### Flux Data

```
Reaction	Sample1	Sample2	Sample3
R00001	12.5	8.5	14.2
R00002	-6.5	13.5	7.2
```

### Sample Classes

```
SampleID	Class
Sample1	Control
Sample2	Treatment
Sample3	Treatment
```

**Note on Metabolic Map**
We provide a default SVG map for the ENGRO2 model. If another model is used, we suggest uploading a custom SVG map.

**File Format Notes:**
- Use **tab-separated** values (TSV) or **comma-separated** (CSV)
- First row must contain column headers
- Sample names must match between flux data and class file
- Class names should not contain spaces

## Statistical Tests

- **t**: Student's t-test (parametric, assumes normality)
- **wilcoxon**: Wilcoxon/Mann-Whitney (non-parametric)
- **ks**: Kolmogorov-Smirnov (distribution-free)

## Comparison Types

- **manyvsmany**: All pairwise class comparisons
- **onevsrest**: Each class vs all others
- **onevsmany**: One reference vs multiple classes

## Output

- `*_map.svg`: Annotated pathway maps
- `comparison_results.tsv`: Statistical results
- `*.log`: Processing log

## Examples

### Basic Comparison

```bash
flux_to_map -input_data fluxes.csv \
            -input_class classes.csv \
            -choice_map ENGRO2 \
            -comparison manyvsmany \
            -pvalue 0.05 \
            -idop results/
```

### With Custom Map

```bash
flux_to_map -input_data fluxes.csv \
            -input_class classes.csv \
            -choice_map Custom \
            -custom_map pathway.svg \
            -comparison onevsrest \
            -test_type wilcoxon \
            -idop results/
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "No matching reactions" | Verify reaction ID consistency |
| "Insufficient data" | Increase sample sizes |

## See Also

- [MAREA](tools/marea)
- [Flux Simulation](tools/flux-simulation)
- [Built-in Models](reference/built-in-models)
