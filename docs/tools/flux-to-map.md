# Flux to Map

Visualize flux distributions on metabolic pathway maps.

## Overview

Flux to Map performs statistical comparison of flux distributions and visualizes results on SVG metabolic maps.

## Galaxy Interface

In Galaxy: **COBRAxy → Flux to Map**

1. Upload flux data and sample class file
2. Select map and configure comparison type
3. Click **Execute**

## Usage

```bash
flux_to_map -input_data fluxes.csv \
            -input_class classes.csv \
            -choice_map ENGRO2 \
            -comparison manyvsmany \
            -pvalue 0.05 \
            -fdecoration true \
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
| Arrows | `-fdecoration` | Show fold-change arrows | false |
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

### With Fold-Change Arrows

```bash
flux_to_map -input_data fluxes.csv \
            -input_class classes.csv \
            -choice_map ENGRO2 \
            -comparison manyvsmany \
            -pvalue 0.01 \
            -fdr true \
            -fdecoration true \
            -idop results/
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "No matching reactions" | Verify reaction ID consistency |
| "Insufficient data" | Increase sample sizes |

## See Also

- [MAREA](marea.md)
- [Flux Simulation](flux-simulation.md)
- [Built-in Models](reference/built-in-models)
