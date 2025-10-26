# MAREA

Metabolic Reaction Enrichment Analysis for statistical comparison and map visualization.

## Overview

MAREA performs statistical enrichment analysis on RAS/RPS data to identify significantly different metabolic reactions between sample groups. It generates enriched pathway maps with color-coded reactions showing statistical significance and fold changes.

## Usage

### Command Line

```bash
marea -td /path/to/COBRAxy \
      -using_RAS true \
      -input_data ras_data.tsv \
      -input_class class_labels.tsv \
      -comparison manyvsmany \
      -test ks \
      -pv 0.05 \
      -fc 1.5 \
      -choice_map ENGRO2 \
      -generate_svg true \
      -idop results/
```

### Galaxy Interface

Select "MAREA" from the COBRAxy tool suite and configure analysis parameters through the web interface.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |

### Data Input Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Use RAS | `-using_RAS, --using_RAS` | Include RAS analysis | true |
| RAS Data | `-input_data, --input_data` | RAS scores TSV file | - |
| RAS Classes | `-input_class, --input_class` | Sample group labels | - |
| Multiple RAS | `-input_datas, --input_datas` | Multiple RAS files (space-separated) | - |
| RAS Names | `-names, --names` | Names for multiple datasets | - |
| Use RPS | `-using_RPS, --using_RPS` | Include RPS analysis | false |
| RPS Data | `-input_data_rps, --input_data_rps` | RPS scores TSV file | - |
| RPS Classes | `-input_class_rps, --input_class_rps` | RPS sample groups | - |

### Statistical Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Comparison Type | `-co, --comparison` | Statistical comparison mode | manyvsmany |
| Statistical Test | `-te, --test` | Statistical test method | ks |
| P-Value Threshold | `-pv, --pValue` | Significance threshold | 0.1 |
| Adjusted P-values | `-adj, --adjusted` | Apply FDR correction | false |
| Fold Change | `-fc, --fChange` | Minimum fold change | 1.5 |
| Net Enrichment | `-ne, --net` | Use net enrichment for RPS | false |
| Analysis Option | `-op, --option` | Analysis mode | datasets |

### Visualization Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Map Choice | `-choice_map, --choice_map` | Built-in metabolic map | - |
| Custom Map | `-custom_map, --custom_map` | Path to custom SVG map | - |
| Generate SVG | `-gs, --generate_svg` | Create SVG output | true |
| Generate PDF | `-gp, --generate_pdf` | Create PDF output | false |
| Generate PNG | `-gpng, --generate_png` | Create PNG output | false |
| Color Map | `-colorm, --color_map` | Color scheme (jet/viridis) | jet |
| Output Directory | `-idop, --output_path` | Results directory | result/ |

### Advanced Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Output Log | `-ol, --out_log` | Log file path | - |
| Control Sample | `-on, --control` | Control group identifier | - |

## Input Formats

### RAS/RPS Data File

Tab-separated format with reactions as rows and samples as columns:

```
Reaction	Sample1	Sample2	Sample3	Sample4
R00001	1.25	0.85	1.42	0.78
R00002	0.65	1.35	0.72	1.28
R00003	2.15	2.05	0.45	0.52
```

### Class Labels File

Sample grouping information:

```
Sample	Class
Sample1	Control
Sample2	Treatment
Sample3	Control
Sample4	Treatment
```

## Comparison Types

### manyvsmany
Compare all possible pairs of groups:
- Group A vs Group B
- Group A vs Group C  
- Group B vs Group C

### onevsrest
Compare each group against all others combined:
- Group A vs (Group B + Group C)
- Group B vs (Group A + Group C)

### onevsmany
Compare one specific group against each other group separately:
- Control vs Treatment1
- Control vs Treatment2

## Statistical Tests

| Test | Description | Use Case |
|------|-------------|----------|
| `ks` | Kolmogorov-Smirnov | Non-parametric, distribution-free |
| `ttest_p` | Paired t-test | Related samples |
| `ttest_ind` | Independent t-test | Unrelated samples |
| `wilcoxon` | Wilcoxon signed-rank | Non-parametric paired |
| `mw` | Mann-Whitney U | Non-parametric independent |
| `DESeq` | DESeq2-style analysis | Count-like data with dispersion |

## Output Files

### Statistical Results
- `comparison_stats.tsv`: P-values, fold changes, and test statistics
- `enriched_reactions.tsv`: Significantly enriched reactions only
- `comparison_summary.txt`: Analysis summary and parameters

### Visualization
- `pathway_map.svg`: Color-coded metabolic map
- `pathway_map.pdf`: PDF version (if requested)
- `pathway_map.png`: PNG version (if requested)
- `legend.svg`: Color scale and significance indicators

## Examples

### Basic RAS Analysis

```bash
# Simple two-group comparison
marea -td /opt/COBRAxy \
      -using_RAS true \
      -input_data ras_scores.tsv \
      -input_class sample_groups.tsv \
      -comparison manyvsmany \
      -test ks \
      -pv 0.05 \
      -choice_map ENGRO2 \
      -idop results/
```

### Combined RAS + RPS Analysis

```bash
# Multi-modal analysis
marea -td /opt/COBRAxy \
      -using_RAS true \
      -input_data ras_scores.tsv \
      -input_class ras_groups.tsv \
      -using_RPS true \
      -input_data_rps rps_scores.tsv \
      -input_class_rps rps_groups.tsv \
      -comparison onevsrest \
      -test DESeq \
      -adj true \
      -fc 2.0 \
      -choice_map HMRcore \
      -generate_pdf true \
      -idop multimodal_results/
```

### Multiple Dataset Analysis

```bash
# Compare multiple experiments
marea -td /opt/COBRAxy \
      -using_RAS true \
      -input_datas exp1_ras.tsv exp2_ras.tsv exp3_ras.tsv \
      -names "Experiment1" "Experiment2" "Experiment3" \
      -comparison onevsmany \
      -test wilcoxon \
      -pv 0.01 \
      -custom_map custom_pathway.svg \
      -idop multi_experiment/
```

## Map Visualization

### Built-in Maps
- **ENGRO2**: Human genome-scale reconstruction
- **HMRcore**: Core human metabolic network  
- **Recon**: Recon3D human model

### Color Coding
- **Red**: Significantly upregulated (high activity)
- **Blue**: Significantly downregulated (low activity)
- **Gray**: Not significant
- **Line Width**: Proportional to fold change magnitude

### Custom Maps
SVG files with reaction elements having IDs matching your data:
```xml
<rect id="R00001" class="reaction" .../>
<path id="R00002" class="reaction" .../>
```

## Quality Control

### Pre-analysis Checks
- Verify sample names match between data and class files
- Check for missing values and outliers
- Ensure adequate sample sizes per group (n ≥ 3 recommended)

### Post-analysis Validation
- Review statistical distribution assumptions
- Check multiple testing correction effects
- Validate biological relevance of enriched pathways

## Tips and Best Practices

### Statistical Considerations
- Use FDR correction (`-adj true`) for multiple comparisons
- Choose appropriate tests based on data distribution
- Consider effect size alongside significance

### Visualization Optimization
- Use high fold change thresholds (>2.0) for cleaner maps
- Export both SVG (editable) and PDF (publication-ready) formats
- Adjust color schemes for colorblind accessibility

## Troubleshooting

### Common Issues

**No significant reactions found**
- Lower p-value threshold (`-pv 0.1`)
- Reduce fold change requirement (`-fc 1.2`)
- Check sample group definitions

**Map rendering errors**
- Verify SVG map file integrity
- Check reaction ID matching between data and map
- Ensure sufficient system memory for large maps

**Statistical test failures**
- Verify data normality for parametric tests
- Check for sufficient sample sizes
- Consider alternative test methods

## Integration

### Upstream Tools
- [RAS Generator](ras-generator.md) - Generate RAS scores
- [RPS Generator](rps-generator.md) - Generate RPS scores

### Downstream Analysis
- [Flux to Map](flux-to-map.md) - Additional visualization options
- [MAREA Cluster](marea-cluster.md) - Sample clustering analysis

## See Also

- [Statistical Tests Documentation](/tutorials/statistical-tests.md)
- [Map Customization Guide](/tutorials/custom-maps.md)
- [Multi-modal Analysis Tutorial](/tutorials/multimodal-analysis.md)