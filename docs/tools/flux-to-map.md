# Flux to Map

Visualize metabolic flux data on pathway maps with statistical analysis and color coding.

## Overview

Flux to Map performs statistical analysis on flux distribution data and generates color-coded metabolic pathway maps. It compares flux values between sample groups and highlights significantly different reactions with appropriate colors and line weights.

## Usage

### Command Line

```bash
flux_to_map -td /path/to/COBRAxy \
            -input_data_fluxes flux_data.tsv \
            -input_class_fluxes sample_groups.tsv \
            -comparison manyvsmany \
            -test ks \
            -pv 0.05 \
            -fc 1.5 \
            -choice_map ENGRO2 \
            -generate_svg true \
            -generate_pdf true \
            -idop flux_maps/
```

### Galaxy Interface

Select "Flux to Map" from the COBRAxy tool suite and configure flux analysis and visualization parameters.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |

### Data Input Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Flux Data | `-idf, --input_data_fluxes` | Flux values TSV file | - |
| Flux Classes | `-icf, --input_class_fluxes` | Sample group labels for fluxes | - |
| Multiple Flux Files | `-idsf, --input_datas_fluxes` | Multiple flux datasets (space-separated) | - |
| Flux Names | `-naf, --names_fluxes` | Names for multiple flux datasets | - |
| Analysis Option | `-op, --option` | Analysis mode (datasets or dataset_class) | - |

### Statistical Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Comparison Type | `-co, --comparison` | Statistical comparison mode | manyvsmany |
| Statistical Test | `-te, --test` | Statistical test method | ks |
| P-Value Threshold | `-pv, --pValue` | Significance threshold | 0.1 |
| Adjusted P-values | `-adj, --adjusted` | Apply FDR correction | false |
| Fold Change | `-fc, --fChange` | Minimum fold change threshold | 1.5 |

### Visualization Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Map Choice | `-mc, --choice_map` | Built-in metabolic map | HMRcore |
| Custom Map | `-cm, --custom_map` | Path to custom SVG map | - |
| Generate SVG | `-gs, --generate_svg` | Create SVG output | true |
| Generate PDF | `-gp, --generate_pdf` | Create PDF output | true |
| Color Map | `-colorm, --color_map` | Color scheme (jet, viridis) | - |
| Output Directory | `-idop, --output_path` | Results directory | result/ |

### Advanced Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Output Log | `-ol, --out_log` | Log file path | - |
| Control Sample | `-on, --control` | Control group identifier | - |

## Input Formats

### Flux Data File

Tab-separated format with reactions as rows and samples as columns:

```
Reaction	Sample1	Sample2	Sample3	Control1	Control2
R00001	15.23	-8.45	22.1	12.8	14.2
R00002	0.0	12.67	-5.3	8.9	7.4
R00003	45.8	38.2	51.7	42.1	39.8
R00004	-12.4	-15.8	-9.2	-11.5	-13.1
```

### Sample Class File

Group assignment for statistical comparisons:

```
Sample	Class
Sample1	Treatment
Sample2	Treatment  
Sample3	Treatment
Control1	Control
Control2	Control
```

### Multiple Dataset Format

When using multiple flux files, provide space-separated paths and corresponding names:

```bash
-idsf "dataset1_flux.tsv dataset2_flux.tsv dataset3_flux.tsv"
-naf "Condition_A Condition_B Condition_C"
```

## Statistical Analysis

### Comparison Types

#### manyvsmany
Compare all possible group pairs:
- Treatment vs Control
- Condition_A vs Condition_B
- Condition_A vs Condition_C
- Condition_B vs Condition_C

#### onevsrest
Compare each group against all others combined:
- Treatment vs (Control + Other)
- Control vs (Treatment + Other)

#### onevsmany
Compare one reference group against each other group:
- Control vs Treatment
- Control vs Condition_A
- Control vs Condition_B

### Statistical Tests

| Test | Description | Best For |
|------|-------------|----------|
| `ks` | Kolmogorov-Smirnov | Non-parametric, distribution-free |
| `ttest_p` | Paired t-test | Related samples, normal distributions |
| `ttest_ind` | Independent t-test | Independent samples, normal distributions |
| `wilcoxon` | Wilcoxon signed-rank | Non-parametric paired comparisons |
| `mw` | Mann-Whitney U | Non-parametric independent comparisons |

### Significance Assessment

Reactions are considered significant when:
1. **P-value** ≤ specified threshold (default: 0.1)
2. **Fold change** ≥ specified threshold (default: 1.5)
3. **FDR correction** (if enabled) maintains significance

## Map Visualization

### Built-in Maps

#### HMRcore (Default)
- **Scope**: Core human metabolic network
- **Reactions**: ~300 essential reactions
- **Coverage**: Central carbon, amino acid, lipid metabolism
- **Use Case**: General overview, publication figures

#### ENGRO2  
- **Scope**: Extended human genome-scale reconstruction
- **Reactions**: ~2,000 reactions
- **Coverage**: Comprehensive metabolic network
- **Use Case**: Detailed analysis, specialized tissues

#### Custom Maps
User-provided SVG files with reaction elements:
```xml
<rect id="R00001" class="reaction" fill="gray" stroke="black"/>
<path id="R00002" class="reaction" fill="gray" stroke="black"/>
```

### Color Coding Scheme

#### Significance Colors
- **Red Gradient**: Significantly upregulated (positive fold change)
- **Blue Gradient**: Significantly downregulated (negative fold change)  
- **Gray**: Not statistically significant
- **White**: No data available

#### Visual Elements
- **Line Width**: Proportional to fold change magnitude
- **Color Intensity**: Proportional to statistical significance (-log10 p-value)
- **Transparency**: Indicates confidence level

### Color Maps

#### Jet (Default)
- High contrast color transitions
- Blue (low) → Green → Yellow → Red (high)
- Good for identifying extreme values

#### Viridis
- Perceptually uniform color scale
- Colorblind-friendly
- Purple (low) → Blue → Green → Yellow (high)

## Output Files

### Statistical Results
- `flux_statistics.tsv`: P-values, fold changes, test statistics for all reactions
- `significant_fluxes.tsv`: Only reactions meeting significance criteria
- `comparison_summary.txt`: Analysis parameters and summary statistics

### Visualizations
- `flux_map.svg`: Interactive color-coded pathway map
- `flux_map.pdf`: High-resolution PDF (if requested)  
- `flux_map.png`: Raster image (if requested)
- `legend.svg`: Color scale and statistical significance legend

### Analysis Files
- `fold_changes.tsv`: Detailed fold change calculations
- `group_statistics.tsv`: Per-group summary statistics
- `comparison_matrix.tsv`: Pairwise comparison results

## Examples

### Basic Flux Comparison

```bash
# Compare treatment vs control fluxes
flux_to_map -td /opt/COBRAxy \
            -idf treatment_vs_control_fluxes.tsv \
            -icf sample_groups.tsv \
            -co manyvsmany \
            -te ks \
            -pv 0.05 \
            -fc 2.0 \
            -mc HMRcore \
            -gs true \
            -gp true \
            -idop flux_comparison/
```

### Multiple Condition Analysis

```bash
# Compare multiple experimental conditions
flux_to_map -td /opt/COBRAxy \
            -idsf "cond1_flux.tsv cond2_flux.tsv cond3_flux.tsv" \
            -naf "Control Treatment1 Treatment2" \
            -co onevsrest \
            -te wilcoxon \
            -adj true \
            -pv 0.01 \
            -fc 1.8 \
            -mc ENGRO2 \
            -colorm viridis \
            -idop multi_condition_flux/
```

### Custom Map Visualization

```bash
# Use tissue-specific custom map
flux_to_map -td /opt/COBRAxy \
            -idf liver_flux_data.tsv \
            -icf liver_conditions.tsv \
            -co manyvsmany \
            -te ttest_ind \
            -pv 0.05 \
            -fc 1.5 \
            -cm maps/liver_specific_map.svg \
            -gs true \
            -gp true \
            -idop liver_flux_analysis/ \
            -ol liver_analysis.log
```

### High-Throughput Analysis

```bash
# Process multiple datasets with stringent criteria
flux_to_map -td /opt/COBRAxy \
            -idsf "exp1.tsv exp2.tsv exp3.tsv exp4.tsv" \
            -naf "Exp1 Exp2 Exp3 Exp4" \
            -co manyvsmany \
            -te ks \
            -adj true \
            -pv 0.001 \
            -fc 3.0 \
            -mc HMRcore \
            -colorm jet \
            -gs true \
            -gp true \
            -idop high_throughput_flux/
```

## Quality Control

### Data Validation

#### Pre-analysis Checks
- Verify flux value distributions (check for outliers)
- Ensure sample names match between data and class files
- Validate reaction coverage across samples
- Check for missing values and their patterns

#### Statistical Validation  
- Assess normality assumptions for parametric tests
- Verify adequate sample sizes per group (n≥3 recommended)
- Check variance homogeneity between groups
- Evaluate multiple testing burden

### Result Interpretation

#### Biological Validation
- Compare results with known pathway activities
- Check for pathway coherence (related reactions should cluster)
- Validate against literature or experimental evidence
- Assess metabolic network connectivity

#### Technical Validation
- Compare results across different statistical tests
- Check sensitivity to parameter changes
- Validate fold change calculations
- Verify map element correspondence

## Tips and Best Practices

### Data Preparation
- **Normalization**: Ensure consistent flux units across samples
- **Filtering**: Remove reactions with excessive missing values (>50%)
- **Outlier Detection**: Identify and handle extreme flux values
- **Batch Effects**: Account for technical variation between experiments

### Statistical Considerations
- Use FDR correction for multiple comparisons (`-adj true`)
- Choose appropriate statistical tests based on data distribution
- Consider effect size (fold change) alongside significance
- Validate results with independent datasets when possible

### Visualization Optimization
- Select appropriate color maps for your audience
- Use high fold change thresholds (>2.0) for cleaner maps
- Export both SVG (editable) and PDF (publication) formats
- Include comprehensive legends and annotations

### Performance Tips
- Use HMRcore for faster processing and clearer visualizations
- Reduce dataset size for initial exploratory analysis
- Process large datasets in batches if memory constrained
- Cache intermediate results for parameter optimization

## Integration Workflow

### Upstream Tools
- [Flux Simulation](flux-simulation.md) - Generate flux distributions for comparison
- [MAREA](marea.md) - Alternative analysis pathway for RAS/RPS data

### Downstream Analysis
- Export results to statistical software (R, Python) for advanced analysis
- Integrate with pathway databases (KEGG, Reactome)
- Combine with other omics data for systems-level insights

### Typical Pipeline

```bash
# 1. Generate flux samples from constrained models
flux_simulation -td /opt/COBRAxy -ms ENGRO2 -in bounds/*.tsv \
                -ni Sample1,Sample2,Control1,Control2 -a CBS \
                -ot mean -idop fluxes/

# 2. Analyze and visualize flux differences
flux_to_map -td /opt/COBRAxy -idf fluxes/mean.csv \
            -icf sample_groups.tsv -co manyvsmany -te ks \
            -mc HMRcore -gs true -gp true -idop flux_maps/

# 3. Further analysis with custom scripts
python analyze_flux_results.py -i flux_maps/ -o final_results/
```

## Troubleshooting

### Common Issues

**No significant reactions found**
- Lower p-value threshold (`-pv 0.2`)
- Reduce fold change requirement (`-fc 1.2`)  
- Check sample group definitions and sizes
- Verify flux data quality and normalization

**Map rendering problems**
- Check SVG map file integrity and format
- Verify reaction ID matching between data and map
- Ensure sufficient system memory for large maps
- Validate XML structure of custom maps

**Statistical test failures**
- Check data distribution assumptions
- Verify sufficient sample sizes per group
- Consider alternative non-parametric tests
- Examine variance patterns between groups

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Map file not found" | Missing/invalid map path | Check file location and format |
| "No matching reactions" | ID mismatch between data and map | Verify reaction naming consistency |
| "Insufficient data" | Too few samples per group | Increase sample sizes or merge groups |
| "Memory allocation failed" | Large dataset/map combination | Reduce data size or increase system memory |

### Performance Issues

**Slow processing**
- Use HMRcore instead of ENGRO2 for faster rendering
- Reduce dataset size for testing
- Process subsets of reactions separately
- Monitor system resource usage

**Large output files**
- Use compressed formats when possible
- Reduce map resolution for preliminary analysis
- Export only essential output formats
- Clean temporary files regularly

## Advanced Usage

### Custom Statistical Functions

Advanced users can implement custom statistical tests by modifying the analysis functions:

```python
def custom_test(group1, group2):
    # Custom statistical test implementation
    statistic, pvalue = your_test_function(group1, group2)
    return statistic, pvalue
```

### Batch Processing Script

Process multiple experiments systematically:

```bash
#!/bin/bash
experiments=("exp1" "exp2" "exp3" "exp4")
for exp in "${experiments[@]}"; do
    flux_to_map -td /opt/COBRAxy \
                -idf "data/${exp}_flux.tsv" \
                -icf "data/${exp}_classes.tsv" \
                -co manyvsmany -te ks -pv 0.05 \
                -mc HMRcore -gs true -gp true \
                -idop "results/${exp}/"
done
```

### Result Aggregation

Combine results across multiple analyses:

```bash
# Merge significant reactions across experiments
python merge_flux_results.py \
    -i results/exp*/significant_fluxes.tsv \
    -o combined_significant_reactions.tsv \
    --method intersection
```

## See Also

- [Flux Simulation](flux-simulation.md) - Generate input flux distributions
- [MAREA](marea.md) - Alternative pathway analysis approach
- [Custom Map Creation Guide](/tutorials/custom-map-creation.md)
- [Statistical Methods Reference](/tutorials/statistical-methods.md)