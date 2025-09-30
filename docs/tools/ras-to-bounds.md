# RAS to Bounds

Apply Reaction Activity Scores (RAS) as constraints to metabolic model bounds.

## Overview

The RAS to Bounds tool integrates RAS values into metabolic model flux bounds, creating sample-specific constrained models for flux sampling. This enables personalized metabolic modeling based on gene expression patterns.

## Usage

### Command Line

```bash
ras_to_bounds -td /path/to/COBRAxy \
              -ms ENGRO2 \
              -ir ras_scores.tsv \
              -rs true \
              -mes allOpen \
              -idop constrained_bounds/
```

### Galaxy Interface

Select "RAS to Bounds" from the COBRAxy tool suite and configure model and constraint parameters.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |
| Model Selector | `-ms, --model_selector` | Built-in model (ENGRO2, Custom) |
| RAS Selector | `-rs, --ras_selector` | Enable RAS constraint application |

### Model Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Model Selector | `-ms, --model_selector` | Built-in model choice | ENGRO2 |
| Custom Model | `-mo, --model` | Path to custom SBML model | - |
| Model Name | `-mn, --model_name` | Custom model filename | - |

### Medium Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Medium Selector | `-mes, --medium_selector` | Medium configuration | allOpen |
| Custom Medium | `-meo, --medium` | Path to custom medium file | - |

### Constraint Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| RAS Input | `-ir, --input_ras` | RAS scores TSV file | - |
| RAS Names | `-rn, --name` | Sample names for RAS data | - |
| Cell Class | `-cc, --cell_class` | Output cell class information | - |

### Output Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Output Path | `-idop, --output_path` | Directory for constrained bounds | ras_to_bounds/ |
| Output Log | `-ol, --out_log` | Log file path | - |

## Input Formats

### RAS Scores File

Tab-separated format with reactions as rows and samples as columns:

```
Reaction	Sample1	Sample2	Sample3	Control1	Control2
R00001	1.25	0.85	1.42	1.05	0.98
R00002	0.65	1.35	0.72	1.15	1.08
R00003	2.15	2.05	0.45	0.95	1.12
```

### Custom Model File (Optional)

SBML format metabolic model:
- XML format (.xml, .sbml)
- Compressed formats supported (.xml.gz, .xml.zip, .xml.bz2)
- Must contain valid reaction, metabolite, and gene definitions

### Custom Medium File (Optional)

Exchange reactions defining growth medium:

```
reaction
EX_glc__D_e
EX_o2_e
EX_pi_e
EX_nh4_e
```

## Algorithm

### Constraint Application

1. **Base Model Loading**: Load specified metabolic model and medium
2. **Bounds Extraction**: Extract original flux bounds for each reaction  
3. **RAS Integration**: For each sample and reaction:
   ```
   if RAS > 1.0:
       new_upper_bound = original_upper_bound * RAS
   if RAS < 1.0:  
       new_lower_bound = original_lower_bound * RAS
   ```
4. **Bounds Output**: Generate sample-specific bounds files

### Scaling Rules

- **RAS > 1**: Upregulated reactions → increased flux capacity
- **RAS < 1**: Downregulated reactions → decreased flux capacity  
- **RAS = 1**: No change from original bounds
- **Missing RAS**: Original bounds retained

## Output Format

### Bounds Files

One TSV file per sample with constrained bounds:

```
# bounds_Sample1.tsv
Reaction	lower_bound	upper_bound
R00001	-1000	1250.5
R00002	-650.2	1000  
R00003	-1000	2150.8
```

### Directory Structure

```
ras_to_bounds/
├── bounds_Sample1.tsv
├── bounds_Sample2.tsv  
├── bounds_Sample3.tsv
├── bounds_Control1.tsv
├── bounds_Control2.tsv
└── constraints_log.txt
```

## Examples

### Basic Usage with Built-in Model

```bash
# Apply RAS constraints to ENGRO2 model
ras_to_bounds -td /opt/COBRAxy \
              -ms ENGRO2 \
              -ir ras_data.tsv \
              -rs true \
              -idop results/bounds/
```

### Custom Model and Medium

```bash
# Use custom model with specific medium
ras_to_bounds -td /opt/COBRAxy \
              -ms Custom \
              -mo models/custom_model.xml \
              -mn custom_model.xml \
              -mes custom \
              -meo media/minimal_medium.tsv \
              -ir patient_ras.tsv \
              -rs true \
              -idop personalized_models/ \
              -ol constraints.log
```

### Multiple Sample Processing

```bash
# Process cohort data with sample classes
ras_to_bounds -td /opt/COBRAxy \
              -ms ENGRO2 \
              -ir cohort_ras_scores.tsv \
              -rn "Patient1,Patient2,Patient3,Healthy1,Healthy2" \
              -rs true \
              -cc sample_classes.tsv \
              -idop cohort_bounds/
```

## Built-in Models

### ENGRO2
- **Species**: Homo sapiens
- **Scope**: Genome-scale reconstruction  
- **Reactions**: ~2,000 reactions
- **Metabolites**: ~1,500 metabolites
- **Use Case**: General human metabolism

### Custom Model Requirements
- Valid SBML format
- Consistent reaction/metabolite naming
- Proper compartment definitions
- Gene-protein-reaction associations

## Medium Configurations

### allOpen (Default)
- All exchange reactions unconstrained
- Maximum metabolic flexibility
- Suitable for exploratory analysis

### Custom Medium
- User-defined nutrient availability
- Tissue-specific conditions
- Disease-specific constraints

## Quality Control

### Pre-processing Checks
- Verify RAS data completeness (recommend >80% reaction coverage)
- Check for extreme RAS values (>10 or <0.1 may indicate issues)
- Validate model consistency and solvability

### Post-processing Validation
- Confirm bounds files generated for all samples
- Check constraint log for warnings
- Test model feasibility with sample bounds

## Tips and Best Practices

### RAS Data Preparation
- **Normalization**: Ensure RAS values are properly normalized (median ~1.0)
- **Filtering**: Remove reactions with consistently missing data
- **Validation**: Check RAS distributions across samples

### Model Selection
- Use ENGRO2 for general human tissue analysis
- Consider custom models for specific organisms or tissues
- Validate model scope matches your biological question

### Medium Configuration  
- Match medium to experimental conditions
- Use minimal medium for growth requirement analysis
- Consider tissue-specific nutrient availability

## Integration Workflow

### Upstream Tools
- [RAS Generator](ras-generator.md) - Generate RAS scores from expression data

### Downstream Tools  
- [Flux Simulation](flux-simulation.md) - Sample fluxes using constrained bounds
- [MAREA](marea.md) - Statistical analysis of constraint effects

### Typical Pipeline

```bash
# 1. Generate RAS from expression data
ras_generator -td /opt/COBRAxy -in expression.tsv -ra ras.tsv

# 2. Apply RAS constraints to model bounds  
ras_to_bounds -td /opt/COBRAxy -ms ENGRO2 -ir ras.tsv -rs true -idop bounds/

# 3. Sample fluxes with constraints
flux_simulation -td /opt/COBRAxy -ms ENGRO2 -in bounds/*.tsv -a CBS -idop fluxes/

# 4. Analyze and visualize results
marea -td /opt/COBRAxy -input_data fluxes/mean.tsv -choice_map ENGRO2 -idop maps/
```

## Troubleshooting

### Common Issues

**No bounds files generated**
- Check RAS file format and sample names
- Verify model loading (check model path/format)
- Ensure sufficient disk space for output

**Model infeasibility after constraints**
- RAS values may be too restrictive
- Consider scaling factor adjustment
- Check medium compatibility with constraints

**Missing reactions in bounds**  
- RAS data may not cover all model reactions
- Original bounds retained for missing reactions
- Consider reaction mapping validation

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Model not found" | Invalid model path | Check model file location |
| "RAS file invalid" | Malformed TSV format | Verify file structure and encoding |
| "Infeasible solution" | Over-constrained model | Relax RAS scaling or medium constraints |

### Performance Issues

**Slow processing**
- Large models may require significant memory
- Consider batch processing for many samples
- Monitor system resource usage

**Memory errors**
- Reduce model size or split processing
- Increase available system memory
- Use more efficient file formats

## Advanced Usage

### Batch Processing Script

```bash
#!/bin/bash
# Process multiple RAS files
for ras_file in ras_data/*.tsv; do
    sample_name=$(basename "$ras_file" .tsv)
    ras_to_bounds -td /opt/COBRAxy \
                  -ms ENGRO2 \
                  -ir "$ras_file" \
                  -rs true \
                  -idop "bounds_$sample_name/"
done
```

### Custom Scaling Functions

For advanced users, RAS scaling can be customized by modifying the constraint application logic in the source code.

## See Also

- [RAS Generator](ras-generator.md) - Generate input RAS data
- [Flux Simulation](flux-simulation.md) - Use constrained bounds for sampling  
- [Model Setting](metabolic-model-setting.md) - Extract model components