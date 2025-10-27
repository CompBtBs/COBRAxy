# RAS to Bounds

Apply Reaction Activity Scores (RAS) as constraints to metabolic model bounds.

## Overview

RAS to Bounds integrates RAS values into metabolic model flux bounds, creating sample-specific constrained models.

## Galaxy Interface

In Galaxy: **COBRAxy → RAS to Bounds**

1. Select model and upload RAS scores file
2. Configure medium and constraint options
3. Click **Execute**

## Usage

```bash
ras_to_bounds -ms ENGRO2 \
              -ir ras_scores.tsv \
              -rs true \
              -mes allOpen \
              -idop constrained_bounds/
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Model Selector | `-ms` | ENGRO2, Recon, or Custom | ENGRO2 |
| RAS Input | `-ir` | RAS scores TSV file | - |
| RAS Selector | `-rs` | Enable RAS constraints | false |
| Medium Selector | `-mes` | Medium configuration | allOpen |
| Save Models | `--save_models` | Save complete models with bounds | false |
| Output Path | `-idop` | Output directory | ras_to_bounds/ |
| Custom Model | `-mo` | Path to custom SBML model | - |
| Custom Medium | `-meo` | Custom medium file | - |

## Input Format

RAS scores file (TSV):

```
Reaction	Sample1	Sample2	Sample3
R00001	1.25	0.85	1.42
R00002	0.65	1.35	0.72
```

**File Format Notes:**
- Use **tab-separated** values (TSV)
- First row must contain column headers (Reaction, Sample names)
- Reaction IDs must match model reaction IDs
- Numeric values for RAS scores

## Output Format

Bounds files for each sample:

```
reaction	lower_bound	upper_bound
R00001	-125.0	125.0
R00002	-65.0	65.0
```

## Output Collections

The tool generates three types of output:

1. **Bounds files** (`ras_to_bounds/`): Individual bound files per sample (TSV format)
2. **Cell classes** (`cell_class`): Sample-to-class mapping file
3. **Complete models** (optional, `saved_models/`): Full tabular models with bounds applied

To save complete models with integrated bounds, set `--save_models true`. This creates ready-to-use model files that can be directly used with Flux Simulation or other downstream tools.

## Examples

### Basic Usage

```bash
ras_to_bounds -ms ENGRO2 \
              -ir ras_scores.tsv \
              -rs true \
              -idop output/
```

### With Custom Model

```bash
ras_to_bounds -ms Custom \
              -mo custom_model.xml \
              -ir ras_scores.tsv \
              -rs true \
              -idop output/
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Model not found" | Check model file path |
| "RAS file invalid" | Verify TSV format |
| "Infeasible solution" | Relax RAS scaling or constraints |

## See Also

- [RAS Generator](ras-generator)
- [Flux Simulation](flux-simulation)
- [Built-in Models](../reference/built-in-models)
