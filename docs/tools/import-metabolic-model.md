# Import Metabolic Model

Import and extract metabolic model components into tabular format.

## Overview

Import Metabolic Model extracts metabolic models from SBML/JSON/MAT/YAML formats into tabular summaries for analysis.

**Input**: Model files or built-in models  
**Output**: Tabular data (CSV/XLSX)

## Galaxy Interface

In Galaxy: **COBRAxy → Import Metabolic Model**

1. Select built-in model or upload custom file
2. Set model name and medium configuration
3. Click **Execute**

## Usage

```bash
# Import built-in model
importMetabolicModel \
  --model ENGRO2 \
  --name ENGRO2 \
  --medium_selector allOpen \
  --out_tabular model_data.csv \
  --out_log extraction.log
```

## Parameters

### Required

| Parameter | Flag | Description |
|-----------|------|-------------|
| Model Name | `--name` | Model identifier |
| Medium Selector | `--medium_selector` | Medium configuration (use `allOpen`) |
| Output Tabular | `--out_tabular` | Output file (CSV/XLSX) |
| Output Log | `--out_log` | Log file |

### Model Selection

| Parameter | Flag | Description |
|-----------|------|-------------|
| Built-in Model | `--model` | ENGRO2 or Recon |
| Custom Model | `--input` | Path to SBML/JSON/MAT/YAML file |

**Note**: Use either `--model` OR `--input`.

### Optional

| Parameter | Flag | Description |
|-----------|------|-------------|
| Custom Medium | `--custom_medium` | CSV file with medium constraints |

## Built-in Models

- **ENGRO2**: ~2,000 reactions (recommended)
- **Recon**: ~10,000 reactions (comprehensive)

See [Built-in Models](../reference/built-in-models.md) for details.

## Output Format

```
Reaction_ID	GPR_Rule	Reaction_Formula	Lower_Bound	Upper_Bound	Objective_Coefficient	Medium_Member	Compartment	Subsystem
R00001	GENE1 or GENE2	A + B -> C + D	-1000.0	1000.0	0.0	FALSE	cytosol	Glycolysis
EX_glc_e	-	glc_e <->	-1000.0	1000.0	0.0	TRUE	extracellular	Exchange
```

## Examples

### Extract Built-in Model

```bash
importMetabolicModel --model ENGRO2 \
                     --name ENGRO2_extraction \
                     --medium_selector allOpen \
                     --out_tabular ENGRO2_data.csv \
                     --out_log ENGRO2_log.txt
```

### Process Custom Model

```bash
importMetabolicModel --input custom_model.xml \
                     --name CustomModel \
                     --medium_selector allOpen \
                     --out_tabular custom_data.csv \
                     --out_log custom_log.txt
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Model file not found" | Check file path |
| "Unsupported format" | Use SBML, JSON, MAT, or YAML |

## See Also

- [Export Metabolic Model](export-metabolic-model.md)
- [RAS Generator](ras-generator.md)
- [RPS Generator](rps-generator.md)
