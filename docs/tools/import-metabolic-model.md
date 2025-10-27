# Import Metabolic Model

Import and extract metabolic model components into tabular format.

## Overview

Import Metabolic Model extracts metabolic models from SBML/JSON/MAT/YAML formats into tabular summary for analysis.

**Input**: Model file or built-in models  
**Output**: Tabular data (CSV/TSV)

## Galaxy Interface

In Galaxy: **COBRAxy → Import Metabolic Model**

1. Select built-in model or upload custom file
2. Set model name and medium configuration
3. Click **Run tool**

## Command-line console

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

### Model Selection

| Parameter | Flag | Description |
|-----------|------|-------------|
| Built-in Model | `--model` | ENGRO2 or Recon |
| Custom Model | `--input` | Path to SBML/JSON/MAT/YAML file |

**Note**: Use either `--model` OR `--input`.


### Required

| Parameter | Flag | Description |
|-----------|------|-------------|
| Model Name | `--name` | Model identifier |
| Medium Selector | `--medium_selector` | Medium configuration (use `allOpen`) |
| Output Tabular | `--out_tabular` | Output file (CSV/XLSX) |
| Output Log | `--out_log` | Log file |

### Optional

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Custom Medium | `--custom_medium` | CSV file with medium constraints | - |
| Gene Format | `--gene_format` | Gene ID conversion: Default, ENSG, HGNC_ID, entrez_id | Default |

## Built-in Models

- **ENGRO2**: ~500 reactions (recommended)
- **Recon**: ~10,000 reactions (genome-wide)

See [Built-in Models](reference/built-in-models) for details.

## Supported Formats

- **Model formats**: SBML (.xml), JSON (.json), MAT (.mat), YAML (.yml)
- **Compression**: .zip, .gz, .bz2 (e.g., `model.xml.gz`)

Compressed files are automatically detected and extracted.

## Output Format

**ENGRO2 model:**
```
ReactionID	Formula	GPR	lower_bound	upper_bound	ObjectiveCoefficient	Pathway_1	Pathway_2	InMedium	TranslationIssues
R00001	A + B -> C + D	GENE1 or GENE2	-1000.0	1000.0	0.0	Glycolysis	Central_Metabolism	FALSE	
EX_glc_e	glc_e <->	-	-1000.0	1000.0	0.0	Exchange	Transport	TRUE	
```

**Other models (Recon):**
```
ReactionID	Formula	GPR	lower_bound	upper_bound	ObjectiveCoefficient	InMedium	TranslationIssues
R00001	A + B -> C + D	GENE1 or GENE2	-1000.0	1000.0	0.0	FALSE	
EX_glc_e	glc_e <->	-	-1000.0	1000.0	0.0	TRUE	
```

**File Format Notes:**
- Output can be **tab-separated** (CSV) or Excel (XLSX)
- Contains all model information in tabular format
- Can be edited and re-imported using Export Metabolic Model

## Understanding Medium Composition

Exchange reactions with `InMedium = TRUE` represent nutrients in the medium:
- **Lower bound**: Uptake rate (negative value, e.g., -10 = uptake 10 mmol/gDW/hr)
- **Upper bound**: Secretion rate (positive value)

Example:
```
EX_glc_e	glc_e <->	-	-10.0	1000.0	0.0	TRUE
```
Glucose uptake: 10 mmol/gDW/hr (lower bound = -10)

More info: [COBRApy Media Documentation](https://cobrapy.readthedocs.io/en/latest/media.html)

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

- [Export Metabolic Model](reference/export-metabolic-model)
- [RAS Generator](tools/ras-generator)
- [RPS Generator](tools/rps-generator)
