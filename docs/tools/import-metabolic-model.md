# Import Metabolic Model

Import and extract metabolic model components into tabular format for analysis and integration.

## Overview

Import Metabolic Model (importMetabolicModel) imports metabolic models from various formats (SBML, JSON, MAT, YAML) and extracts key components into comprehensive tabular summaries. This tool processes built-in or custom models, applies medium constraints, handles gene nomenclature conversion, and outputs structured data for downstream analysis.

**Input**: SBML/JSON/MAT/YAML model files  
**Output**: Tabular model data (CSV/XLSX)

## Usage

### Command Line

```bash
# Import built-in model
importMetabolicModel \
  --model ENGRO2 \
  --name ENGRO2 \
  --medium_selector allOpen \
  --out_tabular model_data.csv \
  --out_log extraction.log
```

### Galaxy Interface

1. Select **Import Metabolic Model** from COBRAxy tools
2. Choose built-in model or upload custom model file
3. Configure extraction parameters and click **Execute**

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Model Name | `--name` | Model identifier for output files |
| Medium Selector | `--medium_selector` | Medium configuration option |
| Output Tabular | `--out_tabular` | Output file path (CSV or XLSX) |
| Output Log | `--out_log` | Log file for processing information |

### Model Selection Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Built-in Model | `--model` | Pre-installed model (ENGRO2, Recon) | - |
| Custom Model | `--input` | Path to custom SBML/JSON model file | - |

**Note**: Provide either `--model` OR `--input`, not both.

### Optional Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Tool Directory | `--tool_dir` | Path to COBRAxy installation directory | Auto-detected |
| Custom Medium | `--custom_medium` | CSV file with medium constraints | - |

## Model Selection

### Built-in Models

Choose from two pre-installed models:
- **ENGRO2**: ~2,000 reactions - recommended for most analyses  
- **Recon**: ~10,000 reactions - comprehensive studies

See [Built-in Models Reference](../reference/built-in-models.md) for detailed specifications.

### Custom Models

Supported formats for custom model import:
- **SBML**: Systems Biology Markup Language (.xml, .sbml)
- **JSON**: COBRApy JSON format (.json)
- **MAT**: MATLAB format (.mat)  
- **YML**: YAML format (.yml, .yaml)
- **Compressed**: All formats support .gz, .zip, .bz2 compression

## Medium Configuration

### allOpen (Default)
- All exchange reactions unconstrained
- Maximum metabolic flexibility
- Suitable for general analysis

### Custom Medium
Users can specify custom medium constraints by providing a CSV file with exchange reaction bounds.

## Output Format

### Tabular Summary File

The output contains comprehensive model information in CSV or XLSX format:

#### Column Structure
```
Reaction_ID	GPR_Rule	Reaction_Formula	Lower_Bound	Upper_Bound	Objective_Coefficient	Medium_Member	Compartment	Subsystem
R00001	GENE1 or GENE2	A + B -> C + D	-1000.0	1000.0	0.0	FALSE	cytosol	Glycolysis
R00002	GENE3 and GENE4	E <-> F	-1000.0	1000.0	0.0	FALSE	mitochondria	TCA_Cycle
EX_glc_e	-	glc_e <->	-1000.0	1000.0	0.0	TRUE	extracellular	Exchange
```

#### Data Fields

| Field | Description | Values |
|-------|-------------|---------|
| Reaction_ID | Unique reaction identifier | String |
| GPR_Rule | Gene-protein-reaction association | Logical expression |
| Reaction_Formula | Stoichiometric equation | Metabolites with coefficients |
| Lower_Bound | Minimum flux constraint | Numeric (typically -1000) |
| Upper_Bound | Maximum flux constraint | Numeric (typically 1000) |
| Objective_Coefficient | Biomass/objective weight | Numeric (0 or 1) |
| Medium_Member | Exchange reaction flag | TRUE/FALSE |
| Compartment | Subcellular location | String (for ENGRO2 only) |
| Subsystem | Metabolic pathway | String |

## Examples

### Extract Built-in Model Data

```bash
# Extract ENGRO2 model with default settings
importMetabolicModel --model ENGRO2 \
                     --name ENGRO2_extraction \
                     --medium_selector allOpen \
                     --out_tabular ENGRO2_data.csv \
                     --out_log ENGRO2_log.txt
```

### Process Custom Model

```bash
# Extract custom SBML model
importMetabolicModel --input /data/custom_model.xml \
                     --name CustomModel \
                     --medium_selector allOpen \
                     --out_tabular custom_model_data.csv \
                     --out_log custom_extraction.log
```

### Batch Processing Multiple Models

```bash
#!/bin/bash
models=("ENGRO2" "Recon")
for model in "${models[@]}"; do
    importMetabolicModel --model "$model" \
                         --name "${model}_extract" \
                         --medium_selector allOpen \
                         --out_tabular "${model}_data.csv" \
                         --out_log "${model}_log.txt"
done
```

## Use Cases

### Model Comparison
Extract multiple models to compare:
- Reaction coverage across different reconstructions  
- Gene-reaction associations
- Pathway representation
- Metabolite compartmentalization

### Data Integration
Prepare model data for:
- Custom analysis pipelines
- Database integration
- Pathway annotation
- Cross-reference mapping

### Quality Control
Validate model properties:
- Check reaction balancing
- Verify gene associations
- Assess network connectivity
- Identify missing annotations

### Custom Analysis
Export structured data for:
- Network analysis (graph theory)
- Machine learning applications
- Statistical modeling
- Comparative genomics

## Integration Workflow

### Downstream Tools

The extracted tabular data serves as input for:

#### COBRAxy Tools
- [RAS Generator](ras-generator.md) - Use extracted GPR rules
- [RPS Generator](rps-generator.md) - Use reaction formulas
- [RAS to Bounds](ras-to-bounds.md) - Use reaction bounds
- [MAREA](marea.md) - Use reaction annotations

#### External Analysis
- **R/Bioconductor**: Import CSV for pathway analysis
- **Python/pandas**: Load data for network analysis  
- **MATLAB**: Process XLSX for modeling
- **Cytoscape**: Network visualization
- **Databases**: Populate reaction databases

### Typical Pipeline

```bash
# 1. Extract model components
importMetabolicModel --model ENGRO2 --name ModelData \
                     --out_tabular model_components.csv

# 2. Use extracted data for RAS analysis
ras_generator -rs Custom \
              -rl model_components.csv \
              -in expression_data.tsv -ra ras_scores.tsv

# 3. Apply constraints and sample fluxes
ras_to_bounds -ms Custom -mo model_components.csv \
              -ir ras_scores.tsv -idop constrained_bounds/

# 4. Visualize results
marea -input_data ras_scores.tsv \
      -choice_map Custom -custom_map custom.svg -idop results/
```

## Quality Control

### Pre-extraction Validation
- Verify model file integrity and format
- Check SBML compliance for custom models
- Validate gene ID formats and coverage
- Confirm medium constraint specifications

### Post-extraction Checks
- **Completeness**: Verify all expected reactions extracted
- **Consistency**: Check stoichiometric balance
- **Annotations**: Validate gene-reaction associations
- **Formatting**: Confirm output file structure

### Data Validation

#### Reaction Balancing
```bash
# Check for unbalanced reactions
awk -F'\t' 'NR>1 && $3 !~ /\<->\|->/ {print $1, $3}' model_data.csv
```

#### Gene Coverage
```bash
# Count reactions with GPR rules  
awk -F'\t' 'NR>1 && $2 != "" {count++} END {print count " reactions with GPR"}' model_data.csv
```

#### Exchange Reactions
```bash
# List medium components
awk -F'\t' 'NR>1 && $7 == "TRUE" {print $1}' model_data.csv
```

## Tips and Best Practices

### Model Selection
- **ENGRO2**: Balanced coverage for human tissue analysis
- **Recon**: Comprehensive analysis requiring computational resources
- **Custom**: Organism-specific or specialized models

### Output Format Optimization
- **CSV**: Lightweight, universal compatibility
- Choose based on downstream analysis requirements

### Performance Considerations
- Large models (Recon) may require substantial memory
- Consider batch processing for multiple extractions

## Troubleshooting

### Common Issues

**Model loading fails**
- Check file format and compression
- Verify SBML/JSON/MAT/YAML validity for custom models
- Ensure sufficient system memory

**Empty output file**
- Model may contain no reactions
- Check model file integrity
- Verify tool directory configuration

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Model file not found" | Invalid file path | Check file location and permissions |
| "Unsupported format" | Invalid model format | Use SBML, JSON, MAT, or YAML |
| "Memory allocation error" | Insufficient system memory | Use smaller model or increase memory |

### Performance Issues

**Slow processing**
- Large models require more time
- Monitor system resource usage

**Memory errors**
- Reduce model size if possible
- Process in smaller batches
- Increase available system memory

**Output file corruption**  
- Check disk space availability
- Verify file write permissions
- Monitor for system interruptions

## Advanced Usage

### Batch Extraction Script

```python
#!/usr/bin/env python3
import subprocess
import sys

models = ['ENGRO2', 'Recon']

for model in models:
    cmd = [
        'importMetabolicModel',
        '--model', model,
        '--name', f'{model}_data',
        '--medium_selector', 'allOpen',
        '--out_tabular', f'{model}.csv',
        '--out_log', f'{model}.log'
    ]
    subprocess.run(cmd, check=True)
```

### Database Integration

Export model data to databases:

```sql
-- Load CSV into PostgreSQL
CREATE TABLE model_reactions (
    reaction_id VARCHAR(50),
    gpr_rule TEXT,
    reaction_formula TEXT,
    lower_bound FLOAT,
    upper_bound FLOAT,
    objective_coefficient FLOAT,
    medium_member BOOLEAN,
    compartment VARCHAR(50),
    subsystem VARCHAR(100)
);

COPY model_reactions FROM 'model_data.csv' WITH CSV HEADER;
```

## See Also

- [Export Metabolic Model](export-metabolic-model.md) - Export tabular data to model formats
- [RAS Generator](ras-generator.md) - Use extracted GPR rules for RAS computation
- [RPS Generator](rps-generator.md) - Use reaction formulas for RPS analysis
- [Custom Model Tutorial](/tutorials/custom-model-integration.md)