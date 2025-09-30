# Metabolic Model Setting

Extract and organize metabolic model components into tabular format for analysis and integration.

## Overview

Metabolic Model Setting (metabolicModel2Tabular) extracts key components from SBML metabolic models and generates comprehensive tabular summaries. This tool processes built-in or custom models, applies medium constraints, handles gene nomenclature conversion, and outputs structured data for downstream analysis.

## Usage

### Command Line

```bash
metabolicModel2Tabular --model ENGRO2 \
                       --name ENGRO2 \
                       --medium_selector allOpen \
                       --gene_format Default \
                       --out_tabular model_data.csv \
                       --out_log extraction.log \
                       --tool_dir /path/to/COBRAxy
```

### Galaxy Interface

Select "Metabolic Model Setting" from the COBRAxy tool suite and configure model extraction parameters.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Model Name | `--name` | Model identifier for output files |
| Medium Selector | `--medium_selector` | Medium configuration option |
| Output Tabular | `--out_tabular` | Output file path (CSV or XLSX) |
| Output Log | `--out_log` | Log file for processing information |
| Tool Directory | `--tool_dir` | COBRAxy installation directory |

### Model Selection Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Built-in Model | `--model` | Pre-installed model (ENGRO2, Recon, HMRcore) | - |
| Custom Model | `--input` | Path to custom SBML/JSON model file | - |

**Note**: Provide either `--model` OR `--input`, not both.

### Optional Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Gene Format | `--gene_format` | Gene ID format conversion | Default |

## Model Selection

### Built-in Models

#### ENGRO2
- **Species**: Homo sapiens
- **Scope**: Genome-scale reconstruction
- **Reactions**: ~2,000 reactions
- **Metabolites**: ~1,500 metabolites  
- **Coverage**: Comprehensive human metabolism

#### Recon  
- **Species**: Homo sapiens
- **Scope**: Recon3D human reconstruction
- **Reactions**: ~10,000+ reactions
- **Metabolites**: ~5,000+ metabolites
- **Coverage**: Most comprehensive human model

#### HMRcore
- **Species**: Homo sapiens  
- **Scope**: Core metabolic network
- **Reactions**: ~300 essential reactions
- **Metabolites**: ~200 core metabolites
- **Coverage**: Central carbon and energy metabolism

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
User can specify custom medium constraints through Galaxy interface or by modifying the tool configuration.

## Gene Format Options

| Format | Description | Example |
|--------|-------------|---------|
| Default | Original model gene IDs | As stored in model |
| ENSNG | Ensembl Gene IDs | ENSG00000139618 |
| HGNC_SYMBOL | HUGO Gene Symbols | BRCA2 |  
| HGNC_ID | HUGO Gene Committee IDs | HGNC:1101 |
| ENTREZ | NCBI Entrez Gene IDs | 675 |

Gene format conversion uses internal mapping tables and may not cover all genes in custom models.

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
metabolicModel2Tabular --model ENGRO2 \
                       --name ENGRO2_extraction \
                       --medium_selector allOpen \
                       --gene_format Default \
                       --out_tabular ENGRO2_data.csv \
                       --out_log ENGRO2_log.txt \
                       --tool_dir /opt/COBRAxy
```

### Process Custom Model

```bash
# Extract custom SBML model with gene conversion
metabolicModel2Tabular --input /data/custom_model.xml \
                       --name CustomModel \
                       --medium_selector allOpen \
                       --gene_format HGNC_SYMBOL \
                       --out_tabular custom_model_data.xlsx \
                       --out_log custom_extraction.log \
                       --tool_dir /opt/COBRAxy
```

### Extract Core Model for Quick Analysis

```bash  
# Extract HMRcore for rapid prototyping
metabolicModel2Tabular --model HMRcore \
                       --name CoreModel \
                       --medium_selector allOpen \
                       --gene_format ENSNG \
                       --out_tabular core_reactions.csv \
                       --out_log core_log.txt \
                       --tool_dir /opt/COBRAxy
```

### Batch Processing Multiple Models

```bash
#!/bin/bash
models=("ENGRO2" "HMRcore" "Recon")
for model in "${models[@]}"; do
    metabolicModel2Tabular --model "$model" \
                           --name "${model}_extract" \
                           --medium_selector allOpen \
                           --gene_format HGNC_SYMBOL \
                           --out_tabular "${model}_data.csv" \
                           --out_log "${model}_log.txt" \
                           --tool_dir /opt/COBRAxy
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
metabolicModel2Tabular --model ENGRO2 --name ModelData \
                       --out_tabular model_components.csv

# 2. Use extracted data for RAS analysis
ras_generator -td /opt/COBRAxy -rs Custom \
              -rl model_components.csv \
              -in expression_data.tsv -ra ras_scores.tsv

# 3. Apply constraints and sample fluxes
ras_to_bounds -td /opt/COBRAxy -ms Custom -mo model_components.csv \
              -ir ras_scores.tsv -idop constrained_bounds/

# 4. Visualize results
marea -td /opt/COBRAxy -input_data ras_scores.tsv \
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
- **HMRcore**: Fast processing for algorithm development  
- **Recon**: Comprehensive analysis requiring computational resources
- **Custom**: Organism-specific or specialized models

### Gene Format Selection
- **Default**: Preserve original model annotations
- **HGNC_SYMBOL**: Human-readable gene names
- **ENSNG**: Stable identifiers for bioinformatics
- **ENTREZ**: Cross-database compatibility

### Output Format Optimization
- **CSV**: Lightweight, universal compatibility
- **XLSX**: Rich formatting, multiple sheets possible
- Choose based on downstream analysis requirements

### Performance Considerations
- Large models (Recon) may require substantial memory
- Gene format conversion adds processing time
- Consider batch processing for multiple extractions

## Troubleshooting

### Common Issues

**Model loading fails**
- Check file format and compression
- Verify SBML validity for custom models
- Ensure sufficient system memory

**Gene format conversion errors**
- Mapping tables may not cover all genes
- Original gene IDs retained when conversion fails
- Check log file for conversion statistics

**Empty output file**
- Model may contain no reactions
- Check model file integrity
- Verify tool directory configuration

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Model file not found" | Invalid file path | Check file location and permissions |
| "Unsupported format" | Invalid model format | Use SBML, JSON, MAT, or YML |
| "Gene mapping failed" | Missing gene conversion data | Use Default format or update mappings |
| "Memory allocation error" | Insufficient system memory | Use smaller model or increase memory |

### Performance Issues

**Slow processing**
- Large models require more time
- Gene conversion adds overhead
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

### Custom Gene Mapping

Advanced users can extend gene format conversion by modifying mapping files in the `local/mappings/` directory.

### Batch Extraction Script

```python
#!/usr/bin env python3
import subprocess
import sys

models = ['ENGRO2', 'HMRcore', 'Recon']
formats = ['Default', 'HGNC_SYMBOL', 'ENSNG']

for model in models:
    for fmt in formats:
        cmd = [
            'metabolicModel2Tabular',
            '--model', model,
            '--name', f'{model}_{fmt}',
            '--medium_selector', 'allOpen',
            '--gene_format', fmt,
            '--out_tabular', f'{model}_{fmt}.csv',
            '--out_log', f'{model}_{fmt}.log',
            '--tool_dir', '/opt/COBRAxy'
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

- [RAS Generator](ras-generator.md) - Use extracted GPR rules for RAS computation
- [RPS Generator](rps-generator.md) - Use reaction formulas for RPS analysis
- [Custom Model Tutorial](../tutorials/custom-model-integration.md)
- [Gene Mapping Reference](../tutorials/gene-id-conversion.md)