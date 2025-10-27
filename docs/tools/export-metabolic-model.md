# Export Metabolic Model

Export tabular data (CSV/TSV) into COBRA metabolic models in various formats.

## Overview

Export Metabolic Model (exportMetabolicModel) converts structured tabular data containing reaction information into fully functional COBRA metabolic models. This tool enables creation of custom models from spreadsheet data and supports multiple output formats including SBML, JSON, MATLAB, and YAML.

**Input**: Tabular model data (CSV/TSV)  
**Output**: SBML/JSON/MAT/YAML model files

## Usage

### Command Line

```bash
# Export tabular data to SBML
exportMetabolicModel \
  --input model_data.csv \
  --format sbml \
  --output custom_model.xml \
  --out_log conversion.log
```

### Galaxy Interface

1. Upload tabular model data to Galaxy
2. Select **Export Metabolic Model** from COBRAxy tools
3. Choose output format and click **Execute**

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Input File | `--input` | Tabular file (CSV/TSV) with model data |
| Output Format | `--format` | Model format (sbml, json, mat, yaml) |
| Output File | `--output` | Output model file path |
| Output Log | `--out_log` | Log file for conversion process |

### Optional Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Tool Directory | `--tool_dir` | Path to COBRAxy installation directory | Auto-detected |

## Input Format

### Tabular Model Data

The input file must contain structured model information with the following columns:

```csv
Reaction_ID,GPR_Rule,Reaction_Formula,Lower_Bound,Upper_Bound,Objective_Coefficient,Medium_Member,Compartment,Subsystem
R00001,GENE1 or GENE2,A + B -> C + D,-1000.0,1000.0,0.0,FALSE,cytosol,Glycolysis
R00002,GENE3 and GENE4,E <-> F,-1000.0,1000.0,0.0,FALSE,mitochondria,TCA_Cycle
EX_glc_e,-,glc_e <->,-1000.0,1000.0,0.0,TRUE,extracellular,Exchange
BIOMASS,GENE5,0.5 A + 0.3 B -> 1 BIOMASS,0.0,1000.0,1.0,FALSE,cytosol,Biomass
```

### Required Columns

| Column | Description | Format |
|--------|-------------|--------|
| **Reaction_ID** | Unique reaction identifier | String |
| **Reaction_Formula** | Stoichiometric equation | Metabolite notation |
| **Lower_Bound** | Minimum flux constraint | Numeric |
| **Upper_Bound** | Maximum flux constraint | Numeric |

### Optional Columns

| Column | Description | Default |
|--------|-------------|---------|
| **GPR_Rule** | Gene-protein-reaction association | Empty string |
| **Objective_Coefficient** | Biomass/objective weight | 0.0 |
| **Medium_Member** | Exchange reaction flag | FALSE |
| **Compartment** | Subcellular location | Empty |
| **Subsystem** | Metabolic pathway | Empty |

## Output Formats

### SBML (Systems Biology Markup Language)
- **Format**: XML-based standard
- **Extension**: `.xml` or `.sbml`
- **Use Case**: Interoperability with other tools
- **Advantages**: Widely supported, standardized

### JSON (JavaScript Object Notation)  
- **Format**: COBRApy native JSON
- **Extension**: `.json`
- **Use Case**: Python/COBRApy workflows
- **Advantages**: Human-readable, lightweight

### MATLAB (.mat)
- **Format**: MATLAB workspace format
- **Extension**: `.mat`
- **Use Case**: MATLAB COBRA Toolbox
- **Advantages**: Direct MATLAB compatibility

### YAML (YAML Ain't Markup Language)
- **Format**: Human-readable data serialization
- **Extension**: `.yml` or `.yaml`  
- **Use Case**: Configuration and documentation
- **Advantages**: Most human-readable format

## Reaction Formula Syntax

Reactions use standard metabolic notation:

```
# Irreversible
A + B -> C + D

# Reversible  
A + B <-> C + D

# With stoichiometry
2 A + 3 B -> 1 C + 4 D

# Compartmentalized
glc_c + atp_c -> g6p_c + adp_c
```

**Compartment suffixes**: `_c` (cytosol), `_m` (mitochondria), `_e` (extracellular)

## GPR Rule Syntax

Gene-Protein-Reaction rules use Boolean logic:

```
# Single gene
GENE1

# Alternative genes (OR)
GENE1 or GENE2

# Required complex (AND)
GENE1 and GENE2

# Nested logic
(GENE1 and GENE2) or GENE3
```

## Model Validation

## GPR Rule Syntax

### Logical Operators
- **AND**: Gene products required together
- **OR**: Alternative gene products
- **Parentheses**: Grouping for complex logic

### Examples
```
# Single gene
GENE1

# Alternative genes (isozymes)
GENE1 or GENE2 or GENE3

# Required genes (complex)
GENE1 and GENE2

# Complex logic
(GENE1 and GENE2) or (GENE3 and GENE4)
```

## Examples

### Create Basic Model

```bash
# Convert simple CSV to SBML model
exportMetabolicModel --input simple_model.csv \
                     --format sbml \
                     --output simple_model.xml \
                     --out_log simple_conversion.log
```

### Multi-format Export

```bash
# Create models in all supported formats
formats=("sbml" "json" "mat" "yaml")
for fmt in "${formats[@]}"; do
    exportMetabolicModel --input comprehensive_model.csv \
                         --format "$fmt" \
                         --output "model.$fmt" \
                         --out_log "conversion_$fmt.log"
done
```

### Custom Model Creation

```bash
# Build tissue-specific model from curated data
exportMetabolicModel --input liver_reactions.tsv \
                     --format sbml \
                     --output liver_model.xml \
                     --out_log liver_model.log
```

### Model Integration Pipeline

```bash
# Extract existing model, modify, and recreate
importMetabolicModel --model ENGRO2 \
                     --out_tabular base_model.csv

# Edit base_model.csv with custom reactions/constraints

# Create modified model
exportMetabolicModel --input modified_model.csv \
                     --format sbml \
                     --output custom_model.xml \
                     --out_log custom_creation.log
```

## Model Validation

### Automatic Checks

The tool performs validation during conversion:
- **Stoichiometric Balance**: Reaction mass balance
- **Metabolite Consistency**: Compartment assignments
- **Bound Validation**: Feasible constraint ranges
- **Objective Function**: Valid biomass reaction

### Post-conversion Validation

```python
import cobra

# Load and validate model
model = cobra.io.read_sbml_model('custom_model.xml')

# Check basic properties
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")

# Test model solvability
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")

# Validate mass balance
unbalanced = cobra.flux_analysis.check_mass_balance(model)
if unbalanced:
    print("Unbalanced reactions found:", unbalanced)
```

## Integration Workflow

### Upstream Data Sources

#### COBRAxy Tools
- [Import Metabolic Model](import-metabolic-model.md) - Extract tabular data for modification

#### External Sources
- **Databases**: KEGG, Reactome, BiGG
- **Literature**: Manually curated reactions
- **Spreadsheets**: User-defined custom models

### Downstream Applications

#### COBRAxy Analysis
- [RAS to Bounds](ras-to-bounds.md) - Apply constraints to custom model
- [Flux Simulation](flux-simulation.md) - Sample fluxes from custom model
- [MAREA](marea.md) - Analyze custom pathways

#### External Tools
- **COBRApy**: Python-based analysis
- **COBRA Toolbox**: MATLAB analysis  
- **OptFlux**: Strain design
- **Escher**: Pathway visualization

### Typical Pipeline

```bash
# 1. Start with existing model data
importMetabolicModel --model ENGRO2 \
                     --out_tabular base_reactions.csv

# 2. Modify/extend the reaction data
# Edit base_reactions.csv to add tissue-specific reactions

# 3. Create custom model
exportMetabolicModel --input modified_reactions.csv \
                     --format sbml \
                     --output tissue_model.xml \
                     --out_log tissue_creation.log

# 4. Validate and use custom model
ras_to_bounds --model Custom --input tissue_model.xml \
              --ras_input tissue_expression.tsv \
              --idop tissue_bounds/

# 5. Perform flux analysis
flux_simulation --model Custom --input tissue_model.xml \
                --bounds tissue_bounds/*.tsv \
                --algorithm CBS --idop tissue_fluxes/
```

## Quality Control

### Input Data Validation

#### Pre-conversion Checks
- **Format Consistency**: Verify column headers and data types
- **Reaction Completeness**: Check for missing required fields  
- **Stoichiometric Validity**: Validate reaction formulas
- **Bound Feasibility**: Ensure lower ≤ upper bounds

#### Common Data Issues
```bash
# Check for missing reaction IDs
awk -F',' 'NR>1 && ($1=="" || $1=="NA") {print "Empty ID in line " NR}' input.csv

# Validate reaction directions  
awk -F',' 'NR>1 && $3 !~ /->|<->/ {print "Invalid formula: " $1 ", " $3}' input.csv

# Check bound consistency
awk -F',' 'NR>1 && $4>$5 {print "Invalid bounds: " $1 ", LB=" $4 " > UB=" $5}' input.csv
```

### Model Quality Assessment

#### Structural Properties
- **Network Connectivity**: Ensure realistic pathway structure
- **Compartmentalization**: Validate transport reactions
- **Exchange Reactions**: Verify medium composition
- **Biomass Function**: Check objective reaction completeness

#### Functional Testing
```python
# Test model functionality
model = cobra.io.read_sbml_model('custom_model.xml')

# Check growth capability
growth = model.optimize().objective_value
print(f"Maximum growth rate: {growth}")

# Flux Variability Analysis
fva_result = cobra.flux_analysis.flux_variability_analysis(model)
blocked_reactions = fva_result[(fva_result.minimum == 0) & (fva_result.maximum == 0)]
print(f"Blocked reactions: {len(blocked_reactions)}")

# Essential gene analysis
essential_genes = cobra.flux_analysis.find_essential_genes(model)
print(f"Essential genes: {len(essential_genes)}")
```

## Tips and Best Practices

### Data Preparation
- **Consistent Naming**: Use systematic metabolite/reaction IDs
- **Compartment Notation**: Follow standard suffixes (_c, _m, _e)  
- **Balanced Reactions**: Verify mass and charge balance
- **Realistic Bounds**: Use physiologically relevant constraints

### Model Design
- **Modular Structure**: Organize reactions by pathway/subsystem
- **Exchange Reactions**: Include all necessary transport processes
- **Biomass Function**: Define appropriate growth objective
- **Gene Associations**: Add GPR rules where available

### Format Selection
- **SBML**: Choose for maximum compatibility and sharing
- **JSON**: Use for COBRApy-specific workflows
- **MATLAB**: Select for COBRA Toolbox integration
- **YAML**: Pick for human-readable documentation

### Performance Optimization
- **Model Size**: Balance comprehensiveness with computational efficiency
- **Reaction Pruning**: Remove unnecessary or blocked reactions
- **Compartmentalization**: Minimize unnecessary compartments
- **Validation**: Test model properties before distribution

## Troubleshooting

### Common Issues

**Conversion fails with format error**
- Check CSV/TSV column headers and data consistency
- Verify reaction formula syntax
- Ensure numeric fields contain valid numbers

**Model is infeasible after conversion**
- Check reaction bounds for conflicts
- Verify exchange reaction setup
- Validate stoichiometric balance

**Missing metabolites or reactions**
- Confirm all required columns present in input
- Check for empty rows or malformed data
- Validate reaction formula parsing

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Input file not found" | Invalid file path | Check file location and permissions |
| "Unknown format" | Invalid output format | Use: sbml, json, mat, or yaml |
| "Formula parsing failed" | Malformed reaction equation | Check reaction formula syntax |
| "Model infeasible" | Conflicting constraints | Review bounds and exchange reactions |

### Performance Issues

**Slow conversion**
- Large input files require more processing time
- Complex GPR rules increase parsing overhead
- Monitor system memory usage

**Memory errors**  
- Reduce model size or split into smaller files
- Increase available system memory
- Use more efficient data structures

**Output file corruption**
- Ensure sufficient disk space
- Check file write permissions
- Verify format-specific requirements

## Advanced Usage

### Batch Model Creation

```python
#!/usr/bin/env python3
import subprocess
import pandas as pd

# Create multiple tissue-specific models
tissues = ['liver', 'muscle', 'brain', 'heart']
base_data = pd.read_csv('base_model.csv')

for tissue in tissues:
    # Modify base data for tissue specificity
    tissue_data = customize_for_tissue(base_data, tissue)
    tissue_data.to_csv(f'{tissue}_model.csv', index=False)
    
    # Convert to SBML
    subprocess.run([
        'exportMetabolicModel',
        '--input', f'{tissue}_model.csv',
        '--format', 'sbml',
        '--output', f'{tissue}_model.xml',
        '--out_log', f'{tissue}_conversion.log'
    ])
```

### Model Merging

Combine multiple tabular files into comprehensive models:

```bash
# Merge core metabolism with tissue-specific pathways
cat core_reactions.csv > combined_model.csv
tail -n +2 tissue_reactions.csv >> combined_model.csv
tail -n +2 disease_reactions.csv >> combined_model.csv

# Create merged model
exportMetabolicModel --input combined_model.csv \
                     --format sbml \
                     --output comprehensive_model.xml
```

### Model Versioning

Track model versions and changes:

```bash
# Version control for model development
git add model_v1.csv
git commit -m "Initial model version"

# Create versioned models
exportMetabolicModel --input model_v1.csv --format sbml \
                     --output model_v1.xml
exportMetabolicModel --input model_v2.csv --format sbml \
                     --output model_v2.xml

# Compare model versions
cobra_compare_models model_v1.xml model_v2.xml
```

## See Also

- [Import Metabolic Model](import-metabolic-model.md) - Extract tabular data from existing models
- [RAS to Bounds](ras-to-bounds.md) - Apply constraints to custom models  
- [Flux Simulation](flux-simulation.md) - Analyze custom models with flux sampling
- [Model Creation Tutorial](/tutorials/custom-model-creation.md)
- [COBRA Model Standards](/tutorials/cobra-model-standards.md)