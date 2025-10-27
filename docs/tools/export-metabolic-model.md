# Export Metabolic Model

Convert tabular data into COBRA metabolic models.

## Overview

Export Metabolic Model converts structured tabular data (CSV/TSV) into functional COBRA models in SBML, JSON, MATLAB, or YAML formats.

**Input**: Tabular model data (CSV/TSV)  
**Output**: SBML/JSON/MAT/YAML model files

## Galaxy Interface

In Galaxy: **COBRAxy → Export Metabolic Model**

1. Upload tabular model data file
2. Select output format (SBML/JSON/MAT/YAML)
3. Click **Execute**

## Usage

```bash
exportMetabolicModel \
  --input model_data.csv \
  --format sbml \
  --output custom_model.xml \
  --out_log conversion.log
```

## Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Input File | `--input` | Tabular file (CSV/TSV) with model data |
| Output Format | `--format` | Model format: sbml, json, mat, yaml |
| Output File | `--output` | Output model file path |
| Output Log | `--out_log` | Log file |

## Input Format

Required columns:

```csv
ReactionID,Formula,GPR,lower_bound,upper_bound,ObjectiveCoefficient,InMedium,TranslationIssues
R00001,A + B -> C + D,GENE1 or GENE2,-1000.0,1000.0,0.0,FALSE,
EX_glc_e,glc_e <->,-,-1000.0,1000.0,0.0,TRUE,
```

**File Format Notes:**
- Use **comma-separated** (CSV) or **tab-separated** (TSV)
- First row must contain column headers
- Required columns: ReactionID, Formula, lower_bound, upper_bound
- Optional columns: GPR, ObjectiveCoefficient, InMedium, Pathway_1, Pathway_2

## Reaction Formula Syntax

```
# Irreversible
A + B -> C + D

# Reversible  
A + B <-> C + D

# With stoichiometry
2 A + 3 B -> 1 C + 4 D
```

## GPR Rule Syntax

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

## Output Formats

- **SBML**: XML standard, maximum compatibility
- **JSON**: COBRApy native format
- **MATLAB**: COBRA Toolbox compatibility
- **YAML**: Human-readable format

## Examples

### Basic Export

```bash
exportMetabolicModel --input model.csv \
                     --format sbml \
                     --output model.xml \
                     --out_log conversion.log
```

### Multi-format Export

```bash
for fmt in sbml json mat yaml; do
    exportMetabolicModel --input model.csv \
                         --format "$fmt" \
                         --output "model.$fmt" \
                         --out_log "log_$fmt.txt"
done
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Formula parsing failed" | Check reaction formula syntax |
| "Model infeasible" | Review bounds and exchange reactions |

## See Also

- [Import Metabolic Model](reference/import-metabolic-model)
- [RAS to Bounds](tools/ras-to-bounds)
- [Flux Simulation](tools/flux-simulation)
