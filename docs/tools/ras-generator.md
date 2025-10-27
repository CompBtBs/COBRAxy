# RAS Generator

Generate Reaction Activity Scores (RAS) from gene expression data and GPR (Gene-Protein-Reaction) rules.

## Overview

The RAS Generator computes metabolic reaction activity by:
1. Mapping gene expression to reactions via GPR rules
2. Applying logical operations (AND/OR) for enzyme complexes
3. Producing activity scores for each reaction in each sample

**Input**: Gene expression data + GPR rules  
**Output**: Reaction activity scores (RAS)

## Usage

### Command Line

```bash
# Basic usage with built-in model
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2
```

### Galaxy Interface

1. Upload gene expression file to Galaxy
2. Select **RAS Generator** from COBRAxy tools
3. Configure parameters and click **Execute**

## Parameters

### Required Parameters

| Parameter | Short | Type | Description |
|-----------|--------|------|-------------|
| `--input` | `-in` | file | Gene expression dataset (TSV format) |
| `--ras_output` | `-ra` | file | Output file for RAS values |
| `--rules_selector` | `-rs` | choice | Built-in model (ENGRO2, Recon) |

### Optional Parameters

| Parameter | Short | Type | Default | Description |
|-----------|--------|------|---------|-------------|
| `--tool_dir` | `-td` | string | Auto-detected | Path to COBRAxy installation directory |
| `--none` | `-n` | boolean | true | Handle missing gene values |
| `--model_upload` | `-rl` | file | - | Custom GPR rules file |
| `--model_upload_name` | `-rn` | string | - | Custom model name |
| `--out_log` | - | file | log.txt | Output log file |

## Input Format

### Gene Expression File
```tsv
Gene_ID	Sample_1	Sample_2	Sample_3	Sample_4
HGNC:5	10.5	11.2	15.7	14.3
HGNC:10	3.2	4.1	8.8	7.9
HGNC:15	7.9	8.2	4.4	5.1
HGNC:25	12.1	13.5	18.2	17.8
```

**Requirements**:
- First column: Gene identifiers (HGNC, Ensembl, Entrez, etc.)
- Subsequent columns: Expression values (numeric)
- Header row with sample names
- Tab-separated format

### Custom GPR Rules File (Optional)
```tsv
Reaction_ID	GPR
R_HEX1	HGNC:4922
R_PGI	HGNC:8906
R_PFK	HGNC:8877 or HGNC:8878
R_ALDOA	HGNC:414 and HGNC:417
```

## Algorithm Details

### GPR Rule Processing

**Gene Mapping**: Each gene in the expression data is mapped to reactions via GPR rules.

**Logical Operations**:
- **OR**: `Gene1 or Gene2` → `expr1 + expr2`
- **AND**: `Gene1 and Gene2` → `min(expr1, expr2)`

**Missing Gene Handling**:
- `-n true`: Ignore missing genes in the GPR rules.
- `-n false`: Missing genes cause reaction score to be NaN

### RAS Computation

**Example**:
```
GPR: (HGNC:5 and HGNC:10) or HGNC:15
Expression: HGNC:5=10.5, HGNC:10=3.2, HGNC:15=7.9
RAS = max(min(10.5, 3.2), 7.9) = max(3.2, 7.9) = 7.9
```

## Output Format

### RAS Values File
```tsv
Reactions	Sample_1	Sample_2	Sample_3	Sample_4
R_HEX1	8.5	9.2	12.1	11.3
R_PGI	7.3	8.1	6.4	7.2
R_PFK	15.2	16.8	20.1	18.9
R_ALDOA	3.2	4.1	4.4	5.1
```

**Format**:
- First column: Reaction identifiers
- Subsequent columns: RAS values for each sample
- Missing values represented as "None"

## Examples

### Basic Usage

```bash
# Generate RAS with built-in ENGRO2 model
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2
```

### Custom Model

```bash
# Use custom GPR rules
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rl custom_rules.tsv \
  -rn "CustomModel"
```

### Strict Missing Gene Handling

```bash
# Treat missing genes as errors
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2 \
  -n false
```

### Batch Processing

```bash
# Process multiple files
for file in data/*.tsv; do
  basename=$(basename "$file" .tsv)
  ras_generator \
    -in "$file" \
    -ra "ras_${basename}.tsv" \
    -rs ENGRO2
done
```

## Algorithm Details

### GPR Rule Processing

**Gene Mapping**: Each gene in the expression data is mapped to reactions via GPR rules.

**Logical Operations**:
- **OR**: `Gene1 or Gene2` → `expr1 + expr2`
- **AND**: `Gene1 and Gene2` → `min(expr1, expr2)`

**Missing Gene Handling**:
- `-n true`: Ignore missing genes (default)
- `-n false`: Missing genes cause reaction score to be NaN

### RAS Computation Example

```
GPR: (HGNC:5 and HGNC:10) or HGNC:15
Expression: HGNC:5=10.5, HGNC:10=3.2, HGNC:15=7.9
RAS = max(min(10.5, 3.2), 7.9) = max(3.2, 7.9) = 7.9
```

## Built-in Models

COBRAxy includes two pre-installed models. See [Built-in Models Reference](../reference/built-in-models.md) for details.

- **ENGRO2** (recommended): ~2,000 reactions, ~500 genes - best for most analyses
- **Recon** (comprehensive): ~10,000 reactions, ~2,000 genes - for genome-wide studies

## Gene ID Formats

Supported gene identifier formats (HGNC ID recommended):

| Format | Example |
|--------|---------|
| HGNC ID | `HGNC:5` |
| HGNC Symbol | `ALDOA` |
| Ensembl | `ENSG00000149925` |
| Entrez | `226` |

## Integration

### Workflow Position

```
Gene Expression → RAS Generator → RAS Values
                                     ↓
                              ┌──────┴──────┐
                              ↓             ↓
                           MAREA      RAS to Bounds
                              ↓             ↓
                        Enriched Maps  Flux Simulation
```

### Downstream Tools

- **[MAREA](marea.md)**: Statistical enrichment analysis
- **[RAS to Bounds](ras-to-bounds.md)**: Flux constraint application
- **[MAREA Cluster](marea-cluster.md)**: Sample clustering

### Preprocessing Recommendations

Before RAS generation:
- **Normalize** expression data (log2, TPM, FPKM, etc.)
- **Filter** low-expression genes (optional)
- **Batch correct** if combining multiple datasets

## Troubleshooting

### Common Issues

**"Gene not found" warnings**
- Check gene ID format matches model (HGNC vs symbols vs Ensembl)
- Verify first column contains gene identifiers
- Use `-n true` to ignore missing genes

**"No computable scores" error**
- Insufficient gene overlap between data and model
- Check gene ID format compatibility
- Try different built-in model (Recon has more genes)

**Empty output file**
- Verify TSV format with tab separators
- Check file paths are correct
- Ensure write permissions for output directory

### Debug Mode

```bash
# Enable detailed logging
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2 \
  --out_log debug.log
```

### Output Validation

```python
import pandas as pd

# Load and inspect RAS output
ras_df = pd.read_csv('ras_output.tsv', sep='\t', index_col=0)

print(f"Shape: {ras_df.shape}")
print(f"Non-null: {ras_df.count().sum()}")
print(f"Range: [{ras_df.min().min():.2f}, {ras_df.max().max():.2f}]")
print(f"Empty reactions: {ras_df.isnull().all(axis=1).sum()}")
```

## Advanced Usage

### Custom GPR Rules

```python
import pandas as pd

# Create custom rules
rules = pd.DataFrame({
    'Reaction_ID': ['R_CUSTOM1', 'R_CUSTOM2', 'R_CUSTOM3'],
    'GPR': [
        'HGNC:5 and HGNC:10',
        'HGNC:15 or HGNC:20',
        '(HGNC:25 and HGNC:30) or HGNC:35'
    ]
})

rules.to_csv('custom_rules.tsv', sep='\t', index=False)
```

### Programmatic Use

```python
from cobraxy import ras_generator

# Process multiple datasets
datasets = ['control.tsv', 'treatment.tsv', 'time_series.tsv']

for dataset in datasets:
    ras_generator.main([
        '-in', dataset,
        '-ra', f'ras_{dataset}',
        '-rs', 'ENGRO2'
    ])
```

## See Also

- [MAREA](marea.md) - Statistical enrichment analysis  
- [RAS to Bounds](ras-to-bounds.md) - Apply RAS constraints
- [MAREA Cluster](marea-cluster.md) - Cluster samples by RAS
- [Import Metabolic Model](import-metabolic-model.md) - Extract custom GPR rules
- [COBRApy GPR Documentation](https://cobrapy.readthedocs.io/en/stable/getting_started.html#gene-protein-reaction-rules)
````