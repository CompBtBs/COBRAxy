# RAS Generator

Generate Reaction Activity Scores (RAS) from gene expression data and GPR (Gene-Protein-Reaction) rules.

## Overview

The RAS Generator computes metabolic reaction activity by:
1. Mapping gene expression to reactions via GPR rules
2. Applying logical operations (AND/OR) for enzyme complexes
3. Producing activity scores for each reaction in each sample

**Input**: Gene expression data + GPR rules  
**Output**: Reaction activity scores (RAS)

## Parameters

### Required Parameters

| Parameter | Short | Type | Description |
|-----------|--------|------|-------------|
| `--input` | `-in` | file | Gene expression dataset (TSV format) |
| `--ras_output` | `-ra` | file | Output file for RAS values |
| `--rules_selector` | `-rs` | choice | Built-in model (ENGRO2, Recon, HMRcore) |

### Optional Parameters

| Parameter | Short | Type | Default | Description |
|-----------|--------|------|---------|-------------|
| `--tool_dir` | `-td` | string | auto-detected | COBRAxy installation directory (automatically detected after pip install) |
| `--none` | `-n` | boolean | true | Handle missing gene values |
| `--model_upload` | `-rl` | file | - | Custom GPR rules file |
| `--model_upload_name` | `-rn` | string | - | Custom model name |
| `--out_log` | - | file | log.txt | Output log file |

> **Note**: After installing COBRAxy via pip, the `--tool_dir` parameter is automatically detected and doesn't need to be specified.

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

## Usage Examples

### Command Line

```bash
# Basic usage with built-in model (after pip install)
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2

# With custom model and strict missing gene handling
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rl custom_rules.tsv \
  -rn "CustomModel" \
  -n false

# Explicitly specify tool directory (only needed if not using pip install)
ras_generator -td /path/to/COBRAxy \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2
```

### Galaxy Usage

1. Upload gene expression file to Galaxy
2. Select **RAS Generator** from COBRAxy tools
3. Configure parameters:
   - **Input dataset**: Your expression file
   - **Rule selector**: ENGRO2 (or other model)
   - **Handle missing genes**: Yes/No
4. Click **Execute**

## Built-in Models

### ENGRO2 (Recommended for most analyses)
- **Scope**: Focused human metabolism
- **Reactions**: ~500
- **Genes**: ~500
- **Use case**: Core metabolic analysis

### Recon (Comprehensive analysis)
- **Scope**: Complete human metabolism  
- **Reactions**: ~10,000
- **Genes**: ~2,000
- **Use case**: Genome-wide metabolic studies

## Gene ID Mapping

COBRAxy supports multiple gene identifier formats:

| Format | Example | Notes |
|--------|---------|--------|
| **HGNC ID** | HGNC:5 | Recommended, most stable |
| **HGNC Symbol** | ALDOA | Human-readable but may change |
| **Ensembl** | ENSG00000149925 | Version-specific |
| **Entrez** | 226 | Numeric identifier |

**Recommendation**: Use HGNC IDs for best compatibility and stability.



## Troubleshooting

### Common Issues

**"Gene not found" warnings**
```
Solution: Check gene ID format matches model expectations
- Verify gene identifiers (HGNC vs symbols vs Ensembl)
- Use gene mapping tools if needed
- Set -n true to handle missing genes
```

**"No computable scores" error**
```
Solution: Insufficient gene overlap between data and model
- Check gene ID format compatibility
- Verify expression file format
- Try different built-in model
```

**Empty output file**
```
Solution: Check input file format and permissions
- Ensure TSV format with proper headers
- Verify file paths are correct
- Check write permissions for output directory
```



### Debug Mode

Enable detailed logging:

```bash
ras_generator -td /path/to/COBRAxy \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2 \
  --out_log detailed_log.txt
```

Check log file for detailed error messages and processing statistics.

## Validation

### Check Output Quality

```python
import pandas as pd

# Read RAS output
ras_df = pd.read_csv('ras_output.tsv', sep='\t', index_col=0)

# Basic statistics
print(f"RAS matrix shape: {ras_df.shape}")
print(f"Non-null values: {ras_df.count().sum()}")
print(f"Value range: {ras_df.min().min():.2f} to {ras_df.max().max():.2f}")

# Check for problematic reactions
null_reactions = ras_df.isnull().all(axis=1).sum()
print(f"Reactions with no data: {null_reactions}")
```


## Integration with Other Tools

### Downstream Analysis

RAS output can be used with:

- **[MAREA](marea.md)**: Statistical enrichment analysis
- **[RAS to Bounds](ras-to-bounds.md)**: Flux constraint application
- **[MAREA Cluster](marea-cluster.md)**: Sample clustering

### Preprocessing Options

Before RAS generation:
- **Normalize** expression data (log2, quantile, etc.)
- **Filter** low-expression genes
- **Batch correct** if multiple datasets

## Advanced Usage

### Custom Model Integration

```python
# Create custom GPR rules
custom_rules = {
    'R_CUSTOM1': 'HGNC:5 and HGNC:10',
    'R_CUSTOM2': 'HGNC:15 or HGNC:20'  
}

# Save as TSV
import pandas as pd
rules_df = pd.DataFrame(list(custom_rules.items()), 
                       columns=['Reaction_ID', 'GPR'])
rules_df.to_csv('custom_rules.tsv', sep='\t', index=False)

# Use with RAS generator
args = ['-rl', 'custom_rules.tsv', '-rn', 'CustomModel']
```

### Batch Processing

```python
# Process multiple expression files
expression_files = ['data1.tsv', 'data2.tsv', 'data3.tsv']

for i, exp_file in enumerate(expression_files):
    output_file = f'ras_output_{i}.tsv'
    
    args = [
        '-td', '/path/to/COBRAxy',
        '-in', exp_file,
        '-ra', output_file,
        '-rs', 'ENGRO2'
    ]
    
    ras_generator.main(args)
    print(f"Processed {exp_file} → {output_file}")
```

## References

- [COBRApy documentation](https://cobrapy.readthedocs.io/) - Underlying metabolic modeling
- [GPR rules format](https://cobrapy.readthedocs.io/en/stable/getting_started.html#gene-protein-reaction-rules) - Standard format specification  
- [HGNC database](https://www.genenames.org/) - Gene nomenclature standards