# RAS Generator

Compute Reaction Activity Scores (RAS) from gene expression data.

## Overview

RAS Generator computes reaction activity scores by evaluating GPR rules with gene expression values.

## Galaxy Interface

In Galaxy: **COBRAxy → RAS Generator**

1. Select built-in model or upload custom GPR rules
2. Upload gene expression data
3. Click **Execute**

## Usage

```bash
ras_generator -rs ENGRO2 \
              -in expression_data.tsv \
              -ra ras_scores.tsv \
              -ol ras_generation.log
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Rules Selector | `-rs` | ENGRO2, Recon, or Custom | ENGRO2 |
| Input Data | `-in` | Gene expression TSV file | - |
| Output RAS | `-ra` | Output RAS scores file | - |
| Output Log | `-ol` | Log file | - |
| Custom Rules | `-rl` | Custom GPR rules file | - |
| Gene Names | `-gn` | Gene ID type | HGNC_Symbol |
| Remove Gene | `-rg` | Remove missing genes | true |

## Input Format

Gene expression file (TSV):

```
Gene	Sample1	Sample2	Sample3
ALDOA	125.5	98.3	142.7
ENO1	85.2	110.4	95.8
PFKM	200.3	185.6	210.1
```

## GPR Rules

- **AND**: All genes required
- **OR**: Any gene sufficient
- Example: `(GENE1 and GENE2) or GENE3`

## Output Format

```
Reaction	Sample1	Sample2	Sample3
R00001	125.5	98.3	142.7
R00002	85.2	110.4	95.8
```

## Examples

### Basic Usage

```bash
ras_generator -rs ENGRO2 \
              -in expression.tsv \
              -ra ras_scores.tsv
```

### Custom Rules

```bash
ras_generator -rs Custom \
              -rl custom_gpr.csv \
              -in expression.tsv \
              -ra ras_scores.tsv
```

### Strict Mode (Keep Missing Genes)

```bash
ras_generator -rs ENGRO2 \
              -in expression.tsv \
              -ra ras_scores.tsv \
              -rg false
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Gene not found" | Check gene ID format |
| "Invalid GPR" | Verify GPR rule syntax |

## See Also

- [RAS to Bounds](ras-to-bounds.md)
- [MAREA](marea.md)
- [Built-in Models](reference/built-in-models)
