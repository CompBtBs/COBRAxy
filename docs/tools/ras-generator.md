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
| Ignore NaN | `--none` | Handle missing gene expression | true |

## Input Format

Gene expression file (TSV):

```
Gene	Sample1	Sample2	Sample3
ALDOA	125.5	98.3	142.7
ENO1	85.2	110.4	95.8
PFKM	200.3	185.6	210.1
```

**File Format Notes:**
- Use **tab-separated** values (TSV)
- First row must contain column headers (Gene, Sample names)
- Gene names must match selected gene ID type
- Numeric values only for expression data

## GPR Rules

- **AND**: All genes required
- **OR**: Any gene sufficient
- Example: `(GENE1 and GENE2) or GENE3`

## NaN Handling

The `--none` parameter controls how missing gene expression values are treated in GPR rules:

**When `--none true` (default):**
- `(GENE1 and NaN)` → evaluated as `GENE1` value
- `(GENE1 or NaN)` → evaluated as `GENE1` value
- Missing genes don't block reaction activity calculation

**When `--none false` (strict mode):**
- `(GENE1 and NaN)` → `NaN` (reaction cannot be evaluated)
- `(GENE1 or NaN)` → `NaN` (reaction cannot be evaluated)
- Any missing gene propagates NaN through the entire GPR expression

**Recommendation**: Use default (`true`) for datasets with incomplete gene coverage.

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

- [RAS to Bounds](tools/ras-to-bounds)
- [MAREA](tools/marea)
- [Built-in Models](reference/built-in-models)
