# RPS Generator

Compute Reaction Presence Scores (RPS) from metabolite abundance data.

## Overview

RPS Generator calculates reaction presence scores based on metabolite availability in reaction formulas.

## Galaxy Interface

In Galaxy: **COBRAxy → RPS Generator**

1. Select built-in model or upload custom reactions
2. Upload metabolite abundance data
3. Click **Execute**

## Usage

```bash
rps_generator -rs ENGRO2 \
              -in metabolite_data.tsv \
              -rps rps_scores.tsv \
              -ol rps_generation.log
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Rules Selector | `-rs` | ENGRO2, Recon, or Custom | ENGRO2 |
| Input Data | `-in` | Metabolite abundance TSV file | - |
| Output RPS | `-rps` | Output RPS scores file | - |
| Output Log | `-ol` | Log file | - |
| Custom Rules | `-rl` | Custom reaction formulas file | - |

## Input Format

Metabolite data file (TSV):

```
Metabolite	Sample1	Sample2	Sample3
glc_c	2.5	1.8	3.2
atp_c	5.2	4.9	5.8
pyr_c	1.5	2.1	1.8
```

**File Format Notes:**
- Use **tab-separated** values (TSV)
- First row must contain column headers (Metabolite, Sample names)
- Metabolite names must include compartment suffix (e.g., _c, _m, _e)
- Numeric values only for abundance data

## Output Format

```
Reaction	Sample1	Sample2	Sample3
R00001	1.25	0.95	1.42
R00002	0.85	1.15	0.92
```

## Examples

### Basic Usage

```bash
rps_generator -rs ENGRO2 \
              -in metabolites.tsv \
              -rps rps_scores.tsv
```

### Custom Reactions

```bash
rps_generator -rs Custom \
              -rl custom_reactions.csv \
              -in metabolites.tsv \
              -rps rps_scores.tsv
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Metabolite not found" | Check metabolite nomenclature |
| "Invalid formula" | Verify reaction formula syntax |

## See Also

- [MAREA](tools/marea)
- [RAS Generator](tools/ras-generator)
- [Built-in Models](reference/built-in-models)
