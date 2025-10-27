# Flux Simulation

Sample flux distributions from constraint-based metabolic models.

## Overview

Flux Simulation generates flux samples using CBS (Constraint-Based Sampling) or OPTGP (MCMC-based) algorithms.

## Galaxy Interface

In Galaxy: **COBRAxy → Flux Simulation**

1. Select model and upload bounds files
2. Choose algorithm (CBS/OPTGP) and sampling parameters
3. Click **Execute**

## Usage

```bash
flux_simulation -ms ENGRO2 \
                -in bounds/*.tsv \
                -ni Sample1,Sample2,Sample3 \
                -a CBS \
                -ns 1000 \
                -idop output/
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Model Selector | `-ms` | ENGRO2, Recon, or Custom | ENGRO2 |
| Input Bounds | `-in` | Bounds files | - |
| Name Input | `-ni` | Sample names (comma-separated) | - |
| Algorithm | `-a` | CBS or OPTGP | CBS |
| Num Samples | `-ns` | Number of samples per batch | 1000 |
| Num Batches | `-nb` | Number of batches | 1 |
| Thinning | `-th` | OPTGP thinning parameter | 100 |
| Output Type | `-ot` | mean, median, quantiles | mean,median |
| Output Path | `-idop` | Output directory | flux_simulation/ |

## Algorithms

### CBS (Constraint-Based Sampling)
- Random objective optimization
- Requires GLPK (recommended) or COBRApy solver
- Suitable for large models

### OPTGP (MCMC Sampling)
- Markov Chain Monte Carlo
- Uniform sampling guarantee
- Requires thinning parameter

## Input Format

Bounds files (TSV):

```
reaction	lower_bound	upper_bound
R00001	-125.0	125.0
R00002	-65.0	65.0
```

## Output

- `mean.csv`: Mean flux values
- `median.csv`: Median flux values
- `quantiles.csv`: Flux quantiles (25%, 75%)
- `*.log`: Processing log

## Examples

### Basic CBS Sampling

```bash
flux_simulation -ms ENGRO2 \
                -in bounds/*.tsv \
                -ni Sample1,Sample2 \
                -a CBS \
                -ns 1000 \
                -idop output/
```

### OPTGP Sampling

```bash
flux_simulation -ms ENGRO2 \
                -in bounds/*.tsv \
                -ni Sample1,Sample2 \
                -a OPTGP \
                -ns 1000 \
                -th 200 \
                -idop output/
```

### Custom Model

```bash
flux_simulation -ms Custom \
                -mo custom_model.xml \
                -in bounds/*.tsv \
                -ni Sample1 \
                -a CBS \
                -ns 2000 \
                -idop output/
```

### Batch Processing

```bash
flux_simulation -ms ENGRO2 \
                -in bounds/*.tsv \
                -ni Sample1,Sample2,Sample3 \
                -a CBS \
                -ns 500 \
                -nb 4 \
                -ot mean,quantiles \
                -idop output/
```

## Algorithm Selection

- **CBS**: Large models (>1000 reactions), GLPK available
- **OPTGP**: Theoretical sampling guarantees needed

## Troubleshooting

| Error | Solution |
|-------|----------|
| "GLPK solver failed" | Install GLPK libraries |
| "Model infeasible" | Check bounds constraints |

## See Also

- [RAS to Bounds](ras-to-bounds.md)
- [Flux to Map](flux-to-map.md)
- [Built-in Models](/reference/built-in-models)
