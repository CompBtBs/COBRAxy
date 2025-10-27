# Flux Simulation

Simulate flux distributions from constraint-based metabolic models using different optimization or sampling strategies.

## Overview

Two types of analysis are available:
- **flux optimization**
- **flux sampling**
For flux optimization, one of the following methods can be performed: parsimonious-FBA, Flux Variability Analysis, Biomass sensitivity analysis (single reaction knock-out)
The objective function, a linear combination of fluxes weighted by specific coefficients, depends on the provided metabolic network.

For flux sampling, one of the following methods can be performed: CBS (Corner-based sampling), OPTGP (Improved Artificial Centering Hit-and-Run sampler).

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
| Input Format | `--model_and_bounds` | Separate files (true) or complete models (false) | true |
| Input Bounds | `-in` | Bounds files | - |
| Name Input | `-ni` | Sample names (comma-separated) | - |
| Algorithm | `-a` | CBS or OPTGP | CBS |
| Num Samples | `-ns` | Number of samples per batch | 1000 |
| Num Batches | `-nb` | Number of batches | 1 |
| Thinning | `-th` | OPTGP thinning parameter | 100 |
| Output Type | `-ot` | mean, median, quantiles, fluxes | mean,median |
| FVA Optimality | `--perc_opt` | Optimality fraction (0.0-1.0) | 0.90 |
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

## Input Modes

The tool supports two different input formats:

### Mode 1: Model + Bounds (default, `--model_and_bounds true`)
Upload one base model + multiple bound files (one per sample/context):
- Base model: Tabular file with reaction structure (from Import Metabolic Model)
- Bounds: Individual TSV files with sample-specific constraints (from RAS to Bounds)
- Use when you have RAS-derived bounds for multiple samples

### Mode 2: Multiple Complete Models (`--model_and_bounds false`)
Upload pre-built model files, each already containing integrated bounds:
- Each file is a complete tabular model with reaction structure + bounds
- Use when models are already prepared with specific constraints
- Useful for comparing different modeling scenarios

## Input Format

Bounds files (TSV):

```
reaction	lower_bound	upper_bound
R00001	-125.0	125.0
R00002	-65.0	65.0
```

**File Format Notes:**
- Use **tab-separated** values (TSV)
- Column headers must be: reaction, lower_bound, upper_bound
- Reaction IDs must match model reaction IDs
- Numeric values for bounds

## Sampling Outputs

The tool can generate different types of output from flux sampling:

| Output Type | Description | Use Case |
|-------------|-------------|----------|
| **mean** | Mean flux across all samples | Central tendency analysis |
| **median** | Median flux across all samples | Robust central tendency |
| **quantiles** | 25th, 50th, 75th percentiles | Distribution spread analysis |
| **fluxes** | Complete flux distributions (all samples, all reactions) | Detailed statistical analysis, custom processing |

**Note**: The `fluxes` output can be very large for many samples. Use summary statistics (mean/median/quantiles) unless you need the complete distribution.

## Optimization Methods

In addition to sampling, the tool can perform constraint-based optimization analyses:

| Method | Description | Output |
|--------|-------------|--------|
| **FVA** | Flux Variability Analysis | Min/max flux ranges for each reaction |
| **pFBA** | Parsimonious FBA | Flux distribution with minimal total flux |
| **sensitivity** | Reaction knockout analysis | Biomass impact of single reaction deletions |

### FVA Optimality Fraction

The `--perc_opt` parameter (default: 0.90) controls the optimality constraint for FVA:
- **1.0**: Only optimal solutions (100% of maximum biomass)
- **0.90**: Allow suboptimal solutions (≥90% of maximum biomass)
- **Lower values**: Explore broader flux ranges

**Recommendation**: Use 0.90-0.95 for more realistic flux ranges; use 1.0 for strict optimality.

## Output

- `mean.csv`: Mean flux values
- `median.csv`: Median flux values
- `quantiles.csv`: Flux quantiles (25%, 50%, 75%)
- `fluxes/`: Complete flux distributions (if requested)
- `fva.csv`: FVA results (if requested)
- `pfba.csv`: pFBA results (if requested)
- `sensitivity.csv`: Knockout sensitivity analysis (if requested)
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

- [RAS to Bounds](tools/ras-to-bounds)
- [Flux to Map](tools/flux-to-map)
- [Built-in Models](reference/built-in-models)
