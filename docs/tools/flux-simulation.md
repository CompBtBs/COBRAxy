# Flux Simulation

Sample metabolic fluxes using constraint-based modeling with CBS or OPTGP algorithms.

## Overview

Flux Simulation performs constraint-based sampling of metabolic flux distributions from constrained models. It supports two sampling algorithms (CBS and OPTGP) and provides comprehensive flux statistics including mean, median, quantiles, pFBA, FVA, and sensitivity analysis.

## Usage

### Command Line

```bash
flux_simulation -td /path/to/COBRAxy \
                -ms ENGRO2 \
                -in bounds1.tsv,bounds2.tsv \
                -ni Sample1,Sample2 \
                -a CBS \
                -ns 1000 \
                -nb 1 \
                -sd 42 \
                -ot mean,median,quantiles \
                -ota pFBA,FVA,sensitivity \
                -idop flux_results/
```

### Galaxy Interface

Select "Flux Simulation" from the COBRAxy tool suite and configure sampling parameters through the web interface.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |
| Input Bounds | `-in, --input` | Comma-separated list of bounds files |
| Sample Names | `-ni, --names` | Comma-separated sample names |
| Algorithm | `-a, --algorithm` | Sampling algorithm (CBS or OPTGP) |
| Number of Samples | `-ns, --n_samples` | Samples per batch |
| Number of Batches | `-nb, --n_batches` | Number of sampling batches |
| Random Seed | `-sd, --seed` | Random seed for reproducibility |
| Output Types | `-ot, --output_type` | Flux statistics to compute |

### Model Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Model Selector | `-ms, --model_selector` | Built-in model (ENGRO2, Custom) | ENGRO2 |
| Custom Model | `-mo, --model` | Path to custom SBML model | - |
| Model Name | `-mn, --model_name` | Custom model filename | - |

### Sampling Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Algorithm | `-a, --algorithm` | CBS or OPTGP | - |
| Thinning | `-th, --thinning` | OPTGP thinning parameter | 100 |
| Samples | `-ns, --n_samples` | Samples per batch | - |
| Batches | `-nb, --n_batches` | Number of batches | - |
| Seed | `-sd, --seed` | Random seed | - |

### Output Parameters

| Parameter | Flag | Description | Options |
|-----------|------|-------------|---------|
| Output Types | `-ot, --output_type` | Flux statistics | mean,median,quantiles,fluxes |
| Analysis Types | `-ota, --output_type_analysis` | Additional analyses | pFBA,FVA,sensitivity |
| Output Path | `-idop, --output_path` | Results directory | flux_simulation/ |
| Output Log | `-ol, --out_log` | Log file path | - |

## Algorithms

### CBS (Constraint-Based Sampling)

**Method**: Random objective function optimization
- Generates random linear combinations of reactions
- Optimizes using LP solver (GLPK preferred, COBRApy fallback)
- Fast and memory-efficient
- Suitable for large models

**Advantages**:
- High performance with GLPK
- Good coverage of solution space
- Robust to model size

### OPTGP (Optimal Growth Perturbation)

**Method**: MCMC-based sampling
- Markov Chain Monte Carlo with growth optimization
- Requires thinning to reduce autocorrelation
- More computationally intensive
- Better theoretical guarantees

**Advantages**:
- Uniform sampling guarantee
- Well-established method
- Good for smaller models

## Input Formats

### Bounds Files

Tab-separated format with reaction bounds:

```
Reaction	lower_bound	upper_bound
R00001	-1000.0	1250.5
R00002	-650.2	1000.0
R00003	0.0	2150.8
```

Multiple bounds files can be processed simultaneously by providing comma-separated paths.

### Custom Model File (Optional)

SBML format metabolic model compatible with COBRApy.

## Output Formats

### Flux Statistics

#### Mean Fluxes (`mean.csv`)
```
Reaction	Sample1	Sample2	Sample3
R00001	15.23	-8.45	22.1
R00002	0.0	12.67	-5.3
R00003	45.8	38.2	51.7
```

#### Median Fluxes (`median.csv`)
```
Reaction	Sample1	Sample2	Sample3
R00001	14.1	-7.8	21.5
R00002	0.0	11.9	-4.8
R00003	44.2	37.1	50.3
```

#### Quantiles (`quantiles.csv`)
```
Reaction	Sample1_q1	Sample1_q2	Sample1_q3	Sample2_q1	...
R00001	10.5	14.1	18.7	-12.3	...
R00002	-2.1	0.0	1.8	8.9	...
R00003	38.9	44.2	49.8	32.1	...
```

### Additional Analyses

#### pFBA (`pFBA.csv`)
Parsimonious Flux Balance Analysis results:
```
Reaction	Sample1	Sample2	Sample3
R00001	12.5	-6.7	19.3
R00002	0.0	8.9	-3.2
R00003	41.2	35.8	47.9
```

#### FVA (`FVA.csv`)
Flux Variability Analysis bounds:
```
Reaction	Sample1_min	Sample1_max	Sample2_min	Sample2_max	...
R00001	-5.2	35.8	-25.3	8.7	...
R00002	-8.9	8.9	0.0	28.4	...
R00003	15.6	78.3	10.2	65.9	...
```

#### Sensitivity (`sensitivity.csv`)
Single reaction deletion effects:
```
Reaction	Sample1	Sample2	Sample3
R00001	0.98	0.95	0.97
R00002	1.0	0.87	1.0
R00003	0.23	0.19	0.31
```

## Examples

### Basic CBS Sampling

```bash
# Simple CBS sampling with statistics
flux_simulation -td /opt/COBRAxy \
                -ms ENGRO2 \
                -in sample1_bounds.tsv,sample2_bounds.tsv \
                -ni Sample1,Sample2 \
                -a CBS \
                -ns 500 \
                -nb 2 \
                -sd 42 \
                -ot mean,median \
                -ota pFBA \
                -idop cbs_results/
```

### Comprehensive OPTGP Analysis

```bash
# Full analysis with OPTGP
flux_simulation -td /opt/COBRAxy \
                -ms ENGRO2 \
                -in bounds/*.tsv \
                -ni Sample1,Sample2,Sample3,Control1,Control2 \
                -a OPTGP \
                -th 200 \
                -ns 1000 \
                -nb 1 \
                -sd 123 \
                -ot mean,median,quantiles,fluxes \
                -ota pFBA,FVA,sensitivity \
                -idop comprehensive_analysis/ \
                -ol sampling.log
```

### Custom Model Sampling

```bash
# Use custom model with CBS
flux_simulation -td /opt/COBRAxy \
                -ms Custom \
                -mo models/tissue_specific.xml \
                -mn tissue_specific.xml \
                -in patient_bounds.tsv \
                -ni PatientA \
                -a CBS \
                -ns 2000 \
                -nb 5 \
                -sd 456 \
                -ot mean,quantiles \
                -ota FVA,sensitivity \
                -idop patient_analysis/
```

### Batch Processing Multiple Conditions

```bash
# Process multiple experimental conditions
flux_simulation -td /opt/COBRAxy \
                -ms ENGRO2 \
                -in ctrl1.tsv,ctrl2.tsv,treat1.tsv,treat2.tsv \
                -ni Control1,Control2,Treatment1,Treatment2 \
                -a CBS \
                -ns 800 \
                -nb 3 \
                -sd 789 \
                -ot mean,median,fluxes \
                -ota pFBA,FVA \
                -idop batch_conditions/
```

## Algorithm Selection Guide

### Choose CBS When:
- Large models (>1000 reactions)
- High sample throughput required  
- GLPK solver available
- Memory constraints present

### Choose OPTGP When:
- Theoretical sampling guarantees needed
- Smaller models (<500 reactions)
- Sufficient computational resources
- Publication-quality sampling required

## Performance Optimization

### CBS Optimization
- Install GLPK and swiglpk for maximum performance
- Increase batch number rather than samples per batch
- Monitor memory usage for large models

### OPTGP Optimization  
- Adjust thinning based on model size (100-500)
- Use parallel processing when available
- Consider warmup period for chain convergence

### General Tips
- Use appropriate sample sizes (500-2000 per condition)
- Balance batches vs samples for memory management
- Set consistent random seeds for reproducibility

## Quality Control

### Convergence Assessment
- Compare statistics across batches
- Check for systematic trends in sampling
- Validate against known flux ranges

### Statistical Validation
- Ensure adequate sample sizes (n≥100 recommended)
- Check for outliers and artifacts
- Validate against experimental flux data when available

### Output Verification
- Confirm mass balance constraints satisfied
- Check thermodynamic consistency
- Verify biological plausibility of results

## Integration Workflow

### Upstream Tools
- [RAS to Bounds](ras-to-bounds.md) - Generate constrained bounds from RAS
- [Import Metabolic Model](import-metabolic-model.md) - Extract model components

### Downstream Tools
- [Flux to Map](flux-to-map.md) - Visualize flux distributions on metabolic maps
- [MAREA](marea.md) - Statistical analysis of flux differences

### Typical Pipeline

```bash
# 1. Generate sample-specific bounds
ras_to_bounds -td /opt/COBRAxy -ms ENGRO2 -ir ras.tsv -idop bounds/

# 2. Sample fluxes from constrained models
flux_simulation -td /opt/COBRAxy -ms ENGRO2 -in bounds/*.tsv \
                -ni Sample1,Sample2,Sample3 -a CBS -ns 1000 \
                -ot mean,quantiles -ota pFBA,FVA -idop fluxes/

# 3. Visualize results on metabolic maps  
flux_to_map -td /opt/COBRAxy -input_data_fluxes fluxes/mean.csv \
            -choice_map ENGRO2 -idop flux_maps/
```

## Troubleshooting

### Common Issues

**CBS sampling fails**
- GLPK installation issues → Install GLPK and swiglpk
- Model infeasibility → Check bounds constraints  
- Memory errors → Reduce samples per batch

**OPTGP convergence problems**
- Poor mixing → Increase thinning parameter
- Slow convergence → Extend sampling time
- Chain stuck → Check model feasibility

**Output files missing**
- Insufficient disk space → Check available storage
- Permission errors → Verify write permissions
- Invalid sample names → Check naming conventions

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "GLPK solver failed" | Missing GLPK/swiglpk | Install GLPK libraries |
| "Model infeasible" | Over-constrained bounds | Relax constraints or check model |
| "Sampling timeout" | Insufficient time/resources | Reduce sample size or increase resources |

### Performance Issues

**Slow sampling**
- Use CBS instead of OPTGP for speed
- Reduce model size if possible  
- Increase system resources

**Memory errors**
- Lower samples per batch
- Process samples sequentially
- Use more efficient data formats

**Disk space issues**
- Monitor output file sizes
- Clean intermediate files
- Use compressed formats when possible

## Advanced Usage

### Custom Sampling Parameters

For fine-tuning sampling behavior, advanced users can modify:
- Objective function generation (CBS)
- MCMC parameters (OPTGP)  
- Convergence criteria
- Output precision and format

### Parallel Processing

```bash
# Split sampling across multiple cores/nodes
for i in {1..4}; do
    flux_simulation -td /opt/COBRAxy -ms ENGRO2 \
                    -in subset_${i}_bounds.tsv \
                    -ni Batch${i} -a CBS -ns 250 \
                    -sd $((42 + i)) -idop batch_${i}/ &
done
wait
```

### Result Aggregation

Combine results from multiple simulation runs:

```bash
# Merge statistics files
python merge_flux_results.py -i batch_*/mean.csv -o combined_mean.csv
```

## See Also

- [RAS to Bounds](ras-to-bounds.md) - Generate input constraints
- [Flux to Map](flux-to-map.md) - Visualize flux results
- [CBS Algorithm Documentation](/tutorials/cbs-algorithm.md)
- [OPTGP Algorithm Documentation](/tutorials/optgp-algorithm.md)