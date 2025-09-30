# Python API Tutorial

Learn how to use COBRAxy tools programmatically in Python scripts and analysis pipelines.

## Overview

This tutorial teaches you to integrate COBRAxy into Python workflows by calling tool main functions directly with parsed arguments.

**Time required**: ~45 minutes  
**Difficulty**: Intermediate  
**Prerequisites**: Basic Python knowledge, COBRAxy installation

## Understanding COBRAxy Architecture

### Tool Structure

Each COBRAxy tool is a Python module with:
- `main(args)` function that accepts argument list
- Command-line argument parsing
- Self-contained execution logic

```python
# General pattern for all tools
import tool_module
tool_module.main(['-arg1', 'value1', '-arg2', 'value2'])
```

### Available Tools

| Python Module | Purpose | Key Arguments |
|---------------|---------|---------------|
| `ras_generator` | Compute reaction activity scores | `-in`, `-ra`, `-rs` |
| `rps_generator` | Compute reaction propensity scores | `-id`, `-rp` |
| `marea` | Statistical pathway analysis | `-input_data`, `-choice_map` |
| `ras_to_bounds` | Apply RAS constraints to model | `-ir`, `-ms`, `-idop` |
| `flux_simulation` | Sample metabolic fluxes | `-ms`, `-in`, `-a`, `-ns` |
| `flux_to_map` | Add flux data to maps | `-if`, `-mp`, `-idop` |

## Setup Your Environment

### Import Required Modules

```python
import sys
import os
from pathlib import Path

# Add COBRAxy to Python path
cobraxy_path = "/path/to/COBRAxy"
sys.path.insert(0, cobraxy_path)

# Import COBRAxy tools
import ras_generator
import rps_generator
import marea
import ras_to_bounds
import flux_simulation  
import flux_to_map
import metabolicModel2Tabular as model_setting
```

### Set Working Directory

```python
# Set up working directory
work_dir = Path("/path/to/analysis")
work_dir.mkdir(exist_ok=True)
os.chdir(work_dir)

# COBRAxy tools expect this parameter
tool_dir = str(Path(cobraxy_path).absolute())
```

## Basic Workflow Example

### Step 1: Prepare Sample Data

```python
import pandas as pd
import numpy as np

# Create sample gene expression data
genes = ['HGNC:5', 'HGNC:10', 'HGNC:15', 'HGNC:25', 'HGNC:30']
samples = ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2']

# Generate random expression values
np.random.seed(42)
data = np.random.lognormal(mean=2, sigma=1, size=(len(genes), len(samples)))

# Create DataFrame
expression_df = pd.DataFrame(data, index=genes, columns=samples)
expression_df.index.name = 'Gene_ID'

# Save to file
expression_file = work_dir / "expression_data.tsv"
expression_df.to_csv(expression_file, sep='\t')
print(f"Created sample data: {expression_file}")
```

### Step 2: Extract Model Information

```python
# Extract model components (optional, for understanding model structure)
model_args = [
    '-td', tool_dir,
    '-ms', 'ENGRO2',  # Use built-in ENGRO2 model
    '-idop', str(work_dir / 'model_info')
]

try:
    model_setting.main(model_args)
    print("✓ Model information extracted")
except Exception as e:
    print(f"Model extraction failed: {e}")
```

### Step 3: Generate RAS Scores

```python
# Generate Reaction Activity Scores
ras_output = work_dir / "ras_scores.tsv"

ras_args = [
    '-td', tool_dir,
    '-in', str(expression_file),
    '-ra', str(ras_output), 
    '-rs', 'ENGRO2',  # Built-in model
    '-n', 'true'  # Handle missing genes
]

try:
    ras_generator.main(ras_args)
    print(f"✓ RAS scores generated: {ras_output}")
except Exception as e:
    print(f"RAS generation failed: {e}")
    raise
```

### Step 4: Generate RPS Scores (Optional)

```python
# Create sample metabolite data
metabolites = ['glucose', 'pyruvate', 'lactate', 'ATP', 'NADH']
met_data = np.random.lognormal(mean=3, sigma=0.5, size=(len(metabolites), len(samples)))

met_df = pd.DataFrame(met_data, index=metabolites, columns=samples)
met_df.index.name = 'Metabolite_ID'

metabolite_file = work_dir / "metabolite_data.tsv"
met_df.to_csv(metabolite_file, sep='\t')

# Generate Reaction Propensity Scores
rps_output = work_dir / "rps_scores.tsv"

rps_args = [
    '-td', tool_dir,
    '-id', str(metabolite_file),
    '-rp', str(rps_output)
]

try:
    rps_generator.main(rps_args)
    print(f"✓ RPS scores generated: {rps_output}")
except Exception as e:
    print(f"RPS generation warning: {e}")
    # RPS generation might fail with sample data - that's OK
```

### Step 5: Statistical Analysis with MAREA

```python
# Create enriched pathway maps
maps_output = work_dir / "pathway_maps"

marea_args = [
    '-td', tool_dir,
    '-using_RAS', 'true',
    '-input_data', str(ras_output),
    '-choice_map', 'ENGRO2',
    '-gs', 'true',  # Gene set analysis
    '-idop', str(maps_output)
]

try:
    marea.main(marea_args)
    print(f"✓ Pathway maps created: {maps_output}")
except Exception as e:
    print(f"MAREA analysis failed: {e}")
```

### Step 6: Flux Simulation Pipeline

```python
# Apply RAS constraints to model
bounds_output = work_dir / "model_bounds"

bounds_args = [
    '-td', tool_dir,
    '-ms', 'ENGRO2',
    '-ir', str(ras_output),
    '-rs', 'true',  # Use RAS values
    '-idop', str(bounds_output)
]

try:
    ras_to_bounds.main(bounds_args)
    print(f"✓ Model constraints applied: {bounds_output}")
except Exception as e:
    print(f"Bounds generation failed: {e}")
    raise

# Sample metabolic fluxes
flux_output = work_dir / "flux_samples"

flux_args = [
    '-td', tool_dir,
    '-ms', 'ENGRO2',
    '-in', str(bounds_output / "*.tsv"),  # Will be expanded by tool
    '-a', 'CBS',  # Sampling algorithm
    '-ns', '1000',  # Number of samples
    '-idop', str(flux_output)
]

try:
    flux_simulation.main(flux_args)
    print(f"✓ Flux samples generated: {flux_output}")
except Exception as e:
    print(f"Flux simulation failed: {e}")
```

### Step 7: Create Final Visualizations

```python
# Add flux data to enriched maps
final_maps = work_dir / "final_visualizations"

# Check if we have both maps and flux data
maps_dir = maps_output
flux_dir = flux_output

if maps_dir.exists() and flux_dir.exists():
    flux_to_map_args = [
        '-td', tool_dir,
        '-if', str(flux_dir / "*.tsv"),
        '-mp', str(maps_dir / "*.svg"),
        '-idop', str(final_maps)
    ]
    
    try:
        flux_to_map.main(flux_to_map_args)
        print(f"✓ Final visualizations created: {final_maps}")
    except Exception as e:
        print(f"Final mapping failed: {e}")
else:
    print("Skipping final visualization - missing input files")
```

## Advanced Usage Patterns

### Error Handling and Validation

```python
def run_cobraxy_tool(tool_module, args, description):
    """Helper function to run COBRAxy tools with error handling."""
    try:
        print(f"Running {description}...")
        tool_module.main(args)
        print(f"✓ {description} completed successfully")
        return True
    except Exception as e:
        print(f"✗ {description} failed: {e}")
        return False

# Usage
success = run_cobraxy_tool(
    ras_generator, 
    ras_args, 
    "RAS generation"
)

if not success:
    print("Pipeline stopped due to error")
    exit(1)
```

### Batch Processing Multiple Datasets

```python
def process_dataset(dataset_path, output_dir):
    """Process a single dataset through COBRAxy pipeline."""
    
    dataset_name = dataset_path.stem
    out_dir = Path(output_dir) / dataset_name
    out_dir.mkdir(exist_ok=True)
    
    # Generate RAS
    ras_file = out_dir / "ras_scores.tsv"
    ras_args = [
        '-td', tool_dir,
        '-in', str(dataset_path),
        '-ra', str(ras_file),
        '-rs', 'ENGRO2'
    ]
    
    if run_cobraxy_tool(ras_generator, ras_args, f"RAS for {dataset_name}"):
        # Continue with MAREA analysis
        maps_dir = out_dir / "maps"
        marea_args = [
            '-td', tool_dir,
            '-using_RAS', 'true',
            '-input_data', str(ras_file),
            '-choice_map', 'ENGRO2',
            '-idop', str(maps_dir)
        ]
        run_cobraxy_tool(marea, marea_args, f"MAREA for {dataset_name}")
    
    return out_dir

# Process multiple datasets
datasets = [
    "/path/to/dataset1.tsv",
    "/path/to/dataset2.tsv", 
    "/path/to/dataset3.tsv"
]

results = []
for dataset in datasets:
    result_dir = process_dataset(Path(dataset), work_dir / "batch_results")
    results.append(result_dir)
    
print(f"Processed {len(results)} datasets")
```

### Custom Analysis Pipelines

```python
class COBRAxyPipeline:
    """Custom COBRAxy analysis pipeline."""
    
    def __init__(self, tool_dir, work_dir):
        self.tool_dir = tool_dir
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(exist_ok=True)
        
    def run_enrichment_analysis(self, expression_file, model='ENGRO2'):
        """Run enrichment-focused analysis."""
        
        # Generate RAS
        ras_file = self.work_dir / "ras_scores.tsv"
        ras_args = ['-td', self.tool_dir, '-in', str(expression_file), 
                   '-ra', str(ras_file), '-rs', model]
        
        if not run_cobraxy_tool(ras_generator, ras_args, "RAS generation"):
            return None
            
        # Run MAREA
        maps_dir = self.work_dir / "enrichment_maps"  
        marea_args = ['-td', self.tool_dir, '-using_RAS', 'true',
                     '-input_data', str(ras_file), '-choice_map', model,
                     '-gs', 'true', '-idop', str(maps_dir)]
                     
        if run_cobraxy_tool(marea, marea_args, "MAREA analysis"):
            return maps_dir
        return None
    
    def run_flux_analysis(self, expression_file, model='ENGRO2', n_samples=1000):
        """Run flux sampling analysis."""
        
        # Generate RAS and apply bounds
        ras_file = self.work_dir / "ras_scores.tsv" 
        bounds_dir = self.work_dir / "bounds"
        flux_dir = self.work_dir / "flux_samples"
        
        # RAS generation
        ras_args = ['-td', self.tool_dir, '-in', str(expression_file),
                   '-ra', str(ras_file), '-rs', model]
        if not run_cobraxy_tool(ras_generator, ras_args, "RAS generation"):
            return None
            
        # Apply bounds
        bounds_args = ['-td', self.tool_dir, '-ms', model, '-ir', str(ras_file),
                      '-rs', 'true', '-idop', str(bounds_dir)]
        if not run_cobraxy_tool(ras_to_bounds, bounds_args, "Bounds application"):
            return None
            
        # Flux sampling
        flux_args = ['-td', self.tool_dir, '-ms', model, 
                    '-in', str(bounds_dir / "*.tsv"),
                    '-a', 'CBS', '-ns', str(n_samples), '-idop', str(flux_dir)]
        
        if run_cobraxy_tool(flux_simulation, flux_args, "Flux simulation"):
            return flux_dir
        return None

# Usage
pipeline = COBRAxyPipeline(tool_dir, work_dir / "custom_analysis")

# Run enrichment analysis
enrichment_results = pipeline.run_enrichment_analysis(expression_file)
if enrichment_results:
    print(f"Enrichment analysis completed: {enrichment_results}")

# Run flux analysis  
flux_results = pipeline.run_flux_analysis(expression_file, n_samples=500)
if flux_results:
    print(f"Flux analysis completed: {flux_results}")
```

## Integration with Data Analysis Libraries

### Pandas Integration

```python
# Read COBRAxy results back into pandas
ras_df = pd.read_csv(ras_output, sep='\t', index_col=0)
print(f"RAS data shape: {ras_df.shape}")
print(f"Sample statistics:\n{ras_df.describe()}")

# Filter highly variable reactions
ras_std = ras_df.std(axis=1)
variable_reactions = ras_std.nlargest(20).index
print(f"Most variable reactions: {list(variable_reactions)}")
```

### Matplotlib Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Visualize RAS distributions
plt.figure(figsize=(12, 8))
sns.heatmap(ras_df.iloc[:50], cmap='RdBu_r', center=0, cbar_kws={'label': 'RAS Score'})
plt.title('Reaction Activity Scores (Top 50 Reactions)')
plt.xlabel('Samples')
plt.ylabel('Reactions')
plt.tight_layout()
plt.savefig(work_dir / 'ras_heatmap.png', dpi=300)
plt.show()
```

## Best Practices

### 1. Environment Management
```python
# Use pathlib for cross-platform compatibility
from pathlib import Path

# Use absolute paths
tool_dir = str(Path(cobraxy_path).absolute())
work_dir = Path("/analysis").absolute()
```

### 2. Error Handling
```python
# Always wrap tool calls in try-except
try:
    ras_generator.main(ras_args)
except Exception as e:
    print(f"RAS generation failed: {e}")
    # Log details, cleanup, or alternative action
```

### 3. Argument Validation
```python
def validate_file_exists(filepath):
    """Validate input file exists."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")
    return str(path.absolute())

# Use before calling tools
expression_file = validate_file_exists(expression_file)
```

### 4. Resource Management
```python
# Monitor memory usage for large datasets
import psutil

def check_memory():
    """Check available memory."""
    memory = psutil.virtual_memory()
    print(f"Memory usage: {memory.percent}%")
    if memory.percent > 90:
        print("Warning: Low memory available")

check_memory()
# Run memory-intensive operations
check_memory()
```

## Troubleshooting

### Common Issues

**Import errors**
```python
# Check if COBRAxy path is correct
import sys
print("Python path includes:")
for p in sys.path:
    print(f"  {p}")
    
# Add COBRAxy path
sys.path.insert(0, "/correct/path/to/COBRAxy")
```

**Tool execution failures**
```python
# Enable verbose output
import logging
logging.basicConfig(level=logging.DEBUG)

# Check working directory
print(f"Current directory: {os.getcwd()}")
print(f"Directory contents: {list(Path('.').iterdir())}")
```

**File path issues**
```python
# Use absolute paths
ras_args = [
    '-td', str(Path(tool_dir).absolute()),
    '-in', str(Path(expression_file).absolute()),
    '-ra', str(Path(ras_output).absolute()),
    '-rs', 'ENGRO2'
]
```

## Next Steps

Now that you can use COBRAxy programmatically:

1. **[Tools Reference](../tools/)** - Detailed parameter documentation  
2. **[Examples](../examples/)** - Real-world analysis scripts
3. **Build custom analysis pipelines** for your research needs
4. **Integrate with workflow managers** like Snakemake or Nextflow

## Resources

- [COBRApy Documentation](https://cobrapy.readthedocs.io/) - Underlying metabolic modeling library
- [Pandas Documentation](https://pandas.pydata.org/) - Data manipulation
- [Matplotlib Gallery](https://matplotlib.org/gallery/) - Visualization examples
- [Python Pathlib](https://docs.python.org/3/library/pathlib.html) - Modern path handling