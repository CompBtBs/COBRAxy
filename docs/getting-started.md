# Getting Started

Welcome to COBRAxy! This guide will help you get up and running with metabolic flux analysis.

## What is COBRAxy?

COBRAxy is a comprehensive toolkit for metabolic flux analysis that bridges the gap between omics data and biological insights. It provides:

- **Data Integration**: Combine gene expression and metabolite data
- **Metabolic Modeling**: Use constraint-based models for flux analysis
- **Visualization**: Generate interactive pathway maps
- **Statistical Analysis**: Perform enrichment and sensitivity analysis

## Core Concepts

### Reaction Activity Scores (RAS)
RAS quantify how active metabolic reactions are based on gene expression data. COBRAxy computes RAS by:
1. Mapping genes to reactions via GPR (Gene-Protein-Reaction) rules
2. Applying logical operations (AND/OR) based on enzyme complexes
3. Producing activity scores for each reaction in each sample

### Reaction Propensity Scores (RPS)
RPS indicate metabolic preferences based on metabolite abundance:
1. Map metabolites to reactions as substrates/products
2. Weight by stoichiometry and frequency
3. Compute propensity scores using log-normalized formulas

### Flux Sampling
Sample feasible flux distributions using:
- **CBS (Coordinate Hit-and-Run with Rounding)**: Fast, uniform sampling
- **OptGP (Optimal Growth Parallel)**: Growth-optimized sampling

## Analysis Workflows

COBRAxy supports two main analysis paths:

### 1. Enrichment Analysis Workflow
```bash
# Generate activity scores
ras_generator → RAS values
rps_generator → RPS values

# Statistical enrichment analysis  
marea → Enriched pathway maps
```

**Use when**: You want to identify significantly altered pathways and create publication-ready maps.

### 2. Flux Simulation Workflow  
```bash
# Apply constraints to model
ras_generator → RAS values
ras_to_bounds → Constrained model

# Sample flux distributions
flux_simulation → Flux samples
flux_to_map → Final visualizations
```

**Use when**: You want to predict metabolic flux distributions and study network-wide changes.

## Your First Analysis

Let's run a basic analysis with sample data:

### Step 1: Prepare Your Data

You'll need:
- **Gene expression data**: TSV file with genes (rows) × samples (columns)
- **Metabolic model**: SBML file or use built-in models (ENGRO2, Recon)
- **Metabolite data** (optional): TSV file with metabolites (rows) × samples (columns)

### Step 2: Generate Activity Scores

```bash
# Generate RAS from expression data
ras_generator \
  -in expression_data.tsv \
  -ra ras_output.tsv \
  -rs ENGRO2
```

### Step 3: Create Pathway Maps

```bash
# Generate enriched pathway maps
marea \
  -using_RAS true \
  -input_data ras_output.tsv \
  -choice_map ENGRO2 \
  -gs true \
  -idop pathway_maps
```

### Step 4: View Results

Your analysis will generate:
- **RAS values**: `ras_output.tsv` - Activity scores for each reaction
- **Statistical maps**: `pathway_maps/` - SVG files with enrichment visualization
- **Log files**: Detailed execution logs for troubleshooting

## Built-in Models

COBRAxy includes ready-to-use metabolic models:

| Model | Organism | Reactions | Genes | Description |
|-------|----------|-----------|-------|-------------|
| **ENGRO2** | Human | ~2,000 | ~500 | Focused human metabolism model |
| **Recon** | Human | ~10,000 | ~2,000 | Comprehensive human metabolism |

Models are stored in the `src/local/` directory and include:
- SBML files
- GPR rules  
- Gene mapping tables
- Pathway templates

## Data Formats

### Gene Expression Format
```tsv
Gene_ID	Sample_1	Sample_2	Sample_3
HGNC:5	12.5	8.3	15.7
HGNC:10	3.2	4.1	2.8
HGNC:15	7.9	11.2	6.4
```

### Metabolite Format
```tsv
Metabolite_ID	Sample_1	Sample_2	Sample_3
glucose	100.5	85.3	120.7
pyruvate	45.2	38.1	52.8
lactate	23.9	41.2	19.4
```

## Next Steps

Now that you understand the basics:

1. **[Quick Start Guide](/quickstart.md)** - Complete walkthrough with example data
2. **[Galaxy Tutorial](/tutorials/galaxy-setup.md)** - Web-based analysis setup
3. **[Tools Reference](/tools/)** - Detailed documentation for each tool
4. **[Examples](/examples/)** - Real-world analysis examples

## Need Help?

- **[Troubleshooting](/troubleshooting.md)** - Common issues and solutions
- **[GitHub Issues](https://github.com/CompBtBs/COBRAxy/issues)** - Report bugs or ask questions
- **[Contributing](/contributing.md)** - Help improve COBRAxy