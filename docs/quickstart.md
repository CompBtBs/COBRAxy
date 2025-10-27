# Quick Start Guide

Get started with COBRAxy! This guide walks you through your first metabolic analysis.

## Step 1: Verify Installation

Test that COBRAxy is working:

```bash
# Check if tools are available
ras_generator --help

# Should display help text without errors
```

## Step 2: Download Sample Data

Create a sample gene expression file:

```bash
# Create sample data
cat > sample_expression.tsv << 'EOF'
Gene_ID	Control_1	Control_2	Treatment_1	Treatment_2
HGNC:5	8.5	9.2	15.7	14.3
HGNC:10	3.2	4.1	8.8	7.9
HGNC:15	7.9	8.2	4.4	5.1
HGNC:25	12.1	13.5	18.2	17.8
HGNC:30	6.3	7.1	11.5	10.8
HGNC:55	14.2	15.8	22.1	21.3
HGNC:80	5.7	6.4	2.8	3.2
HGNC:100	9.8	10.5	16.7	15.9
EOF
```

## Step 3: Generate Activity Scores

Compute Reaction Activity Scores (RAS) from your gene expression:

```bash
# Generate RAS scores using built-in ENGRO2 model
ras_generator \
  -in sample_expression.tsv \
  -ra ras_scores.tsv \
  -rs ENGRO2

# Check output
head ras_scores.tsv
```

**Expected output**:
```tsv
Reactions	Control_1	Control_2	Treatment_1	Treatment_2
R_HEX1	8.5	9.2	15.7	14.3
R_PGI	7.9	8.2	4.4	5.1
...
```

## Step 4: Create Pathway Visualizations

Generate enriched pathway maps with statistical analysis:

```bash
# Create pathway maps with statistical analysis
marea \
  -using_RAS true \
  -input_data ras_scores.tsv \
  -choice_map ENGRO2 \
  -gs true \
  -idop pathway_maps

# Check results
ls pathway_maps/
```

**Expected output**: SVG files with colored pathway maps showing metabolic changes.

## Step 5: View Results

Open the generated pathway maps:

```bash
# Open SVG files in your browser or image viewer
# Files will be in pathway_maps/ directory
firefox pathway_maps/*.svg  # Linux
open pathway_maps/*.svg     # macOS  
```

## What Just Happened?

1. **RAS Generation**: Mapped gene expression to metabolic reactions using GPR rules
2. **Statistical Analysis**: Identified significantly altered pathways between conditions
3. **Visualization**: Created colored pathway maps highlighting metabolic changes

## Next Steps

### Learn More About the Analysis

- **[Understanding RAS](tools/ras-generator)** - How activity scores are computed
- **[MAREA Analysis](tools/marea)** - Statistical enrichment methods  
- **[Data Flow](getting-started.md#analysis-workflows)** - Complete workflow overview

### Try Advanced Features

- **[Flux Sampling](tutorials/workflow.md#flux-simulation-workflow)** - Predict metabolic flux distributions
- **[Galaxy Interface](tutorials/galaxy-setup)** - Web-based analysis

### Use Your Own Data

- **[Data Formats](tutorials/data-formats)** - Prepare your expression data
- **[Troubleshooting](troubleshooting)** - Common issues and solutions

## Complete Example Pipeline

Here's the full command sequence for reference:

```bash
# Set up
cd /path/to/analysis/

# Generate sample data (or use your own)
cat > expression.tsv << 'EOF'
[your gene expression data]
EOF

# Run analysis pipeline
ras_generator -in expression.tsv -ra ras.tsv -rs ENGRO2
marea -using_RAS true -input_data ras.tsv -choice_map ENGRO2 -gs true -idop maps

# View results
ls maps/*.svg
```

## Getting Help

If something doesn't work:

1. **Check Prerequisites**: Ensure COBRAxy is properly installed
2. **Verify File Format**: Make sure your data is tab-separated TSV
3. **Review Logs**: Look for error messages in the terminal output
4. **Consult Guides**: [Troubleshooting](troubleshooting) and [Installation](installation)