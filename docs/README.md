# COBRAxy Documentation

> A Python-based command-line suite for metabolic flux analysis and visualization, with Galaxy integration.

COBRAxy enables the integration of transcriptomics data with COBRA-based metabolic models, offering a comprehensive framework for studying metabolism in both health and disease. With COBRAxy, users can load and enrich metabolic models by incorporating transcriptomic data and adjusting the model's medium conditions.

## Overview

COBRAxy enables constraint-based modeling and sampling techniques, allowing users to compute metabolic
flux distributions for multiple biological samples. The tool also enables the integration of medium
composition information to refine flux predictions. Additionally, COBRAxy provides a user-friendly interface
for visualizing significant flux differences between populations on an enriched metabolic map. This
extension provides a comprehensive and accessible framework for advanced metabolic analysis, enabling
researchers without extensive programming expertise to explore complex metabolic processes

By offering an intuitive and accessible platform for multi-omics integration and metabolic analysis, COBRAxy meets the growing need for tools that help researchers explore complex metabolic processes with ease.

## Key Features

- **Galaxy Tools** - Web-based analysis with intuitive interface
- **Reaction Activity Scores (RAS)** - Compute metabolic activity from gene expression data
- **Reaction Propensity Scores (RPS)** - Infer metabolic preferences from metabolite abundance
- **Flux computation** - Compute metabolic flux distributions using different optimization or sampling algorithms
- **Statistical Analysis** - Perform statistically significant flux differences between groups of samples and report on an enriched metabolic map
- **Built-in Models** - Ready-to-use models including ENGRO2 and Recon3D

## Quick Navigation

### [Installation](installation)
Install COBRAxy and get it running on your system

### [Tutorials](tutorials/)
Step-by-step guides for Galaxy usage

### [Tools Documentation](tools/)
Complete reference for all COBRAxy tools

## Data Flow

COBRAxy supports three main analysis workflows:

1. **RAS-based Enrichment Analysis**: 
   Gene Expression → RAS Generation → MAREA → Enriched Pathway Maps

2. **Flux Sampling Analysis**: 
   Gene Expression → RAS Generation → RAS to Bounds → Flux Simulation → Flux to Map → Flux-enriched Maps

3. **RAS+RPS Combined Enrichment**: 
   Gene Expression + Metabolite Data → RAS + RPS Generation → MAREA → Multi-omics Enriched Maps

## Community & Support

- **Documentation**: Complete guides and references
- **Issues**: Report bugs and request features on [GitHub](https://github.com/CompBtBs/COBRAxy/issues)
- **Contributing**: Help improve COBRAxy

**Ready to start?** Follow the [Installation Guide](installation) to get COBRAxy up and running!