# Tutorials

Learn COBRAxy through hands-on tutorials for web-based analysis.

To set up Galaxy and start using it for web-based analyses, see the [Galaxy Setup](tutorials/galaxy-setup)

## Available Workflows

This is a collection of GALAXY workflows illustrating different applications of the tool.
The general repository is at the following link: [Galaxy workflows](http://marea4galaxy.cloud.ba.infn.it/galaxy/workflows/list_published). 

To use a workflow, click the "Import" button, and it will be added to your personal workflow page.

| Tutorial | Description | More details |
|----------|-------------|-------------|
|[Flux Enrichment Analysis - separated datasets](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=a64417ff266b740e) | Creation of maps of the fluxes differently expressed between two conditions. One gene expression dataset different for each condition. | [Flux Enrichment Analysis (Sampling Mean) — Separated Datasets](#flux-enrichment-analysis-sampling-mean--separated-datasets) |
| [Flux Enrichment Analysis (sampling mean) - separated datasets](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=16e792953f5b45db) |  Creation of maps of the fluxes differently expressed between two conditions. One gene expression dataset different for each condition. ||
| [Flux clustering (sampling mean) + Flux Enrichment Analys](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=c851ab275e52f8af) | Creation of maps of the fluxes, using one dataset differently expressed for each condition and its sample group specification||
| [Flux Enrichment Analysis (pFBA) - separated datasets](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=bf0806da5b28c6d9) | Creation of maps of the fluxes differently expressed between two conditions. One gene expression dataset different for each condition. ||
| [Flux clustering (pFBA) + Flux Enrichment Analysis](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=be0a27b9edd0db03) | Creation of maps of the fluxes, using one dataset differently expressed for each condition and its sample group specification ||
| [RAS clustering + Reaction Enrichment Analysis](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=81991b32733a4fc4) | Creation of RAS maps, one single expression gene dataset and its sample group specification ||
| [Reaction Enrichment Analysis - unified datasets](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=0d16186aaff7cbfd) |Creation of RAS maps starting from an expression dataset and its corresponding classes. One gene expression dataset as input and its classes to compare. ||
| [Reaction Enrichment Analysis - separated datasets](http://marea4galaxy.cloud.ba.infn.it/galaxy/published/workflow?id=290670ee50ab85f0) | Creation of RAS maps using the tool MaREA. Confrontation of two datasets that must be different from one another. ||

A more detailed description of the tools is available on the corresponding GALAXY page.

### Flux Enrichment Analysis (Sampling Mean) — Separated Datasets

#### Goal
Perform a **statistical analysis of fluxes** from two different datasets derived from flux simulations using the **sampling mean** method.

#### Scenario
You have **one gene expression dataset per condition** (e.g., *Cancer vs Normal*).

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform a **gene expression file** into a **RAS dataset**.  
- This step must be applied **individually for each dataset**.

#### 3. RAS to Bounds
- Use the **two RAS datasets** (one per condition) as input.  
- Generate the corresponding **flux bounds**.

#### 4. Flux Simulation
- Use the **output from the RAS to Bounds** step as input.  
- Select **sampling mean (mean)** as the simulation method.

#### 5. Metabolic Flux Enrichment Analysis
- Create a **map of significant differences** between fluxes from the two datasets.  
- Use the **flux simulation output** together with the **RASToBounds results** to identify enriched pathways or reactions.


### Flux Clustering (Sampling Mean) + Flux Enrichment Analysis

#### Goal
Creation of **flux maps** from two different datasets and **clustering** based on flux simulations using the **sampling mean** method.

#### Scenario
You have **one gene expression dataset**.

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform a **gene expression file** into a **RAS dataset**.  
- This step must be applied **for each dataset**.

#### 3. RAS to Bounds
- Use **two different RAS datasets** as input.  
- Generate the corresponding **flux bounds**.

#### 4. Flux Simulation
- Use the **output from the RAS to Bounds** step as input.  
- Select **sampling mean (mean)** as the simulation method.

#### 5. Cluster Analysis
- Perform **clustering** on the **flux dataset** obtained from the simulation.  
- Identify patterns or groups within the flux profiles.

#### 6. Metabolic Flux Enrichment Analysis
- Create **flux maps** showing **significant differences** between clusters.  
- Use:
  - The **clusters** as the *sample group specification*.  
  - The **mean of each sample** from flux sampling as the *input flux data*.  


### Flux Enrichment Analysis (pFBA) — Separated Datasets

#### Goal
Perform a **statistical analysis of fluxes** from two different datasets obtained from flux simulations using **pFBA** (parsimonious Flux Balance Analysis).

#### Scenario
You have **one gene expression dataset per condition** (e.g., *Cancer vs Normal*).

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform a **gene expression file** into a **RAS dataset**.  
- This step must be applied **individually for each dataset**.

#### 3. RAS to Bounds
- Use the **two RAS datasets** (one per condition) as input.  
- Generate the corresponding **flux bounds**.

#### 4. Flux Simulation
- Use the **output from the RAS to Bounds** step as input.  
- Select **pFBA** as the simulation method.

#### 5. Metabolic Flux Enrichment Analysis
- Perform **analysis and visualization** of **significant differences** between fluxes of the two groups (e.g., *Normal* vs *Cancer*).  
- Use the **flux simulation output** together with the **RASToBounds results** to identify enriched or altered metabolic pathways.  


### Flux Clustering (pFBA) + Flux Enrichment Analysis

#### Goal
Perform a **statistical analysis of fluxes** from two datasets using **clusters derived from flux simulations** with **pFBA** (parsimonious Flux Balance Analysis).

#### Scenario
You have **two gene expression datasets**.

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform each **gene expression file** into a **RAS dataset**.  
- This step must be applied **for each dataset**.

#### 3. RAS to Bounds
- Use the **two RAS datasets** as input.  
- Generate the corresponding **flux bounds**.

#### 4. Flux Simulation
- Use the **output from the RAS to Bounds** step as input.  
- Select **pFBA** as the simulation method.

#### 5. Cluster Analysis
- Perform **clustering** on the **flux dataset** obtained from the pFBA simulation.  
- Identify clusters or groups within the flux profiles.

#### 6. Metabolic Flux Enrichment Analysis
- Perform **analysis and visualization** of **significant differences** between fluxes of different clusters.  
- Use:
  - The **clusters** as the *sample group specification*.  
  - The **output from the pFBA flux simulation** as the *input flux data*.  
- Optionally, specify **p-value** and **fold change** thresholds to refine the analysis.  


### RAS Clustering + Reaction Enrichment Analysis

#### Goal
Perform **RAS statistical analysis** using the **MaREA** tool.  
Compare **RAS clusters** obtained from a **single gene expression dataset**.

#### Scenario
You have **one gene expression dataset**.

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform the **gene expression file** into a **RAS dataset**.

#### 3. Cluster Analysis
- Perform **clustering** on the **RAS dataset**.  
- Identify distinct clusters or groups within the data.

#### 4. Metabolic Reaction Enrichment Analysis (MaREA)
- Perform **analysis and visualization** of **significant differences** between RAS values of different clusters.  
- Use **MaREA** to detect enriched reactions and metabolic changes.  
- Optionally, specify **p-value** and **fold change** thresholds to refine the analysis.  

### Reaction Enrichment Analysis — Unified Datasets

#### Goal
Perform **RAS statistical analysis** using the **MaREA** tool.  
Compare **groups within the same gene expression dataset**.

#### Scenario
You have **one gene expression dataset** along with the **corresponding class labels** (e.g., *Normal* vs *Cancer*).

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform the **gene expression file** into a **RAS dataset**.

#### 3. Metabolic Reaction Enrichment Analysis (MaREA)
- Perform **analysis and visualization** of **significant differences** between RAS values of different groups (e.g., *Normal* vs *Cancer*).  
- The **classes** are provided as input and used for **sample group specification** in the tool.  
- Optionally, specify **p-value** and **fold change** thresholds to refine the analysis.  

### Reaction Enrichment Analysis — Separated Datasets

#### Goal
Perform **RAS statistical analysis** using the **MaREA** tool with **different gene expression datasets**.

#### Scenario
You have **one gene expression dataset per condition** (e.g., *Cancer* vs *Normal*).

#### 1. Import Metabolic Model
- Load the metabolic model.  
- Define the **medium** and the **gene nomenclature format**.

#### 2. Expression to RAS
- Transform each **gene expression file** into a **RAS dataset**.  
- This step must be applied **individually for each dataset**.

#### 3. Metabolic Reaction Enrichment Analysis (MaREA)
- Perform **analysis and visualization** of **significant differences** between RAS values from two different datasets.  
- In this scenario, the **two RAS datasets** are provided as **separate inputs**.  
- Optionally, specify **p-value** and **fold change** thresholds to refine the analysis.  


## Tutorial Data

Download example datasets used in tutorials:

```bash
# Download tutorial data
wget https://github.com/CompBtBs/COBRAxy/blob/main/data_tutorial/data_tutorial.zip
unzip tutorial_data.zip
```

The tutorial data includes Sample gene expression datasets (Cancer.txt and Normal.txt)

## Getting Help

If you encounter issues during tutorials:

1. Check the specific tutorial's troubleshooting section
2. Refer to the main [Troubleshooting Guide](troubleshooting)

## Contributing

Found an error or want to improve a tutorial? 

- [Edit on GitHub](https://github.com/CompBtBs/COBRAxy/tree/main/docs/tutorials)
- [Report issues](https://github.com/CompBtBs/COBRAxy/issues)
- Suggest new tutorial topics

Ready to start? Pick your first tutorial above! 🚀
