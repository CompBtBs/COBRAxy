# MAREA Cluster

Cluster analysis for metabolic data (RAS/RPS scores, flux distributions).

## Overview

MAREA Cluster performs unsupervised clustering on metabolic data using K-means, DBSCAN, or hierarchical algorithms.

## Galaxy Interface

In Galaxy: **COBRAxy → MAREA Cluster**

1. Upload metabolic data file
2. Select clustering algorithm and parameters
3. Click **Execute**

## Usage

```bash
marea_cluster -in metabolic_data.tsv \
              -cy kmeans \
              -sc true \
              -k1 2 \
              -k2 10 \
              -idop output/
```

## Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Input Data | `-in` | Metabolic data TSV file | - |
| Algorithm | `-cy` | kmeans, dbscan, hierarchy | kmeans |
| Scaling | `-sc` | Scale data | false |
| K Min | `-k1` | Minimum clusters (K-means/hierarchy) | 2 |
| K Max | `-k2` | Maximum clusters (K-means/hierarchy) | 10 |
| Epsilon | `-ep` | DBSCAN radius | 0.5 |
| Min Samples | `-ms` | DBSCAN minimum samples | 5 |
| Elbow Plot | `-el` | Generate elbow plot | false |
| Silhouette | `-si` | Compute silhouette scores | false |
| Output Path | `-idop` | Output directory | marea_cluster/ |

## Input Format

```
Reaction	Sample1	Sample2	Sample3
R00001	1.25	0.85	1.42
R00002	0.65	1.35	0.72
```

**File Format Notes:**
- Use **tab-separated** values (TSV) or **comma-separated** (CSV)
- First row must contain column headers (Reaction, Sample names)
- Numeric values only for metabolic data
- Missing values should be avoided or handled before clustering

## Algorithms

- **K-means**: Fast, requires number of clusters
- **DBSCAN**: Density-based, handles noise and irregular shapes
- **Hierarchical**: Tree-based, good for small datasets

## Output

- `clusters.tsv`: Sample assignments
- `silhouette_scores.tsv`: Cluster quality metrics
- `elbow_plot.svg`: Optimal K visualization (K-means)
- `*.log`: Processing log

## Examples

### K-means Clustering

```bash
marea_cluster -in ras_data.tsv \
              -cy kmeans \
              -sc true \
              -k1 2 \
              -k2 10 \
              -el true \
              -si true \
              -idop kmeans_results/
```

### DBSCAN Clustering

```bash
marea_cluster -in flux_data.tsv \
              -cy dbscan \
              -sc true \
              -ep 0.5 \
              -ms 5 \
              -idop dbscan_results/
```

### Hierarchical Clustering

```bash
marea_cluster -in rps_data.tsv \
              -cy hierarchy \
              -sc true \
              -k1 2 \
              -k2 6 \
              -idop hierarchy_results/
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| "Convergence failed" | Increase max iterations or check data |
| "No clusters found" | Adjust DBSCAN parameters |

## See Also

- [MAREA](tools/marea)
- [RAS Generator](tools/ras-generator)
- [Flux Simulation](tools/flux-simulation)
