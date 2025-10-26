# MAREA Cluster

Perform clustering analysis on metabolic data to identify sample groups and patterns.

## Overview

MAREA Cluster performs unsupervised clustering analysis on RAS, RPS, or flux data to identify natural groupings among samples. It supports multiple clustering algorithms (K-means, DBSCAN, Hierarchical) with optional data scaling and validation metrics including elbow plots and silhouette analysis.

## Usage

### Command Line

```bash
marea_cluster -td /path/to/COBRAxy \
              -in metabolic_data.tsv \
              -cy kmeans \
              -sc true \
              -k1 2 \
              -k2 8 \
              -el true \
              -si true \
              -idop clustering_results/ \
              -ol cluster.log
```

### Galaxy Interface

Select "MAREA Cluster" from the COBRAxy tool suite and configure clustering parameters through the web interface.

## Parameters

### Required Parameters

| Parameter | Flag | Description |
|-----------|------|-------------|
| Tool Directory | `-td, --tool_dir` | Path to COBRAxy installation directory |
| Input Data | `-in, --input` | Metabolic data file (TSV format) |

### Clustering Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Cluster Type | `-cy, --cluster_type` | Clustering algorithm | kmeans |
| Data Scaling | `-sc, --scaling` | Apply data normalization | true |
| Minimum K | `-k1, --k_min` | Minimum number of clusters | 2 |
| Maximum K | `-k2, --k_max` | Maximum number of clusters | 7 |

### Analysis Options

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Elbow Plot | `-el, --elbow` | Generate elbow plot for K-means | false |
| Silhouette Analysis | `-si, --silhouette` | Generate silhouette plots | false |

### DBSCAN Specific Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Min Samples | `-ms, --min_samples` | Minimum samples per cluster | - |
| Epsilon | `-ep, --eps` | Maximum distance between samples | - |

### Output Parameters

| Parameter | Flag | Description | Default |
|-----------|------|-------------|---------|
| Output Path | `-idop, --output_path` | Results directory | clustering/ |
| Output Log | `-ol, --out_log` | Log file path | - |
| Best Cluster | `-bc, --best_cluster` | Best clustering result file | - |

## Clustering Algorithms

### K-means
**Method**: Partitional clustering using centroids
- Assumes spherical clusters
- Requires pre-specified number of clusters (k)
- Fast and scalable
- Works well with normalized data

**Best for**:
- Well-separated, compact clusters
- Large datasets
- When cluster number is approximately known

### DBSCAN
**Method**: Density-based clustering  
- Identifies clusters of varying shapes
- Automatically determines cluster number
- Robust to outliers and noise
- Requires epsilon and min_samples parameters

**Best for**:
- Irregular cluster shapes
- Datasets with noise/outliers
- Unknown number of clusters

### Hierarchical
**Method**: Agglomerative clustering with dendrograms
- Creates tree-like cluster hierarchy
- No need to specify cluster number initially  
- Deterministic results
- Provides multiple resolution levels

**Best for**:
- Small to medium datasets
- When cluster hierarchy is important
- Exploratory analysis

## Input Format

### Metabolic Data File

Tab-separated format with samples as rows and reactions/metabolites as columns:

```
Sample	R00001	R00002	R00003	R00004	...
Sample1	1.25	0.85	1.42	0.78	...
Sample2	0.65	1.35	0.72	1.28	...
Sample3	2.15	2.05	0.45	0.52	...
Control1	1.05	0.98	1.15	1.08	...
Control2	0.95	1.12	0.88	0.92	...
```

**Requirements**:
- First column: sample identifiers
- Subsequent columns: feature values (RAS, RPS, fluxes)
- Missing values: use 0 or leave empty
- Numeric data only (excluding sample names)

## Data Preprocessing

### Scaling Options

#### Standard Scaling (Recommended)
- Mean centering and unit variance scaling
- Formula: `(x - mean) / std`
- Ensures equal feature contribution
- Required for distance-based algorithms

#### No Scaling
- Use original data values
- May be appropriate for already normalized data
- Risk of feature dominance by high-magnitude variables

### Feature Selection

Consider preprocessing steps:
- Remove low-variance features
- Apply dimensionality reduction (PCA)
- Select most variable reactions/metabolites
- Handle missing data appropriately

## Output Files

### Cluster Assignments

#### Best Clustering Result (`best_clusters.tsv`)
```
Sample	Cluster	Silhouette_Score
Sample1	1	0.73
Sample2	1	0.68  
Sample3	2	0.81
Control1	0	0.59
Control2	0	0.62
```

#### All K Results (`clustering_results_k{n}.tsv`)
Individual files for each tested cluster number.

### Validation Metrics

#### Elbow Plot (`elbow_plot.png`)
- X-axis: Number of clusters (k)
- Y-axis: Within-cluster sum of squares (WCSS)
- Identifies optimal k at the "elbow" point

#### Silhouette Plots (`silhouette_k{n}.png`)
- Individual sample silhouette scores
- Average silhouette width per cluster
- Overall clustering quality assessment

### Summary Statistics

#### Clustering Summary (`clustering_summary.txt`)
```
Algorithm: kmeans
Scaling: true
Optimal K: 3
Best Silhouette Score: 0.72
Number of Samples: 20
Feature Dimensions: 150
```

#### Cluster Characteristics (`cluster_stats.tsv`)
```
Cluster	Size	Centroid_R00001	Centroid_R00002	Avg_Silhouette
0	8	0.95	1.12	0.68
1	7	1.35	0.82	0.74
2	5	0.65	1.55	0.69
```

## Examples

### Basic K-means Clustering

```bash
# Simple K-means with elbow analysis
marea_cluster -td /opt/COBRAxy \
              -in ras_data.tsv \
              -cy kmeans \
              -sc true \
              -k1 2 \
              -k2 10 \
              -el true \
              -si true \
              -idop kmeans_results/ \
              -ol kmeans.log
```

### DBSCAN Analysis

```bash
# Density-based clustering with custom parameters
marea_cluster -td /opt/COBRAxy \
              -in flux_samples.tsv \
              -cy dbscan \
              -sc true \
              -ms 5 \
              -ep 0.5 \
              -idop dbscan_results/ \
              -bc best_dbscan_clusters.tsv \
              -ol dbscan.log
```

### Hierarchical Clustering

```bash
# Hierarchical clustering for small dataset
marea_cluster -td /opt/COBRAxy \
              -in rps_scores.tsv \
              -cy hierarchy \
              -sc true \
              -k1 2 \
              -k2 6 \
              -si true \
              -idop hierarchical_results/ \
              -ol hierarchy.log
```

### Comprehensive Clustering Analysis

```bash
# Compare multiple algorithms
algorithms=("kmeans" "dbscan" "hierarchy")
for alg in "${algorithms[@]}"; do
    marea_cluster -td /opt/COBRAxy \
                  -in metabolomics_data.tsv \
                  -cy "$alg" \
                  -sc true \
                  -k1 2 \
                  -k2 8 \
                  -el true \
                  -si true \
                  -idop "${alg}_clustering/" \
                  -ol "${alg}_cluster.log"
done
```

## Parameter Optimization

### K-means Optimization

#### Elbow Method
1. Run K-means for k = 2 to k_max
2. Plot WCSS vs k
3. Identify "elbow" point where improvement diminishes
4. Select k at elbow as optimal

#### Silhouette Analysis
1. Compute silhouette scores for each k
2. Select k with highest average silhouette score
3. Validate with silhouette plots
4. Ensure clusters are well-separated

### DBSCAN Parameter Tuning

#### Epsilon (eps) Selection
- Use k-distance plot to identify knee point
- Start with eps = average distance to k-th nearest neighbor
- Adjust based on cluster quality metrics

#### Min Samples Selection
- Rule of thumb: min_samples ≥ dimensionality + 1
- Higher values create denser clusters
- Lower values may increase noise sensitivity

### Hierarchical Clustering

#### Linkage Method
- Ward: Minimizes within-cluster variance
- Complete: Maximum distance between clusters
- Average: Mean distance between clusters
- Single: Minimum distance (prone to chaining)

## Quality Assessment

### Internal Validation Metrics

#### Silhouette Score
- Range: [-1, 1]
- >0.7: Strong clustering
- 0.5-0.7: Reasonable clustering
- <0.5: Weak clustering

#### Calinski-Harabasz Index
- Higher values indicate better clustering
- Ratio of between-cluster to within-cluster variance

#### Davies-Bouldin Index
- Lower values indicate better clustering
- Average similarity between clusters

### External Validation

When ground truth labels available:
- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
- Homogeneity and Completeness scores

## Biological Interpretation

### Cluster Characterization

#### Metabolic Pathway Analysis
- Identify enriched pathways per cluster
- Compare metabolic profiles between clusters
- Relate clusters to biological conditions

#### Sample Annotation
- Map clusters to experimental conditions
- Identify batch effects or confounders
- Validate with independent datasets

#### Feature Importance
- Determine reactions/metabolites driving clustering
- Analyze cluster centroids for biological insights
- Connect to known metabolic phenotypes

## Integration Workflow

### Upstream Data Sources

#### COBRAxy Tools
- [RAS Generator](ras-generator.md) - Cluster based on reaction activities
- [RPS Generator](rps-generator.md) - Cluster based on reaction propensities  
- [Flux Simulation](flux-simulation.md) - Cluster flux distributions

#### External Data
- Gene expression matrices
- Metabolomics datasets
- Clinical metadata

### Downstream Analysis

#### Supervised Learning
Use cluster labels for:
- Classification model training
- Biomarker discovery
- Outcome prediction

#### Differential Analysis
- Compare clusters with [MAREA](marea.md)
- Identify cluster-specific metabolic signatures
- Pathway enrichment analysis

### Typical Pipeline

```bash
# 1. Generate metabolic scores
ras_generator -td /opt/COBRAxy -in expression.tsv -ra ras.tsv

# 2. Perform clustering analysis
marea_cluster -td /opt/COBRAxy -in ras.tsv -cy kmeans \
              -sc true -k1 2 -k2 8 -el true -si true \
              -idop clusters/ -bc best_clusters.tsv

# 3. Analyze cluster differences
marea -td /opt/COBRAxy -input_data ras.tsv \
      -input_class best_clusters.tsv -comparison manyvsmany \
      -test ks -choice_map ENGRO2 -idop cluster_analysis/
```

## Tips and Best Practices

### Data Preparation
- **Normalization**: Always scale features for distance-based methods
- **Dimensionality**: Consider PCA for high-dimensional data (>1000 features)
- **Missing Values**: Handle appropriately (imputation or removal)
- **Outliers**: Identify and consider removal for K-means

### Algorithm Selection
- **K-means**: Start here for most applications
- **DBSCAN**: Use when clusters have irregular shapes or noise present
- **Hierarchical**: Choose for small datasets or when hierarchy matters

### Parameter Selection
- **Start Simple**: Begin with default parameters
- **Use Validation**: Always employ silhouette analysis
- **Cross-Validate**: Test stability across parameter ranges
- **Biological Validation**: Ensure clusters make biological sense

### Result Interpretation
- **Multiple Algorithms**: Compare results across methods
- **Stability Assessment**: Check clustering reproducibility
- **Biological Context**: Integrate with known sample characteristics
- **Statistical Testing**: Validate cluster differences formally

## Troubleshooting

### Common Issues

**Poor clustering quality**
- Check data scaling and normalization
- Assess feature selection and dimensionality
- Try different algorithms or parameters
- Evaluate data structure with PCA/t-SNE

**Algorithm doesn't converge**
- Increase iteration limits for K-means
- Adjust epsilon/min_samples for DBSCAN
- Check for numerical stability issues
- Verify input data format

**Memory or performance issues**
- Reduce dataset size or dimensionality
- Use sampling for large datasets
- Consider approximate algorithms
- Monitor system resources

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Convergence failed" | K-means iteration limit | Increase max iterations or check data |
| "No clusters found" | DBSCAN parameters too strict | Reduce eps or min_samples |
| "Memory allocation error" | Dataset too large | Reduce size or increase memory |
| "Invalid silhouette score" | Single cluster found | Adjust parameters or algorithm |

### Performance Optimization

**Large Datasets**
- Use mini-batch K-means for speed
- Sample data for parameter optimization
- Employ dimensionality reduction
- Consider distributed computing

**High-Dimensional Data**
- Apply feature selection
- Use PCA preprocessing
- Consider specialized algorithms
- Validate results carefully

## Advanced Usage

### Custom Distance Metrics

For specialized applications, modify distance calculations:

```python
# Custom distance function for metabolic data
def metabolic_distance(x, y):
    # Implement pathway-aware distance metric
    return custom_distance_value
```

### Ensemble Clustering

Combine multiple clustering results:

```bash
# Run multiple algorithms and combine
for method in kmeans dbscan hierarchy; do
    marea_cluster -cy $method -in data.tsv -idop ${method}_results/
done

# Consensus clustering (requires custom script)
python consensus_clustering.py -i *_results/best_clusters.tsv -o consensus.tsv
```

### Interactive Analysis

Generate interactive plots for exploration:

```python
import plotly.express as px
import pandas as pd

# Load clustering results  
results = pd.read_csv('best_clusters.tsv', sep='\t')
data = pd.read_csv('metabolic_data.tsv', sep='\t')

# Interactive scatter plot
fig = px.scatter(data, x='PC1', y='PC2', color=results['Cluster'])
fig.show()
```

## See Also

- [MAREA](marea.md) - Statistical analysis of cluster differences
- [RAS Generator](ras-generator.md) - Generate clustering input data
- [Flux Simulation](flux-simulation.md) - Alternative clustering data source
- [Clustering Tutorial](/tutorials/clustering-analysis.md)
- [Validation Methods Reference](/tutorials/cluster-validation.md)