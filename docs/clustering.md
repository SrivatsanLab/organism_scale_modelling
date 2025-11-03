# Clustering Results - Usage Guide

## Overview

The clustering parameter sweep tests 26 different clustering configurations across 5 methods **on 50D PCA space** (which captures 91.8% of variance from the 1152D ESM embeddings):

- **HDBSCAN** (8 configs): Density-based, finds clusters of varying density
- **K-means** (5 configs): Partition-based, requires specifying k
- **MiniBatch K-means** (3 configs): Faster variant of K-means for quick testing
- **Leiden** (4 configs): Graph-based community detection
- **DBSCAN** (3 configs): Density-based with fixed eps
- **Gaussian Mixture** (3 configs): Probabilistic clustering

**Important**: Clustering is performed on the **50D PCA space**, not the 2D UMAP coordinates. This preserves much more information from the original 1152D ESM embeddings for better clustering quality.

## Files Generated

```
results/clustering/
├── clusters_hdbscan_pca_minclust50.npz
├── clusters_hdbscan_pca_minclust100.npz
├── clusters_hdbscan_pca_minclust200.npz
├── clusters_hdbscan_pca_minclust500.npz
├── clusters_kmeans_pca_k50.npz
├── clusters_leiden_pca_res1.0.npz
├── clusters_minibatch_pca_k50.npz
└── ... (26 total)
```

Each `.npz` file contains:
- `cluster_labels`: Array of cluster assignments
- `gene_ids`: Gene identifiers
- `genome_ids`: Genome identifiers
- `method`: Clustering method used
- `params`: Parameters as JSON
- `metrics`: Clustering metrics (n_clusters, n_noise, etc.)

## Using in Jupyter Notebook

### 1. Load and Compare All Clustering Results

```python
import sys
sys.path.insert(0, '/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling')

from scripts.embeddings.load_clustering_results import (
    get_clustering_summary,
    load_cluster_labels,
    load_cluster_data,
    plot_clustering_comparison
)

# Get summary of all clustering results
summary = get_clustering_summary()
print(summary)

# Filter to see only HDBSCAN results
hdbscan_results = summary[summary['method'] == 'hdbscan']
print(hdbscan_results[['suffix', 'n_clusters', 'pct_noise', 'median_cluster_size']])
```

### 2. Load UMAP and Specific Clustering

```python
import numpy as np

# Load UMAP coordinates
umap_data = np.load('results/umap/umap_n15_subsample1000000.npz', allow_pickle=True)
umap_coords = umap_data['umap_embedding']
gene_ids = umap_data['gene_ids']
genome_ids = umap_data['genome_ids']

# Load cluster labels
cluster_labels = load_cluster_labels('hdbscan_pca_minclust500')

print(f"UMAP shape: {umap_coords.shape}")
print(f"Cluster labels shape: {cluster_labels.shape}")
print(f"Number of clusters: {len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)}")
```

### 3. Load COG Annotations

```python
sys.path.insert(0, 'scripts/functional_annotation')
from load_cog_metadata import load_all_cog_metadata

# Load COG metadata
cog_df = load_all_cog_metadata('results/metadata/cog_annotations')

# Create a lookup dict for faster access
cog_lookup = dict(zip(cog_df['gene_id'], cog_df['cog_category']))

# Map to your subsampled genes
cog_categories = np.array([cog_lookup.get(gid, 'No COG') for gid in gene_ids])

print(f"COG categories loaded: {len(cog_categories)}")
print(f"Annotation rate: {(cog_categories != 'No COG').sum() / len(cog_categories):.1%}")
```

### 4. Plot UMAP Colored by COG and Clusters

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# Plot 1: Color by COG category
ax = axes[0]
cog_colors = {
    'J': '#e41a1c', 'K': '#377eb8', 'L': '#4daf4a', 'C': '#fc8d62',
    'E': '#e78ac3', 'G': '#8da0cb', 'M': '#a65628', 'S': '#fb8072',
    'R': '#bebada', 'No COG': '#d3d3d3'
}

for cog, color in cog_colors.items():
    mask = cog_categories == cog
    if mask.sum() > 0:
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                  c=color, s=0.5, alpha=0.3, label=cog, rasterized=True)

ax.set_title('UMAP colored by COG Category', fontsize=14)
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=10)
ax.grid(True, alpha=0.2)

# Plot 2: Color by cluster
ax = axes[1]

# Plot noise points
mask_noise = cluster_labels == -1
if mask_noise.sum() > 0:
    ax.scatter(umap_coords[mask_noise, 0], umap_coords[mask_noise, 1],
              c='lightgray', s=0.1, alpha=0.3, label='Noise', rasterized=True)

# Plot clusters
mask_clustered = ~mask_noise
if mask_clustered.sum() > 0:
    ax.scatter(umap_coords[mask_clustered, 0], umap_coords[mask_clustered, 1],
              c=cluster_labels[mask_clustered], s=0.5, alpha=0.5,
              cmap='tab20', rasterized=True)

n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
ax.set_title(f'UMAP colored by Cluster ({n_clusters} clusters)', fontsize=14)
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('umap_cog_and_clusters.png', dpi=300, bbox_inches='tight')
plt.show()
```

### 5. Compare Multiple Clustering Results Side-by-Side

```python
# Load multiple clustering results
clustering_results = {
    'HDBSCAN min=100': load_cluster_labels('hdbscan_pca_minclust100'),
    'HDBSCAN min=500': load_cluster_labels('hdbscan_pca_minclust500'),
    'HDBSCAN min=1000': load_cluster_labels('hdbscan_pca_minclust1000'),
    'K-means k=50': load_cluster_labels('kmeans_pca_k50'),
    'K-means k=100': load_cluster_labels('kmeans_pca_k100'),
    'Leiden res=1.0': load_cluster_labels('leiden_pca_res1.0'),
}

# Plot comparison
fig = plot_clustering_comparison(umap_coords, clustering_results, figsize=(20, 12))
plt.savefig('clustering_comparison.png', dpi=300, bbox_inches='tight')
plt.show()
```

### 6. Analyze COG Distribution Per Cluster

```python
import pandas as pd

# Create DataFrame
df = pd.DataFrame({
    'gene_id': gene_ids,
    'genome_id': genome_ids,
    'umap_1': umap_coords[:, 0],
    'umap_2': umap_coords[:, 1],
    'cluster': cluster_labels,
    'cog_category': cog_categories
})

# Get COG distribution per cluster
cog_by_cluster = df[df['cluster'] != -1].groupby(['cluster', 'cog_category']).size().unstack(fill_value=0)

# Normalize to percentages
cog_by_cluster_pct = cog_by_cluster.div(cog_by_cluster.sum(axis=1), axis=0) * 100

print("COG Category Distribution by Cluster (%):")
print(cog_by_cluster_pct.round(1))

# Find clusters enriched for specific COGs
print("\nClusters enriched for Translation (J):")
translation_enriched = cog_by_cluster_pct[cog_by_cluster_pct['J'] > 10].sort_values('J', ascending=False)
print(translation_enriched[['J']])
```

### 7. Export Cluster Assignments for Further Analysis

```python
# Save cluster assignments with metadata
output_df = df.copy()
output_df.to_csv('results/umap_with_clusters_and_cogs.csv', index=False)

print(f"Saved {len(output_df):,} genes with cluster and COG annotations")
```

## Choosing the Best Clustering

Consider these factors:

1. **Number of clusters**: Should be biologically meaningful (10-500 for functional categories)
2. **Noise percentage**: <30% noise is good for HDBSCAN
3. **Cluster size distribution**: Avoid clustering dominated by 1-2 huge clusters
4. **Visual inspection**: Clusters should correspond to UMAP structure
5. **COG enrichment**: Good clusters should show COG category enrichment

Recommended starting points:
- **HDBSCAN**: `hdbscan_pca_minclust500` (balanced sensitivity, density-based)
- **K-means**: `kmeans_pca_k50` or `kmeans_pca_k100` (fast, partition-based)
- **MiniBatch K-means**: `minibatch_pca_k100` (fastest for quick exploration)
- **Leiden**: `leiden_pca_res1.0` (graph-based communities)

## Troubleshooting

**Q: Clustering file not found?**
```python
# Check which clustering results exist
import os
files = os.listdir('results/clustering')
print([f.replace('clusters_', '').replace('.npz', '') for f in files if f.endswith('.npz')])
```

**Q: How to see all available clustering configurations?**
```python
summary = get_clustering_summary()
print(summary[['suffix', 'method', 'n_clusters', 'pct_noise']].to_string())
```

**Q: Jobs still running?**
```bash
squeue -u $USER -n cluster_sweep
```
