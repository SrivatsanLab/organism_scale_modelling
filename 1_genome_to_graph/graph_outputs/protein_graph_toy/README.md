# Toy Protein Graph (GitHub-Friendly Version)

This is a compressed version of the protein cluster embeddings suitable for GitHub sharing.

## Filtering Criteria

- **Sequence identity**: ≥ 90% (MMseqs high-stringency clustering)
- **Cluster size**: ≥ 10 proteins
- **Dimensionality reduction**: PCA from 1152D → 50D

## Files

### 1. cluster_embeddings_pca50.npz
**PCA-reduced cluster embeddings** (50-dimensional)
- Shape: (365211, 50)
- File size: 66.9 MB
- dtype: float32

```python
import numpy as np
data = np.load('cluster_embeddings_pca50.npz')
embeddings = data['embeddings']  # Shape: (365211, 50)
cluster_ids = data['cluster_ids']
```

### 2. cluster_index.csv
**Cluster metadata**
- Number of clusters: 365,211
- Columns: cluster_id, genome_id, cluster_size

```python
import pandas as pd
index = pd.read_csv('cluster_index.csv')
```

### 3. pca_model.npz
**PCA transformation parameters** (for reconstruction or analysis)
- Variance explained: 0.909

```python
data = np.load('pca_model.npz')
components = data['components']  # Shape: (50, 1152)
mean = data['mean']  # Shape: (1152,)
```

### 4. gene_gene_matrix.npz
**Gene-gene (cluster-cluster) edge matrix** from operon predictions
- Shape: (365211, 365211) sparse CSR matrix
- File size: 15 MB
- Edges: Gene clusters predicted to be in operons (adjacent on genome, same strand)

```python
import scipy.sparse as sp
matrix = sp.load_npz('gene_gene_matrix.npz')
# matrix[i,j] = 1 if clusters i and j are predicted to be in same operon
```

### 5. genome_gene_matrix.npz
**Genome-gene bipartite graph edge matrix**
- Shape: (7664, 365211) sparse CSR matrix
- File size: 35 MB
- Edges: Which gene clusters are present in which genomes

```python
import scipy.sparse as sp
matrix = sp.load_npz('genome_gene_matrix.npz')
# matrix[i,j] = 1 if genome i contains gene cluster j
```

### 6. metadata/
**Cluster metadata and documentation**
- `cluster_metadata_simple.csv`: Parsed cluster IDs (genome, contig, gene number)
- `README.md`: Detailed documentation on cluster naming and COG annotations

## Statistics

- **Total clusters**: 365,211
- **Total genomes**: 7,664
- **Gene-gene edges**: ~19M (operon predictions)
- **Genome-gene edges**: ~19M (bipartite graph)
- **Mean cluster size**: 55.8
- **Median cluster size**: 20
- **Embedding dimension**: 50 (reduced from 1152)
- **Total file size**: ~135 MB

## Comparison to Full Version

The full protein graph (in `graph_outputs/protein_graph/`) will contain:
- All ~388k clusters (not just high-confidence)
- Full 1152D ESM embeddings (not PCA-reduced)
- Total size: ~1.7 GB

This toy version is suitable for:
- Quick prototyping and testing
- Sharing on GitHub (< 100 MB)
- Educational purposes
- Method development

For production use, use the full version from `graph_outputs/protein_graph/` once generated.
