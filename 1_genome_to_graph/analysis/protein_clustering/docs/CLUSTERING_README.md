# Clustering Scripts

Dimensionality reduction, clustering, and evaluation for protein embeddings.

## Pipeline

```
PCA (1152D → 50D) → UMAP (50D → 2D) → Leiden Clustering → Evaluation
```

## Components

### 1. Dimensionality Reduction

#### PCA Preprocessing
- `add_pca_model_to_cache.py` - Add PCA model to cache for full dataset processing

#### UMAP Embedding
- `compute_umap_array.py` - Compute UMAP embeddings (GPU-accelerated)
- `compute_umap_cogonly.py` - UMAP for COG-annotated genes only
- `cluster_umap_array.py` - Alternative UMAP with clustering

**Submission scripts**:
- `submit_umap_array.sh` - Array job for multiple n_neighbors values
- `submit_umap_cogonly_array.sh` - COG-only UMAP
- `submit_umap_full_array.sh` - Full parameter sweep
- `submit_umap_rapids_container.sh` - RAPIDS GPU container

### 2. Clustering

#### Leiden Algorithm
- `cluster_leiden_comprehensive.py` - Comprehensive Leiden clustering
- `evaluate_leiden_sweep.py` - Evaluate Leiden parameter sweep

**Submission scripts**:
- `submit_leiden_sweep.sh` - Parameter sweep (multiple resolutions)
- `submit_leiden_high_res.sh` - High resolution clustering
- `submit_leiden_res1500.sh` - Specific resolution 1500
- `submit_clustering_sweep.sh` - General clustering sweep
- `submit_clustering_high_res.sh` - High resolution sweep

### 3. Evaluation

#### Quality Metrics
- `evaluate_clustering_quality.py` - ARI, AMI, silhouette, Davies-Bouldin
- `submit_quality_evaluation.sh` - SLURM submission

**Metrics computed**:
- **ARI**: Adjusted Rand Index vs COG annotations
- **AMI**: Adjusted Mutual Information vs COG
- **Silhouette**: Cluster separation quality
- **Davies-Bouldin**: Cluster compactness

#### Stability Analysis
- `evaluate_clustering_stability_efficient.py` - Memory-efficient stability testing
- `submit_stability_efficient.sh` - Array job testing 50 configurations

**Tests**: Reproducibility of clusters across independent random subsamples.

### 4. Full Dataset Assignment

- `assign_full_dataset_to_clusters.py` - Assign all 29M genes to clusters
- `submit_full_dataset_assignment.sh` - SLURM submission

**Output**:
- Gene → cluster assignments for all 29M genes
- Cluster means and standard deviations (for each PCA dimension)
- Cluster sizes and centroids

### 5. Sampling Strategies

- `diversity_subsample.py` - Diversity-based subsampling (vs random)
- `compare_sampling_strategies.py` - Compare sampling methods
- `submit_diversity_subsample.sh` - SLURM submission

**Strategies**:
- Cluster-based uniform sampling
- MaxMin (farthest-point) sampling
- Hybrid approach

## Key Results

### Current Best Configuration
- **Resolution**: 1500
- **n_neighbors**: 15
- **Dataset**: COG-only (798K genes)
- **Clusters**: 23,591
- **Homogeneity**: 61.7% weighted

### Stability Testing (In Progress)
- Testing 50 configurations: 5 resolutions × 5 n_neighbors × 2 datasets
- Job array 41488384 currently running
- Expected completion: ~6-8 hours

## Usage Examples

### Run UMAP
```bash
sbatch submit_umap_array.sh
```

### Run Leiden Clustering
```bash
sbatch submit_leiden_sweep.sh
```

### Evaluate Stability
```bash
sbatch submit_stability_efficient.sh
```

### Assign Full Dataset
```bash
sbatch submit_full_dataset_assignment.sh
```

## Output Locations

- **UMAP results**: `../../results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/`
- **Clustering results**: `../../results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/`
- **Stability results**: `../../results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability/`
- **Full assignments**: `../../results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/full_dataset_assignments_*.npz`

## See Also

- **Visualization scripts**: `../visualization/`
- **Documentation**: `../../docs/evaluation.md`
- **Current status**: `../../CURRENT_STATUS.md`
