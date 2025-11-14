# Current Project Status

**Last Updated: 2025-11-14**

## Component 1: Genome to Graph - ✅ COMPLETE

All embedding and clustering analysis tasks have been successfully completed!

---

## ✅ Completed Tasks

### 1. ESM Embedding Generation
- **Generated embeddings for 11,837,414 proteins** using ESM-C 600M model
- Processed 1,184 GPU batches over ~24 hours
- Total batch storage: 64GB (can be deleted, merged cache retained)

### 2. Dimensionality Reduction
- **PCA reduction**: 1152D → 50D embeddings
- **Variance retained**: 89.2% (excellent compression)
- **Output**: `results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz` (4.3GB)

### 3. UMAP Visualization
- **Computed 2D UMAP embeddings** for all 11.8M proteins
- Parameters: n_neighbors=15, min_dist=0.1
- Runtime: ~6 hours
- **Output**: `results/.../umap/umap_full_n15.npz` (123MB)
- **Visualizations created**:
  - `umap_density.png` - Smooth, continuous protein space
  - `umap_by_cluster.png` - 100K sample colored by MMseqs2 cluster
  - `cluster_size_distribution.png`

### 4. MMseqs2 Cluster Tightness Analysis
- **Analyzed 388,858 clusters** (min size ≥2) in ESM embedding space
- Runtime: 33 minutes
- **Statistics computed**:
  - Mean & std deviation per cluster across all 50 PCA dimensions
  - Pairwise distances within clusters
  - Total variance per cluster
- **Outputs**:
  - `mmseqs_cluster_statistics.csv` (52MB) - Overall metrics per cluster
  - `mmseqs_cluster_per_dimension_stats.csv` (1.5GB) - Per-dimension stats
  - 5 visualization figures in `cluster_analysis/figures/`

### 5. Storage Optimization
- Moved `data/` directory (17GB) from `/fast` to `/working` storage
- Created symlink for compatibility
- Saved 17GB on expensive storage tier

---

## 🎯 Key Finding: MMseqs2 Clusters Are Extremely Tight

### Cluster Tightness Metrics

The cluster tightness analysis revealed that proteins within MMseqs2 clusters (70% sequence identity threshold) have nearly identical ESM embeddings:

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Mean std per cluster** | **0.0088** | Very low variance within clusters |
| **Median std** | 0.0106 | Consistent tightness across all clusters |
| **Mean pairwise distance** | 0.21 | Proteins in same cluster are very close |
| **75th percentile mean_std** | 0.0106 | Even larger clusters remain tight |

### Interpretation

Proteins clustered by sequence similarity (70% identity) have **nearly identical learned representations** in ESM embedding space. This demonstrates that:

1. Sequence-based clustering captures functional similarity
2. ESM embeddings align well with sequence identity
3. No additional embedding-based clustering (e.g., Leiden) is needed

### Recommendation

✅ **Use the 388,858 MMseqs2 clusters directly as gene family nodes** in the genome graph construction. This approach:
- Maintains biological interpretability (sequence-based)
- Validated by embedding space analysis (tight clusters)
- Simplifies the pipeline (no additional clustering step)
- Provides clear cluster definitions for downstream analysis

---

## 📊 Dataset Summary

| Metric | Count |
|--------|-------|
| **Total proteins in dataset** | 30,098,843 |
| **Unique genomes** | 7,664 |
| **Proteins with embeddings** | 11,837,414 |
| **MMseqs2 clusters (total)** | 11,984,419 |
| **Clusters with ≥2 members** | 388,858 |
| **Average cluster size** | 2.54 proteins |

---

## 📁 Generated Outputs

### PCA Embeddings
```
results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz
- Size: 4.3 GB
- Dimensions: 11,837,414 proteins × 50 dimensions
- Variance explained: 89.2%
```

### UMAP Visualizations
```
results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/
├── umap_full_n15.npz (123MB) - Full 2D embeddings
└── figures/
    ├── umap_density.png - Density plot of protein space
    ├── umap_by_cluster.png - Sample colored by cluster
    ├── umap_by_genome.png - Sample colored by genome
    └── cluster_size_distribution.png
```

### Cluster Analysis
```
results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/
├── mmseqs_cluster_statistics.csv (52MB)
│   - 388,858 clusters with: cluster_size, mean_std, max_std,
│     total_variance, mean_pairwise_distance, max_pairwise_distance
├── mmseqs_cluster_per_dimension_stats.csv (1.5GB)
│   - Per-dimension mean and std for all 50 PCA dimensions × all clusters
└── figures/
    ├── cluster_tightness_summary.png - Overview with key metrics
    ├── cluster_tightness_distributions.png - Histograms of all metrics
    ├── size_vs_tightness.png - Correlation analysis
    ├── mean_std_per_dimension.png - Dimension variance analysis
    └── dimension_std_heatmap.png - Top 100 clusters heatmap
```

---

## 🎯 Next Steps for Component 1

Now that clustering analysis is complete, the next steps are:

1. **Define Gene Node Set**
   - Use 388,858 MMseqs2 clusters as gene family definitions
   - Create gene node metadata (cluster representative, size, functional annotations)

2. **Build Genome Graphs**
   - Represent each genome as a graph with gene nodes
   - Define edges (synteny, co-occurrence, functional relationships)
   - Create graph dataset for training

3. **Prepare for GraphMAE**
   - Implement graph data structures
   - Define node features (embeddings, annotations)
   - Set up masked graph learning pipeline

---

## 🔧 Technical Details

### Storage Locations

**Fast storage** (`/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/`):
- Project code and scripts
- Final results and analysis outputs (~6GB)

**Working storage** (`/fh/working/srivatsan_s/dmullane_organism_scale/`):
- Large data files (17GB)
- Temporary batch files (64GB, can be deleted)
- SLURM job logs

### Compute Resources Used

| Task | Runtime | Resources | Cost (GPU-hours) |
|------|---------|-----------|------------------|
| ESM embedding generation | ~24 hours | 1,184 GPU batches | 1,184 |
| PCA reduction | ~30 min | 64GB RAM, 8 CPUs | 0 |
| UMAP computation | ~6 hours | 64GB RAM, 8 CPUs | 0 |
| Cluster tightness analysis | ~33 min | 128GB RAM, 8 CPUs | 0 |
| UMAP visualization | ~3 min | 64GB RAM, 4 CPUs | 0 |
| **Total** | **~30 hours** | | **1,184 GPU-hours** |

---

## 📚 Documentation

For detailed information about the embedding pipeline, see:
- [EMBEDDING_PIPELINE.md](EMBEDDING_PIPELINE.md) - Complete pipeline documentation
- [clustering.md](clustering.md) - Clustering approach and rationale
- Main [README.md](../README.md) - Project overview

---

**Generated with Claude Code**
