# Embeddings Scripts

This directory contains scripts for the RefSeq genome-scale protein clustering pipeline.

## Documentation

**📚 Complete documentation has been moved to the `docs/` directory:**

- **[docs/README.md](../../docs/README.md)** - Master pipeline overview
- **[docs/evaluation.md](../../docs/evaluation.md)** - Multi-metric evaluation framework
- **[docs/umap_array.md](../../docs/umap_array.md)** - UMAP parallelization guide
- **[docs/clustering.md](../../docs/clustering.md)** - Clustering usage examples

## Quick Reference

### Dimensionality Reduction
- `compute_umap_array.py` - UMAP with PCA preprocessing
- `submit_umap_array.sh` - SLURM array for parallel UMAP
- `compute_umap_cogonly.py` - UMAP for COG-annotated genes only

### Clustering
- `cluster_leiden_comprehensive.py` - Leiden clustering on PCA space
- `submit_leiden_sweep.sh` - Parameter sweep (res 5-500)
- `submit_leiden_high_res.sh` - High resolutions (750, 1000)
- `submit_leiden_res1500.sh` - Resolution 1500 testing

### Evaluation
- `evaluate_leiden_sweep.py` - COG homogeneity evaluation
- `evaluate_clustering_quality.py` - ARI, AMI, silhouette, DB index
- `evaluate_clustering_stability.py` - Co-clustering across subsamples
- `compare_clustering_metrics.py` - Comprehensive comparison

### Visualization
- `plot_leiden_clustering.py` - Single clustering plot (clusters + COG)
- `submit_plot_leiden_array.sh` - Array job for all plots
- `visualize_clustering_stability.py` - Stability analysis visualizations

### Diversity Sampling
- `diversity_subsample.py` - Three strategies for diverse 29M → 1M sampling
- `compare_sampling_strategies.py` - Evaluate sampling quality
- `submit_diversity_subsample.sh` - SLURM submission

## Usage

For detailed usage instructions, see the documentation in `docs/`.

For script-specific help:
```bash
python <script_name>.py --help
```
