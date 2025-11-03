# Documentation Index

Complete documentation for the organism-scale modeling project.

## 📘 Main Documentation

### [README.md](README.md)
**Master Pipeline Overview** - RefSeq genome-scale protein clustering pipeline
- Complete workflow from genomes → clustering → evaluation
- Current status and key results
- Quick start guides
- File organization and troubleshooting

## 📗 Pipeline Components

### Data Preparation
- **[README_gene_prediction.md](README_gene_prediction.md)** - Gene prediction with Prodigal
- **[README_embeddings.md](README_embeddings.md)** - ESM-C embedding generation

### Dimensionality Reduction & Visualization
- **[umap_array.md](umap_array.md)** - UMAP parallelization system
  - Job arrays for testing n_neighbors
  - PCA caching strategy
  - GPU vs CPU performance

### Clustering
- **[clustering.md](clustering.md)** - Clustering usage examples
  - Loading results in Jupyter
  - COG annotation integration
  - Visualization templates

### Evaluation
- **[evaluation.md](evaluation.md)** - Multi-metric evaluation framework
  - **Why it matters**: Beyond simple homogeneity
  - ARI/AMI quality metrics
  - Stability analysis (co-clustering across subsamples)
  - Comprehensive comparison workflow
  - Interpretation guidelines

## 🧬 Alternative Models

### [README_glm2.md](README_glm2.md)
Genomic Language Model (gLM2) integration
- DNA sequence modeling
- Comparison with ESM protein embeddings
- Jacobian analysis

### [gLM2_tokenization_guide.md](gLM2_tokenization_guide.md)
Technical guide for gLM2 tokenization

## 🗂️ Quick Navigation

### Getting Started
1. Start with **[README.md](README.md)** for the big picture
2. Check **[evaluation.md](evaluation.md)** to understand metrics
3. Refer to component docs for specific tasks

### For Analysis
- **Loading results**: See [clustering.md](clustering.md)
- **Understanding metrics**: See [evaluation.md](evaluation.md)
- **UMAP parameters**: See [umap_array.md](umap_array.md)

### For Development
- **Running jobs**: Check individual component docs
- **Script reference**: See `scripts/embeddings/README.md`
- **Troubleshooting**: See main [README.md](README.md#troubleshooting)

## 📊 Current Status (2025-11-03)

### ✅ Completed
- 7,664 RefSeq genomes processed
- 29M proteins embedded with ESM-C
- 129 Leiden clustering configurations
- 645 UMAP visualization plots
- COG functional annotations

### 🔄 In Progress
- Quality evaluation (ARI, AMI, silhouette, Davies-Bouldin)
- Stability analysis (co-clustering across subsamples)

### 📋 Next
- Comprehensive metric comparison
- Stability visualization
- Optimal parameter identification

## 🔗 Related Resources

- **Scripts**: `scripts/embeddings/`
- **Results**: `results/clustering/`, `results/umap/`
- **Notebooks**: `notebooks/clustering_sweep_analysis.ipynb`
- **Logs**: `logs/`
