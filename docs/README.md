# RefSeq Genome-Scale Protein Clustering Pipeline

Comprehensive workflow for analyzing 7,664 RefSeq bacterial genomes (~29M proteins) using ESM-C embeddings, dimensionality reduction, and functional clustering.

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Current Status](#current-status)
3. [Complete Workflow](#complete-workflow)
4. [Key Results](#key-results)
5. [Documentation Index](#documentation-index)
6. [Quick Start Guide](#quick-start-guide)

---

## Pipeline Overview

```
RefSeq Genomes (7,664)
    ↓
Gene Prediction (Prodigal)
    ↓
ESM-C Embeddings (600M model, 1152D)
    ↓ 29M proteins
PCA (50D, 91.8% variance)
    ↓ subsample: 1M proteins
UMAP (2D visualization, n_neighbors: 15/25/50/100/200)
    ↓
Leiden Clustering (resolution: 5-1500)
    ↓
COG Functional Annotation
    ↓
Comprehensive Evaluation (ARI, AMI, Stability)
```

---

## Current Status

### ✅ Completed

1. **Gene Prediction**: All 7,664 genomes processed with Prodigal → ~29M proteins
2. **ESM-C Embeddings**: All proteins embedded (1152D) using ESM-C 600M model
3. **PCA Preprocessing**: 50D PCA (91.8% variance) for efficient clustering
4. **UMAP Embeddings**:
   - 1M gene subsample
   - 5 n_neighbors values: 15, 25, 50, 100, 200
   - COG-only variants (798K annotated genes)
5. **COG Annotation**: DIAMOND search against COG database for functional labels
6. **Leiden Clustering Sweep**:
   - 129 configurations total
   - Resolutions: 5, 10, 20, 50, 100, 200, 300, 400, 500, 750, 1000, 1500
   - n_neighbors: 15, 25, 50, 100, 200
   - Both all genes and COG-only variants
7. **Visualization**:
   - 645 UMAP plots (clustering × UMAP combinations)
   - COG category annotations on all plots

### 🔄 In Progress

1. **Quality Evaluation** (Job 41274292):
   - Computing ARI, AMI, silhouette, Davies-Bouldin for all 129 clusterings
   - Runtime: ~4 hours
   - Output: `results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/quality_metrics_comprehensive.csv`

2. **Stability Evaluation** (Job 41274295):
   - Testing clustering robustness across independent subsamples
   - 4 configurations: res1500/1000/750 with n=15
   - Runtime: ~12 hours
   - Output: `results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability/stability_*.npz`

### 📋 Next Steps

1. Run visualization and comparison scripts after evaluations complete:
   ```bash
   python scripts/embeddings/visualize_clustering_stability.py
   python scripts/embeddings/compare_clustering_metrics.py
   ```

2. Identify optimal clustering parameters balancing:
   - Biological relevance (ARI/AMI vs COG)
   - Robustness (stability across subsamples)
   - Geometric quality (silhouette, Davies-Bouldin)

3. Scale best configuration to full 29M protein dataset

---

## Complete Workflow

### Phase 1: Data Preparation (Completed)

#### 1.1 Gene Prediction
```bash
# Predict genes for all RefSeq genomes
scripts/preprocessing/predict_genes_refseq.sh
```
- Input: 7,664 RefSeq genome FASTA files
- Output: ~29M protein sequences
- Tool: Prodigal (optimized for bacteria)

#### 1.2 ESM-C Embedding Generation
```bash
# Generate embeddings for all proteins
scripts/embeddings/generate_esm_embeddings.sh
```
- Input: Protein FASTA files
- Output: 1152D embeddings per protein
- Model: ESM-C 600M (evolutionary scale modeling)
- Storage: ~500GB total

### Phase 2: Dimensionality Reduction (Completed)

#### 2.1 PCA Preprocessing
```bash
# Computed automatically by UMAP/clustering scripts
# Output: results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz
```
- 1152D → 50D (preserves 91.8% variance)
- Enables efficient clustering on large dataset
- Cached to avoid recomputation

#### 2.2 UMAP Visualization
```bash
# Array job for multiple n_neighbors values
sbatch scripts/embeddings/submit_umap_array.sh
```
- Subsample: 1M genes (representative)
- Parameters tested: n_neighbors ∈ {15, 25, 50, 100, 200}
- Output: 2D coordinates for visualization
- GPU-accelerated: ~3 min per configuration
- See: [umap_array.md](umap_array.md)

### Phase 3: Functional Annotation (Completed)

#### 3.1 COG Annotation via DIAMOND
```bash
# BLAST-like search against COG database
scripts/functional_annotation/run_diamond_cog.sh
```
- Tool: DIAMOND (fast protein alignment)
- Database: COG-2020 (Clusters of Orthologous Genes)
- Coverage: 798K / 1M genes (79.8%)
- 24 functional categories (J, K, L, C, E, G, M, ...)

### Phase 4: Clustering Parameter Sweep (Completed)

#### 4.1 Leiden Clustering Sweep
```bash
# Comprehensive parameter sweep
sbatch scripts/embeddings/submit_leiden_sweep.sh      # res 5-500
sbatch scripts/embeddings/submit_leiden_high_res.sh    # res 750,1000
sbatch scripts/embeddings/submit_leiden_res1500.sh     # res 1500
```

**Parameters tested**:
- Resolutions: 5, 10, 20, 50, 100, 200, 300, 400, 500, 750, 1000, 1500
- n_neighbors: 15, 25, 50, 100, 200
- Variants: All genes vs COG-annotated only
- **Total**: 129 configurations

**Algorithm**: Leiden (graph-based community detection)
- Operates on 50D PCA space
- Constructs k-NN graph (cosine similarity)
- Optimizes modularity with resolution parameter
- Higher resolution → more, smaller clusters

#### 4.2 Visualization Generation
```bash
# Generate UMAP plots for all clustering results
sbatch scripts/embeddings/submit_plot_leiden_array.sh
```
- Creates 645 plots (129 clusterings × 5 UMAPs)
- Each plot: clusters (left) + COG categories (right)
- Output: `results/plots/umap_n*/leiden_*.png`

### Phase 5: Clustering Evaluation (In Progress)

#### 5.1 Quality Metrics
```bash
sbatch scripts/embeddings/submit_quality_evaluation.sh
```

**Metrics computed**:
1. **Adjusted Rand Index (ARI)**: Similarity to COG annotations (0-1)
2. **Adjusted Mutual Information (AMI)**: Information-theoretic similarity (0-1)
3. **Silhouette Score**: Cluster separation quality (-1 to 1)
4. **Davies-Bouldin Index**: Cluster compactness (lower better)
5. **Size Distribution**: Gini coefficient, percentiles

**Why these matter**: Tests if clusters align with known biology, not just internal consistency.

See: [evaluation.md](evaluation.md)

#### 5.2 Stability Analysis
```bash
sbatch scripts/embeddings/submit_stability_evaluation.sh
```

**Approach**:
1. Generate 10 independent 100K subsamples
2. Cluster each with same parameters
3. Track gene pairs appearing in multiple subsamples
4. Compute **co-clustering rate**: (times together) / (times co-occurred)

**Why this matters**: Robust, meaningful clusters should persist across random samples. Prevents overfitting to specific subsample.

**Key insight** (from user): Simply increasing resolution trivially increases homogeneity by making smaller clusters. Stability tests true robustness.

#### 5.3 Comprehensive Comparison
```bash
# After evaluations complete
python scripts/embeddings/compare_clustering_metrics.py
```

Produces:
- `all_metrics_combined.csv`: All metrics merged
- Rankings by composite scores (quality-focused, stability-focused, balanced)
- Plots showing metrics vs resolution
- Identification of optimal parameter sweet spot

#### 5.4 Stability Visualization
```bash
python scripts/embeddings/visualize_clustering_stability.py
```

Creates:
- Co-clustering rate distributions
- Stability comparison across configs
- Network visualization of stable gene pairs
- COG agreement analysis
- Lists of stable/unstable pairs

---

## Key Results

### Best Clustering (Homogeneity-Based)

**Configuration**: Resolution 1500, n_neighbors=15, COG-only

**Performance**:
- **61.7% weighted COG homogeneity** (highest achieved)
- 23,591 clusters from 798K annotated genes
- 28% of clusters have ≥80% homogeneity
- 6.6% perfect clusters (100% homogeneity)
- ~270K genes in highly homogeneous clusters

**Trend**: Homogeneity continues increasing with resolution, but diminishing returns suggest approaching plateau.

**Caveat**: High homogeneity may be artifact of smaller clusters. Awaiting quality/stability metrics to confirm.

### Dataset Statistics

- **Total genomes**: 7,664 RefSeq bacterial genomes
- **Total proteins**: ~29 million
- **Subsample for UMAP/clustering**: 1 million (representative)
- **COG annotation rate**: 79.8% (798K / 1M)
- **PCA variance captured**: 91.8% (1152D → 50D)
- **Embedding dimension**: 1152D (ESM-C 600M)

### COG Category Distribution

Top categories in dataset:
- **J**: Translation, ribosomal structure
- **K**: Transcription
- **L**: Replication, recombination, repair
- **C**: Energy production and conversion
- **E**: Amino acid transport and metabolism
- **G**: Carbohydrate transport and metabolism

---

## Documentation Index

### Core Workflows
- **[README.md](README.md)** (this file): Master overview
- **[umap_array.md](umap_array.md)**: UMAP parallelization system
- **[clustering.md](clustering.md)**: Old clustering guide (pre-Leiden sweep)
- **[evaluation.md](evaluation.md)**: New evaluation framework

### Script Categories

#### Embedding Generation
- `generate_esm_embeddings.sh`: Batch ESM-C embedding generation
- `predict_genes_refseq.sh`: Prodigal gene prediction

#### Dimensionality Reduction
- `compute_umap_array.py`: UMAP with PCA preprocessing
- `submit_umap_array.sh`: SLURM array for parallel UMAP
- `compute_umap_cogonly.py`: UMAP for COG-annotated genes only

#### Clustering
- `cluster_leiden_comprehensive.py`: Leiden clustering on PCA space
- `submit_leiden_sweep.sh`: Parameter sweep (res 5-500)
- `submit_leiden_high_res.sh`: High resolutions (750, 1000)
- `submit_leiden_res1500.sh`: Resolution 1500 testing

#### Evaluation
- `evaluate_leiden_sweep.py`: COG homogeneity evaluation
- `evaluate_clustering_quality.py`: ARI, AMI, silhouette, DB index
- `evaluate_clustering_stability.py`: Co-clustering across subsamples
- `compare_clustering_metrics.py`: Comprehensive comparison

#### Visualization
- `plot_leiden_clustering.py`: Single clustering plot (clusters + COG)
- `submit_plot_leiden_array.sh`: Array job for all plots
- `visualize_clustering_stability.py`: Stability analysis visualizations

#### Diversity Sampling (Ready, Not Run)
- `diversity_subsample.py`: Three strategies for diverse 29M → 1M sampling
- `compare_sampling_strategies.py`: Evaluate sampling quality
- `submit_diversity_subsample.sh`: SLURM submission

#### Functional Annotation
- `run_diamond_cog.sh`: DIAMOND search against COG database
- COG database files in `data/cog/`

---

## Quick Start Guide

### For Analysis (Using Existing Results)

#### 1. Load UMAP and Clustering in Jupyter

```python
import numpy as np
import matplotlib.pyplot as plt

# Load UMAP (1M genes, n_neighbors=15)
umap_data = np.load('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample1000000.npz', allow_pickle=True)
umap_coords = umap_data['umap_embedding']
gene_ids = umap_data['gene_ids']

# Load best clustering (res=1500, n=15, COG-only)
clust_data = np.load('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/clusters_leiden_res1500_nn15_cogonly.npz',
                     allow_pickle=True)
cluster_labels = clust_data['labels']

# Plot
plt.figure(figsize=(12, 10))
plt.scatter(umap_coords[:, 0], umap_coords[:, 1],
           c=cluster_labels, s=1, alpha=0.5, cmap='tab20')
plt.title('Leiden Clustering (res=1500, n=15, COG-only)')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.colorbar(label='Cluster')
plt.show()
```

#### 2. Check Evaluation Results

```bash
# View homogeneity results
cat results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/leiden_sweep_summary.csv | column -t -s','

# When quality eval completes:
cat results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/quality_metrics_comprehensive.csv | column -t -s','

# Compare metrics
python scripts/embeddings/compare_clustering_metrics.py
```

#### 3. Explore Visualizations

```bash
# Browse UMAP plots
ls results/plots/umap_n15/leiden_*.png

# View specific clustering
open results/plots/umap_n15/leiden_res1500_nn15_cogonly.png

# After stability completes:
python scripts/embeddings/visualize_clustering_stability.py
ls results/plots/stability_analysis/
```

### For New Analyses

#### Run Clustering on Different Parameters

```bash
# Edit parameters in submit_leiden_sweep.sh, then:
sbatch scripts/embeddings/submit_leiden_sweep.sh
```

#### Test Different UMAP Settings

```bash
# Edit N_NEIGHBORS_VALUES in submit_umap_array.sh, then:
sbatch scripts/embeddings/submit_umap_array.sh
```

#### Compute Diversity-Based Subsample

```bash
# For scaling to 29M genes with diverse sampling:
sbatch scripts/embeddings/submit_diversity_subsample.sh
```

---

## File Organization

```
organism_scale_modelling/
├── data/
│   ├── refseq_genomes/           # 7,664 genome FASTA files
│   ├── refseq_proteins/          # Predicted protein sequences
│   ├── refseq_esm_embeddings/    # ESM-C embeddings (1152D)
│   ├── refseq_cog_annotations/   # DIAMOND COG search results
│   └── cog/                      # COG database files
│
├── results/
│   ├── umap/
│   │   ├── pca_cache.npz         # PCA preprocessing (50D)
│   │   ├── umap_n15_*.npz        # UMAP embeddings
│   │   └── ...
│   ├── clustering/
│   │   ├── clusters_leiden_*.npz # Clustering results (129 files)
│   │   ├── leiden_sweep_summary.csv
│   │   ├── quality_metrics_comprehensive.csv  # In progress
│   │   ├── all_metrics_combined.csv           # After comparison
│   │   └── stability/            # Stability results
│   └── plots/
│       ├── umap_n15/             # 129 plots for n=15
│       ├── umap_n25/             # 129 plots for n=25
│       ├── ...
│       ├── metric_comparison/    # After comparison
│       └── stability_analysis/   # After visualization
│
├── scripts/
│   └── embeddings/
│       ├── README.md             # This file
│       ├── evaluation.md
│       ├── umap_array.md
│       ├── clustering.md
│       └── *.py, *.sh            # Pipeline scripts
│
├── notebooks/
│   └── clustering_sweep_analysis.ipynb  # Interactive analysis
│
└── logs/                         # SLURM output/error logs
```

---

## Computing Resources

### Completed Jobs

- **Gene prediction**: 7,664 × 5 min = ~27 CPU-hours
- **ESM embeddings**: 7,664 × 30 min = ~160 GPU-hours (A100)
- **UMAP array**: 5 jobs × 3 min = 15 GPU-minutes (parallelized)
- **Leiden sweep**: 129 jobs × 10 min = ~22 CPU-hours (parallelized)
- **Plotting array**: 645 jobs × 2 min = ~22 CPU-hours (parallelized)

### In Progress

- **Quality evaluation**: 1 job × 4 hours = 4 CPU-hours
- **Stability evaluation**: 4 jobs × 12 hours = 48 CPU-hours (parallelized)

### Typical Resource Requirements

- **UMAP (GPU)**: 8 CPUs, 64GB RAM, 1 GPU, ~3 min
- **UMAP (CPU)**: 8 CPUs, 64GB RAM, ~45 min
- **Leiden clustering**: 8 CPUs, 64GB RAM, ~10 min
- **Quality evaluation**: 8 CPUs, 64GB RAM, ~4 hours
- **Stability evaluation**: 16 CPUs, 128GB RAM, ~12 hours

---

## Key Findings & Insights

### 1. PCA Preprocessing is Critical

- Direct UMAP on 1152D: very slow, memory-intensive
- 50D PCA (91.8% variance): 10× faster, minimal quality loss
- Enables clustering on larger datasets

### 2. Leiden Outperforms Other Methods

Early testing showed:
- **HDBSCAN**: Too many noise points (~50%)
- **K-means**: Requires pre-specifying k, less flexible
- **Leiden**: Balanced, good COG enrichment, flexible resolution

### 3. COG-Only vs All Genes

COG-only clustering consistently shows:
- **Higher homogeneity** (fewer unannotated genes to confuse signal)
- **Better biological interpretability**
- **Trade-off**: Excludes 20% of genes (novel functions?)

### 4. Resolution Selection Challenge

**Problem**: Homogeneity increases linearly with resolution
- res=500: 54.9% homogeneity, 10.9K clusters
- res=1000: 59.2% homogeneity, 17.8K clusters
- res=1500: 61.7% homogeneity, 23.6K clusters

**User insight**: This is trivial - smaller clusters are more homogeneous!

**Solution**: Multi-metric evaluation (ARI, AMI, stability) to find true optimum

### 5. Evaluation Framework Necessity

Need metrics that assess:
- **Biological relevance**: Do clusters match known functions? (ARI/AMI)
- **Geometric quality**: Are clusters well-separated? (silhouette/DB)
- **Robustness**: Do clusters persist across samples? (stability)
- **Balance**: Composite scores prevent gaming single metric

---

## Troubleshooting

### Jobs Not Running

```bash
# Check status
squeue -u $USER

# Check partition availability
sinfo -p campus-new

# View detailed job info
scontrol show job JOBID
```

### Out of Memory

```bash
# Increase memory in submit script
#SBATCH --mem=256G

# Or reduce subsample size
--subsample 500000
```

### Missing Files

```bash
# Check if PCA cache exists
ls results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz

# Regenerate if missing
python scripts/embeddings/compute_umap_array.py --save-pca results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz ...
```

### GPU Issues

```bash
# Check GPU availability
sinfo -p campus-new --Format=Partition,Gres,Nodes,NodeList

# Or disable GPU
# Remove --use-gpu flag and --gres=gpu:1
```

---

## Citation & Acknowledgments

### Tools Used

- **ESM** (Evolutionary Scale Modeling): Lin et al., 2023
- **UMAP**: McInnes et al., 2018
- **Leiden**: Traag et al., 2019
- **DIAMOND**: Buchfink et al., 2015
- **COG**: Galperin et al., 2021
- **Prodigal**: Hyatt et al., 2010

### Compute Resources

- Fred Hutch Scientific Computing
- SLURM job scheduling
- NVIDIA A100 GPUs

---

## Contact & Support

For questions or issues:
1. Check documentation in `scripts/embeddings/`
2. Review SLURM logs in `logs/`
3. Consult individual script help: `python script.py --help`

---

## Version History

### Current (2025-11-03)

- Comprehensive Leiden parameter sweep (129 configs)
- Multi-metric evaluation framework (quality + stability)
- Visualization pipeline for all results
- Resolution 1500 testing (best homogeneity so far)

### Previous

- UMAP array system with GPU acceleration
- COG functional annotation
- Basic clustering methods (HDBSCAN, K-means)
- PCA preprocessing implementation
