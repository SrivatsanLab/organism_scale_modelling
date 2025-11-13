# ESM Embedding Generation Pipeline

## Overview

This pipeline generates ESM-C embeddings for all proteins in MMseqs2 filtered clusters (70% sequence identity, 10+ members), enabling full-scale clustering analysis and UMAP visualization.

## Dataset Summary

- **Total proteins**: 30,098,843 (across 7,664 genomes)
- **Filtered clusters** (70% ID, 10+ members): 388,858 clusters
- **Proteins in filtered clusters**: 12,246,048 proteins
- **Existing embeddings**: 408,634 proteins (3.34%)
- **Proteins needing embeddings**: 11,837,414 proteins (96.66%)

## Pipeline Components

### 1. Protein Identification
**Script**: `scripts/analysis/identify_proteins_needing_embeddings.py`

Identifies proteins in filtered clusters that don't have embeddings yet.

**Output**: `results/clustering/filtered_0p7/proteins_needing_embeddings.txt` (11.8M protein IDs)

### 2. Sequence Extraction
**Script**: `scripts/analysis/extract_sequences_for_embedding.py`

Extracts FASTA sequences for proteins needing embeddings from the full dataset.

**Input**: `data/all_proteins.faa` (30M proteins, 13 GB)

**Output**: `data/proteins_for_embedding.faa` (~12 GB, 11.8M proteins)

**Runtime**: ~20-30 minutes

### 3. Batch Embedding Generation
**Script**: `scripts/analysis/batch_generate_embeddings.py`
**Submission**: `scripts/analysis/submit_batch_embeddings.sh`

Generates ESM-C embeddings in batches using GPU array jobs.

**Configuration**:
- Batch size: 10,000 proteins per job
- Total batches: 1,184 jobs
- Concurrent jobs: 50 (max)
- Resources per job: 1 GPU, 64GB RAM, 8 CPUs, 12 hours
- Output location: `/fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches/`

**Estimated time**: 24-48 hours (depending on queue)

**Storage**: ~50-100 GB (batch files)

**Submit**:
```bash
sbatch scripts/analysis/submit_batch_embeddings.sh
```

**Monitor**:
```bash
# Check running jobs
squeue -u $USER | grep batch_embeddings

# Count completed batches
ls /fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches/*.npz | wc -l

# Check logs
tail -f /fh/working/srivatsan_s/dmullane_organism_scale/logs/batch_embeddings_*.out
```

### 4. Merge Embeddings
**Script**: `scripts/analysis/merge_embedding_batches.py`

Merges all batch embeddings and combines with existing PCA cache.

**Process**:
1. Loads existing 1M embeddings (50D) from PCA cache
2. Loads all 1,184 batch files with new embeddings
3. Applies PCA to reduce new embeddings to 50D
4. Combines into single file

**Output**: `results/umap/pca_cache_full.npz` (~5-10 GB, 12.2M proteins)

**Runtime**: ~1-2 hours

**Execute**:
```bash
/home/dmullane/micromamba/envs/esm3_env/bin/python \
    scripts/analysis/merge_embedding_batches.py
```

### 5. Full-Scale Clustering Analysis
**Script**: `scripts/analysis/process_filtered_clusters_fast.py`

Recomputes cluster statistics and samples proteins for UMAP using full embedding coverage.

**Output**:
- Cluster statistics (mean/SD per dimension)
- Sampled proteins for UMAP (~3.9M proteins from 388K clusters)

**Execute**:
```bash
/home/dmullane/micromamba/envs/esm3_env/bin/python \
    scripts/analysis/process_filtered_clusters_fast.py \
    --mmseqs-dir data/mmseqs_seqid_0p7 \
    --pca-cache results/umap/pca_cache_full.npz \
    --min-size 10 \
    --n-sample 10 \
    --output-dir results/clustering/filtered_0p7_full
```

### 6. UMAP Visualization
**Script**: `scripts/analysis/run_umap_and_visualize.py`

Generates 2D UMAP embeddings and visualizations for all 388K clusters.

**Output**:
- `umap_by_mmseqs_cluster.png` - Colored by 388K clusters
- `umap_by_cog_category.png` - Colored by functional categories
- `umap_combined.png` - Side-by-side comparison

**Runtime**: ~30-60 minutes (for ~3.9M proteins)

## Resource Summary

| Stage | Runtime | Storage | GPU Required |
|-------|---------|---------|--------------|
| Sequence extraction | 20-30 min | ~12 GB | No |
| Batch embedding generation | 24-48 hrs | ~50-100 GB | Yes (1,184 GPU-hours) |
| Merge embeddings | 1-2 hrs | ~10 GB | No |
| Clustering analysis | 10-20 min | ~1 GB | No |
| UMAP visualization | 30-60 min | <1 GB | No |

**Total storage**: ~75-125 GB (working directory)
**Total time**: ~26-51 hours (mostly embedding generation)

## Directory Structure

```
/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/
├── data/
│   ├── all_proteins.faa                    # 30M proteins (13 GB)
│   └── proteins_for_embedding.faa          # 11.8M proteins (12 GB)
├── results/
│   ├── umap/
│   │   ├── pca_cache.npz                   # Original 1M embeddings
│   │   └── pca_cache_full.npz              # Full 12.2M embeddings
│   └── clustering/
│       ├── filtered_0p7/                    # Current results (3.7K clusters)
│       └── filtered_0p7_full/               # Full results (388K clusters)
└── /fh/working/srivatsan_s/dmullane_organism_scale/
    ├── embeddings/
    │   └── batches/                         # 1,184 batch files (~50-100 GB)
    └── logs/                                # SLURM job logs
```

## Troubleshooting

**If embedding generation jobs fail:**
1. Check logs: `tail /fh/working/srivatsan_s/dmullane_organism_scale/logs/batch_embeddings_*.err`
2. Identify failed batches: Compare expected (1,184) vs actual batch files
3. Resubmit failed batches using specific array indices:
   ```bash
   sbatch --array=<failed_indices> scripts/analysis/submit_batch_embeddings.sh
   ```

**If running out of storage:**
- Working directory (`/fh/working/srivatsan_s/`) is designed for large temporary files
- Batch files can be deleted after merging
- Final merged cache (`pca_cache_full.npz`) is the only file needed for downstream analysis

**If jobs are queued for too long:**
- Reduce concurrent job limit in submission script (default: 50)
- Check cluster status: `sinfo -p gpu`
- Consider using different GPU partition if available

## Next Steps After Completion

Once the full embedding cache is generated:

1. **Rerun clustering evaluation** with full coverage
2. **Generate COG annotations** for the 12.2M proteins
3. **Create comprehensive visualizations** showing:
   - All 388K clusters in UMAP space
   - Functional category distributions
   - Cluster quality metrics
4. **Compare** with previous results (3.7K clusters with embeddings)

## Status Update - 2025-11-12

### Completed Tasks

✅ **All 1,184 embedding batches completed successfully**
- Total proteins: 11,837,414
- Total size: 64 GB (batch files)
- All batches verified and complete

✅ **Embeddings merged into PCA cache**
- Output: `results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz`
- Size: 4.3 GB
- Dimensions: 11.8M proteins × 50D
- PCA variance explained: **89.2%** (excellent retention)

✅ **Storage optimization**
- Moved `data/` directory (17GB) from `/fast` to `/working` storage
- Created symlink to maintain compatibility
- Saved 17GB on expensive storage

### Currently Running

🔄 **UMAP computation** (Job 41771479)
- Computing 2D UMAP embeddings for all 11.8M proteins
- Parameters: n_neighbors=15, min_dist=0.1
- Input: 50D PCA embeddings (4.3GB)
- Estimated completion: ~1-2 hours
- Output: `results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_full_n15.npz`

🔄 **MMseqs2 cluster tightness analysis** (Job 41772028)
- Analyzing embedding space tightness of 70% sequence identity clusters
- Processing 30M cluster assignments
- Computing per-cluster statistics:
  - Mean & std dev across all 50 PCA dimensions
  - Pairwise distances within clusters
  - Overall cluster variance metrics
- Outputs:
  - `mmseqs_cluster_statistics.csv` - Overall metrics per cluster
  - `mmseqs_cluster_per_dimension_stats.csv` - Per-dimension stats
  - Visualization plots in `cluster_analysis/figures/`

### Purpose of Cluster Tightness Analysis

This analysis determines whether MMseqs2 clusters (70% sequence identity) are already tight enough in ESM embedding space that we don't need additional Leiden clustering. If clusters show low variance in embedding space, we can use MMseqs2 clusters directly for downstream analysis.

### Next Steps

Once current jobs complete:
1. Review cluster tightness results to decide on clustering strategy
2. Generate UMAP visualizations colored by genome and cluster
3. If needed, perform Leiden clustering on embeddings
4. Create final comprehensive visualizations

## Generated By

Pipeline created: 2025-11-10
Last updated: 2025-11-12
Generated with Claude Code
