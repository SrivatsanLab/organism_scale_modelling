# UMAP Array Job System

Parallelized UMAP computation using SLURM job arrays to test multiple n_neighbors values simultaneously.

## Overview

This system allows you to efficiently compute UMAP embeddings with different parameters in parallel using SLURM job arrays. Each job in the array:
1. Loads gene embeddings (or uses cached PCA)
2. Computes UMAP with a specific n_neighbors value
3. Saves results independently

## Files

- **`compute_umap_array.py`** - Python script that computes UMAP for one parameter set
- **`submit_umap_array.sh`** - SLURM submission script for job array
- **`analyze_umap_results.py`** - Analysis script to compare results

## Quick Start

### 1. Configure Parameters

Edit `submit_umap_array.sh` to set your parameters:

```bash
N_NEIGHBORS_VALUES=(15 25 50 100 200)  # Values to test
N_PCS=50                                # PCA components
SUBSAMPLE=100000                        # Genes to use (or empty for all)
USE_GPU="--use-gpu"                     # Enable GPU
```

SLURM settings:
```bash
#SBATCH --array=0-4           # Array size matches N_NEIGHBORS_VALUES length
#SBATCH --cpus-per-task=8     # CPU cores per job
#SBATCH --mem=64G             # Memory per job
#SBATCH --time=4:00:00        # Time limit
#SBATCH --gres=gpu:1          # Request GPU
```

### 2. Submit Job Array

```bash
cd /home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling
sbatch scripts/embeddings/submit_umap_array.sh
```

### 3. Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View output logs
tail -f logs/umap_*.out

# Check for errors
tail -f logs/umap_*.err
```

### 4. Analyze Results

After jobs complete:

```bash
python scripts/embeddings/analyze_umap_results.py results/umap/
```

This creates:
- `umap_comparison.png` - Side-by-side UMAP plots
- `umap_metrics.csv` - Quantitative comparison
- `umap_metrics.png` - Metrics visualization

## Workflow Details

### Job Array Execution

```
Job 0 (n_neighbors=15):
  ├─ Load embeddings
  ├─ Compute PCA (save cache)
  └─ Compute UMAP → umap_n15_subsample100000.npz

Job 1 (n_neighbors=25):
  ├─ Wait for PCA cache
  ├─ Load PCA from cache
  └─ Compute UMAP → umap_n25_subsample100000.npz

Job 2 (n_neighbors=50):
  ├─ Wait for PCA cache
  ├─ Load PCA from cache
  └─ Compute UMAP → umap_n50_subsample100000.npz

... (parallel execution)
```

### PCA Caching

- **Job 0** computes PCA and saves to `results/umap/pca_cache.npz`
- **Jobs 1-N** wait for the cache, then load it
- This avoids redundant PCA computation across jobs

### Output Files

Each job creates:
```
results/umap/
├── pca_cache.npz                       # Shared PCA results
├── umap_n15_subsample100000.npz       # UMAP with n_neighbors=15
├── umap_n25_subsample100000.npz       # UMAP with n_neighbors=25
├── umap_n50_subsample100000.npz       # UMAP with n_neighbors=50
├── umap_n100_subsample100000.npz      # UMAP with n_neighbors=100
└── umap_n200_subsample100000.npz      # UMAP with n_neighbors=200
```

## Command-Line Usage

You can also run jobs manually:

```bash
# Compute UMAP with n_neighbors=15
python scripts/embeddings/compute_umap_array.py \
    --embeddings-dir data/esm_embeddings \
    --output results/umap/umap_n15.npz \
    --n-neighbors 15 \
    --n-pcs 50 \
    --subsample 100000 \
    --use-gpu

# Use cached PCA
python scripts/embeddings/compute_umap_array.py \
    --load-cached-pca results/umap/pca_cache.npz \
    --output results/umap/umap_n25.npz \
    --n-neighbors 25 \
    --use-gpu
```

## Options

### `compute_umap_array.py`

```
--embeddings-dir    Directory with embedding files
--output            Output file (.npz)
--n-neighbors       UMAP n_neighbors parameter (required)
--n-pcs             Number of PCA components (default: 50)
--min-dist          UMAP min_dist (default: 0.1)
--subsample         Number of genes to sample (optional)
--seed              Random seed (default: 42)
--use-gpu           Enable GPU acceleration
--load-cached-pca   Load pre-computed PCA from file
--save-pca          Save PCA to file
```

### `analyze_umap_results.py`

```
results_dir         Directory containing UMAP .npz files
--output            Output plot filename (default: umap_comparison.png)
```

## Performance

### With GPU (recommended)

- **PCA**: ~5 minutes (100k genes)
- **UMAP**: ~2-3 minutes per n_neighbors value
- **Total**: ~20 minutes for 5 values (parallel)

### CPU only

- **PCA**: ~10 minutes (100k genes)
- **UMAP**: ~30-60 minutes per n_neighbors value
- **Total**: ~5 hours for 5 values (parallel)

## Tips

1. **Test first**: Run with `--subsample 10000` to test quickly
2. **Use GPU**: Much faster for UMAP computation
3. **Monitor memory**: Adjust `--mem` based on your subsample size
4. **Adjust time limits**: Longer for larger datasets or CPU-only
5. **Check logs**: Look at `logs/umap_*.out` for progress

## Troubleshooting

**Jobs pending?**
```bash
squeue -u $USER
# Check REASON column for issues
```

**Out of memory?**
```bash
# Increase memory in submit script:
#SBATCH --mem=128G
```

**GPU not found?**
```bash
# Check GPU availability:
sinfo -p campus-new
# Or remove --use-gpu flag
```

**PCA cache timeout?**
```bash
# Job 0 might have failed
# Check logs/umap_*_0.err
# May need to manually create cache first
```

## Advanced Usage

### Run Specific n_neighbors Values

Modify the array in `submit_umap_array.sh`:

```bash
# Test only n_neighbors=15,30,60
N_NEIGHBORS_VALUES=(15 30 60)
#SBATCH --array=0-2
```

### Process All Genes (No Subsampling)

```bash
SUBSAMPLE=""  # Empty = use all genes
```

Warning: Requires ~200GB+ RAM for full dataset!

### Different PCA Components

```bash
N_PCS=100  # More PCs = captures more variance but slower
```

### Custom UMAP Parameters

Edit in `submit_umap_array.sh`:
```bash
MIN_DIST=0.5  # Larger = more spread out
```

## Output Data Format

Each `.npz` file contains:

```python
import numpy as np

data = np.load('umap_n15_subsample100000.npz', allow_pickle=True)

# Arrays:
umap_embedding = data['umap_embedding']  # (n_genes, 2) - UMAP coordinates
gene_ids = data['gene_ids']              # (n_genes,) - Gene identifiers
genome_ids = data['genome_ids']          # (n_genes,) - Genome identifiers

# Metadata:
n_neighbors = data['n_neighbors']        # UMAP parameter
min_dist = data['min_dist']              # UMAP parameter
n_pcs = data['n_pcs']                    # Number of PCs used
seed = data['seed']                      # Random seed
```

## Comparison Workflow

```bash
# 1. Submit jobs
sbatch scripts/embeddings/submit_umap_array.sh

# 2. Wait for completion (check with squeue)
watch -n 30 squeue -u $USER

# 3. Analyze results
python scripts/embeddings/analyze_umap_results.py results/umap/

# 4. View comparison
open umap_comparison.png
open results/umap/umap_metrics.png
```

## Integration with Notebook

Load results in Jupyter:

```python
import numpy as np
import matplotlib.pyplot as plt

# Load specific result
data = np.load('results/umap/umap_n15_subsample100000.npz', allow_pickle=True)
umap_emb = data['umap_embedding']

# Plot
plt.figure(figsize=(10, 8))
plt.scatter(umap_emb[:, 0], umap_emb[:, 1], s=1, alpha=0.3)
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title(f'UMAP (n_neighbors={data["n_neighbors"]})')
plt.show()
```
