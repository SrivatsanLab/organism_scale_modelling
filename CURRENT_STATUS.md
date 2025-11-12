# Current Clustering Pipeline Status

**Updated**: 2025-11-07 (MSA-Based Approach Test)
**Status**: 🔄 Testing MSA + balanced sampling approach

---

## 🔬 MSA-Based Clustering Test (In Progress)

### Background: ESM Stability Results
**Previous approach**: Random sampling → ESM embeddings → Leiden clustering
**Result**: ❌ **Poor stability across all configs** (ARI 0.13-0.32)
- 47/50 configurations tested (Job Array 41488384)
- Highest ARI: 0.32 (res=300, nn=15)
- All configs below "moderate" threshold (ARI < 0.4)
- **Conclusion**: Random sampling introduces too much bias OR ESM-C doesn't capture needed structure

### New Approach: Sequence Clustering → Balanced Sampling → ESM Clustering
**Hypothesis**: Balanced sampling from sequence-based protein families will improve stability

### Current Job: MMseqs2 Test (Job 41552086)
**Status**: 🔄 Running
**Script**: `scripts/analysis/submit_mmseqs_test.sh`
**Runtime**: ~30 minutes
**Expected completion**: ~1 hour

**Test Data**:
- 50 genomes downloaded from NCBI (160,840 proteins)
- File: `data/test_proteins_50genomes.faa` (61 MB)
- Coverage: 5,526 genes overlap with PCA cache (0.55% of 1M subsample)

**MMseqs2 Parameters**:
- Min sequence identity: 0.5 (50%)
- Coverage: 0.8 (80%)
- Cluster mode: 0 (greedy)
- Expected output: ~5,000-50,000 sequence clusters

### Test Pipeline Steps

1. ✅ **Download test genomes** - Complete (160K proteins from 50 genomes)
2. 🔄 **MMseqs2 sequence clustering** - Running (Job 41552086)
3. ⏸️ **Balanced sampling** - Pending (after MMseqs2 completes)
4. ⏸️ **ESM clustering on balanced sample** - Pending
5. ⏸️ **Stability evaluation** - Pending (compare to random sampling baseline)

---

## Why This Might Work Better

### Problem with Random Sampling
- Random sampling may over-represent abundant protein families
- Under-represents rare but functionally important proteins
- Creates unstable clusters due to sampling bias

### Balanced Sampling Solution
1. **Sequence clustering first**: Group proteins by homology (MMseqs2)
2. **Balanced sampling**: Sample proportionally from each homology group
3. **ESM clustering second**: Use ESM embeddings on balanced sample
4. **Expected benefit**: Better representation → more stable functional clusters

### Success Criteria
- **ARI > 0.4**: At least "moderate" stability (vs 0.13-0.32 for random)
- **Improvement > 50%**: Clear benefit over random sampling
- **Gene stability > 0.7**: Most genes consistently assigned

---

## Next Steps (After MMseqs2 Completes)

### 1. Check MMseqs2 Results (~30 minutes from now)

```bash
# Check job status
tail logs/mmseqs_test_41552086.out

# View cluster summary
head data/mmseqs_test/cluster_summary.csv
```

**Expected output**:
- Number of sequence clusters (5,000-50,000)
- Cluster size distribution
- Average proteins per cluster

### 2. Create Balanced Sample

```bash
python scripts/analysis/balanced_sample_from_clusters.py \
    --clusters data/mmseqs_test/cluster_summary_assignments.csv \
    --pca-cache results/umap/pca_cache.npz \
    --output data/balanced_sample_gene_ids.txt \
    --n-samples 5000 \
    --strategy sqrt
```

**Expected**:
- Sample ~5,000 genes balanced across sequence clusters
- Genes selected from PCA cache (have ESM embeddings)

### 3. Run Clustering on Balanced Sample

Adapt existing Leiden clustering scripts to use balanced sample instead of random sample

### 4. Evaluate Stability

Run stability evaluation and compare ARI to random sampling baseline (0.13-0.32)

---

## Expected Outcomes

### Scenario 1: Balanced Sampling Works (ARI > 0.4)
✅ **Scale to full dataset**
1. Download all 7,664 genomes (or use full 29M gene sequences if available)
2. Run MMseqs2 on full 29M proteins
3. Create balanced 1M sample from sequence clusters
4. Run full clustering pipeline on balanced sample
5. Assign remaining 28M genes to clusters

### Scenario 2: Marginal Improvement (ARI 0.3-0.4)
⚠️ **Consider hybrid approach**
- Use sequence clusters as coarse grouping
- Within each sequence cluster, use ESM for fine-grained functional clustering
- May need to adjust MMseqs2 stringency (higher/lower identity threshold)

### Scenario 3: No Improvement (ARI < 0.3)
❌ **Problem is with ESM-C embeddings, not sampling**
- Consider gLM2 embeddings (include genomic context: 5' UTR, promoters)
- User's coworker had encouraging results with gLM2
- Would need to generate gLM2 embeddings for all 29M genes

---

## Key Metrics to Watch

### 1. **ARI (Adjusted Rand Index)**
Measures clustering agreement between independent subsamples
- **> 0.8**: Excellent - clusters highly reproducible
- **0.6-0.8**: Good - clusters generally stable
- **0.4-0.6**: Moderate - some variability
- **< 0.4**: Poor - not reproducible

### 2. **Gene Stability Mean**
Fraction of genes assigned to same cluster across samples
- **> 0.9**: Most genes have consistent assignments
- **0.7-0.9**: Good consistency
- **< 0.7**: Many genes have uncertain placements

### 3. **% Genes with Stability > 0.9**
What fraction of genes are highly stable?
- **> 70%**: Most genes confidently assigned
- **50-70%**: Moderate confidence
- **< 50%**: Many uncertain assignments

### 4. **Cluster Count Stability (Std)**
How much does # clusters vary across subsamples?
- **Low std**: Consistent cluster counts
- **High std**: Variable (less stable)

---

## Monitoring Progress

```bash
# Check how many jobs are running
squeue -u dmullane | grep 41488384 | wc -l

# Watch a specific job
tail -f logs/stability_eff_41488384_0.out  # res=300, nn=15, COG-only
tail -f logs/stability_eff_41488384_24.out # res=750, nn=15, COG-only
tail -f logs/stability_eff_41488384_48.out # res=1500, nn=15, COG-only

# Check which ones have finished
ls results/clustering/stability/stability_eff_*.npz | wc -l

# When all done (should have 50 files)
ls results/clustering/stability/ | wc -l
```

---

## Understanding Low ARI vs COG (Not a Problem!)

You correctly noted: **24 COG categories** but **10,000-25,000 clusters**

**This means**:
- Each COG category → hundreds of sub-clusters (ortholog groups)
- Example: COG category "J" (translation) → L1, L2, L3, ..., L12, S1, S2, ... (100+ ribosomal proteins)
- **ARI vs COG will be low** - this is expected!
- Clusters represent **orthology within functional categories**, not categories themselves

**What matters for your use case**:
1. ✅ **Homogeneity**: 61.7% (clusters are functionally coherent)
2. 🔄 **Stability**: Testing now (reproducibility across samples)
3. ✅ **Genes without COG**: Fine to cluster with known genes

**Hierarchy**:
```
COG Category (e.g., "J - Translation")
  └─ Your Cluster 1234 (e.g., "50S ribosomal protein L7/L12")
      ├─ Gene 1 (E. coli rplL)
      ├─ Gene 2 (S. aureus rplL)
      └─ Gene 3 (Unknown species, unannotated but orthologous)
```

---

## Next Steps After Results

### 1. Generate Summary
```bash
python scripts/embeddings/summarize_stability_results.py
```

### 2. Create Heatmaps (Optional)
```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load all results
stability_files = sorted(Path('results/clustering/stability').glob('*.npz'))

# Build matrix: resolutions × n_neighbors
resolutions = [300, 500, 750, 1000, 1500]
n_neighbors = [15, 25, 50, 100, 200]
ari_matrix = np.zeros((len(resolutions), len(n_neighbors)))

for f in stability_files:
    data = np.load(f)
    if data['cog_only']:  # Focus on COG-only
        res_idx = resolutions.index(int(data['resolution']))
        nn_idx = n_neighbors.index(int(data['n_neighbors']))
        ari_matrix[res_idx, nn_idx] = data['ari_mean']

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(ari_matrix, cmap='RdYlGn', vmin=0, vmax=1)
plt.colorbar(label='ARI')
plt.xticks(range(len(n_neighbors)), n_neighbors)
plt.yticks(range(len(resolutions)), resolutions)
plt.xlabel('n_neighbors')
plt.ylabel('Resolution')
plt.title('Clustering Stability Heatmap (COG-only)')
plt.savefig('results/plots/stability_heatmap.png', dpi=150)
```

### 3. Choose Best Configuration
Based on stability results and your downstream needs:
- **If you need maximum resolution**: Use highest stable resolution
- **If you prioritize stability**: Use resolution with highest ARI
- **If you need balanced**: Use composite score

### 4. Run Full Assignment
Once you know the optimal parameters:
```bash
# Edit submit_full_dataset_assignment.sh to use optimal clustering
# Change: CLUSTERING_FILE to best configuration
sbatch scripts/embeddings/submit_full_dataset_assignment.sh
```

---

## Files & Directories

### Scripts
- ✅ `evaluate_clustering_stability_efficient.py` - Stability evaluation
- ✅ `submit_stability_efficient.sh` - 50-task array job
- ✅ `summarize_stability_results.py` - Results summary
- ✅ `assign_full_dataset_to_clusters.py` - Full 29M assignment

### Output
- `results/clustering/stability/` - Will contain 50 .npz files
- `results/clustering/stability_summary.csv` - Summary table
- `results/clustering/full_dataset_assignments_*.npz` - Final assignments

### Logs
- `logs/stability_eff_41488384_*.out` - Progress logs
- `logs/stability_eff_41488384_*.err` - Error logs

---

## Resource Usage

**Current jobs**:
- 50 tasks × 128GB memory × 16 CPUs
- ~150-300 node-hours total
- All running in parallel (cluster permitting)

**Memory**: No OOM failures - using efficient cluster-level metrics

**I/O**: Moderate - loading PCA cache repeatedly

**Time**: Most should finish in 3-6 hours

---

## Why This Comprehensive Sweep Matters

Testing the full parameter space lets you:

1. **Understand trade-offs**: Resolution vs stability, n_neighbors vs stability
2. **Avoid overfitting**: Ensure clusters generalize beyond your specific 1M sample
3. **Optimize confidence**: Choose parameters that maximize both granularity and reproducibility
4. **Publication ready**: Can show robustness analysis with 50 configurations
5. **Future-proof**: If you re-run with different samples, you know expected stability

**This is the right way to validate clustering** - not just "it looks good" but "it's reproducible and robust."

---

## Summary

- ✅ **50 stability evaluations** running in parallel
- ✅ **Comprehensive parameter sweep**: 5 resolutions × 5 n_neighbors × 2 datasets
- ✅ **Memory efficient**: No OOM failures
- ✅ **Full coverage**: Tests all existing clustering configurations
- ⏸️ **Full dataset assignment**: Pending, will use optimal parameters
- 📊 **Results expected**: ~6-8 hours

**Next action**: Wait for jobs to finish, then run summary script!
