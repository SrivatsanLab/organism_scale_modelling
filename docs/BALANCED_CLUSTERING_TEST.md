# Balanced Clustering Stability Test

**Status**: Running (Job 41552397)
**Started**: 2025-11-07
**Expected duration**: 2-4 hours

## Hypothesis

**Random sampling introduces bias** → Poor clustering stability (ARI 0.13-0.32)

**Balanced sampling** (based on sequence similarity) → Better stability (target: ARI > 0.4)

## Test Design

### 1. Sequence-Based Sampling
- MMseqs2 clustered 169,238 proteins (50 genomes) by sequence similarity
- Result: 120,199 sequence clusters
- **Balanced sample**: 5,333 genes across 5,308 sequence families
- Sampling strategy: sqrt-proportional (balances rare and common families)

### 2. Stability Evaluation
- **Input**: 5,333 genes with ESM-C embeddings (50D PCA)
- **Method**: Create 10 random subsamples (4,000 genes each)
- **Clustering**: Leiden algorithm (res=750, nn=15)
- **Metrics**:
  - **ARI** (Adjusted Rand Index): Agreement between subsample clusterings
  - **Gene stability**: % of genes consistently assigned to same cluster
  - **Cluster count stability**: Variance in number of clusters

### 3. Comparison Baseline

**Random sampling results** (from previous 47 configurations):
- Best ARI: 0.32 (res=300, nn=15, COG-only)
- Typical ARI: 0.13-0.25
- **All configs below "moderate" threshold (ARI < 0.4)**

## Expected Outcomes

### ✅ Success (ARI > 0.4)
Balanced sampling improves stability → Scale to full dataset:
1. Run MMseqs2 on all 29M proteins (or use full Prodigal proteins)
2. Create balanced 1M sample from sequence families
3. Run full clustering pipeline
4. Assign remaining 28M genes to clusters

### ⚠️ Marginal (ARI 0.32-0.4)
Modest improvement → Consider hybrid approach:
- Use sequence clusters as coarse grouping
- Within each, use ESM for fine-grained functional clustering
- Adjust MMseqs2 stringency (test different identity thresholds)

### ❌ No Improvement (ARI < 0.32)
Problem is ESM-C embeddings, not sampling → Try gLM2:
- User's coworker had "encouraging results" with gLM2
- gLM2 includes genomic context (5' UTR, promoters, intergenic)
- Generate gLM2 embeddings for test set
- Compare stability

## Why This Might Work

### Problem with Random Sampling
1. **Oversampling abundant families**: Common proteins (housekeeping genes) dominate
2. **Undersampling rare families**: Functionally important but rare proteins missed
3. **Unstable boundaries**: Cluster boundaries vary wildly across subsamples
4. **High variance**: Different samples → very different clusterings

### Balanced Sampling Solution
1. **Representative coverage**: Each protein family contributes proportionally
2. **Stable rare proteins**: Even small families get sampled
3. **Consistent boundaries**: Cluster boundaries more reproducible
4. **Lower variance**: Similar representation → similar clusterings

## Technical Details

### Balanced Sample Statistics
- Total genes: 5,333
- Sequence clusters represented: 5,308
- Average genes per cluster: 1.0
- Max genes per cluster: 2
- Coverage: All clusters from 50-genome test set

### Clustering Parameters
- **Resolution**: 750 (targets ~13,500 clusters on full 1M sample)
- **n_neighbors**: 15 (local neighborhood structure)
- **Subsamples**: 10 (sufficient for stability estimation)
- **Subsample size**: 4,000 (75% of balanced sample)

### Computational Resources
- Memory: 64GB
- CPUs: 16
- Time: 2-4 hours (10 subsamples × k-NN + Leiden)

## Files Created

```
scripts/analysis/
├── evaluate_balanced_clustering.py    # Stability evaluation script
└── submit_balanced_clustering_test.sh # SLURM submission

data/
├── balanced_sample_gene_ids.txt       # 5,333 gene IDs
├── mmseqs_prodigal_test/              # Sequence clustering results
│   ├── clusters.tsv
│   ├── cluster_summary.csv
│   └── cluster_summary_assignments.csv
└── prodigal_test_50genomes.faa        # 169K proteins (50 genomes)

results/clustering/
└── balanced_stability_test.npz        # Output (pending)
```

## How to Monitor

```bash
# Check job status
squeue -u dmullane | grep balanced

# Watch progress
tail -f logs/balanced_clustering_41552397.out

# Check when done
ls -lh results/clustering/balanced_stability_test.npz

# View results
python -c "
import numpy as np
data = np.load('results/clustering/balanced_stability_test.npz')
print('ARI:', data['ari_mean'], '±', data['ari_std'])
print('Gene stability:', data['gene_stability_mean'])
print('N clusters:', data['n_clusters_mean'], '±', data['n_clusters_std'])
"
```

## Interpretation Guide

### ARI (Adjusted Rand Index)
- **1.0**: Perfect agreement (identical clusterings)
- **> 0.8**: Excellent stability
- **0.6-0.8**: Good stability
- **0.4-0.6**: Moderate stability ⚠️
- **< 0.4**: Poor stability ❌

### Gene Stability
- **> 0.9**: Most genes consistently assigned
- **0.7-0.9**: Good consistency
- **< 0.7**: High uncertainty

### Comparison to Random Sampling
| Metric | Random (Best) | Random (Typical) | Balanced (Target) |
|--------|---------------|------------------|-------------------|
| ARI | 0.32 | 0.13-0.25 | > 0.4 |
| Gene Stability | 0.25-0.35 | 0.20-0.30 | > 0.7 |

## Next Steps

Will be determined by results (expected in 2-4 hours):

1. **If successful** → Scale to full dataset with balanced sampling
2. **If marginal** → Investigate hybrid approaches
3. **If no improvement** → Move to gLM2 embeddings

---

**Job ID**: 41552397
**Log**: `logs/balanced_clustering_41552397.out`
**Output**: `results/clustering/balanced_stability_test.npz`
