# MSA-Based Clustering Approach: Test Implementation

## Motivation

The initial ESM-C embedding clustering showed **poor stability** (ARI 0.13-0.32) across all tested configurations. This suggests that either:
1. ESM-C embeddings don't capture the functional structure we need, OR
2. Random sampling introduces bias that destabilizes clustering

## New Approach: MSA + Balanced Sampling

Instead of random sampling → ESM embeddings → clustering, we're testing:

**Sequence similarity clustering → Balanced sampling → ESM embeddings → Functional clustering**

### Rationale

1. **Address sampling bias**: By clustering proteins by sequence similarity first, we ensure our sample represents the full diversity of protein families
2. **Leverage complementary information**:
   - Sequence similarity captures homology and evolutionary relationships
   - ESM embeddings capture structural and functional properties
3. **Proven approach**: This mirrors traditional protein family classification (Pfam, COG, etc.)

## Implementation

### Phase 1: Proof of Concept (In Progress)

#### Step 1: Download Test Data ✓
- Downloaded protein sequences for 50 genomes from NCBI (160,840 proteins)
- Coverage: 5,526 genes overlap with our PCA cache (0.55% of 1M subsample)
- File: `data/test_proteins_50genomes.faa` (61 MB)

#### Step 2: Sequence Clustering with MMseqs2 (Running)
- **Job**: 41552086 (`logs/mmseqs_test_41552086.out`)
- **Tool**: MMseqs2 (ultra-fast sequence clustering, 100-1000x faster than BLAST)
- **Parameters**:
  - Min sequence identity: 0.5 (50%)
  - Min coverage: 0.8 (80%)
  - Cluster mode: 0 (greedy set cover)
  - Threads: 8
- **Output**: `data/mmseqs_test/clusters.tsv`

#### Step 3: Balanced Sampling
- **Script**: `scripts/analysis/balanced_sample_from_clusters.py`
- **Strategies**:
  - `equal`: Equal number from each cluster
  - `proportional`: Proportional to cluster size
  - `sqrt`: Proportional to sqrt(cluster size) - **recommended**
- **Target**: Sample genes that (a) are in sequence clusters, (b) have ESM embeddings in PCA cache

#### Step 4: ESM Clustering on Balanced Sample
- Run Leiden clustering on the balanced sample
- Compare stability metrics to random sampling approach
- Evaluate if balanced sampling improves ARI scores

## Expected Outcomes

### Success Criteria
- **Sequence clusters**: Expect 5,000-50,000 clusters from 160K proteins (depending on stringency)
- **Balanced sample**: Successfully sample ~10,000-100,000 genes with better family representation
- **Improved stability**: ARI > 0.4 (vs 0.13-0.32 for random sampling)

### If This Works
1. Scale to full dataset:
   - Download all 7,664 genomes (requires ~200GB storage for protein FASTAs)
   - Run MMseqs2 on full 29M proteins (may need to split into batches)
   - Create balanced 1M sample from sequence clusters
   - Run full clustering pipeline on balanced sample

2. Use sequence clusters as ground truth:
   - Compare ESM clustering to sequence clustering
   - Evaluate if ESM adds value beyond sequence similarity

### If This Doesn't Work
Consider gLM2 embeddings as planned - they include genomic context (5' UTR, promoters) which may capture regulatory information that ESM-C misses.

## Scripts Created

### Data Acquisition
- `scripts/analysis/download_test_proteins.py`: Download protein FASTAs from NCBI using Entrez API
- `scripts/analysis/submit_mmseqs_test.sh`: SLURM script for MMseqs2 clustering

### Analysis
- `scripts/analysis/cluster_proteins_mmseqs.py`: Python wrapper for MMseqs2
- `scripts/analysis/analyze_mmseqs_clusters.py`: Analyze clustering results (size distribution, statistics)
- `scripts/analysis/balanced_sample_from_clusters.py`: Create balanced samples from sequence clusters

## Files and Locations

```
data/
├── test_proteins.faa              # 5 genomes, 8.4K proteins (initial test)
├── test_proteins_50genomes.faa    # 50 genomes, 160K proteins
├── download_log.txt               # Download log
└── mmseqs_test/                   # MMseqs2 clustering results (pending)
    ├── clusters.tsv               # Cluster assignments
    └── cluster_summary.csv        # Size distribution

results/umap/
└── pca_cache.npz                  # 1M genes × 50D PCA embeddings
                                   # (5,526 genes overlap with test set)
```

## Next Steps (After MMseqs2 Completes)

1. **Check MMseqs2 results**:
   ```bash
   tail logs/mmseqs_test_41552086.out
   head data/mmseqs_test/cluster_summary.csv
   ```

2. **Create balanced sample**:
   ```bash
   python scripts/analysis/balanced_sample_from_clusters.py \
       --clusters data/mmseqs_test/cluster_summary_assignments.csv \
       --pca-cache results/umap/pca_cache.npz \
       --output data/balanced_sample_gene_ids.txt \
       --n-samples 5000 \
       --strategy sqrt
   ```

3. **Run clustering on balanced sample**: Adapt existing Leiden clustering scripts to use the balanced sample instead of random sampling

4. **Evaluate stability**: Run stability evaluation (5-10 subsamples) and compare ARI to random sampling baseline

## Technical Details

### Gene ID Matching
The NCBI FASTA headers have complex formats:
```
>lcl|NZ_CP191850.1_prot_WP_123456.1_1069 [gene=...] [protein=...] [protein_id=...] [location=...]
```

We extract the contig_genenum format (`NZ_CP191850.1_1069`) to match our PCA cache gene IDs.

### MMseqs2 Parameters
- **Min sequence identity (0.5)**: Proteins must share ≥50% sequence identity to cluster together
  - Lower = broader clusters (more like protein superfamilies)
  - Higher = narrower clusters (more like orthologs)

- **Coverage (0.8)**: At least 80% of shorter sequence must align
  - Prevents partial domain matches
  - Ensures full-length homology

- **Cluster mode (0 - greedy)**: Fast, produces non-overlapping clusters
  - Alternative: Connected component (mode 1) - slower but may find better clusters

## Computational Resources

- MMseqs2 clustering: ~30 min, 32GB RAM, 8 CPUs
- Balanced sampling: ~5 min, 16GB RAM
- Total test time: ~1 hour (including queue wait)

---

**Status**: MMseqs2 clustering job running (Job 41552086)
**Date**: 2025-11-07
