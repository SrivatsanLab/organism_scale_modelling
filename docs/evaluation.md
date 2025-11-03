# Clustering Evaluation Framework

Comprehensive evaluation system for assessing Leiden clustering quality beyond simple homogeneity metrics.

## Problem

Simply maximizing COG homogeneity is trivial - higher resolution creates smaller, more homogeneous clusters. We need metrics that assess **true quality** and **robustness**.

## Solution: Multi-Metric Evaluation

### 1. Quality Metrics (`evaluate_clustering_quality.py`)

Evaluates each clustering using:

- **Adjusted Rand Index (ARI)**: Similarity to COG annotations (0-1, higher better)
- **Adjusted Mutual Information (AMI)**: Information-theoretic similarity to COG (0-1, higher better)
- **Silhouette Score**: Cluster separation quality (-1 to 1, higher better)
- **Davies-Bouldin Index**: Cluster compactness (lower better)
- **Size Distribution**: Gini coefficient, percentiles, mean/median

**Key insight**: Tests if clusters align with known biology (COG), not just internal consistency.

**Usage**:
```bash
# Single clustering
python evaluate_clustering_quality.py \
    --clustering results/clustering/clusters_leiden_res1500_nn15_cogonly.npz \
    --pca results/embeddings/refseq_esm_pca.npz \
    --output results/clustering/quality_single.csv

# All clusterings (batch mode)
sbatch submit_quality_evaluation.sh
# Output: results/clustering/quality_metrics_comprehensive.csv
```

### 2. Stability Analysis (`evaluate_clustering_stability.py`)

Tests clustering robustness using **co-clustering across independent subsamples**.

**Approach**:
1. Generate N independent random subsamples (e.g., 10 × 100K genes)
2. Cluster each subsample with same parameters
3. For gene pairs appearing in multiple subsamples:
   - Track co-occurrence frequency
   - Track co-clustering frequency (same cluster)
4. Compute **co-clustering rate** = (times together) / (times co-occurred)

**Key insight**: Robust, meaningful clusters should persist across different random samples.

**High stability (rate ~1.0)**: Genes consistently cluster together → strong signal
**Low stability (rate ~0.0)**: Random clustering → weak signal
**Random baseline (rate ~0.5)**: No structure

**Usage**:
```bash
# Single configuration
python evaluate_clustering_stability.py \
    --pca results/embeddings/refseq_esm_pca.npz \
    --n-subsamples 10 \
    --subsample-size 100000 \
    --resolution 1500 \
    --n-neighbors 15 \
    --cog-only \
    --output results/clustering/stability/stability_res1500_nn15_cogonly.npz

# Multiple configurations (job array)
sbatch submit_stability_evaluation.sh
# Tests: res1500, res1000, res750 with n=15
```

### 3. Stability Visualization (`visualize_clustering_stability.py`)

Creates comprehensive visualizations:

1. **Distribution plots**: Histogram of co-clustering rates
2. **Comparison plots**: Stability metrics across configurations
3. **Hexbin plots**: Stability vs co-occurrence frequency
4. **COG analysis**: Stability for same vs different COG pairs
5. **Network plots**: Graph of highly stable gene pairs
6. **Summary report**: Text file with detailed statistics

**Usage**:
```bash
# After stability jobs complete
python visualize_clustering_stability.py \
    --stability-dir results/clustering/stability \
    --output-dir results/plots/stability_analysis \
    --stable-threshold 0.9 \
    --unstable-threshold 0.3

# Creates:
# - coclustering_rate_distributions.png
# - stability_comparison.png
# - stability_vs_cooccurrence.png
# - stability_by_cog_*.png (for COG-only configs)
# - stable_network_*.png
# - stable_pairs_*.txt (lists of gene pairs)
# - stability_summary_report.txt
```

### 4. Comprehensive Comparison (`compare_clustering_metrics.py`)

Integrates all metrics to identify best clustering:

1. Loads and merges:
   - Homogeneity metrics (existing)
   - Quality metrics (ARI, AMI, silhouette, etc.)
   - Stability metrics (co-clustering rates)

2. Normalizes all to [0, 1] scale

3. Computes **composite scores**:
   - **Quality-focused**: ARI/AMI weighted 30% each
   - **Stability-focused**: Co-clustering rate weighted 40%
   - **Balanced**: Equal weights across all metrics

4. Creates rankings and visualizations

**Usage**:
```bash
# After all evaluations complete
python compare_clustering_metrics.py

# Creates:
# - results/clustering/all_metrics_combined.csv
# - results/plots/metric_comparison/metrics_vs_resolution.png
# - results/plots/metric_comparison/composite_score_vs_resolution.png
```

## Workflow

### Step 1: Quality Evaluation (~4 hours)
```bash
sbatch submit_quality_evaluation.sh
# Wait for Job 41235321 to complete
```

### Step 2: Stability Evaluation (~12 hours)
```bash
sbatch submit_stability_evaluation.sh
# Wait for Job 41235333 (4 tasks) to complete
```

### Step 3: Visualize Stability
```bash
python visualize_clustering_stability.py
```

### Step 4: Comprehensive Comparison
```bash
python compare_clustering_metrics.py
```

### Step 5: Interpret Results

Look for configurations that:
- **High ARI/AMI**: Align with biological function (COG)
- **High stability**: Robust across subsamples (not random)
- **Good silhouette**: Well-separated in embedding space
- **Moderate size Gini**: Not too uniform, not too skewed
- **High composite score**: Balance all metrics

## Expected Insights

After completion, you can answer:

1. **Does resolution 1500 have better quality or just smaller clusters?**
   - Compare ARI at different resolutions
   - Check if stability decreases (overfitting signal)

2. **Which parameters give stable, reproducible clusters?**
   - Compare co-clustering rates
   - Identify sweet spot balancing quality and stability

3. **Are clusters biologically meaningful?**
   - Check ARI/AMI vs COG
   - Examine if stable pairs share COG categories

4. **What's the optimal configuration?**
   - Check composite_balanced ranking
   - Trade-off curves in metric comparison plots

## Output Files

### Quality Evaluation
- `results/clustering/quality_metrics_comprehensive.csv`: All quality metrics for all clusterings

### Stability Evaluation
- `results/clustering/stability/stability_res*_nn*_*.npz`: Raw stability data
  - Contains: gene_pairs, cooccur_counts, cocluster_counts, coclustering_rates

### Stability Visualization
- `results/plots/stability_analysis/`:
  - Distributions, comparisons, networks
  - Lists of stable/unstable gene pairs
  - Summary report

### Comprehensive Comparison
- `results/clustering/all_metrics_combined.csv`: All metrics merged
- `results/plots/metric_comparison/`:
  - Metrics vs resolution
  - Composite scores

## Interpretation Guide

### ARI/AMI Scores
- **0.0-0.2**: Poor alignment with COG (random)
- **0.2-0.4**: Weak alignment
- **0.4-0.6**: Moderate alignment
- **0.6-0.8**: Good alignment
- **0.8-1.0**: Excellent alignment (rare, COG not perfect ground truth)

### Co-clustering Rates
- **<0.3**: Very unstable (basically random)
- **0.3-0.5**: Unstable (near random baseline)
- **0.5-0.7**: Moderate stability
- **0.7-0.9**: Stable
- **>0.9**: Very stable (strong consistent signal)

### Silhouette Scores
- **<0**: Poor separation (overlapping clusters)
- **0-0.25**: Weak separation
- **0.25-0.5**: Moderate separation
- **0.5-0.75**: Good separation
- **>0.75**: Excellent separation

### Davies-Bouldin Index
- **Lower is better** (measures compactness)
- Typical range: 0.5-3.0
- Compare relative values across configs

## Notes

- Quality evaluation uses 50K sample for silhouette (too slow on full data)
- Stability uses 100K subsamples (balance between coverage and speed)
- COG comparison only works for COG-annotated genes (~798K of 1M)
- Network visualization subsamples if >500 genes or >5000 edges
- All metrics normalized to [0, 1] for composite scores

## Advantages Over Homogeneity Alone

1. **Prevents trivial optimization**: Can't just increase resolution
2. **Tests biological relevance**: ARI/AMI vs known annotations
3. **Assesses robustness**: Stability across samples
4. **Multi-dimensional view**: No single metric is perfect
5. **Enables trade-off analysis**: Balance different objectives

## Questions?

See individual script help:
```bash
python evaluate_clustering_quality.py --help
python evaluate_clustering_stability.py --help
python visualize_clustering_stability.py --help
python compare_clustering_metrics.py --help
```
