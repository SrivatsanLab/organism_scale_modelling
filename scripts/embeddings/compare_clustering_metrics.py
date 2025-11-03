#!/usr/bin/env python
"""
Compare clustering results across all evaluation metrics.

Integrates:
1. COG homogeneity (existing metric)
2. ARI/AMI vs COG (quality metrics)
3. Silhouette/Davies-Bouldin (embedding quality)
4. Cluster size distribution
5. Stability (co-clustering rate across subsamples)

Produces rankings and visualizations to identify best clustering parameters.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


def load_homogeneity_metrics():
    """Load existing homogeneity metrics."""
    csv_file = Path('results/clustering/leiden_sweep_summary.csv')
    if not csv_file.exists():
        print(f"Warning: {csv_file} not found")
        return None

    df = pd.read_csv(csv_file)
    print(f"Loaded homogeneity metrics for {len(df)} clusterings")
    return df


def load_quality_metrics():
    """Load comprehensive quality metrics."""
    csv_file = Path('results/clustering/quality_metrics_comprehensive.csv')
    if not csv_file.exists():
        print(f"Warning: {csv_file} not found")
        return None

    df = pd.read_csv(csv_file)
    print(f"Loaded quality metrics for {len(df)} clusterings")
    return df


def load_stability_metrics():
    """Load stability metrics from multiple configurations."""
    stability_dir = Path('results/clustering/stability')
    if not stability_dir.exists():
        print(f"Warning: {stability_dir} not found")
        return None

    stability_files = list(stability_dir.glob('stability_*.npz'))
    if len(stability_files) == 0:
        print("Warning: No stability files found")
        return None

    results = []
    for npz_file in stability_files:
        data = np.load(npz_file, allow_pickle=True)

        result = {
            'resolution': float(data['resolution']),
            'n_neighbors': int(data['n_neighbors']),
            'cog_only': bool(data['cog_only']),
            'stability_mean': float(data['mean_coclustering_rate']),
            'stability_median': float(data['median_coclustering_rate']),
            'stability_std': float(data['std_coclustering_rate']),
            'stability_n_pairs': int(data['n_pairs'])
        }
        results.append(result)

    df = pd.DataFrame(results)
    print(f"Loaded stability metrics for {len(df)} configurations")
    return df


def merge_all_metrics():
    """Merge all metric sources into single dataframe."""
    print("\nMerging all metrics...")

    # Load all sources
    df_homog = load_homogeneity_metrics()
    df_quality = load_quality_metrics()
    df_stability = load_stability_metrics()

    if df_homog is None:
        print("Error: Cannot proceed without homogeneity metrics")
        return None

    # Start with homogeneity
    df = df_homog.copy()

    # Merge quality metrics
    if df_quality is not None:
        # Match on key columns
        df = df.merge(
            df_quality,
            on=['file', 'resolution', 'n_neighbors', 'cog_only'],
            how='left',
            suffixes=('', '_quality')
        )

    # Merge stability metrics
    if df_stability is not None:
        # Match on resolution, n_neighbors, cog_only
        df = df.merge(
            df_stability,
            on=['resolution', 'n_neighbors', 'cog_only'],
            how='left',
            suffixes=('', '_stability')
        )

    print(f"\nMerged dataframe: {len(df)} rows, {len(df.columns)} columns")
    return df


def compute_composite_scores(df):
    """
    Compute composite quality scores by combining multiple metrics.

    Higher is better for all scores.
    """
    print("\nComputing composite scores...")

    # Normalize metrics to [0, 1] scale
    def normalize(series, higher_is_better=True):
        """Normalize to [0, 1] with 0=worst, 1=best."""
        if series.isna().all():
            return series
        min_val = series.min()
        max_val = series.max()
        if max_val == min_val:
            return pd.Series(0.5, index=series.index)
        normalized = (series - min_val) / (max_val - min_val)
        if not higher_is_better:
            normalized = 1 - normalized
        return normalized

    # Select metrics to include
    scores = pd.DataFrame(index=df.index)

    # 1. COG homogeneity (higher is better)
    if 'weighted_homogeneity' in df.columns:
        scores['homog_score'] = normalize(df['weighted_homogeneity'])

    # 2. ARI vs COG (higher is better)
    if 'ari_vs_cog' in df.columns:
        scores['ari_score'] = normalize(df['ari_vs_cog'])

    # 3. AMI vs COG (higher is better)
    if 'ami_vs_cog' in df.columns:
        scores['ami_score'] = normalize(df['ami_vs_cog'])

    # 4. Silhouette (higher is better)
    if 'silhouette_score_sampled' in df.columns:
        scores['silhouette_score'] = normalize(df['silhouette_score_sampled'])

    # 5. Davies-Bouldin (LOWER is better)
    if 'davies_bouldin_score' in df.columns:
        scores['db_score'] = normalize(df['davies_bouldin_score'], higher_is_better=False)

    # 6. Cluster size distribution (prefer moderate Gini, not too high or low)
    if 'size_gini' in df.columns:
        # Target Gini around 0.3-0.5 (some inequality but not extreme)
        target_gini = 0.4
        gini_penalty = (df['size_gini'] - target_gini).abs()
        scores['gini_score'] = normalize(gini_penalty, higher_is_better=False)

    # 7. Stability (higher is better)
    if 'stability_mean' in df.columns:
        scores['stability_score'] = normalize(df['stability_mean'])

    # Compute composite scores with different weightings

    # Quality-focused (emphasize ARI/AMI)
    quality_weights = {
        'homog_score': 0.2,
        'ari_score': 0.3,
        'ami_score': 0.3,
        'silhouette_score': 0.1,
        'db_score': 0.1
    }

    # Stability-focused (emphasize consistency)
    stability_weights = {
        'homog_score': 0.2,
        'ari_score': 0.2,
        'stability_score': 0.4,
        'silhouette_score': 0.1,
        'db_score': 0.1
    }

    # Balanced (equal weights)
    available_scores = [col for col in scores.columns if not scores[col].isna().all()]
    balanced_weight = 1.0 / len(available_scores) if available_scores else 0

    df['composite_quality'] = sum(
        scores[metric] * quality_weights.get(metric, 0)
        for metric in available_scores
    )

    df['composite_stability'] = sum(
        scores[metric] * stability_weights.get(metric, 0)
        for metric in available_scores
    )

    df['composite_balanced'] = sum(
        scores[metric] * balanced_weight
        for metric in available_scores
    )

    print(f"  Computed {len(available_scores)} normalized scores")
    print(f"  Metrics used: {', '.join(available_scores)}")

    return df


def print_rankings(df, metric='composite_balanced', n=15):
    """Print top N clusterings by specified metric."""
    print(f"\n{'='*100}")
    print(f"Top {n} Clusterings by {metric}")
    print(f"{'='*100}\n")

    # Filter to valid values
    df_valid = df[df[metric].notna()].copy()
    df_sorted = df_valid.sort_values(metric, ascending=False)

    # Print rankings
    for i, (idx, row) in enumerate(df_sorted.head(n).iterrows()):
        print(f"{i+1}. {row['file']}")
        print(f"   Resolution: {row['resolution']:.0f}, N neighbors: {row['n_neighbors']}, "
              f"Version: {'COG-only' if row['cog_only'] else 'All genes'}")
        print(f"   {metric}: {row[metric]:.4f}")

        # Show key individual metrics
        if 'weighted_homogeneity' in row:
            print(f"   Weighted homogeneity: {row['weighted_homogeneity']:.4f}")
        if 'ari_vs_cog' in row and not pd.isna(row['ari_vs_cog']):
            print(f"   ARI vs COG: {row['ari_vs_cog']:.4f}")
        if 'ami_vs_cog' in row and not pd.isna(row['ami_vs_cog']):
            print(f"   AMI vs COG: {row['ami_vs_cog']:.4f}")
        if 'stability_mean' in row and not pd.isna(row['stability_mean']):
            print(f"   Stability: {row['stability_mean']:.4f}")
        print()


def plot_metric_comparison(df, output_dir='results/plots/metric_comparison'):
    """Create visualization comparing metrics across parameters."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nCreating comparison plots in {output_dir}...")

    # Filter to evaluated clusterings
    df_eval = df[df['composite_balanced'].notna()].copy()

    if len(df_eval) == 0:
        print("  No evaluated clusterings to plot")
        return

    # 1. Resolution vs multiple metrics
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Clustering Metrics vs Resolution', fontsize=16, y=0.995)

    metrics_to_plot = [
        ('weighted_homogeneity', 'COG Homogeneity'),
        ('ari_vs_cog', 'ARI vs COG'),
        ('ami_vs_cog', 'AMI vs COG'),
        ('silhouette_score_sampled', 'Silhouette Score'),
        ('davies_bouldin_score', 'Davies-Bouldin (lower=better)'),
        ('stability_mean', 'Stability (Co-clustering Rate)')
    ]

    for ax, (metric, title) in zip(axes.flat, metrics_to_plot):
        if metric not in df_eval.columns or df_eval[metric].isna().all():
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')
            ax.set_title(title)
            continue

        # Plot separately for COG-only vs all genes
        for cog_only, label, marker in [(True, 'COG-only', 'o'), (False, 'All genes', 's')]:
            subset = df_eval[df_eval['cog_only'] == cog_only]
            if len(subset) > 0:
                # Group by resolution, show mean + std
                grouped = subset.groupby('resolution')[metric].agg(['mean', 'std'])
                ax.errorbar(grouped.index, grouped['mean'], yerr=grouped['std'],
                           marker=marker, label=label, capsize=3, alpha=0.7)

        ax.set_xlabel('Resolution')
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'metrics_vs_resolution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: metrics_vs_resolution.png")

    # 2. Composite scores comparison
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    for cog_only, label, marker in [(True, 'COG-only', 'o'), (False, 'All genes', 's')]:
        subset = df_eval[df_eval['cog_only'] == cog_only]
        if len(subset) > 0:
            grouped = subset.groupby('resolution')['composite_balanced'].agg(['mean', 'std'])
            ax.errorbar(grouped.index, grouped['mean'], yerr=grouped['std'],
                       marker=marker, label=label, capsize=3, alpha=0.7, linewidth=2)

    ax.set_xlabel('Resolution', fontsize=12)
    ax.set_ylabel('Composite Balanced Score', fontsize=12)
    ax.set_title('Overall Clustering Quality vs Resolution', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'composite_score_vs_resolution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: composite_score_vs_resolution.png")


def main():
    print("="*100)
    print("Comprehensive Clustering Metrics Comparison")
    print("="*100)

    # Merge all metrics
    df = merge_all_metrics()
    if df is None:
        print("Error: Could not load metrics")
        return

    # Compute composite scores
    df = compute_composite_scores(df)

    # Save combined metrics
    output_file = Path('results/clustering/all_metrics_combined.csv')
    df.to_csv(output_file, index=False)
    print(f"\nSaved combined metrics to {output_file}")
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {len(df.columns)}")

    # Print rankings
    print_rankings(df, metric='composite_balanced', n=15)
    print_rankings(df, metric='composite_quality', n=10)
    print_rankings(df, metric='composite_stability', n=10)

    # Create plots
    plot_metric_comparison(df)

    print("\n" + "="*100)
    print("Done!")
    print("="*100)


if __name__ == '__main__':
    main()
