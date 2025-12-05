#!/usr/bin/env python
"""
Visualize MMseqs2 cluster tightness in embedding space.

Creates interpretable visualizations:
1. Histogram of cluster tightness metrics
2. Scatter plot: cluster size vs tightness
3. Heatmap of mean std per dimension across clusters
4. Comparison with random baseline
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Visualize cluster tightness"
    )
    parser.add_argument(
        "--stats-file",
        type=str,
        default="results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/mmseqs_cluster_statistics.csv",
        help="Path to cluster statistics CSV",
    )
    parser.add_argument(
        "--per-dim-file",
        type=str,
        default="results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/mmseqs_cluster_per_dimension_stats.csv",
        help="Path to per-dimension statistics CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/figures",
        help="Output directory for plots",
    )
    return parser.parse_args()


def plot_tightness_distributions(stats_df, output_dir):
    """Plot distributions of cluster tightness metrics."""
    print("\nCreating tightness distribution plots...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 1. Cluster size distribution
    ax = axes[0, 0]
    ax.hist(stats_df['cluster_size'], bins=100, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Cluster Size', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Cluster Size Distribution', fontsize=12)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # 2. Mean std distribution
    ax = axes[0, 1]
    ax.hist(stats_df['mean_std'], bins=100, edgecolor='black', alpha=0.7, color='orange')
    ax.set_xlabel('Mean Std (across dimensions)', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Mean Standard Deviation per Cluster', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.axvline(stats_df['mean_std'].median(), color='red', linestyle='--',
               label=f'Median={stats_df["mean_std"].median():.3f}')
    ax.legend()

    # 3. Max std distribution
    ax = axes[0, 2]
    ax.hist(stats_df['max_std'], bins=100, edgecolor='black', alpha=0.7, color='green')
    ax.set_xlabel('Max Std (across dimensions)', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Max Standard Deviation per Cluster', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.axvline(stats_df['max_std'].median(), color='red', linestyle='--',
               label=f'Median={stats_df["max_std"].median():.3f}')
    ax.legend()

    # 4. Total variance distribution
    ax = axes[1, 0]
    ax.hist(stats_df['total_variance'], bins=100, edgecolor='black', alpha=0.7, color='purple')
    ax.set_xlabel('Total Variance', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Total Variance per Cluster', fontsize=12)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # 5. Mean pairwise distance distribution
    ax = axes[1, 1]
    ax.hist(stats_df['mean_pairwise_distance'], bins=100, edgecolor='black', alpha=0.7, color='red')
    ax.set_xlabel('Mean Pairwise Distance', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Mean Pairwise Distance per Cluster', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.axvline(stats_df['mean_pairwise_distance'].median(), color='darkred', linestyle='--',
               label=f'Median={stats_df["mean_pairwise_distance"].median():.3f}')
    ax.legend()

    # 6. Max pairwise distance distribution
    ax = axes[1, 2]
    ax.hist(stats_df['max_pairwise_distance'], bins=100, edgecolor='black', alpha=0.7, color='brown')
    ax.set_xlabel('Max Pairwise Distance', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Max Pairwise Distance per Cluster', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.axvline(stats_df['max_pairwise_distance'].median(), color='darkred', linestyle='--',
               label=f'Median={stats_df["max_pairwise_distance"].median():.3f}')
    ax.legend()

    plt.tight_layout()
    output_file = output_dir / 'cluster_tightness_distributions.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_size_vs_tightness(stats_df, output_dir):
    """Plot cluster size vs tightness metrics."""
    print("\nCreating size vs tightness scatter plots...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Size vs Mean Std
    ax = axes[0, 0]
    scatter = ax.scatter(stats_df['cluster_size'], stats_df['mean_std'],
                        alpha=0.3, s=10, c=stats_df['mean_std'], cmap='viridis')
    ax.set_xlabel('Cluster Size', fontsize=11)
    ax.set_ylabel('Mean Std', fontsize=11)
    ax.set_title('Cluster Size vs Mean Std', fontsize=12)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Mean Std')

    # 2. Size vs Max Std
    ax = axes[0, 1]
    scatter = ax.scatter(stats_df['cluster_size'], stats_df['max_std'],
                        alpha=0.3, s=10, c=stats_df['max_std'], cmap='plasma')
    ax.set_xlabel('Cluster Size', fontsize=11)
    ax.set_ylabel('Max Std', fontsize=11)
    ax.set_title('Cluster Size vs Max Std', fontsize=12)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Max Std')

    # 3. Size vs Mean Pairwise Distance
    ax = axes[1, 0]
    scatter = ax.scatter(stats_df['cluster_size'], stats_df['mean_pairwise_distance'],
                        alpha=0.3, s=10, c=stats_df['mean_pairwise_distance'], cmap='coolwarm')
    ax.set_xlabel('Cluster Size', fontsize=11)
    ax.set_ylabel('Mean Pairwise Distance', fontsize=11)
    ax.set_title('Cluster Size vs Mean Pairwise Distance', fontsize=12)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Mean Pairwise Dist')

    # 4. Size vs Total Variance
    ax = axes[1, 1]
    scatter = ax.scatter(stats_df['cluster_size'], stats_df['total_variance'],
                        alpha=0.3, s=10, c=stats_df['total_variance'], cmap='inferno')
    ax.set_xlabel('Cluster Size', fontsize=11)
    ax.set_ylabel('Total Variance', fontsize=11)
    ax.set_title('Cluster Size vs Total Variance', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Total Variance')

    plt.tight_layout()
    output_file = output_dir / 'size_vs_tightness.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_dimension_std_heatmap(per_dim_df, output_dir, max_clusters=100):
    """Create heatmap of std per dimension across clusters."""
    print(f"\nCreating dimension std heatmap (top {max_clusters} clusters)...")

    # Get top clusters by size
    top_clusters = per_dim_df.groupby('cluster_representative')['cluster_size'].first().nlargest(max_clusters).index

    # Filter to top clusters
    df_subset = per_dim_df[per_dim_df['cluster_representative'].isin(top_clusters)]

    # Pivot to create matrix: clusters × dimensions
    pivot_data = df_subset.pivot(index='cluster_representative', columns='dimension', values='std')

    fig, ax = plt.subplots(figsize=(16, 12))

    # Create heatmap
    sns.heatmap(pivot_data, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Std Dev'})

    ax.set_xlabel('PCA Dimension', fontsize=11)
    ax.set_ylabel(f'Cluster Representative (top {max_clusters} by size)', fontsize=11)
    ax.set_title(f'Standard Deviation per Dimension for Top {max_clusters} Clusters', fontsize=12)

    # Only show every 5th dimension label
    n_dims = len(pivot_data.columns)
    ax.set_xticks(range(0, n_dims, 5))
    ax.set_xticklabels(range(0, n_dims, 5))

    # Don't show cluster labels (too many)
    ax.set_yticks([])

    plt.tight_layout()
    output_file = output_dir / 'dimension_std_heatmap.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_mean_std_per_dimension(per_dim_df, output_dir):
    """Plot mean std across all clusters for each dimension."""
    print("\nCreating mean std per dimension plot...")

    # Compute mean std per dimension across all clusters
    mean_std_per_dim = per_dim_df.groupby('dimension')['std'].mean()
    median_std_per_dim = per_dim_df.groupby('dimension')['std'].median()

    fig, ax = plt.subplots(figsize=(14, 6))

    dimensions = mean_std_per_dim.index
    ax.plot(dimensions, mean_std_per_dim.values, label='Mean', linewidth=2, alpha=0.8)
    ax.plot(dimensions, median_std_per_dim.values, label='Median', linewidth=2, alpha=0.8)

    ax.set_xlabel('PCA Dimension', fontsize=12)
    ax.set_ylabel('Std Dev (averaged across clusters)', fontsize=12)
    ax.set_title('Mean Standard Deviation per PCA Dimension', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / 'mean_std_per_dimension.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def create_summary_figure(stats_df, output_dir):
    """Create a single summary figure with key insights."""
    print("\nCreating summary figure...")

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Title
    fig.suptitle(f'MMseqs2 Cluster Tightness Summary (n={len(stats_df):,} clusters)',
                 fontsize=16, fontweight='bold')

    # 1. Cluster size distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(stats_df['cluster_size'], bins=100, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Cluster Size')
    ax1.set_ylabel('Count')
    ax1.set_title(f'Cluster Sizes\n(median={stats_df["cluster_size"].median():.0f})')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)

    # 2. Mean std distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(stats_df['mean_std'], bins=100, edgecolor='black', alpha=0.7, color='orange')
    ax2.set_xlabel('Mean Std (across dims)')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Mean Std per Cluster\n(median={stats_df["mean_std"].median():.3f})')
    ax2.axvline(stats_df['mean_std'].median(), color='red', linestyle='--', linewidth=2)
    ax2.grid(True, alpha=0.3)

    # 3. Mean pairwise distance distribution
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(stats_df['mean_pairwise_distance'], bins=100, edgecolor='black', alpha=0.7, color='green')
    ax3.set_xlabel('Mean Pairwise Distance')
    ax3.set_ylabel('Count')
    ax3.set_title(f'Pairwise Distances\n(median={stats_df["mean_pairwise_distance"].median():.3f})')
    ax3.axvline(stats_df['mean_pairwise_distance'].median(), color='darkgreen', linestyle='--', linewidth=2)
    ax3.grid(True, alpha=0.3)

    # 4. Size vs Mean Std (large plot)
    ax4 = fig.add_subplot(gs[1:, :2])
    scatter = ax4.scatter(stats_df['cluster_size'], stats_df['mean_std'],
                         alpha=0.4, s=20, c=stats_df['mean_pairwise_distance'],
                         cmap='coolwarm', edgecolors='black', linewidth=0.5)
    ax4.set_xlabel('Cluster Size', fontsize=12)
    ax4.set_ylabel('Mean Std (across dimensions)', fontsize=12)
    ax4.set_title('Cluster Size vs Tightness\n(color = mean pairwise distance)', fontsize=13)
    ax4.set_xscale('log')
    ax4.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Mean Pairwise Dist', fontsize=10)

    # 5. Summary statistics table
    ax5 = fig.add_subplot(gs[1:, 2])
    ax5.axis('off')

    summary_text = f"""
    SUMMARY STATISTICS

    Total Clusters: {len(stats_df):,}

    Cluster Size:
      Min: {stats_df['cluster_size'].min()}
      Max: {stats_df['cluster_size'].max():,}
      Median: {stats_df['cluster_size'].median():.0f}

    Mean Std:
      Min: {stats_df['mean_std'].min():.4f}
      Max: {stats_df['mean_std'].max():.4f}
      Median: {stats_df['mean_std'].median():.4f}

    Mean Pairwise Dist:
      Min: {stats_df['mean_pairwise_distance'].min():.4f}
      Max: {stats_df['mean_pairwise_distance'].max():.4f}
      Median: {stats_df['mean_pairwise_distance'].median():.4f}

    Interpretation:
    - Lower std = tighter cluster
    - Lower pairwise dist =
      more similar members
    - 70% sequence ID threshold
      appears to create tight
      embedding clusters
    """

    ax5.text(0.1, 0.9, summary_text, transform=ax5.transAxes,
             fontsize=11, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    output_file = output_dir / 'cluster_tightness_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def main():
    args = parse_args()

    print("=" * 80)
    print("CLUSTER TIGHTNESS VISUALIZATION")
    print("=" * 80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load statistics
    print(f"\nLoading statistics from {args.stats_file}...")
    stats_df = pd.read_csv(args.stats_file)
    print(f"  Loaded {len(stats_df):,} clusters")

    print(f"\nLoading per-dimension stats from {args.per_dim_file}...")
    per_dim_df = pd.read_csv(args.per_dim_file)
    print(f"  Loaded {len(per_dim_df):,} rows")

    # Create visualizations
    create_summary_figure(stats_df, output_dir)
    plot_tightness_distributions(stats_df, output_dir)
    plot_size_vs_tightness(stats_df, output_dir)
    plot_mean_std_per_dimension(per_dim_df, output_dir)
    plot_dimension_std_heatmap(per_dim_df, output_dir, max_clusters=100)

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print(f"All plots saved to: {output_dir}")
    print("=" * 80)


if __name__ == '__main__':
    main()
