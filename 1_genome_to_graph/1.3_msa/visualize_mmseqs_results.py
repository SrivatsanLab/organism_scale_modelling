#!/usr/bin/env python
"""
Visualize MMseqs2 clustering results across different sequence identity thresholds.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 10


def load_cluster_data(cluster_summary_path, seq_id):
    """Load and process cluster data."""
    df = pd.read_csv(cluster_summary_path)
    return df, seq_id


def plot_cluster_size_distributions(data_dict, output_path):
    """Plot cluster size distributions for different sequence identities."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    seq_ids = sorted(data_dict.keys())

    for idx, seq_id in enumerate(seq_ids):
        df = data_dict[seq_id]
        ax = axes[idx]

        # Histogram of cluster sizes (log scale)
        sizes = df['size'].values

        # Create bins on log scale
        bins = [1, 2, 5, 10, 20, 50, 100, 500, 1000, 5000, 10000, 50000]

        ax.hist(sizes, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Cluster Size (proteins)')
        ax.set_ylabel('Number of Clusters')
        ax.set_title(f'Seq ID = {seq_id:.1f}\n{len(df):,} clusters, mean size = {sizes.mean():.1f}')
        ax.grid(True, alpha=0.3)

        # Add text with key statistics
        n_singletons = (sizes == 1).sum()
        pct_singletons = n_singletons / len(sizes) * 100
        ax.text(0.98, 0.97, f'Singletons: {pct_singletons:.1f}%',
                transform=ax.transAxes, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Hide extra subplot
    axes[-1].axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved cluster size distributions to {output_path}")
    plt.close()


def plot_protein_distribution(comparison_df, output_path):
    """Plot how proteins are distributed across cluster size categories."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Prepare data
    categories = ['Singleton', 'Small\n(2-9)', 'Medium\n(10-99)', 'Large\n(100+)']

    seq_ids = comparison_df['seq_id'].values
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(seq_ids)))

    # Plot 1: Stacked bar chart
    x = np.arange(len(seq_ids))
    width = 0.6

    bottom = np.zeros(len(seq_ids))

    cat_cols = ['pct_proteins_singleton', 'pct_proteins_small',
                'pct_proteins_medium', 'pct_proteins_large']
    cat_colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']

    for cat_idx, (col, cat_name) in enumerate(zip(cat_cols, categories)):
        values = comparison_df[col].values
        ax1.bar(x, values, width, label=cat_name, bottom=bottom, color=cat_colors[cat_idx])
        bottom += values

    ax1.set_xlabel('Sequence Identity Threshold')
    ax1.set_ylabel('Percentage of Proteins')
    ax1.set_title('Protein Distribution by Cluster Size Category')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{sid:.1f}' for sid in seq_ids])
    ax1.legend(loc='upper left')
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3, axis='y')

    # Plot 2: Line plot showing key metrics
    ax2.plot(seq_ids, comparison_df['pct_proteins_singleton'],
             'o-', linewidth=2, markersize=8, label='In Singletons', color='#d62728')
    ax2.plot(seq_ids, comparison_df['pct_proteins_large'],
             's-', linewidth=2, markersize=8, label='In Large Clusters (100+)', color='#1f77b4')
    ax2.plot(seq_ids, 100 - comparison_df['pct_proteins_singleton'],
             '^-', linewidth=2, markersize=8, label='In Multi-member Clusters', color='#2ca02c')

    ax2.set_xlabel('Sequence Identity Threshold')
    ax2.set_ylabel('Percentage of Proteins')
    ax2.set_title('Key Metrics vs Sequence Identity')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 100)

    # Highlight optimal region
    ax2.axvspan(0.25, 0.35, alpha=0.1, color='green', label='Optimal range')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved protein distribution plot to {output_path}")
    plt.close()


def plot_cluster_statistics(comparison_df, output_path):
    """Plot cluster statistics across parameters."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    seq_ids = comparison_df['seq_id'].values

    # Plot 1: Total clusters
    ax = axes[0, 0]
    ax.plot(seq_ids, comparison_df['total_clusters'] / 1e6, 'o-', linewidth=2, markersize=8, color='steelblue')
    ax.set_xlabel('Sequence Identity Threshold')
    ax.set_ylabel('Total Clusters (millions)')
    ax.set_title('Total Number of Clusters')
    ax.grid(True, alpha=0.3)

    # Plot 2: Mean cluster size
    ax = axes[0, 1]
    ax.plot(seq_ids, comparison_df['mean_size'], 's-', linewidth=2, markersize=8, color='coral')
    ax.set_xlabel('Sequence Identity Threshold')
    ax.set_ylabel('Mean Cluster Size (proteins)')
    ax.set_title('Average Cluster Size')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='Size = 10')
    ax.legend()

    # Plot 3: Number of useful clusters (10+)
    ax = axes[1, 0]
    useful_clusters = comparison_df['n_medium_10-99'] + comparison_df['n_large_100+']
    ax.bar(seq_ids, useful_clusters / 1000, color='seagreen', alpha=0.7, width=0.08)
    ax.set_xlabel('Sequence Identity Threshold')
    ax.set_ylabel('Clusters with 10+ members (thousands)')
    ax.set_title('Number of "Useful" Clusters for Sampling')
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 4: Singleton percentage
    ax = axes[1, 1]
    ax.plot(seq_ids, comparison_df['pct_singletons'], 'o-', linewidth=2, markersize=8,
            color='#d62728', label='% clusters that are singletons')
    ax.plot(seq_ids, comparison_df['pct_proteins_singleton'], 's-', linewidth=2, markersize=8,
            color='#ff7f0e', label='% proteins in singletons')
    ax.set_xlabel('Sequence Identity Threshold')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Singleton Analysis')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved cluster statistics to {output_path}")
    plt.close()


def plot_cumulative_distributions(data_dict, output_path):
    """Plot cumulative distribution of cluster sizes."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    seq_ids = sorted(data_dict.keys())
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(seq_ids)))

    for seq_id, color in zip(seq_ids, colors):
        df = data_dict[seq_id]
        sizes = np.sort(df['size'].values)

        # Cumulative distribution of cluster sizes
        cumulative = np.arange(1, len(sizes) + 1) / len(sizes) * 100
        ax1.plot(sizes, cumulative, label=f'{seq_id:.1f}', linewidth=2, color=color)

    ax1.set_xscale('log')
    ax1.set_xlabel('Cluster Size (proteins)')
    ax1.set_ylabel('Cumulative % of Clusters')
    ax1.set_title('Cumulative Distribution: Cluster Sizes')
    ax1.legend(title='Seq ID', loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=10, color='red', linestyle='--', alpha=0.5, label='Size = 10')

    # Cumulative distribution of proteins
    for seq_id, color in zip(seq_ids, colors):
        df = data_dict[seq_id]

        # Sort by size and compute cumulative proteins
        df_sorted = df.sort_values('size', ascending=True)
        cumulative_proteins = df_sorted['size'].cumsum() / df_sorted['size'].sum() * 100
        cumulative_clusters = np.arange(1, len(df_sorted) + 1) / len(df_sorted) * 100

        ax2.plot(cumulative_clusters, cumulative_proteins, label=f'{seq_id:.1f}', linewidth=2, color=color)

    ax2.set_xlabel('Cumulative % of Clusters (sorted by size)')
    ax2.set_ylabel('Cumulative % of Proteins')
    ax2.set_title('Protein Coverage by Cluster Rank')
    ax2.legend(title='Seq ID', loc='lower right')
    ax2.grid(True, alpha=0.3)
    ax2.plot([0, 100], [0, 100], 'k--', alpha=0.3, label='Uniform')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved cumulative distributions to {output_path}")
    plt.close()


def main():
    print("=" * 80)
    print("VISUALIZING MMseqs2 CLUSTERING RESULTS")
    print("=" * 80)
    print()

    # Load all clustering results
    runs = [
        ('data/mmseqs_seqid_0p3/cluster_summary.csv', 0.3),
        ('data/mmseqs_seqid_0p4/cluster_summary.csv', 0.4),
        ('data/mmseqs_seqid_0p5/cluster_summary.csv', 0.5),
        ('data/mmseqs_seqid_0p6/cluster_summary.csv', 0.6),
        ('data/mmseqs_seqid_0p7/cluster_summary.csv', 0.7),
    ]

    data_dict = {}
    for path, seq_id in runs:
        if Path(path).exists():
            print(f"Loading {seq_id:.1f} identity clustering...")
            df, sid = load_cluster_data(path, seq_id)
            data_dict[sid] = df

    # Load comparison data
    comparison_df = pd.read_csv('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/mmseqs_parameter_comparison.csv')

    # Create output directory
    output_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/figures')
    output_dir.mkdir(parents=True, exist_ok=True)

    print()
    print("Generating plots...")
    print()

    # Generate plots
    plot_cluster_size_distributions(data_dict, output_dir / 'cluster_size_distributions.png')
    plot_protein_distribution(comparison_df, output_dir / 'protein_distribution.png')
    plot_cluster_statistics(comparison_df, output_dir / 'cluster_statistics.png')
    plot_cumulative_distributions(data_dict, output_dir / 'cumulative_distributions.png')

    print()
    print("=" * 80)
    print("All plots saved to:", output_dir)
    print("=" * 80)
    print()
    print("Plots created:")
    print(f"  1. {output_dir / 'cluster_size_distributions.png'}")
    print(f"  2. {output_dir / 'protein_distribution.png'}")
    print(f"  3. {output_dir / 'cluster_statistics.png'}")
    print(f"  4. {output_dir / 'cumulative_distributions.png'}")
    print()


if __name__ == '__main__':
    main()
