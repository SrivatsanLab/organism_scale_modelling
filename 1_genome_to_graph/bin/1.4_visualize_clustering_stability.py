#!/usr/bin/env python
"""
Visualize clustering stability results.

Creates comprehensive visualizations of co-clustering patterns:
1. Distribution of co-clustering rates
2. Stable vs unstable gene pairs
3. Comparison across configurations
4. Network visualization of highly stable pairs
5. Relationship between stability and COG annotations

Usage:
    python visualize_clustering_stability.py \
        --stability-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability \
        --output-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/stability_analysis
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import sys


def load_stability_result(npz_file):
    """Load a single stability result file."""
    data = np.load(npz_file, allow_pickle=True)

    # Extract metadata
    metadata = {
        'file': npz_file.name,
        'resolution': float(data['resolution']),
        'n_neighbors': int(data['n_neighbors']),
        'cog_only': bool(data['cog_only']),
        'n_subsamples': int(data['n_subsamples']),
        'subsample_size': int(data['subsample_size']),
        'stability_mean': float(data['mean_coclustering_rate']),
        'stability_median': float(data['median_coclustering_rate']),
        'stability_std': float(data['std_coclustering_rate']),
        'n_pairs': int(data['n_pairs'])
    }

    # Extract pair-level data
    gene_pairs = data['gene_pairs']
    cooccur_counts = data['cooccur_counts']
    cocluster_counts = data['cocluster_counts']

    # Compute co-clustering rates
    coclustering_rates = cocluster_counts / cooccur_counts

    return metadata, gene_pairs, cooccur_counts, cocluster_counts, coclustering_rates


def plot_coclustering_distribution(stability_results, output_dir):
    """
    Plot distribution of co-clustering rates for each configuration.
    """
    print("Creating co-clustering rate distributions...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Distribution of Gene Pair Co-clustering Rates Across Subsamples',
                 fontsize=14, y=0.995)

    for i, (metadata, gene_pairs, cooccur, cocluster, rates) in enumerate(stability_results):
        if i >= 4:
            break

        ax = axes.flat[i]

        # Create histogram
        ax.hist(rates, bins=50, alpha=0.7, color='steelblue', edgecolor='black')

        # Add statistics
        mean_rate = metadata['stability_mean']
        median_rate = metadata['stability_median']

        ax.axvline(mean_rate, color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_rate:.3f}')
        ax.axvline(median_rate, color='orange', linestyle='--', linewidth=2,
                   label=f'Median: {median_rate:.3f}')

        # Labels
        version = 'COG-only' if metadata['cog_only'] else 'All genes'
        title = f"Res={metadata['resolution']:.0f}, n={metadata['n_neighbors']}, {version}"
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Co-clustering Rate', fontsize=10)
        ax.set_ylabel('Number of Gene Pairs', fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Add text with pair count
        ax.text(0.02, 0.98, f'N pairs: {len(rates):,}',
                transform=ax.transAxes, verticalalignment='top',
                fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_file = output_dir / 'coclustering_rate_distributions.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def plot_stability_comparison(stability_results, output_dir):
    """
    Compare stability metrics across configurations.
    """
    print("Creating stability comparison plot...")

    # Extract summary statistics
    summaries = []
    for metadata, _, _, _, rates in stability_results:
        # Compute additional percentiles
        percentiles = np.percentile(rates, [10, 25, 50, 75, 90, 95])

        summaries.append({
            'config': f"Res={metadata['resolution']:.0f}\nn={metadata['n_neighbors']}\n" +
                     ('COG' if metadata['cog_only'] else 'All'),
            'resolution': metadata['resolution'],
            'mean': metadata['stability_mean'],
            'median': metadata['stability_median'],
            'std': metadata['stability_std'],
            'p10': percentiles[0],
            'p25': percentiles[1],
            'p50': percentiles[2],
            'p75': percentiles[3],
            'p90': percentiles[4],
            'p95': percentiles[5],
        })

    df = pd.DataFrame(summaries)

    # Create box plot style visualization
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    x = np.arange(len(df))

    # Plot error bars showing different percentiles
    ax.errorbar(x, df['median'],
                yerr=[df['median'] - df['p25'], df['p75'] - df['median']],
                fmt='o', markersize=8, capsize=5, capthick=2,
                label='Median ± IQR', color='steelblue', linewidth=2)

    ax.errorbar(x + 0.15, df['mean'],
                yerr=df['std'],
                fmt='s', markersize=6, capsize=4, capthick=1.5,
                label='Mean ± Std', color='coral', linewidth=1.5, alpha=0.7)

    # Reference line at 0.5 (random expectation)
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=1.5,
               label='Random baseline (0.5)', alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(df['config'], fontsize=9)
    ax.set_ylabel('Co-clustering Rate', fontsize=12)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.set_title('Clustering Stability Comparison Across Configurations', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1)

    plt.tight_layout()
    output_file = output_dir / 'stability_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def identify_stable_unstable_pairs(stability_results, output_dir,
                                   stable_threshold=0.9, unstable_threshold=0.3):
    """
    Identify consistently stable and unstable gene pairs.
    """
    print(f"Identifying stable (>{stable_threshold}) and unstable (<{unstable_threshold}) pairs...")

    for metadata, gene_pairs, cooccur, cocluster, rates in stability_results:
        # Categorize pairs
        stable_mask = rates >= stable_threshold
        unstable_mask = rates <= unstable_threshold

        n_stable = stable_mask.sum()
        n_unstable = unstable_mask.sum()
        n_total = len(rates)

        version = 'cogonly' if metadata['cog_only'] else 'all'
        config_name = f"res{metadata['resolution']:.0f}_nn{metadata['n_neighbors']}_{version}"

        print(f"\n  Config: {config_name}")
        print(f"    Stable pairs (≥{stable_threshold}): {n_stable:,} ({100*n_stable/n_total:.1f}%)")
        print(f"    Unstable pairs (≤{unstable_threshold}): {n_unstable:,} ({100*n_unstable/n_total:.1f}%)")
        print(f"    Intermediate: {n_total - n_stable - n_unstable:,} ({100*(n_total - n_stable - n_unstable)/n_total:.1f}%)")

        # Save lists of stable/unstable pairs
        stable_pairs_file = output_dir / f'stable_pairs_{config_name}.txt'
        with open(stable_pairs_file, 'w') as f:
            f.write(f"# Stable gene pairs (co-clustering rate ≥ {stable_threshold})\n")
            f.write(f"# Config: {config_name}\n")
            f.write(f"# Total: {n_stable:,} pairs\n")
            f.write(f"#\n")
            f.write("gene1\tgene2\tcooccur\tcocluster\trate\n")

            stable_indices = np.where(stable_mask)[0]
            for idx in stable_indices:
                g1, g2 = gene_pairs[idx]
                f.write(f"{g1}\t{g2}\t{cooccur[idx]}\t{cocluster[idx]}\t{rates[idx]:.4f}\n")

        print(f"    Saved: {stable_pairs_file.name}")


def plot_stability_vs_cooccurrence(stability_results, output_dir):
    """
    Plot relationship between co-occurrence count and stability.

    Hypothesis: Pairs that co-occur more often should show clearer stability pattern.
    """
    print("Creating stability vs co-occurrence plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Co-clustering Rate vs Co-occurrence Frequency', fontsize=14, y=0.995)

    for i, (metadata, gene_pairs, cooccur, cocluster, rates) in enumerate(stability_results):
        if i >= 4:
            break

        ax = axes.flat[i]

        # Create 2D histogram / hexbin plot
        hexbin = ax.hexbin(cooccur, rates, gridsize=30, cmap='viridis',
                          mincnt=1, linewidths=0.2)

        # Add colorbar
        cb = plt.colorbar(hexbin, ax=ax)
        cb.set_label('Number of pairs', fontsize=9)

        # Add mean rate per co-occurrence count
        unique_cooccur = np.unique(cooccur)
        mean_rates = [rates[cooccur == c].mean() for c in unique_cooccur]
        ax.plot(unique_cooccur, mean_rates, 'r-', linewidth=2, alpha=0.7,
                label='Mean rate')

        # Labels
        version = 'COG-only' if metadata['cog_only'] else 'All genes'
        title = f"Res={metadata['resolution']:.0f}, n={metadata['n_neighbors']}, {version}"
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Co-occurrence Count', fontsize=10)
        ax.set_ylabel('Co-clustering Rate', fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)

    plt.tight_layout()
    output_file = output_dir / 'stability_vs_cooccurrence.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def analyze_stability_by_cog(stability_results, output_dir):
    """
    Analyze if gene pairs with same COG category are more stable.

    Only works for COG-only configurations.
    """
    print("Analyzing stability by COG category agreement...")

    from evaluate_clustering_quality import load_cog_annotations

    for metadata, gene_pairs, cooccur, cocluster, rates in stability_results:
        if not metadata['cog_only']:
            print(f"  Skipping non-COG config: {metadata['file']}")
            continue

        version = 'cogonly'
        config_name = f"res{metadata['resolution']:.0f}_nn{metadata['n_neighbors']}_{version}"

        print(f"\n  Config: {config_name}")

        # Load COG annotations
        # Need to get gene_ids and genome_ids from PCA file
        pca_file = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz')
        pca_data = np.load(pca_file, allow_pickle=True)
        gene_ids = pca_data['gene_ids']
        genome_ids = pca_data['genome_ids']

        print("    Loading COG annotations...")
        cog_lookup = load_cog_annotations(gene_ids, genome_ids)

        # Categorize pairs by COG agreement
        same_cog_rates = []
        diff_cog_rates = []

        for i, (g1, g2) in enumerate(gene_pairs):
            if g1 in cog_lookup and g2 in cog_lookup:
                if cog_lookup[g1] == cog_lookup[g2]:
                    same_cog_rates.append(rates[i])
                else:
                    diff_cog_rates.append(rates[i])

        if len(same_cog_rates) == 0 or len(diff_cog_rates) == 0:
            print("    Insufficient data for COG comparison")
            continue

        same_cog_rates = np.array(same_cog_rates)
        diff_cog_rates = np.array(diff_cog_rates)

        print(f"    Pairs with same COG: {len(same_cog_rates):,}")
        print(f"      Mean rate: {same_cog_rates.mean():.4f}")
        print(f"      Median rate: {np.median(same_cog_rates):.4f}")
        print(f"    Pairs with different COG: {len(diff_cog_rates):,}")
        print(f"      Mean rate: {diff_cog_rates.mean():.4f}")
        print(f"      Median rate: {np.median(diff_cog_rates):.4f}")

        # Create comparison plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        # Side-by-side histograms
        bins = np.linspace(0, 1, 51)
        ax.hist(same_cog_rates, bins=bins, alpha=0.6, label='Same COG category',
                color='green', edgecolor='black')
        ax.hist(diff_cog_rates, bins=bins, alpha=0.6, label='Different COG category',
                color='red', edgecolor='black')

        # Add mean lines
        ax.axvline(same_cog_rates.mean(), color='darkgreen', linestyle='--',
                   linewidth=2, label=f'Same COG mean: {same_cog_rates.mean():.3f}')
        ax.axvline(diff_cog_rates.mean(), color='darkred', linestyle='--',
                   linewidth=2, label=f'Diff COG mean: {diff_cog_rates.mean():.3f}')

        ax.set_xlabel('Co-clustering Rate', fontsize=12)
        ax.set_ylabel('Number of Gene Pairs', fontsize=12)
        ax.set_title(f'Stability by COG Agreement - {config_name}', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        output_file = output_dir / f'stability_by_cog_{config_name}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {output_file.name}")


def plot_gene_stability_network(stability_results, output_dir,
                                min_rate=0.9, max_genes=500, max_edges=5000):
    """
    Create network visualization of highly stable gene pairs.

    Shows genes as nodes, edges for high co-clustering rates.
    Only visualizes a subset to keep it readable.
    """
    print(f"Creating network visualization of stable pairs (rate ≥ {min_rate})...")

    try:
        import networkx as nx
    except ImportError:
        print("  Warning: networkx not available, skipping network plots")
        return

    for metadata, gene_pairs, cooccur, cocluster, rates in stability_results[:2]:  # Only top 2 configs
        # Filter to highly stable pairs
        stable_mask = rates >= min_rate
        stable_pairs = gene_pairs[stable_mask]
        stable_rates = rates[stable_mask]

        if len(stable_pairs) == 0:
            print(f"  No pairs with rate ≥ {min_rate} for {metadata['file']}")
            continue

        version = 'cogonly' if metadata['cog_only'] else 'all'
        config_name = f"res{metadata['resolution']:.0f}_nn{metadata['n_neighbors']}_{version}"

        print(f"\n  Config: {config_name}")
        print(f"    Stable pairs: {len(stable_pairs):,}")

        # Build network
        G = nx.Graph()
        for (g1, g2), rate in zip(stable_pairs, stable_rates):
            G.add_edge(g1, g2, weight=rate)

        print(f"    Network: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

        # If too large, extract largest connected component
        if G.number_of_nodes() > max_genes or G.number_of_edges() > max_edges:
            components = list(nx.connected_components(G))
            components.sort(key=len, reverse=True)

            # Take largest component that's not too big
            G_sub = None
            for comp in components:
                G_comp = G.subgraph(comp).copy()
                if G_comp.number_of_nodes() <= max_genes and G_comp.number_of_edges() <= max_edges:
                    G_sub = G_comp
                    break

            if G_sub is None:
                # Just take top max_genes nodes by degree
                degrees = dict(G.degree())
                top_nodes = sorted(degrees.keys(), key=lambda x: degrees[x], reverse=True)[:max_genes]
                G_sub = G.subgraph(top_nodes).copy()

            G = G_sub
            print(f"    Subsampled to: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

        # Layout
        print("    Computing layout...")
        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)

        # Create plot
        fig, ax = plt.subplots(1, 1, figsize=(14, 12))

        # Draw edges
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5, edge_color=edge_weights,
                              edge_cmap=plt.cm.YlOrRd, edge_vmin=min_rate, edge_vmax=1.0,
                              ax=ax)

        # Draw nodes sized by degree
        degrees = dict(G.degree())
        node_sizes = [20 + 10 * degrees[node] for node in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='steelblue',
                              alpha=0.7, ax=ax)

        # Add title
        ax.set_title(f'Network of Stable Gene Pairs (rate ≥ {min_rate}) - {config_name}',
                    fontsize=13)
        ax.axis('off')

        # Add colorbar for edge weights
        sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd,
                                   norm=plt.Normalize(vmin=min_rate, vmax=1.0))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Co-clustering Rate', fontsize=11)

        plt.tight_layout()
        output_file = output_dir / f'stable_network_{config_name}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {output_file.name}")


def create_summary_report(stability_results, output_dir):
    """
    Create a summary report of all stability analyses.
    """
    print("\nCreating summary report...")

    report_file = output_dir / 'stability_summary_report.txt'

    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("Clustering Stability Analysis Summary\n")
        f.write("="*80 + "\n\n")

        for metadata, gene_pairs, cooccur, cocluster, rates in stability_results:
            version = 'COG-only' if metadata['cog_only'] else 'All genes'

            f.write(f"\nConfiguration: Resolution={metadata['resolution']:.0f}, "
                   f"n_neighbors={metadata['n_neighbors']}, {version}\n")
            f.write("-"*80 + "\n")

            f.write(f"Parameters:\n")
            f.write(f"  N subsamples: {metadata['n_subsamples']}\n")
            f.write(f"  Subsample size: {metadata['subsample_size']:,}\n")
            f.write(f"  Gene pairs analyzed: {metadata['n_pairs']:,}\n\n")

            f.write(f"Co-clustering rate statistics:\n")
            f.write(f"  Mean: {metadata['stability_mean']:.4f}\n")
            f.write(f"  Median: {metadata['stability_median']:.4f}\n")
            f.write(f"  Std: {metadata['stability_std']:.4f}\n")
            f.write(f"  Min: {rates.min():.4f}\n")
            f.write(f"  Max: {rates.max():.4f}\n\n")

            # Percentiles
            percentiles = np.percentile(rates, [10, 25, 50, 75, 90, 95, 99])
            f.write(f"Percentiles:\n")
            f.write(f"  10th: {percentiles[0]:.4f}\n")
            f.write(f"  25th: {percentiles[1]:.4f}\n")
            f.write(f"  50th: {percentiles[2]:.4f}\n")
            f.write(f"  75th: {percentiles[3]:.4f}\n")
            f.write(f"  90th: {percentiles[4]:.4f}\n")
            f.write(f"  95th: {percentiles[5]:.4f}\n")
            f.write(f"  99th: {percentiles[6]:.4f}\n\n")

            # Categorize pairs
            n_very_stable = (rates >= 0.9).sum()
            n_stable = ((rates >= 0.7) & (rates < 0.9)).sum()
            n_moderate = ((rates >= 0.5) & (rates < 0.7)).sum()
            n_unstable = ((rates >= 0.3) & (rates < 0.5)).sum()
            n_very_unstable = (rates < 0.3).sum()
            n_total = len(rates)

            f.write(f"Pair stability categories:\n")
            f.write(f"  Very stable (≥0.9): {n_very_stable:,} ({100*n_very_stable/n_total:.1f}%)\n")
            f.write(f"  Stable (0.7-0.9): {n_stable:,} ({100*n_stable/n_total:.1f}%)\n")
            f.write(f"  Moderate (0.5-0.7): {n_moderate:,} ({100*n_moderate/n_total:.1f}%)\n")
            f.write(f"  Unstable (0.3-0.5): {n_unstable:,} ({100*n_unstable/n_total:.1f}%)\n")
            f.write(f"  Very unstable (<0.3): {n_very_unstable:,} ({100*n_very_unstable/n_total:.1f}%)\n")
            f.write("\n" + "="*80 + "\n")

    print(f"  Saved: {report_file.name}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize clustering stability results'
    )
    parser.add_argument('--stability-dir', type=str,
                       default='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability',
                       help='Directory containing stability .npz files')
    parser.add_argument('--output-dir', type=str,
                       default='results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/stability_analysis',
                       help='Output directory for plots')
    parser.add_argument('--stable-threshold', type=float, default=0.9,
                       help='Threshold for identifying stable pairs')
    parser.add_argument('--unstable-threshold', type=float, default=0.3,
                       help='Threshold for identifying unstable pairs')

    args = parser.parse_args()

    stability_dir = Path(args.stability_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("Clustering Stability Visualization")
    print("="*80)
    print(f"Stability directory: {stability_dir}")
    print(f"Output directory: {output_dir}")
    print()

    # Load all stability results
    stability_files = sorted(stability_dir.glob('stability_*.npz'))

    if len(stability_files) == 0:
        print(f"Error: No stability files found in {stability_dir}")
        return

    print(f"Found {len(stability_files)} stability result files\n")

    stability_results = []
    for npz_file in stability_files:
        print(f"Loading: {npz_file.name}")
        result = load_stability_result(npz_file)
        stability_results.append(result)

    print()

    # Create visualizations
    plot_coclustering_distribution(stability_results, output_dir)
    plot_stability_comparison(stability_results, output_dir)
    plot_stability_vs_cooccurrence(stability_results, output_dir)

    identify_stable_unstable_pairs(stability_results, output_dir,
                                   args.stable_threshold, args.unstable_threshold)

    analyze_stability_by_cog(stability_results, output_dir)

    plot_gene_stability_network(stability_results, output_dir)

    create_summary_report(stability_results, output_dir)

    print("\n" + "="*80)
    print("Visualization complete!")
    print("="*80)
    print(f"Output directory: {output_dir}")


if __name__ == '__main__':
    main()
