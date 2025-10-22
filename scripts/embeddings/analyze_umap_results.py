#!/usr/bin/env python
"""
Analyze and compare UMAP results from multiple n_neighbors values.

Usage:
    python analyze_umap_results.py results/umap/
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze UMAP results from job array"
    )
    parser.add_argument(
        "results_dir",
        type=str,
        help="Directory containing UMAP result files",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="umap_comparison.png",
        help="Output plot filename",
    )
    return parser.parse_args()


def load_umap_results(results_dir):
    """Load all UMAP result files."""
    results_dir = Path(results_dir)
    result_files = sorted(results_dir.glob("umap_n*.npz"))

    results = []

    for result_file in result_files:
        print(f"Loading {result_file.name}...")

        data = np.load(result_file, allow_pickle=True)

        results.append({
            'n_neighbors': int(data['n_neighbors']),
            'umap_embedding': data['umap_embedding'],
            'gene_ids': data['gene_ids'],
            'genome_ids': data['genome_ids'],
            'n_pcs': int(data['n_pcs']),
            'min_dist': float(data['min_dist']),
            'filename': result_file.name
        })

    print(f"\nLoaded {len(results)} UMAP results")
    return sorted(results, key=lambda x: x['n_neighbors'])


def plot_umap_comparison(results, output_file="umap_comparison.png"):
    """Create comparison plot of all UMAP results."""

    n_results = len(results)
    ncols = min(3, n_results)
    nrows = (n_results + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))

    if n_results == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, result in enumerate(results):
        ax = axes[idx]

        umap_emb = result['umap_embedding']
        n_neighbors = result['n_neighbors']

        # Plot
        ax.scatter(
            umap_emb[:, 0],
            umap_emb[:, 1],
            c='gray',
            alpha=0.3,
            s=1,
            rasterized=True
        )

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title(f'n_neighbors = {n_neighbors}\n{len(umap_emb):,} genes')
        ax.grid(True, alpha=0.3)

        # Add coordinate ranges
        x_range = umap_emb[:, 0].max() - umap_emb[:, 0].min()
        y_range = umap_emb[:, 1].max() - umap_emb[:, 1].min()
        ax.text(
            0.02, 0.98,
            f'Range: {x_range:.1f} × {y_range:.1f}',
            transform=ax.transAxes,
            fontsize=8,
            va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

    # Hide extra subplots
    for idx in range(n_results, len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison plot to: {output_file}")


def compute_spread_metrics(results):
    """Compute spread metrics for each UMAP result."""

    metrics = []

    for result in results:
        umap_emb = result['umap_embedding']

        # Compute spread metrics
        x_range = umap_emb[:, 0].max() - umap_emb[:, 0].min()
        y_range = umap_emb[:, 1].max() - umap_emb[:, 1].min()
        x_std = umap_emb[:, 0].std()
        y_std = umap_emb[:, 1].std()

        # Compute average nearest neighbor distance (sample)
        from sklearn.neighbors import NearestNeighbors
        sample_size = min(10000, len(umap_emb))
        sample_idx = np.random.choice(len(umap_emb), sample_size, replace=False)
        sample = umap_emb[sample_idx]

        nbrs = NearestNeighbors(n_neighbors=2).fit(sample)
        distances, _ = nbrs.kneighbors(sample)
        avg_nn_dist = distances[:, 1].mean()

        metrics.append({
            'n_neighbors': result['n_neighbors'],
            'x_range': x_range,
            'y_range': y_range,
            'x_std': x_std,
            'y_std': y_std,
            'avg_nn_distance': avg_nn_dist,
            'total_variance': x_std**2 + y_std**2
        })

    return pd.DataFrame(metrics)


def main():
    args = parse_args()

    print("=" * 70)
    print("UMAP Results Analysis")
    print("=" * 70)

    # Load results
    results = load_umap_results(args.results_dir)

    if not results:
        print("No UMAP results found!")
        return

    # Create comparison plot
    plot_umap_comparison(results, args.output)

    # Compute metrics
    print("\nComputing spread metrics...")
    metrics_df = compute_spread_metrics(results)

    print("\nSpread Metrics:")
    print(metrics_df.to_string(index=False))

    # Save metrics
    metrics_file = Path(args.results_dir) / "umap_metrics.csv"
    metrics_df.to_csv(metrics_file, index=False)
    print(f"\nSaved metrics to: {metrics_file}")

    # Plot metrics
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Plot 1: Range vs n_neighbors
    ax = axes[0]
    ax.plot(metrics_df['n_neighbors'], metrics_df['x_range'], 'o-', label='X range')
    ax.plot(metrics_df['n_neighbors'], metrics_df['y_range'], 's-', label='Y range')
    ax.set_xlabel('n_neighbors')
    ax.set_ylabel('Coordinate Range')
    ax.set_title('UMAP Spread vs n_neighbors')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Average NN distance
    ax = axes[1]
    ax.plot(metrics_df['n_neighbors'], metrics_df['avg_nn_distance'], 'o-')
    ax.set_xlabel('n_neighbors')
    ax.set_ylabel('Avg Nearest Neighbor Distance')
    ax.set_title('Local Density vs n_neighbors')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    metrics_plot = Path(args.results_dir) / "umap_metrics.png"
    plt.savefig(metrics_plot, dpi=150, bbox_inches='tight')
    print(f"Saved metrics plot to: {metrics_plot}")

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
