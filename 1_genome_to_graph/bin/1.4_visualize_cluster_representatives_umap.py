#!/usr/bin/env python
"""
Create UMAP visualization using cluster representative coordinates.

Much faster than computing cluster means - uses the UMAP coordinates
of the cluster representative protein for each cluster.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Visualize UMAP with cluster representatives"
    )
    parser.add_argument(
        "--umap-file",
        type=str,
        required=True,
        help="Path to UMAP embeddings .npz file",
    )
    parser.add_argument(
        "--cluster-stats",
        type=str,
        required=True,
        help="Path to cluster statistics CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for plots",
    )
    return parser.parse_args()


def load_data(umap_file, cluster_stats_file):
    """Load UMAP embeddings and cluster stats."""
    print("Loading data...")

    # Load UMAP embeddings
    print(f"  Loading UMAP from {umap_file}...")
    umap_data = np.load(umap_file, allow_pickle=True)
    umap_coords = umap_data['umap_embedding']
    gene_ids = umap_data['gene_ids']
    genome_ids = umap_data['genome_ids']

    # Reconstruct full gene IDs
    full_gene_ids = np.array([f"{genome}_{gene}" for genome, gene in zip(genome_ids, gene_ids)])

    print(f"    Loaded {len(full_gene_ids):,} UMAP coordinates")

    # Create gene_id -> UMAP index mapping
    gene_to_umap_idx = {gene: idx for idx, gene in enumerate(full_gene_ids)}

    # Load cluster statistics
    print(f"  Loading cluster stats from {cluster_stats_file}...")
    cluster_stats = pd.read_csv(cluster_stats_file)
    print(f"    Loaded {len(cluster_stats):,} clusters")

    # Get UMAP coordinates for cluster representatives
    print("  Mapping cluster representatives to UMAP coordinates...")
    cluster_coords = []
    cluster_sizes = []
    cluster_tightness = []
    cluster_reps_found = []

    for _, row in cluster_stats.iterrows():
        rep = row['cluster_representative']
        if rep in gene_to_umap_idx:
            idx = gene_to_umap_idx[rep]
            cluster_coords.append(umap_coords[idx])
            cluster_sizes.append(row['cluster_size'])
            cluster_tightness.append(row['mean_std'])
            cluster_reps_found.append(rep)

    cluster_coords = np.array(cluster_coords)
    cluster_sizes = np.array(cluster_sizes)
    cluster_tightness = np.array(cluster_tightness)

    print(f"    Found UMAP coordinates for {len(cluster_coords):,} / {len(cluster_stats):,} cluster representatives")

    return cluster_coords, cluster_sizes, cluster_tightness


def plot_cluster_umap(cluster_coords, cluster_sizes, cluster_tightness, output_dir):
    """Create UMAP plots with one point per cluster."""
    print("\nCreating cluster UMAP visualizations...")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Plot 1: Colored by cluster size
    print("  Creating size-colored plot...")
    fig, ax = plt.subplots(figsize=(14, 10))

    scatter = ax.scatter(
        cluster_coords[:, 0],
        cluster_coords[:, 1],
        c=np.log10(cluster_sizes),
        cmap='viridis',
        s=20,
        alpha=0.6,
        edgecolors='none'
    )

    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title(f'Cluster Representatives UMAP (n={len(cluster_coords):,} clusters, colored by log10(size))', fontsize=16)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('log10(Cluster Size)', fontsize=12)

    # Add statistics box
    stats_text = f'Total clusters: {len(cluster_coords):,}\n'
    stats_text += f'Mean size: {cluster_sizes.mean():.1f}\n'
    stats_text += f'Median size: {np.median(cluster_sizes):.0f}\n'
    stats_text += f'Size range: {cluster_sizes.min()}-{cluster_sizes.max():,}'

    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    output_file = output_dir / 'cluster_reps_umap_by_size.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    Saved {output_file}")
    plt.close()

    # Plot 2: Colored by cluster tightness (variance in ESM embeddings)
    print("  Creating tightness-colored plot...")
    fig, ax = plt.subplots(figsize=(14, 10))

    scatter = ax.scatter(
        cluster_coords[:, 0],
        cluster_coords[:, 1],
        c=cluster_tightness,
        cmap='RdYlGn_r',  # Red = high variance, Green = tight
        s=20,
        alpha=0.6,
        edgecolors='none',
        vmin=0,
        vmax=0.02
    )

    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title(f'Cluster Representatives UMAP (n={len(cluster_coords):,} clusters, colored by ESM variance)', fontsize=16)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Mean Std in ESM Embeddings (lower = tighter)', fontsize=12)

    # Add statistics box
    stats_text = f'Mean variance: {cluster_tightness.mean():.4f}\n'
    stats_text += f'Median: {np.median(cluster_tightness):.4f}\n'
    stats_text += f'75th percentile: {np.percentile(cluster_tightness, 75):.4f}\n'
    stats_text += f'Range: {cluster_tightness.min():.4f}-{cluster_tightness.max():.4f}'

    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    output_file = output_dir / 'cluster_reps_umap_by_tightness.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    Saved {output_file}")
    plt.close()

    # Plot 3: Size vs tightness (point size = cluster size, color = variance)
    print("  Creating combined size/tightness plot...")
    fig, ax = plt.subplots(figsize=(14, 10))

    # Normalize sizes for visualization (10-100 point size range)
    size_normalized = 10 + 90 * (cluster_sizes - cluster_sizes.min()) / (cluster_sizes.max() - cluster_sizes.min())

    scatter = ax.scatter(
        cluster_coords[:, 0],
        cluster_coords[:, 1],
        c=cluster_tightness,
        s=size_normalized,
        cmap='RdYlGn_r',
        alpha=0.5,
        edgecolors='black',
        linewidths=0.3,
        vmin=0,
        vmax=0.02
    )

    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title(f'Cluster Representatives UMAP (n={len(cluster_coords):,}, size=cluster size, color=ESM variance)', fontsize=16)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Mean Std in ESM Embeddings (lower = tighter)', fontsize=12)

    plt.tight_layout()
    output_file = output_dir / 'cluster_reps_umap_size_and_tightness.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    Saved {output_file}")
    plt.close()


def main():
    args = parse_args()

    print("=" * 80)
    print("CLUSTER REPRESENTATIVE UMAP VISUALIZATION")
    print("=" * 80)
    print(f"UMAP file: {args.umap_file}")
    print(f"Cluster stats: {args.cluster_stats}")
    print(f"Output directory: {args.output_dir}")
    print("=" * 80)

    # Load data
    cluster_coords, cluster_sizes, cluster_tightness = load_data(
        args.umap_file, args.cluster_stats
    )

    # Create visualizations
    plot_cluster_umap(cluster_coords, cluster_sizes, cluster_tightness, args.output_dir)

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print(f"Created 3 cluster-level UMAP plots in: {args.output_dir}")
    print("  - cluster_reps_umap_by_size.png")
    print("  - cluster_reps_umap_by_tightness.png (colored by ESM embedding variance)")
    print("  - cluster_reps_umap_size_and_tightness.png")
    print("=" * 80)


if __name__ == '__main__':
    main()
