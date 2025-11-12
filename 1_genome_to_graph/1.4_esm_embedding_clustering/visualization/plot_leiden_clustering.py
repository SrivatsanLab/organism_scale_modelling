#!/usr/bin/env python
"""
Plot a single Leiden clustering result with COG annotations.

Usage:
    python plot_leiden_clustering.py --clustering <file.npz> --umap-n <N> --output <output.png>
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
sys.path.insert(0, '/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling')

from scripts.embeddings.plot_all_clusterings import (
    load_cog_annotations_for_genes,
    COG_COLORS
)

def load_umap_embedding(n_neighbors, subsample=1000000):
    """Load UMAP embedding."""
    umap_file = f'results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n{n_neighbors}_subsample{subsample}.npz'
    print(f"Loading UMAP n={n_neighbors}...")
    data = np.load(umap_file, allow_pickle=True)
    return data['umap_embedding'], data['gene_ids'], data['genome_ids']


def align_clustering_to_umap(cluster_labels, cluster_gene_ids, umap_gene_ids):
    """Align clustering labels to UMAP gene order."""
    gene_to_cluster = dict(zip(cluster_gene_ids, cluster_labels))
    aligned_labels = np.array([gene_to_cluster.get(gid, -1) for gid in umap_gene_ids])
    return aligned_labels


def map_cogs_to_genes(gene_ids, cog_lookup):
    """Map COG categories to gene IDs."""
    cog_categories = np.array([cog_lookup.get(gid, 'No COG') for gid in gene_ids])
    n_annotated = (cog_categories != 'No COG').sum()
    pct_annotated = 100 * n_annotated / len(cog_categories)
    print(f"  COG annotation rate: {pct_annotated:.1f}% ({n_annotated:,}/{len(cog_categories):,})")
    return cog_categories


def plot_clustering_with_cog(umap_coords, cluster_labels, cog_categories,
                              clustering_name, n_neighbors, output_file):
    """Create a figure with 2 subplots: clusters and COG annotations."""
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))

    # === LEFT PLOT: Clusters ===
    ax = axes[0]

    mask_noise = cluster_labels == -1
    mask_clustered = ~mask_noise

    if np.any(mask_noise):
        ax.scatter(
            umap_coords[mask_noise, 0],
            umap_coords[mask_noise, 1],
            c='lightgray',
            s=0.1,
            alpha=0.3,
            rasterized=True,
            label='Noise'
        )

    if np.any(mask_clustered):
        scatter = ax.scatter(
            umap_coords[mask_clustered, 0],
            umap_coords[mask_clustered, 1],
            c=cluster_labels[mask_clustered],
            s=0.5,
            alpha=0.5,
            cmap='tab20',
            rasterized=True
        )

    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = np.sum(mask_noise)
    pct_noise = 100 * n_noise / len(cluster_labels)

    ax.set_title(
        f'{clustering_name}\n{n_clusters:,} clusters, {n_noise:,} noise ({pct_noise:.1f}%)',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel(f'UMAP 1 (n_neighbors={n_neighbors})', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.grid(True, alpha=0.2)

    # === RIGHT PLOT: COG Categories ===
    ax = axes[1]

    for cog, color in COG_COLORS.items():
        mask = cog_categories == cog
        if mask.sum() > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=color,
                s=0.5,
                alpha=0.3,
                label=cog,
                rasterized=True
            )

    n_annotated = (cog_categories != 'No COG').sum()
    pct_annotated = 100 * n_annotated / len(cog_categories)

    ax.set_title(
        f'COG Categories\n{pct_annotated:.1f}% annotated ({n_annotated:,}/{len(cog_categories):,})',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel(f'UMAP 1 (n_neighbors={n_neighbors})', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.grid(True, alpha=0.2)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
              fontsize=8, markerscale=3)

    plt.tight_layout()

    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Plot Leiden clustering with COG annotations')
    parser.add_argument('--clustering', required=True, help='Clustering .npz file')
    parser.add_argument('--umap-n', type=int, required=True, help='UMAP n_neighbors')
    parser.add_argument('--output', required=True, help='Output PNG file')

    args = parser.parse_args()

    print("=" * 80)
    print(f"Plotting: {Path(args.clustering).name}")
    print(f"UMAP n_neighbors: {args.umap_n}")
    print("=" * 80)
    print()

    # Load clustering
    clustering_data = np.load(args.clustering, allow_pickle=True)
    cluster_labels = clustering_data['labels']
    cluster_gene_ids = clustering_data['gene_ids']
    n_clusters = clustering_data['n_clusters']

    print(f"Clustering: {n_clusters:,} clusters, {len(cluster_gene_ids):,} genes")
    print()

    # Load UMAP
    umap_coords, umap_gene_ids, umap_genome_ids = load_umap_embedding(args.umap_n)
    print(f"UMAP: {len(umap_gene_ids):,} genes")
    print()

    # Align clustering to UMAP
    aligned_labels = align_clustering_to_umap(cluster_labels, cluster_gene_ids, umap_gene_ids)
    print(f"Aligned: {(aligned_labels != -1).sum():,} genes in clusters")
    print()

    # Load COG annotations
    print("Loading COG annotations...")
    cog_lookup = load_cog_annotations_for_genes(umap_gene_ids, umap_genome_ids)
    cog_categories = map_cogs_to_genes(umap_gene_ids, cog_lookup)
    print()

    # Create plot
    clustering_name = Path(args.clustering).stem.replace('clusters_', '')
    plot_clustering_with_cog(
        umap_coords=umap_coords,
        cluster_labels=aligned_labels,
        cog_categories=cog_categories,
        clustering_name=clustering_name,
        n_neighbors=args.umap_n,
        output_file=Path(args.output)
    )

    print()
    print("=" * 80)
    print("✓ Done!")
    print("=" * 80)


if __name__ == '__main__':
    main()
