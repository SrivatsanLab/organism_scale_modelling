#!/usr/bin/env python
"""
Run UMAP on sampled proteins and create visualizations colored by COG and MMseqs cluster.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import umap
from tqdm import tqdm


def load_cog_annotations(gene_ids):
    """Load COG annotations for genes."""
    # Try to find COG annotations
    cog_file = Path('data/cog_annotations.csv')

    if not cog_file.exists():
        print(f"  COG annotations not found at {cog_file}")
        print(f"  Will create placeholder annotations")
        # Create placeholder
        return pd.DataFrame({
            'gene_id': gene_ids,
            'cog_category': ['Unknown'] * len(gene_ids)
        })

    print(f"  Loading COG annotations from {cog_file}")
    cog_df = pd.read_csv(cog_file)
    return cog_df


def main():
    print("=" * 80)
    print("UMAP AND VISUALIZATION")
    print("=" * 80)
    print()

    # Paths
    input_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7')
    output_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7/figures')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load sampled data
    print("Loading sampled data...")
    data = np.load(input_dir / 'umap_sample.npz', allow_pickle=True)

    gene_ids = data['gene_ids']
    cluster_reps = data['cluster_reps']
    embeddings = data['embeddings']

    print(f"  Loaded {len(gene_ids):,} proteins")
    print(f"  From {len(np.unique(cluster_reps)):,} clusters")
    print(f"  Embedding dimensions: {embeddings.shape[1]}")

    # Load COG annotations
    print("\nLoading COG annotations...")
    cog_df = load_cog_annotations(gene_ids)

    # Create gene_id -> COG mapping
    gene_to_cog = dict(zip(cog_df['gene_id'], cog_df['cog_category']))

    # Assign COG categories to sampled genes
    cog_categories = [gene_to_cog.get(gene, 'Unknown') for gene in gene_ids]

    print(f"  COG categories assigned")
    print(f"  Unique COG categories: {len(set(cog_categories))}")

    # Run UMAP
    print("\nRunning UMAP...")
    print("  This may take a few minutes...")

    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        metric='euclidean',
        random_state=42,
        verbose=True
    )

    umap_embedding = reducer.fit_transform(embeddings)

    print(f"  UMAP complete: {umap_embedding.shape}")

    # Save UMAP coordinates
    print("\nSaving UMAP coordinates...")
    np.savez(
        input_dir / 'umap_coordinates.npz',
        gene_ids=gene_ids,
        cluster_reps=cluster_reps,
        cog_categories=cog_categories,
        umap_coords=umap_embedding
    )

    # Create visualizations
    print("\nCreating visualizations...")

    # Figure 1: Colored by MMseqs cluster
    print("  Plot 1: Colored by MMseqs cluster...")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Get unique clusters
    unique_clusters = np.unique(cluster_reps)
    n_clusters = len(unique_clusters)

    print(f"    Plotting {n_clusters:,} clusters")

    # Create cluster ID mapping
    cluster_to_id = {cluster: idx for idx, cluster in enumerate(unique_clusters)}
    cluster_ids = [cluster_to_id[cluster] for cluster in cluster_reps]

    # Plot with colormap
    scatter = ax.scatter(
        umap_embedding[:, 0],
        umap_embedding[:, 1],
        c=cluster_ids,
        cmap='tab20',
        s=1,
        alpha=0.5,
        rasterized=True
    )

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title(f'UMAP colored by MMseqs2 Cluster (70% seq ID, {n_clusters:,} clusters)')

    plt.tight_layout()
    plt.savefig(output_dir / 'umap_by_mmseqs_cluster.png', dpi=300, bbox_inches='tight')
    print(f"    Saved to {output_dir / 'umap_by_mmseqs_cluster.png'}")
    plt.close()

    # Figure 2: Colored by COG category
    print("  Plot 2: Colored by COG category...")

    fig, ax = plt.subplots(figsize=(14, 10))

    # Get unique COG categories
    unique_cogs = sorted(set(cog_categories))
    n_cogs = len(unique_cogs)

    print(f"    Plotting {n_cogs} COG categories")

    # Create COG color mapping
    cog_palette = sns.color_palette("husl", n_cogs)
    cog_to_color = {cog: color for cog, color in zip(unique_cogs, cog_palette)}
    colors = [cog_to_color[cog] for cog in cog_categories]

    # Plot
    ax.scatter(
        umap_embedding[:, 0],
        umap_embedding[:, 1],
        c=colors,
        s=1,
        alpha=0.5,
        rasterized=True
    )

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title(f'UMAP colored by COG Category ({n_cogs} categories)')

    # Add legend (only if not too many categories)
    if n_cogs <= 25:
        # Create legend elements
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=cog_to_color[cog], label=cog) for cog in unique_cogs]
        ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
                 fontsize=8, frameon=True)

    plt.tight_layout()
    plt.savefig(output_dir / 'umap_by_cog_category.png', dpi=300, bbox_inches='tight')
    print(f"    Saved to {output_dir / 'umap_by_cog_category.png'}")
    plt.close()

    # Figure 3: Combined view (2x1 subplots)
    print("  Plot 3: Combined view...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))

    # Subplot 1: MMseqs clusters
    ax1.scatter(
        umap_embedding[:, 0],
        umap_embedding[:, 1],
        c=cluster_ids,
        cmap='tab20',
        s=1,
        alpha=0.5,
        rasterized=True
    )
    ax1.set_xlabel('UMAP 1')
    ax1.set_ylabel('UMAP 2')
    ax1.set_title(f'MMseqs2 Clusters (70% ID, n={n_clusters:,})')

    # Subplot 2: COG categories
    ax2.scatter(
        umap_embedding[:, 0],
        umap_embedding[:, 1],
        c=colors,
        s=1,
        alpha=0.5,
        rasterized=True
    )
    ax2.set_xlabel('UMAP 1')
    ax2.set_ylabel('UMAP 2')
    ax2.set_title(f'COG Categories (n={n_cogs})')

    plt.tight_layout()
    plt.savefig(output_dir / 'umap_combined.png', dpi=300, bbox_inches='tight')
    print(f"    Saved to {output_dir / 'umap_combined.png'}")
    plt.close()

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total proteins visualized: {len(gene_ids):,}")
    print(f"MMseqs2 clusters: {n_clusters:,}")
    print(f"COG categories: {n_cogs}")
    print(f"\nOutput files saved to: {output_dir}")
    print(f"  1. umap_by_mmseqs_cluster.png")
    print(f"  2. umap_by_cog_category.png")
    print(f"  3. umap_combined.png")
    print(f"\nUMAP coordinates saved to: {input_dir / 'umap_coordinates.npz'}")
    print("=" * 80)


if __name__ == '__main__':
    main()
