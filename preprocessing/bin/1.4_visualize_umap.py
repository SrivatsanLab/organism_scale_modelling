#!/usr/bin/env python
"""
Create visualizations from UMAP embeddings.

This script creates various plots to explore the UMAP embedding space:
1. Density plot of overall distribution
2. Colored by genome
3. Colored by MMseqs2 cluster (representative)
4. Random sample for clearer visualization
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
import argparse


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Visualize UMAP embeddings"
    )
    parser.add_argument(
        "--umap-file",
        type=str,
        required=True,
        help="Path to UMAP embeddings .npz file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for plots",
    )
    parser.add_argument(
        "--mmseqs-clusters",
        type=str,
        default="results/1_genome_to_graph/1.3_msa/mmseqs_seqid_0p7/clusters.tsv",
        help="Path to MMseqs2 cluster assignments",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=100000,
        help="Number of points to plot in sample visualizations (default: 100k)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    return parser.parse_args()


def load_umap_embeddings(umap_file):
    """Load UMAP embeddings from file."""
    print(f"Loading UMAP embeddings from {umap_file}...")
    data = np.load(umap_file, allow_pickle=True)

    umap_coords = data['umap_embedding']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Loaded {len(gene_ids):,} proteins")
    print(f"  UMAP shape: {umap_coords.shape}")
    print(f"  Unique genomes: {len(np.unique(genome_ids)):,}")

    return umap_coords, gene_ids, genome_ids


def load_mmseqs_clusters(cluster_file):
    """Load MMseqs2 cluster assignments."""
    print(f"\nLoading MMseqs2 clusters from {cluster_file}...")

    # Read cluster file (representative \t member)
    clusters = pd.read_csv(cluster_file, sep='\t', names=['representative', 'member'])

    print(f"  Loaded {len(clusters):,} cluster assignments")
    print(f"  Unique representatives: {clusters['representative'].nunique():,}")

    # Create gene -> representative mapping
    gene_to_rep = dict(zip(clusters['member'], clusters['representative']))

    return gene_to_rep


def plot_density(umap_coords, output_dir):
    """Create density plot of UMAP embedding."""
    print("\nCreating density plot...")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Create hexbin plot for density
    hexbin = ax.hexbin(
        umap_coords[:, 0],
        umap_coords[:, 1],
        gridsize=100,
        cmap='viridis',
        mincnt=1,
        bins='log'
    )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'UMAP Embedding Density (n={len(umap_coords):,} proteins)', fontsize=14)

    # Add colorbar
    cbar = plt.colorbar(hexbin, ax=ax)
    cbar.set_label('log10(count)', fontsize=12)

    plt.tight_layout()
    output_file = output_dir / 'umap_density.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_sample_by_genome(umap_coords, genome_ids, output_dir, sample_size, seed):
    """Plot a random sample colored by genome."""
    print(f"\nCreating genome-colored plot (sample={sample_size:,})...")

    # Random sample
    np.random.seed(seed)
    n_total = len(umap_coords)
    if sample_size < n_total:
        indices = np.random.choice(n_total, sample_size, replace=False)
        coords_sample = umap_coords[indices]
        genomes_sample = genome_ids[indices]
    else:
        coords_sample = umap_coords
        genomes_sample = genome_ids

    # Get unique genomes and assign colors
    unique_genomes = np.unique(genomes_sample)
    n_genomes = len(unique_genomes)

    print(f"  Sampled {len(coords_sample):,} points from {n_genomes:,} genomes")

    # Create genome -> color mapping
    genome_to_id = {genome: idx for idx, genome in enumerate(unique_genomes)}
    genome_colors = [genome_to_id[g] for g in genomes_sample]

    fig, ax = plt.subplots(figsize=(14, 10))

    scatter = ax.scatter(
        coords_sample[:, 0],
        coords_sample[:, 1],
        c=genome_colors,
        cmap='tab20',
        s=1,
        alpha=0.5,
        rasterized=True
    )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'UMAP colored by Genome ({n_genomes:,} genomes, {len(coords_sample):,} proteins)', fontsize=14)

    plt.tight_layout()
    output_file = output_dir / 'umap_by_genome.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_sample_by_cluster(umap_coords, gene_ids, gene_to_rep, output_dir, sample_size, seed):
    """Plot a random sample colored by MMseqs2 cluster."""
    print(f"\nCreating cluster-colored plot (sample={sample_size:,})...")

    # Assign cluster representatives
    print("  Assigning cluster representatives...")
    cluster_reps = []
    for gene in tqdm(gene_ids, desc="  Mapping"):
        rep = gene_to_rep.get(gene, gene)  # If not in cluster file, gene is its own rep
        cluster_reps.append(rep)
    cluster_reps = np.array(cluster_reps)

    # Random sample
    np.random.seed(seed)
    n_total = len(umap_coords)
    if sample_size < n_total:
        indices = np.random.choice(n_total, sample_size, replace=False)
        coords_sample = umap_coords[indices]
        reps_sample = cluster_reps[indices]
    else:
        coords_sample = umap_coords
        reps_sample = cluster_reps

    # Get unique clusters
    unique_reps = np.unique(reps_sample)
    n_clusters = len(unique_reps)

    print(f"  Sampled {len(coords_sample):,} points from {n_clusters:,} clusters")

    # Create cluster -> color mapping
    rep_to_id = {rep: idx for idx, rep in enumerate(unique_reps)}
    cluster_colors = [rep_to_id[r] for r in reps_sample]

    fig, ax = plt.subplots(figsize=(14, 10))

    scatter = ax.scatter(
        coords_sample[:, 0],
        coords_sample[:, 1],
        c=cluster_colors,
        cmap='tab20',
        s=1,
        alpha=0.5,
        rasterized=True
    )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'UMAP colored by MMseqs2 Cluster (70% ID, {n_clusters:,} clusters, {len(coords_sample):,} proteins)', fontsize=14)

    plt.tight_layout()
    output_file = output_dir / 'umap_by_cluster.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def plot_cluster_stats(gene_ids, gene_to_rep, output_dir):
    """Plot cluster size distribution."""
    print("\nCreating cluster statistics plot...")

    # Assign cluster representatives
    print("  Counting cluster sizes...")
    cluster_reps = [gene_to_rep.get(gene, gene) for gene in gene_ids]

    # Count cluster sizes
    from collections import Counter
    cluster_sizes = Counter(cluster_reps)
    sizes = list(cluster_sizes.values())

    print(f"  Total clusters: {len(cluster_sizes):,}")
    print(f"  Mean cluster size: {np.mean(sizes):.1f}")
    print(f"  Median cluster size: {np.median(sizes):.1f}")
    print(f"  Max cluster size: {np.max(sizes):,}")

    # Create histogram
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Linear scale
    ax1.hist(sizes, bins=100, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Cluster Size', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Cluster Size Distribution (n={len(cluster_sizes):,} clusters)', fontsize=14)
    ax1.grid(True, alpha=0.3)

    # Log scale
    ax2.hist(sizes, bins=100, edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Cluster Size', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Cluster Size Distribution (log scale)', fontsize=14)
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / 'cluster_size_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved to {output_file}")
    plt.close()


def main():
    args = parse_args()

    print("=" * 80)
    print("UMAP VISUALIZATION")
    print("=" * 80)
    print(f"UMAP file: {args.umap_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Sample size: {args.sample_size:,}")
    print("=" * 80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load UMAP embeddings
    umap_coords, gene_ids, genome_ids = load_umap_embeddings(args.umap_file)

    # Load MMseqs2 clusters
    gene_to_rep = load_mmseqs_clusters(args.mmseqs_clusters)

    # Create visualizations
    plot_density(umap_coords, output_dir)
    plot_sample_by_genome(umap_coords, genome_ids, output_dir, args.sample_size, args.seed)
    plot_sample_by_cluster(umap_coords, gene_ids, gene_to_rep, output_dir, args.sample_size, args.seed)
    plot_cluster_stats(gene_ids, gene_to_rep, output_dir)

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print(f"All plots saved to: {output_dir}")
    print("=" * 80)


if __name__ == '__main__':
    main()
