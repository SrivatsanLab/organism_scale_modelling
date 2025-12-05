#!/usr/bin/env python
"""
Analyze the tightness of MMseqs2 clusters in ESM embedding space.

For each MMseqs2 cluster, compute:
- Mean embedding across all cluster members
- Standard deviation for each dimension
- Overall cluster variance/tightness metrics

This helps determine if MMseqs2 clusters are already tight enough
in embedding space that we don't need additional Leiden clustering.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import argparse
from collections import defaultdict


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze MMseqs2 cluster tightness in embedding space"
    )
    parser.add_argument(
        "--pca-cache",
        type=str,
        default="results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz",
        help="Path to PCA cache with embeddings",
    )
    parser.add_argument(
        "--mmseqs-clusters",
        type=str,
        default="data/mmseqs_prodigal_test/clusters.tsv",
        help="Path to MMseqs2 cluster assignments",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis",
        help="Output directory",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=2,
        help="Minimum cluster size to analyze (default: 2)",
    )
    return parser.parse_args()


def load_pca_embeddings(pca_file):
    """Load PCA embeddings."""
    print(f"Loading PCA embeddings from {pca_file}...")
    data = np.load(pca_file, allow_pickle=True)

    embeddings_pca = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    # Reconstruct full gene IDs by prepending the genome ID prefix
    # The PCA cache splits "GCF_000006985.1_NC_002932.3_2" into:
    #   genome_ids: "GCF"
    #   gene_ids: "000006985.1_NC_002932.3_2"
    # We need to reconstruct: "GCF_000006985.1_NC_002932.3_2"
    full_gene_ids = np.array([f"{genome}_{gene}" for genome, gene in zip(genome_ids, gene_ids)])

    print(f"  Loaded {len(full_gene_ids):,} proteins")
    print(f"  Embedding dimensions: {embeddings_pca.shape[1]}")
    print(f"  Sample gene IDs: {full_gene_ids[:3]}")

    return embeddings_pca, full_gene_ids


def load_mmseqs_clusters(cluster_file):
    """Load MMseqs2 cluster assignments."""
    print(f"\nLoading MMseqs2 clusters from {cluster_file}...")

    # Read cluster file (representative \t member)
    clusters = pd.read_csv(cluster_file, sep='\t', names=['representative', 'member'])

    print(f"  Loaded {len(clusters):,} assignments")
    print(f"  Unique clusters: {clusters['representative'].nunique():,}")

    return clusters


def compute_cluster_statistics(embeddings, gene_ids, clusters_df, min_cluster_size=2):
    """
    Compute mean and std for each cluster.

    Returns:
        cluster_stats: DataFrame with cluster statistics
        cluster_embeddings: Dict mapping cluster rep -> mean embedding
    """
    print(f"\nComputing cluster statistics (min size={min_cluster_size})...")

    # Create gene_id -> embedding index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_ids)}

    # Group by cluster representative
    cluster_groups = clusters_df.groupby('representative')

    cluster_stats_list = []
    cluster_mean_embeddings = {}

    n_dims = embeddings.shape[1]

    for rep, group in tqdm(cluster_groups, desc="Processing clusters"):
        cluster_size = len(group)

        if cluster_size < min_cluster_size:
            continue

        # Get embeddings for this cluster
        member_indices = []
        for member in group['member']:
            if member in gene_to_idx:
                member_indices.append(gene_to_idx[member])

        if len(member_indices) == 0:
            continue

        cluster_embeddings = embeddings[member_indices]

        # Compute statistics
        mean_emb = np.mean(cluster_embeddings, axis=0)
        std_emb = np.std(cluster_embeddings, axis=0)

        # Overall metrics
        mean_std = np.mean(std_emb)  # Average std across dimensions
        max_std = np.max(std_emb)    # Max std across dimensions

        # Total variance
        total_variance = np.sum(np.var(cluster_embeddings, axis=0))

        # Average pairwise distance
        from scipy.spatial.distance import pdist
        if len(member_indices) > 1:
            pairwise_dists = pdist(cluster_embeddings, metric='euclidean')
            mean_pairwise_dist = np.mean(pairwise_dists)
            max_pairwise_dist = np.max(pairwise_dists)
        else:
            mean_pairwise_dist = 0.0
            max_pairwise_dist = 0.0

        # Store stats
        cluster_stats_list.append({
            'cluster_representative': rep,
            'cluster_size': cluster_size,
            'mean_std': mean_std,
            'max_std': max_std,
            'total_variance': total_variance,
            'mean_pairwise_distance': mean_pairwise_dist,
            'max_pairwise_distance': max_pairwise_dist,
        })

        cluster_mean_embeddings[rep] = mean_emb

    cluster_stats = pd.DataFrame(cluster_stats_list)

    print(f"  Processed {len(cluster_stats):,} clusters")
    print(f"  Size range: {cluster_stats['cluster_size'].min()}-{cluster_stats['cluster_size'].max()}")

    return cluster_stats, cluster_mean_embeddings


def save_cluster_statistics(cluster_stats, output_dir):
    """Save cluster statistics to CSV."""
    output_file = output_dir / 'mmseqs_cluster_statistics.csv'
    print(f"\nSaving cluster statistics to {output_file}...")
    cluster_stats.to_csv(output_file, index=False)
    print(f"  Saved {len(cluster_stats):,} clusters")


def save_detailed_per_dimension_stats(embeddings, gene_ids, clusters_df, output_dir, min_cluster_size=2):
    """
    Save detailed per-dimension statistics for each cluster.
    This creates a CSV with: cluster_rep, dimension, mean, std
    """
    print(f"\nComputing per-dimension statistics...")

    # Create gene_id -> embedding index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_ids)}

    # Group by cluster representative
    cluster_groups = clusters_df.groupby('representative')

    n_dims = embeddings.shape[1]

    results = []

    for rep, group in tqdm(cluster_groups, desc="Computing per-dimension stats"):
        cluster_size = len(group)

        if cluster_size < min_cluster_size:
            continue

        # Get embeddings for this cluster
        member_indices = []
        for member in group['member']:
            if member in gene_to_idx:
                member_indices.append(gene_to_idx[member])

        if len(member_indices) == 0:
            continue

        cluster_embeddings = embeddings[member_indices]

        # Compute mean and std for each dimension
        mean_per_dim = np.mean(cluster_embeddings, axis=0)
        std_per_dim = np.std(cluster_embeddings, axis=0)

        for dim_idx in range(n_dims):
            results.append({
                'cluster_representative': rep,
                'cluster_size': cluster_size,
                'dimension': dim_idx,
                'mean': mean_per_dim[dim_idx],
                'std': std_per_dim[dim_idx]
            })

    df = pd.DataFrame(results)

    output_file = output_dir / 'mmseqs_cluster_per_dimension_stats.csv'
    print(f"  Saving to {output_file}...")
    df.to_csv(output_file, index=False)
    print(f"  Saved {len(df):,} rows ({len(df['cluster_representative'].unique()):,} clusters × {n_dims} dimensions)")


def print_summary(cluster_stats):
    """Print summary statistics."""
    print("\n" + "=" * 80)
    print("CLUSTER TIGHTNESS SUMMARY")
    print("=" * 80)

    print(f"\nTotal clusters analyzed: {len(cluster_stats):,}")

    print(f"\nCluster sizes:")
    print(f"  Min: {cluster_stats['cluster_size'].min()}")
    print(f"  Max: {cluster_stats['cluster_size'].max():,}")
    print(f"  Mean: {cluster_stats['cluster_size'].mean():.1f}")
    print(f"  Median: {cluster_stats['cluster_size'].median():.1f}")

    print(f"\nMean standard deviation (across dimensions):")
    print(f"  Min: {cluster_stats['mean_std'].min():.4f}")
    print(f"  Max: {cluster_stats['mean_std'].max():.4f}")
    print(f"  Mean: {cluster_stats['mean_std'].mean():.4f}")
    print(f"  Median: {cluster_stats['mean_std'].median():.4f}")

    print(f"\nMax standard deviation (across dimensions):")
    print(f"  Min: {cluster_stats['max_std'].min():.4f}")
    print(f"  Max: {cluster_stats['max_std'].max():.4f}")
    print(f"  Mean: {cluster_stats['max_std'].mean():.4f}")
    print(f"  Median: {cluster_stats['max_std'].median():.4f}")

    print(f"\nMean pairwise distance:")
    print(f"  Min: {cluster_stats['mean_pairwise_distance'].min():.4f}")
    print(f"  Max: {cluster_stats['mean_pairwise_distance'].max():.4f}")
    print(f"  Mean: {cluster_stats['mean_pairwise_distance'].mean():.4f}")
    print(f"  Median: {cluster_stats['mean_pairwise_distance'].median():.4f}")

    print(f"\nMax pairwise distance:")
    print(f"  Min: {cluster_stats['max_pairwise_distance'].min():.4f}")
    print(f"  Max: {cluster_stats['max_pairwise_distance'].max():.4f}")
    print(f"  Mean: {cluster_stats['max_pairwise_distance'].mean():.4f}")
    print(f"  Median: {cluster_stats['max_pairwise_distance'].median():.4f}")

    # Quartiles
    print(f"\nMean std quartiles:")
    for q in [0.25, 0.5, 0.75]:
        val = cluster_stats['mean_std'].quantile(q)
        print(f"  {int(q*100)}th percentile: {val:.4f}")

    print("=" * 80)


def main():
    args = parse_args()

    print("=" * 80)
    print("MMSEQS2 CLUSTER TIGHTNESS ANALYSIS")
    print("=" * 80)
    print(f"PCA cache: {args.pca_cache}")
    print(f"MMseqs2 clusters: {args.mmseqs_clusters}")
    print(f"Output directory: {args.output_dir}")
    print(f"Min cluster size: {args.min_cluster_size}")
    print("=" * 80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    embeddings, gene_ids = load_pca_embeddings(args.pca_cache)
    clusters_df = load_mmseqs_clusters(args.mmseqs_clusters)

    # Compute cluster statistics
    cluster_stats, cluster_mean_embeddings = compute_cluster_statistics(
        embeddings, gene_ids, clusters_df, args.min_cluster_size
    )

    # Save results
    save_cluster_statistics(cluster_stats, output_dir)
    save_detailed_per_dimension_stats(embeddings, gene_ids, clusters_df, output_dir, args.min_cluster_size)

    # Print summary
    print_summary(cluster_stats)

    print(f"\n✓ Analysis complete!")
    print(f"  Outputs saved to: {output_dir}")
    print(f"    - mmseqs_cluster_statistics.csv")
    print(f"    - mmseqs_cluster_per_dimension_stats.csv")


if __name__ == '__main__':
    main()
