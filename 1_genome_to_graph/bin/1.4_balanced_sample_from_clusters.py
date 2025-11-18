#!/usr/bin/env python
"""
Create a balanced sample from MMseqs2 sequence clusters.

This script samples genes from each sequence cluster to create a balanced
dataset for ESM embedding-based clustering, avoiding sampling bias.

Usage:
    python balanced_sample_from_clusters.py \
        --clusters data/mmseqs_test/cluster_summary_assignments.csv \
        --pca-cache results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
        --output data/balanced_sample_gene_ids.txt \
        --n-samples 100000
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def extract_gene_id_from_fasta_header(header):
    """
    Extract gene ID from header - handles both FASTA headers and clean gene IDs.

    For NCBI format:
        >lcl|NZ_CP191850.1_prot_WP_123456.1_1069 [...] -> NZ_CP191850.1_1069

    For Prodigal/clean format:
        NC_002932.3_1 # 2 # 1126 [...] -> NC_002932.3_1
        NC_002932.3_1 -> NC_002932.3_1
    """
    # Remove leading '>' if present
    header = header.lstrip('>')

    # Split on whitespace and take first part
    parts = header.split()
    if not parts:
        return None

    gene_id = parts[0]

    # Remove 'lcl|' prefix if present (NCBI format)
    if gene_id.startswith('lcl|'):
        gene_id = gene_id[4:]

    # Extract contig_genenum if _prot_ present (NCBI format)
    if '_prot_' in gene_id:
        gene_id = gene_id.split('_prot_')[0]

    # Otherwise, assume it's already in correct format (Prodigal)
    return gene_id


def match_clusters_to_pca_genes(cluster_assignments, pca_gene_ids):
    """
    Match MMseqs2 cluster assignments to genes in PCA cache.

    Args:
        cluster_assignments: DataFrame with columns [representative, member]
        pca_gene_ids: Array of gene IDs from PCA cache

    Returns:
        DataFrame with matched gene IDs and cluster assignments
    """
    print("Matching cluster assignments to PCA cache genes...")

    # Extract gene IDs - member IDs are already clean from MMseqs2 for Prodigal proteins
    # but apply extraction function in case they have FASTA header artifacts
    cluster_assignments['gene_id'] = cluster_assignments['member'].apply(
        extract_gene_id_from_fasta_header
    )

    # Remove entries where gene_id extraction failed
    valid_assignments = cluster_assignments[cluster_assignments['gene_id'].notna()].copy()

    print(f"  Processed gene IDs: {len(valid_assignments)}/{len(cluster_assignments)}")

    # Match to PCA cache
    pca_gene_set = set(pca_gene_ids)
    valid_assignments['in_pca'] = valid_assignments['gene_id'].isin(pca_gene_set)

    matched = valid_assignments[valid_assignments['in_pca']]

    print(f"  Matched {len(matched)} genes to PCA cache ({len(matched)/len(pca_gene_ids)*100:.1f}% of PCA genes)")
    print(f"  Matched genes span {matched['representative'].nunique()} sequence clusters")

    return matched


def balanced_sample(matched_assignments, n_samples, strategy='proportional', min_per_cluster=1):
    """
    Create a balanced sample from sequence clusters.

    Args:
        matched_assignments: DataFrame with columns [representative, gene_id, in_pca]
        n_samples: Total number of genes to sample
        strategy: Sampling strategy:
            - 'proportional': Sample proportional to cluster size
            - 'equal': Sample equal number from each cluster
            - 'sqrt': Sample proportional to sqrt of cluster size
        min_per_cluster: Minimum samples per cluster

    Returns:
        Array of sampled gene IDs
    """
    print(f"\nCreating balanced sample using '{strategy}' strategy...")
    print(f"  Target samples: {n_samples:,}")
    print(f"  Min per cluster: {min_per_cluster}")

    # Count cluster sizes
    cluster_sizes = matched_assignments.groupby('representative').size()
    n_clusters = len(cluster_sizes)

    print(f"  Number of clusters: {n_clusters:,}")

    # Determine samples per cluster
    if strategy == 'equal':
        # Equal samples from each cluster
        samples_per_cluster = pd.Series(
            n_samples // n_clusters,
            index=cluster_sizes.index
        )

    elif strategy == 'proportional':
        # Proportional to cluster size
        proportions = cluster_sizes / cluster_sizes.sum()
        samples_per_cluster = (proportions * n_samples).round().astype(int)

    elif strategy == 'sqrt':
        # Proportional to sqrt of cluster size
        sqrt_sizes = np.sqrt(cluster_sizes)
        proportions = sqrt_sizes / sqrt_sizes.sum()
        samples_per_cluster = (proportions * n_samples).round().astype(int)

    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    # Ensure minimum per cluster
    samples_per_cluster = samples_per_cluster.clip(lower=min_per_cluster)

    # Cap at cluster size
    samples_per_cluster = pd.Series({
        cluster: min(n_samples, cluster_sizes[cluster])
        for cluster, n_samples in samples_per_cluster.items()
    })

    print(f"  Samples per cluster - mean: {samples_per_cluster.mean():.1f}, "
          f"median: {samples_per_cluster.median():.1f}, "
          f"max: {samples_per_cluster.max()}")

    # Sample from each cluster
    sampled_genes = []

    for cluster_id, n in samples_per_cluster.items():
        cluster_genes = matched_assignments[
            matched_assignments['representative'] == cluster_id
        ]['gene_id'].values

        # Sample without replacement (or take all if n >= cluster size)
        if n >= len(cluster_genes):
            sampled = cluster_genes
        else:
            sampled = np.random.choice(cluster_genes, size=n, replace=False)

        sampled_genes.extend(sampled)

    sampled_genes = np.array(sampled_genes)

    print(f"\nSampled {len(sampled_genes):,} genes from {n_clusters:,} clusters")
    print(f"  Average per cluster: {len(sampled_genes) / n_clusters:.1f}")

    return sampled_genes


def main():
    parser = argparse.ArgumentParser(description='Create balanced sample from sequence clusters')
    parser.add_argument('--clusters', type=str, required=True,
                        help='MMseqs2 cluster assignments CSV')
    parser.add_argument('--pca-cache', type=str, default='results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz',
                        help='PCA cache with gene IDs')
    parser.add_argument('--output', type=str, required=True,
                        help='Output file with sampled gene IDs')
    parser.add_argument('--n-samples', type=int, default=100000,
                        help='Number of genes to sample (default: 100,000)')
    parser.add_argument('--strategy', type=str, default='sqrt',
                        choices=['equal', 'proportional', 'sqrt'],
                        help='Sampling strategy (default: sqrt)')
    parser.add_argument('--min-per-cluster', type=int, default=1,
                        help='Minimum samples per cluster (default: 1)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')

    args = parser.parse_args()

    # Set random seed
    np.random.seed(args.seed)

    # Load data
    print(f"Loading PCA cache from {args.pca_cache}")
    pca_data = np.load(args.pca_cache, allow_pickle=True)
    pca_gene_ids = pca_data['gene_ids']
    print(f"  PCA cache has {len(pca_gene_ids):,} genes")

    print(f"\nLoading cluster assignments from {args.clusters}")
    cluster_assignments = pd.read_csv(args.clusters)
    print(f"  Loaded {len(cluster_assignments):,} protein sequences")
    print(f"  Spanning {cluster_assignments['representative'].nunique():,} clusters")

    # Match clusters to PCA genes
    matched = match_clusters_to_pca_genes(cluster_assignments, pca_gene_ids)

    if len(matched) == 0:
        print("\nERROR: No genes could be matched between clusters and PCA cache!")
        print("This might be due to:")
        print("  1. Different gene ID formats between FASTA headers and PCA cache")
        print("  2. Non-overlapping genome sets")
        print("  3. Gene ID extraction issues")
        return 1

    # Create balanced sample
    sampled_genes = balanced_sample(
        matched,
        args.n_samples,
        strategy=args.strategy,
        min_per_cluster=args.min_per_cluster
    )

    # Save to file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        for gene_id in sampled_genes:
            f.write(f"{gene_id}\n")

    print(f"\nSaved {len(sampled_genes):,} sampled gene IDs to {output_path}")

    # Print sampling statistics
    matched_sampled = matched[matched['gene_id'].isin(sampled_genes)]
    cluster_counts = matched_sampled['representative'].value_counts()

    print("\nSampling statistics:")
    print(f"  Clusters represented: {len(cluster_counts):,}")
    print(f"  Genes per cluster - mean: {cluster_counts.mean():.1f}, "
          f"median: {cluster_counts.median():.1f}")
    print(f"  Smallest cluster sampled: {cluster_counts.min()} genes")
    print(f"  Largest cluster sampled: {cluster_counts.max()} genes")

    return 0


if __name__ == '__main__':
    exit(main())
