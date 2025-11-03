#!/usr/bin/env python
"""
Evaluate clustering stability across multiple independent subsamples.

Key idea: If clustering is robust and meaningful, pairs of genes that co-occur
in multiple subsamples should consistently be assigned to the same cluster.

Approach:
1. Generate N independent random subsamples of size M from full dataset
2. Cluster each subsample with the same parameters
3. For all gene pairs that co-occur in multiple subsamples:
   - Compute co-clustering rate (fraction of times they're in same cluster)
4. Aggregate statistics:
   - Mean co-clustering rate (higher = more stable)
   - Distribution of co-clustering rates
   - Comparison across different parameters

Usage:
    python evaluate_clustering_stability.py \
        --pca results/umap/pca_cache.npz \
        --n-subsamples 10 \
        --subsample-size 100000 \
        --resolution 1500 \
        --n-neighbors 15 \
        --cog-only \
        --output results/clustering/stability_res1500_nn15.npz
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from collections import defaultdict
from tqdm import tqdm
import igraph as ig
import leidenalg


def load_pca_data(pca_file, cog_only=False):
    """Load PCA embeddings, optionally filtering to COG-annotated genes."""
    data = np.load(pca_file, allow_pickle=True)
    embeddings = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    if cog_only:
        print("Filtering to COG-annotated genes...")
        from evaluate_clustering_quality import load_cog_annotations
        cog_lookup = load_cog_annotations(gene_ids, genome_ids)
        annotated_mask = np.array([gid in cog_lookup for gid in gene_ids])

        embeddings = embeddings[annotated_mask]
        gene_ids = gene_ids[annotated_mask]
        genome_ids = genome_ids[annotated_mask]

        print(f"  Filtered to {len(gene_ids):,} COG-annotated genes")

    return embeddings, gene_ids, genome_ids


def cluster_leiden(embeddings, n_neighbors=15, resolution=1.0):
    """
    Perform Leiden clustering on embeddings.

    Args:
        embeddings: (N, D) array of embeddings
        n_neighbors: Number of neighbors for kNN graph
        resolution: Leiden resolution parameter

    Returns:
        Array of cluster labels
    """
    from sklearn.neighbors import NearestNeighbors

    # Build kNN graph
    print(f"    Building kNN graph (k={n_neighbors})...")
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric='cosine', n_jobs=-1)
    nbrs.fit(embeddings)
    distances, indices = nbrs.kneighbors(embeddings)

    # Create edges
    edges = []
    weights = []
    for i in range(len(embeddings)):
        for j, dist in zip(indices[i][1:], distances[i][1:]):  # Skip self
            if i < j:  # Avoid duplicates
                edges.append((i, j))
                # Convert distance to similarity (inverse)
                weights.append(1.0 / (1.0 + dist))

    # Create igraph
    print(f"    Creating graph ({len(edges):,} edges)...")
    g = ig.Graph(n=len(embeddings), edges=edges, directed=False)
    g.es['weight'] = weights

    # Leiden clustering
    print(f"    Running Leiden (resolution={resolution})...")
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        n_iterations=-1,
        seed=42
    )

    labels = np.array(partition.membership)
    print(f"    Found {len(set(labels)):,} clusters")

    return labels


def generate_subsamples(n_total, n_subsamples, subsample_size, seed=42):
    """
    Generate multiple independent random subsamples.

    Returns:
        List of index arrays, each of length subsample_size
    """
    print(f"Generating {n_subsamples} independent subsamples of size {subsample_size:,}...")

    rng = np.random.RandomState(seed)
    subsamples = []

    for i in range(n_subsamples):
        # Random sample without replacement
        indices = rng.choice(n_total, size=subsample_size, replace=False)
        subsamples.append(indices)

    return subsamples


def compute_coclustering_matrix(labels_list, gene_ids_list):
    """
    Compute co-clustering statistics across multiple clusterings.

    For each pair of genes that co-occur in multiple subsamples,
    compute the fraction of times they're in the same cluster.

    Args:
        labels_list: List of cluster label arrays
        gene_ids_list: List of gene ID arrays

    Returns:
        dict mapping (gene1, gene2) → (n_cooccur, n_cocluster)
    """
    print("Computing co-clustering statistics...")

    # Track co-occurrences and co-clusterings
    pair_stats = defaultdict(lambda: {'cooccur': 0, 'cocluster': 0})

    for labels, gene_ids in tqdm(zip(labels_list, gene_ids_list),
                                  total=len(labels_list),
                                  desc="Processing subsamples"):
        # Build gene → cluster mapping
        gene_to_cluster = dict(zip(gene_ids, labels))

        # For each pair in this subsample
        genes = list(gene_ids)
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                g1, g2 = genes[i], genes[j]
                # Canonicalize pair order
                pair = tuple(sorted([g1, g2]))

                pair_stats[pair]['cooccur'] += 1

                # Check if same cluster
                if gene_to_cluster[g1] == gene_to_cluster[g2]:
                    pair_stats[pair]['cocluster'] += 1

    return pair_stats


def analyze_stability(pair_stats, min_cooccur=2):
    """
    Analyze co-clustering stability.

    Args:
        pair_stats: Dict from compute_coclustering_matrix
        min_cooccur: Minimum co-occurrences to include pair

    Returns:
        dict with stability metrics
    """
    print(f"Analyzing stability (min_cooccur={min_cooccur})...")

    # Filter to pairs with sufficient co-occurrences
    filtered_pairs = {
        pair: stats for pair, stats in pair_stats.items()
        if stats['cooccur'] >= min_cooccur
    }

    if len(filtered_pairs) == 0:
        print("  Warning: No pairs with sufficient co-occurrences!")
        return {
            'n_pairs': 0,
            'mean_coclustering_rate': 0.0,
            'median_coclustering_rate': 0.0,
            'std_coclustering_rate': 0.0
        }

    # Compute co-clustering rates
    coclustering_rates = []
    for stats in filtered_pairs.values():
        rate = stats['cocluster'] / stats['cooccur']
        coclustering_rates.append(rate)

    coclustering_rates = np.array(coclustering_rates)

    # Compute statistics
    mean_rate = np.mean(coclustering_rates)
    median_rate = np.median(coclustering_rates)
    std_rate = np.std(coclustering_rates)

    # Percentiles
    percentiles = {}
    for p in [10, 25, 50, 75, 90, 95, 99]:
        percentiles[f'p{p}'] = np.percentile(coclustering_rates, p)

    print(f"  Pairs analyzed: {len(filtered_pairs):,}")
    print(f"  Mean co-clustering rate: {mean_rate:.4f}")
    print(f"  Median co-clustering rate: {median_rate:.4f}")
    print(f"  Std co-clustering rate: {std_rate:.4f}")

    return {
        'n_pairs': len(filtered_pairs),
        'mean_coclustering_rate': float(mean_rate),
        'median_coclustering_rate': float(median_rate),
        'std_coclustering_rate': float(std_rate),
        **{k: float(v) for k, v in percentiles.items()}
    }


def evaluate_stability(pca_file, n_subsamples, subsample_size,
                       resolution, n_neighbors, cog_only=False,
                       min_cooccur=2):
    """
    Full stability evaluation pipeline.

    Returns:
        dict with stability metrics and detailed pair statistics
    """
    print(f"\n{'='*80}")
    print("Clustering Stability Evaluation")
    print(f"{'='*80}\n")
    print(f"Parameters:")
    print(f"  Resolution: {resolution}")
    print(f"  N neighbors: {n_neighbors}")
    print(f"  COG-only: {cog_only}")
    print(f"  N subsamples: {n_subsamples}")
    print(f"  Subsample size: {subsample_size:,}")
    print()

    # Load data
    embeddings, gene_ids, genome_ids = load_pca_data(pca_file, cog_only)
    n_total = len(gene_ids)
    print(f"Total genes: {n_total:,}\n")

    # Generate subsamples
    subsample_indices = generate_subsamples(
        n_total, n_subsamples, subsample_size
    )

    # Cluster each subsample
    print("\nClustering subsamples...")
    labels_list = []
    gene_ids_list = []

    for i, indices in enumerate(subsample_indices):
        print(f"\nSubsample {i+1}/{n_subsamples}:")
        subsample_embeddings = embeddings[indices]
        subsample_gene_ids = gene_ids[indices]

        labels = cluster_leiden(subsample_embeddings, n_neighbors, resolution)

        labels_list.append(labels)
        gene_ids_list.append(subsample_gene_ids)

    # Compute co-clustering statistics
    print("\n")
    pair_stats = compute_coclustering_matrix(labels_list, gene_ids_list)

    # Analyze stability
    print("\n")
    stability_metrics = analyze_stability(pair_stats, min_cooccur)

    # Add metadata
    results = {
        'resolution': resolution,
        'n_neighbors': n_neighbors,
        'cog_only': cog_only,
        'n_subsamples': n_subsamples,
        'subsample_size': subsample_size,
        'n_total_genes': n_total,
        'min_cooccur': min_cooccur,
        **stability_metrics
    }

    return results, pair_stats


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate clustering stability across subsamples'
    )
    parser.add_argument('--pca', type=str, required=True,
                        help='Path to PCA embeddings .npz file')
    parser.add_argument('--n-subsamples', type=int, default=10,
                        help='Number of independent subsamples')
    parser.add_argument('--subsample-size', type=int, default=100000,
                        help='Size of each subsample')
    parser.add_argument('--resolution', type=float, default=1500.0,
                        help='Leiden resolution parameter')
    parser.add_argument('--n-neighbors', type=int, default=15,
                        help='Number of neighbors for kNN graph')
    parser.add_argument('--cog-only', action='store_true',
                        help='Use only COG-annotated genes')
    parser.add_argument('--min-cooccur', type=int, default=2,
                        help='Minimum co-occurrences to analyze pair')
    parser.add_argument('--output', type=str, required=True,
                        help='Output .npz file for results')

    args = parser.parse_args()

    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Run evaluation
    results, pair_stats = evaluate_stability(
        pca_file=args.pca,
        n_subsamples=args.n_subsamples,
        subsample_size=args.subsample_size,
        resolution=args.resolution,
        n_neighbors=args.n_neighbors,
        cog_only=args.cog_only,
        min_cooccur=args.min_cooccur
    )

    # Save results
    print(f"\n{'='*80}")
    print(f"Saving results to {output_file}")
    print(f"{'='*80}\n")

    # Convert pair_stats to saveable format
    pair_data = {
        'gene_pairs': np.array([list(pair) for pair in pair_stats.keys()]),
        'cooccur_counts': np.array([stats['cooccur'] for stats in pair_stats.values()]),
        'cocluster_counts': np.array([stats['cocluster'] for stats in pair_stats.values()])
    }

    np.savez_compressed(
        output_file,
        **results,
        **pair_data
    )

    # Print summary
    print("Summary metrics:")
    for key, value in results.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")


if __name__ == '__main__':
    main()
