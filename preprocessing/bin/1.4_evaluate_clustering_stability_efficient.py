#!/usr/bin/env python
"""
Memory-efficient clustering stability evaluation.

Instead of tracking all pairwise co-clustering (5 billion pairs for 100K genes),
we use cluster-level metrics:
1. Adjusted Rand Index (ARI) between repeated clusterings
2. Cluster centroid stability across subsamples
3. Per-cluster stability: fraction of genes staying in same functional group

This is tractable and directly measures what we care about: do the same
clusters emerge from different random samples?

Usage:
    python evaluate_clustering_stability_efficient.py \
        --pca results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
        --n-subsamples 10 \
        --subsample-size 100000 \
        --resolution 1500 \
        --n-neighbors 15 \
        --cog-only \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability/stability_res1500_nn15_cogonly.npz
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from collections import defaultdict, Counter
from tqdm import tqdm
import igraph as ig
import leidenalg
import sys


def load_pca_data(pca_file, cog_only=False):
    """Load PCA embeddings, optionally filtering to COG-annotated genes."""
    data = np.load(pca_file, allow_pickle=True)
    embeddings = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    if cog_only:
        print("Filtering to COG-annotated genes...")
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))
        from scripts.embeddings.evaluate_clustering_quality import load_cog_annotations
        cog_lookup = load_cog_annotations(gene_ids, genome_ids)
        annotated_mask = np.array([gid in cog_lookup for gid in gene_ids])

        embeddings = embeddings[annotated_mask]
        gene_ids = gene_ids[annotated_mask]
        genome_ids = genome_ids[annotated_mask]

        print(f"  Filtered to {len(gene_ids):,} COG-annotated genes")

    return embeddings, gene_ids, genome_ids


def cluster_leiden(embeddings, n_neighbors=15, resolution=1.0, seed=None):
    """
    Perform Leiden clustering on embeddings.

    Args:
        embeddings: (N, D) array of embeddings
        n_neighbors: Number of neighbors for kNN graph
        resolution: Leiden resolution parameter
        seed: Random seed for reproducibility

    Returns:
        Array of cluster labels
    """
    from sklearn.neighbors import NearestNeighbors

    # Build kNN graph
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
                weights.append(1.0 / (1.0 + dist))

    # Create igraph
    g = ig.Graph(n=len(embeddings), edges=edges, directed=False)
    g.es['weight'] = weights

    # Leiden clustering
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        n_iterations=-1,
        seed=seed if seed is not None else 42
    )

    labels = np.array(partition.membership)
    return labels


def generate_subsamples(n_total, n_subsamples, subsample_size, seed=42):
    """Generate multiple independent random subsamples."""
    print(f"Generating {n_subsamples} independent subsamples of size {subsample_size:,}...")

    rng = np.random.RandomState(seed)
    subsamples = []

    for i in range(n_subsamples):
        indices = rng.choice(n_total, size=subsample_size, replace=False)
        subsamples.append(indices)

    return subsamples


def compute_ari_stability(labels_list, gene_ids_list):
    """
    Compute pairwise ARI between all clusterings.

    High ARI = clusterings are similar = stable.
    Only compares genes that appear in both subsamples.

    Returns:
        List of pairwise ARI scores
    """
    print("Computing pairwise ARI between clusterings...")

    ari_scores = []
    n = len(labels_list)

    for i in tqdm(range(n), desc="Computing ARIs"):
        for j in range(i + 1, n):
            # Find common genes
            genes_i = set(gene_ids_list[i])
            genes_j = set(gene_ids_list[j])
            common_genes = genes_i & genes_j

            if len(common_genes) < 100:  # Need reasonable overlap
                continue

            # Get labels for common genes
            gene_to_label_i = dict(zip(gene_ids_list[i], labels_list[i]))
            gene_to_label_j = dict(zip(gene_ids_list[j], labels_list[j]))

            labels_i = [gene_to_label_i[g] for g in common_genes]
            labels_j = [gene_to_label_j[g] for g in common_genes]

            # Compute ARI
            ari = adjusted_rand_score(labels_i, labels_j)
            ari_scores.append(ari)

    return ari_scores


def compute_cluster_membership_stability(labels_list, gene_ids_list):
    """
    For each gene that appears in multiple subsamples, compute:
    - How many different clusters it was assigned to
    - Stability = 1 / (number of unique clusters)

    Returns:
        Dict mapping gene_id -> stability score
    """
    print("Computing per-gene cluster membership stability...")

    gene_clusters = defaultdict(list)

    # Collect cluster assignments for each gene
    for labels, gene_ids in zip(labels_list, gene_ids_list):
        for gene_id, cluster in zip(gene_ids, labels):
            gene_clusters[gene_id].append(cluster)

    # Compute stability
    gene_stability = {}
    for gene_id, clusters in gene_clusters.items():
        n_appearances = len(clusters)
        n_unique_clusters = len(set(clusters))
        # Stability = 1 if always same cluster, 1/n if always different
        stability = 1.0 / n_unique_clusters
        gene_stability[gene_id] = {
            'stability': stability,
            'n_appearances': n_appearances,
            'n_unique_clusters': n_unique_clusters
        }

    return gene_stability


def compute_cluster_size_stability(labels_list):
    """
    Check if the number and sizes of clusters are stable across subsamples.

    Returns:
        Dict with cluster count and size statistics
    """
    print("Computing cluster size stability...")

    n_clusters_list = []
    mean_sizes = []
    median_sizes = []

    for labels in labels_list:
        cluster_counts = Counter(labels[labels >= 0])
        sizes = list(cluster_counts.values())

        n_clusters_list.append(len(sizes))
        mean_sizes.append(np.mean(sizes))
        median_sizes.append(np.median(sizes))

    return {
        'n_clusters_mean': np.mean(n_clusters_list),
        'n_clusters_std': np.std(n_clusters_list),
        'n_clusters_min': np.min(n_clusters_list),
        'n_clusters_max': np.max(n_clusters_list),
        'mean_size_mean': np.mean(mean_sizes),
        'mean_size_std': np.std(mean_sizes),
        'median_size_mean': np.mean(median_sizes),
        'median_size_std': np.std(median_sizes)
    }


def evaluate_stability(pca_file, n_subsamples, subsample_size,
                       resolution, n_neighbors, cog_only=False):
    """
    Efficient stability evaluation using cluster-level metrics.

    Returns:
        dict with stability metrics
    """
    print(f"\n{'='*80}")
    print("Clustering Stability Evaluation (Efficient)")
    print(f"{'='*80}\n")
    print(f"Parameters:")
    print(f"  Resolution: {resolution}")
    print(f"  N neighbors: {n_neighbors}")
    print(f"  COG-only: {cog_only}")
    print(f"  N subsamples: {n_subsamples}")
    print(f"  Subsample size: {subsample_size:,}")
    print()
    sys.stdout.flush()

    # Load data
    print("Loading PCA data...")
    sys.stdout.flush()
    embeddings, gene_ids, genome_ids = load_pca_data(pca_file, cog_only)
    n_total = len(gene_ids)
    print(f"Total genes: {n_total:,}\n")
    sys.stdout.flush()

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

        print(f"    Building kNN graph (k={n_neighbors})...")
        print(f"    Creating graph and running Leiden (resolution={resolution})...")
        labels = cluster_leiden(subsample_embeddings, n_neighbors, resolution, seed=42+i)
        print(f"    Found {len(set(labels)):,} clusters")

        labels_list.append(labels)
        gene_ids_list.append(subsample_gene_ids)
        sys.stdout.flush()

    # Compute stability metrics
    print("\n" + "="*80)
    print("Computing stability metrics...")
    print("="*80 + "\n")

    # 1. Pairwise ARI between clusterings
    ari_scores = compute_ari_stability(labels_list, gene_ids_list)

    # 2. Per-gene cluster membership stability
    gene_stability = compute_cluster_membership_stability(labels_list, gene_ids_list)

    # 3. Cluster size stability
    size_stability = compute_cluster_size_stability(labels_list)

    # Aggregate results
    print("\n" + "="*80)
    print("Summary Statistics")
    print("="*80)

    print(f"\n1. Pairwise ARI (clustering agreement):")
    print(f"   Mean ARI: {np.mean(ari_scores):.4f}")
    print(f"   Std ARI: {np.std(ari_scores):.4f}")
    print(f"   Min ARI: {np.min(ari_scores):.4f}")
    print(f"   Max ARI: {np.max(ari_scores):.4f}")

    stabilities = [stats['stability'] for stats in gene_stability.values()]
    appearances = [stats['n_appearances'] for stats in gene_stability.values()]

    print(f"\n2. Per-gene membership stability:")
    print(f"   Genes tracked: {len(gene_stability):,}")
    print(f"   Mean stability: {np.mean(stabilities):.4f}")
    print(f"   Median stability: {np.median(stabilities):.4f}")
    print(f"   Std stability: {np.std(stabilities):.4f}")
    print(f"   Genes with stability >= 0.9: {np.sum(np.array(stabilities) >= 0.9):,} "
          f"({100*np.sum(np.array(stabilities) >= 0.9)/len(stabilities):.1f}%)")

    print(f"\n3. Cluster count/size stability:")
    print(f"   Mean # clusters: {size_stability['n_clusters_mean']:.1f} ± {size_stability['n_clusters_std']:.1f}")
    print(f"   Range # clusters: {size_stability['n_clusters_min']} - {size_stability['n_clusters_max']}")
    print(f"   Mean cluster size: {size_stability['mean_size_mean']:.1f} ± {size_stability['mean_size_std']:.1f}")

    # Package results
    results = {
        'resolution': resolution,
        'n_neighbors': n_neighbors,
        'cog_only': cog_only,
        'n_subsamples': n_subsamples,
        'subsample_size': subsample_size,
        'n_total_genes': n_total,
        # ARI metrics
        'ari_mean': float(np.mean(ari_scores)),
        'ari_std': float(np.std(ari_scores)),
        'ari_min': float(np.min(ari_scores)),
        'ari_max': float(np.max(ari_scores)),
        'ari_median': float(np.median(ari_scores)),
        # Gene stability metrics
        'gene_stability_mean': float(np.mean(stabilities)),
        'gene_stability_median': float(np.median(stabilities)),
        'gene_stability_std': float(np.std(stabilities)),
        'pct_genes_stable_90': float(100*np.sum(np.array(stabilities) >= 0.9)/len(stabilities)),
        'pct_genes_stable_80': float(100*np.sum(np.array(stabilities) >= 0.8)/len(stabilities)),
        'pct_genes_stable_70': float(100*np.sum(np.array(stabilities) >= 0.7)/len(stabilities)),
        'n_genes_tracked': len(gene_stability),
        # Cluster size stability
        **size_stability
    }

    # Also save per-gene stability for detailed analysis
    gene_stability_data = {
        'gene_ids': np.array(list(gene_stability.keys())),
        'stabilities': np.array([gene_stability[g]['stability'] for g in gene_stability.keys()]),
        'n_appearances': np.array([gene_stability[g]['n_appearances'] for g in gene_stability.keys()]),
        'n_unique_clusters': np.array([gene_stability[g]['n_unique_clusters'] for g in gene_stability.keys()])
    }

    return results, gene_stability_data, ari_scores


def main():
    parser = argparse.ArgumentParser(
        description='Efficient clustering stability evaluation'
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
    parser.add_argument('--output', type=str, required=True,
                        help='Output .npz file for results')

    args = parser.parse_args()

    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Run evaluation
    results, gene_stability_data, ari_scores = evaluate_stability(
        pca_file=args.pca,
        n_subsamples=args.n_subsamples,
        subsample_size=args.subsample_size,
        resolution=args.resolution,
        n_neighbors=args.n_neighbors,
        cog_only=args.cog_only
    )

    # Save results
    print(f"\n{'='*80}")
    print(f"Saving results to {output_file}")
    print(f"{'='*80}\n")

    np.savez_compressed(
        output_file,
        **results,
        **gene_stability_data,
        ari_scores=np.array(ari_scores)
    )

    print("Done!")


if __name__ == '__main__':
    main()
