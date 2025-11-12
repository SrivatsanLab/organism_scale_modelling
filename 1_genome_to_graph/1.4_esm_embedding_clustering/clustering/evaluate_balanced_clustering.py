#!/usr/bin/env python
"""
Evaluate Leiden clustering stability on balanced sample from sequence clusters.

This tests if balanced sampling (based on sequence similarity) improves
clustering stability compared to random sampling.

Usage:
    python evaluate_balanced_clustering.py \
        --gene-ids data/balanced_sample_gene_ids.txt \
        --resolution 750 \
        --n-neighbors 15 \
        --n-subsamples 10 \
        --subsample-size 4000 \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/balanced_stability.npz
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import adjusted_rand_score
import igraph as ig
import leidenalg
from tqdm import tqdm


def load_balanced_sample(gene_ids_file, pca_cache_file):
    """Load PCA embeddings for balanced sample genes."""
    print(f"Loading balanced sample from {gene_ids_file}")

    # Load gene IDs from balanced sample
    with open(gene_ids_file, 'r') as f:
        sample_gene_ids = set(line.strip() for line in f)

    print(f"  Balanced sample: {len(sample_gene_ids):,} genes")

    # Load PCA cache
    print(f"Loading PCA cache from {pca_cache_file}")
    data = np.load(pca_cache_file, allow_pickle=True)

    all_gene_ids = data['gene_ids']
    all_embeddings = data['embeddings_pca']

    # Filter to balanced sample
    mask = np.isin(all_gene_ids, list(sample_gene_ids))

    sample_embeddings = all_embeddings[mask]
    sample_genes = all_gene_ids[mask]

    print(f"  Loaded {len(sample_genes):,} genes with embeddings")

    return sample_embeddings, sample_genes


def build_knn_graph(embeddings, n_neighbors=15):
    """Build k-nearest neighbors graph."""
    from sklearn.neighbors import NearestNeighbors

    print(f"Building k-NN graph (k={n_neighbors})...")
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto',
                            metric='euclidean', n_jobs=-1)
    nbrs.fit(embeddings)

    distances, indices = nbrs.kneighbors(embeddings)

    # Build igraph
    edges = []
    weights = []

    for i in range(len(embeddings)):
        for j, neighbor_idx in enumerate(indices[i]):
            if neighbor_idx != i:  # Skip self-loops
                edges.append((i, neighbor_idx))
                # Use inverse distance as weight
                dist = distances[i][j]
                weight = 1.0 / (dist + 1e-10)
                weights.append(weight)

    g = ig.Graph(n=len(embeddings), edges=edges, directed=False)
    g.es['weight'] = weights

    # Simplify (remove duplicate edges)
    g.simplify(combine_edges='max')

    print(f"  Graph: {g.vcount()} nodes, {g.ecount()} edges")

    return g


def run_leiden_clustering(graph, resolution=1.0, seed=42):
    """Run Leiden clustering on graph."""
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed,
        n_iterations=-1
    )

    labels = np.array(partition.membership)
    return labels


def evaluate_stability(embeddings, gene_ids, resolution, n_neighbors,
                      n_subsamples=10, subsample_size=4000, seed=42):
    """
    Evaluate clustering stability via subsampling.

    Creates multiple random subsamples, clusters each, and measures:
    1. Pairwise ARI between subsamples
    2. Per-gene stability (fraction assigned to same cluster)
    3. Cluster count stability
    """
    np.random.seed(seed)

    n_genes = len(gene_ids)
    print(f"\nEvaluating stability:")
    print(f"  Dataset: {n_genes:,} genes (balanced sample)")
    print(f"  Resolution: {resolution}")
    print(f"  n_neighbors: {n_neighbors}")
    print(f"  Subsamples: {n_subsamples}")
    print(f"  Subsample size: {subsample_size:,}")
    print()

    # Create subsamples
    all_labels = []
    all_indices = []
    n_clusters_list = []

    for i in tqdm(range(n_subsamples), desc="Processing subsamples"):
        # Random subsample
        subsample_idx = np.random.choice(n_genes, size=subsample_size, replace=False)
        subsample_idx = np.sort(subsample_idx)  # Keep order for consistency

        sub_embeddings = embeddings[subsample_idx]

        # Build k-NN graph
        graph = build_knn_graph(sub_embeddings, n_neighbors=n_neighbors)

        # Run Leiden clustering
        labels = run_leiden_clustering(graph, resolution=resolution, seed=seed+i)

        all_labels.append(labels)
        all_indices.append(subsample_idx)
        n_clusters_list.append(len(np.unique(labels)))

    print(f"\nComputing stability metrics...")

    # 1. Pairwise ARI between subsamples
    ari_scores = []
    for i in range(n_subsamples):
        for j in range(i+1, n_subsamples):
            # Find genes present in both subsamples
            indices_i = set(all_indices[i])
            indices_j = set(all_indices[j])
            common_indices = indices_i & indices_j

            if len(common_indices) < 100:
                continue  # Skip if too few common genes

            # Get labels for common genes
            common_list = sorted(common_indices)

            # Map back to subsample positions
            labels_i = []
            labels_j = []
            for gene_idx in common_list:
                pos_i = np.where(all_indices[i] == gene_idx)[0][0]
                pos_j = np.where(all_indices[j] == gene_idx)[0][0]
                labels_i.append(all_labels[i][pos_i])
                labels_j.append(all_labels[j][pos_j])

            ari = adjusted_rand_score(labels_i, labels_j)
            ari_scores.append(ari)

    ari_mean = np.mean(ari_scores)
    ari_std = np.std(ari_scores)
    ari_min = np.min(ari_scores)
    ari_max = np.max(ari_scores)

    print(f"  ARI: {ari_mean:.3f} ± {ari_std:.3f} (min={ari_min:.3f}, max={ari_max:.3f})")

    # 2. Per-gene stability
    # Track which genes appear in which subsamples and their cluster assignments
    gene_cluster_assignments = {i: [] for i in range(n_genes)}

    for subsample_idx, labels in zip(all_indices, all_labels):
        for pos, gene_idx in enumerate(subsample_idx):
            gene_cluster_assignments[gene_idx].append(labels[pos])

    # Compute stability for each gene (fraction of times assigned to modal cluster)
    gene_stabilities = []
    genes_tracked = []

    for gene_idx in range(n_genes):
        assignments = gene_cluster_assignments[gene_idx]
        if len(assignments) >= 2:  # Need at least 2 occurrences
            # Find modal cluster
            unique, counts = np.unique(assignments, return_counts=True)
            modal_count = counts.max()
            stability = modal_count / len(assignments)

            gene_stabilities.append(stability)
            genes_tracked.append(gene_idx)

    gene_stabilities = np.array(gene_stabilities)

    gene_stability_mean = gene_stabilities.mean()
    gene_stability_median = np.median(gene_stabilities)
    pct_stable_90 = (gene_stabilities > 0.9).sum() / len(gene_stabilities) * 100
    pct_stable_80 = (gene_stabilities > 0.8).sum() / len(gene_stabilities) * 100
    pct_stable_70 = (gene_stabilities > 0.7).sum() / len(gene_stabilities) * 100

    print(f"  Gene stability: {gene_stability_mean:.3f} (median={gene_stability_median:.3f})")
    print(f"  Genes tracked: {len(genes_tracked):,}")
    print(f"  % genes with stability > 0.9: {pct_stable_90:.1f}%")
    print(f"  % genes with stability > 0.8: {pct_stable_80:.1f}%")
    print(f"  % genes with stability > 0.7: {pct_stable_70:.1f}%")

    # 3. Cluster count stability
    n_clusters_mean = np.mean(n_clusters_list)
    n_clusters_std = np.std(n_clusters_list)
    n_clusters_min = np.min(n_clusters_list)
    n_clusters_max = np.max(n_clusters_list)

    print(f"  N clusters: {n_clusters_mean:.0f} ± {n_clusters_std:.0f} "
          f"(min={n_clusters_min}, max={n_clusters_max})")

    return {
        'resolution': resolution,
        'n_neighbors': n_neighbors,
        'n_genes': n_genes,
        'n_subsamples': n_subsamples,
        'subsample_size': subsample_size,
        'ari_mean': ari_mean,
        'ari_std': ari_std,
        'ari_min': ari_min,
        'ari_max': ari_max,
        'gene_stability_mean': gene_stability_mean,
        'gene_stability_median': gene_stability_median,
        'pct_genes_stable_90': pct_stable_90,
        'pct_genes_stable_80': pct_stable_80,
        'pct_genes_stable_70': pct_stable_70,
        'n_genes_tracked': len(genes_tracked),
        'n_clusters_mean': n_clusters_mean,
        'n_clusters_std': n_clusters_std,
        'n_clusters_min': n_clusters_min,
        'n_clusters_max': n_clusters_max,
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate clustering stability on balanced sample')
    parser.add_argument('--gene-ids', type=str, required=True,
                        help='File with balanced sample gene IDs (one per line)')
    parser.add_argument('--pca-cache', type=str, default='results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz',
                        help='PCA cache file')
    parser.add_argument('--resolution', type=float, default=750,
                        help='Leiden resolution parameter')
    parser.add_argument('--n-neighbors', type=int, default=15,
                        help='Number of nearest neighbors for k-NN graph')
    parser.add_argument('--n-subsamples', type=int, default=10,
                        help='Number of subsamples for stability evaluation')
    parser.add_argument('--subsample-size', type=int, default=4000,
                        help='Size of each subsample')
    parser.add_argument('--output', type=str, required=True,
                        help='Output NPZ file for results')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')

    args = parser.parse_args()

    print("=" * 80)
    print("BALANCED SAMPLE CLUSTERING STABILITY EVALUATION")
    print("=" * 80)
    print()

    # Load balanced sample
    embeddings, gene_ids = load_balanced_sample(args.gene_ids, args.pca_cache)

    # Evaluate stability
    results = evaluate_stability(
        embeddings, gene_ids,
        resolution=args.resolution,
        n_neighbors=args.n_neighbors,
        n_subsamples=args.n_subsamples,
        subsample_size=args.subsample_size,
        seed=args.seed
    )

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez(output_path, **results)

    print(f"\nResults saved to {output_path}")
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Balanced sample: {results['n_genes']:,} genes")
    print(f"ARI: {results['ari_mean']:.3f} ± {results['ari_std']:.3f}")
    print(f"Gene stability: {results['gene_stability_mean']:.3f}")
    print(f"N clusters: {results['n_clusters_mean']:.0f} ± {results['n_clusters_std']:.0f}")
    print()
    print("Compare to random sampling baseline (from previous results):")
    print("  Best ARI: ~0.32 (res=300, nn=15)")
    print("  Typical ARI: 0.13-0.25")
    print()

    if results['ari_mean'] > 0.4:
        print("✅ SUCCESS! Balanced sampling improves stability (ARI > 0.4)")
    elif results['ari_mean'] > 0.32:
        print("⚠️  MARGINAL: Balanced sampling shows modest improvement")
    else:
        print("❌ NO IMPROVEMENT: ARI still below random sampling baseline")

    print("=" * 80)


if __name__ == '__main__':
    main()
