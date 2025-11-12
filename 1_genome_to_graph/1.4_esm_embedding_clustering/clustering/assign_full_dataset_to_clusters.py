#!/usr/bin/env python
"""
Assign all 29M genes to clusters based on subsample clustering.

Strategy:
1. Load best clustering from subsample (e.g., res1500, nn15, cogonly)
2. Compute cluster centroids in PCA space
3. Process all 29M gene embeddings in batches
4. Assign each gene to nearest cluster centroid
5. Compute cluster statistics (mean, SD per dimension)

This scales to full dataset without re-clustering.

Usage:
    python assign_full_dataset_to_clusters.py \
        --clustering results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/clusters_leiden_res1500_nn15_cogonly.npz \
        --embedding-dir data/refseq_esm_embeddings \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/full_dataset_assignments_res1500_nn15.npz \
        --batch-size 100000
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import sys


def load_clustering(clustering_file):
    """Load clustering results from subsample."""
    print(f"Loading clustering from {clustering_file.name}...")
    data = np.load(clustering_file, allow_pickle=True)

    labels = data['labels']
    gene_ids = data['gene_ids']

    # Extract metadata
    resolution = float(data.get('resolution', 0))
    n_neighbors = int(data.get('n_neighbors', 0))
    n_clusters = len(set(labels[labels >= 0]))

    print(f"  Resolution: {resolution}")
    print(f"  N neighbors: {n_neighbors}")
    print(f"  N genes in subsample: {len(gene_ids):,}")
    print(f"  N clusters: {n_clusters:,}")

    return labels, gene_ids, {'resolution': resolution, 'n_neighbors': n_neighbors}


def compute_cluster_centroids(embeddings_pca, labels, gene_ids):
    """
    Compute centroid of each cluster in PCA space.

    Returns:
        cluster_centroids: (n_clusters, pca_dim) array
        cluster_ids: Array of cluster IDs
    """
    print("Computing cluster centroids...")

    unique_clusters = sorted(set(labels[labels >= 0]))
    n_clusters = len(unique_clusters)
    pca_dim = embeddings_pca.shape[1]

    cluster_centroids = np.zeros((n_clusters, pca_dim))
    cluster_ids = np.array(unique_clusters)

    for i, cluster_id in enumerate(tqdm(unique_clusters, desc="Computing centroids")):
        mask = labels == cluster_id
        cluster_centroids[i] = embeddings_pca[mask].mean(axis=0)

    print(f"  Computed {n_clusters:,} centroids in {pca_dim}D space")

    return cluster_centroids, cluster_ids


def fit_pca_model(embedding_dir, n_genes_sample=1000000, n_components=50, seed=42):
    """
    Fit PCA model from a sample of genes.

    Returns:
        pca_components, pca_mean
    """
    print(f"Fitting PCA model from sample of {n_genes_sample:,} genes...")

    embedding_files = sorted(Path(embedding_dir).glob('*_embeddings.npz'))

    # Collect sample
    all_embeddings = []
    all_gene_ids = []

    from tqdm import tqdm
    for emb_file in tqdm(embedding_files, desc="Loading embeddings for PCA"):
        data = np.load(emb_file, allow_pickle=True)
        all_embeddings.append(data['embeddings'])
        all_gene_ids.extend(data['gene_ids'])

        if len(all_gene_ids) >= n_genes_sample:
            break

    all_embeddings = np.vstack(all_embeddings)[:n_genes_sample]
    print(f"  Loaded {len(all_embeddings):,} genes, dim={all_embeddings.shape[1]}")

    # Fit PCA
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components, random_state=seed)
    print(f"  Fitting PCA to {n_components} components...")
    pca.fit(all_embeddings)

    print(f"  Explained variance: {pca.explained_variance_ratio_.sum():.4f}")

    return pca.components_, pca.mean_


def process_embeddings_batch(embedding_files, pca_components, pca_mean,
                            cluster_centroids, cluster_ids, batch_size=100000):
    """
    Process all embedding files in batches and assign to clusters.

    Returns:
        all_gene_ids: Array of all gene IDs
        all_cluster_assignments: Array of cluster assignments
        all_distances: Array of distances to assigned cluster
    """
    print(f"\nProcessing {len(embedding_files)} embedding files...")

    all_gene_ids = []
    all_genome_ids = []
    all_cluster_assignments = []
    all_distances = []

    for emb_file in tqdm(embedding_files, desc="Processing genomes"):
        try:
            # Load embeddings
            data = np.load(emb_file, allow_pickle=True)
            embeddings = data['embeddings']  # (n_genes, 1152)
            gene_ids = data['gene_ids']
            genome_id = emb_file.stem.replace('_embeddings', '')

            # Apply PCA transform
            embeddings_centered = embeddings - pca_mean
            embeddings_pca = embeddings_centered @ pca_components.T

            # Assign to clusters (find nearest centroid)
            # Use cosine similarity
            similarities = cosine_similarity(embeddings_pca, cluster_centroids)
            nearest_cluster_idx = np.argmax(similarities, axis=1)
            max_similarities = np.max(similarities, axis=1)

            # Convert to cluster IDs
            assigned_clusters = cluster_ids[nearest_cluster_idx]

            # Store results
            all_gene_ids.extend(gene_ids)
            all_genome_ids.extend([genome_id] * len(gene_ids))
            all_cluster_assignments.extend(assigned_clusters)
            all_distances.extend(1.0 - max_similarities)  # Convert similarity to distance

        except Exception as e:
            print(f"\nError processing {emb_file.name}: {e}")
            continue

    print(f"\nTotal genes processed: {len(all_gene_ids):,}")

    return (np.array(all_gene_ids),
            np.array(all_genome_ids),
            np.array(all_cluster_assignments),
            np.array(all_distances))


def compute_cluster_statistics(embedding_files, pca_components, pca_mean,
                               cluster_assignments, gene_ids, cluster_ids):
    """
    Compute mean and standard deviation for each cluster across all dimensions.

    This requires a second pass through the data to accumulate statistics.
    Uses Welford's online algorithm for numerical stability.

    Returns:
        cluster_means: (n_clusters, pca_dim) array
        cluster_stds: (n_clusters, pca_dim) array
        cluster_sizes: (n_clusters,) array
    """
    print("\nComputing cluster statistics (mean, SD per dimension)...")

    n_clusters = len(cluster_ids)
    pca_dim = pca_components.shape[0]

    # Initialize accumulators
    cluster_counts = np.zeros(n_clusters, dtype=np.int64)
    cluster_means = np.zeros((n_clusters, pca_dim))
    cluster_m2 = np.zeros((n_clusters, pca_dim))  # For variance computation

    # Map cluster IDs to indices
    cluster_id_to_idx = {cid: i for i, cid in enumerate(cluster_ids)}

    # Create gene_id to cluster mapping
    gene_to_cluster = dict(zip(gene_ids, cluster_assignments))

    print("  First pass: computing means and variances...")
    for emb_file in tqdm(embedding_files, desc="Processing genomes"):
        try:
            data = np.load(emb_file, allow_pickle=True)
            embeddings = data['embeddings']
            file_gene_ids = data['gene_ids']

            # Apply PCA
            embeddings_centered = embeddings - pca_mean
            embeddings_pca = embeddings_centered @ pca_components.T

            # Accumulate statistics per cluster using Welford's algorithm
            for i, gene_id in enumerate(file_gene_ids):
                if gene_id not in gene_to_cluster:
                    continue

                cluster_id = gene_to_cluster[gene_id]
                if cluster_id not in cluster_id_to_idx:
                    continue

                cluster_idx = cluster_id_to_idx[cluster_id]
                x = embeddings_pca[i]

                cluster_counts[cluster_idx] += 1
                delta = x - cluster_means[cluster_idx]
                cluster_means[cluster_idx] += delta / cluster_counts[cluster_idx]
                delta2 = x - cluster_means[cluster_idx]
                cluster_m2[cluster_idx] += delta * delta2

        except Exception as e:
            print(f"\nError processing {emb_file.name}: {e}")
            continue

    # Compute standard deviations
    cluster_stds = np.zeros((n_clusters, pca_dim))
    for i in range(n_clusters):
        if cluster_counts[i] > 1:
            cluster_stds[i] = np.sqrt(cluster_m2[i] / (cluster_counts[i] - 1))

    print(f"\nCluster statistics computed:")
    print(f"  Clusters: {n_clusters:,}")
    print(f"  Dimensions: {pca_dim}")
    print(f"  Size range: {cluster_counts.min():,} - {cluster_counts.max():,} genes")
    print(f"  Mean size: {cluster_counts.mean():.1f} genes")

    return cluster_means, cluster_stds, cluster_counts


def main():
    parser = argparse.ArgumentParser(
        description='Assign full 29M dataset to clusters from subsample clustering'
    )
    parser.add_argument('--clustering', type=str, required=True,
                        help='Path to subsample clustering .npz file')
    parser.add_argument('--pca', type=str,
                        default='results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz',
                        help='Path to PCA cache with model parameters')
    parser.add_argument('--embedding-dir', type=str,
                        default='data/refseq_esm_embeddings',
                        help='Directory containing all ESM embedding files')
    parser.add_argument('--output', type=str, required=True,
                        help='Output .npz file for full dataset assignments')
    parser.add_argument('--batch-size', type=int, default=100000,
                        help='Batch size for processing embeddings')
    parser.add_argument('--compute-statistics', action='store_true',
                        help='Compute cluster mean/SD statistics (requires second pass)')

    args = parser.parse_args()

    # Setup paths
    clustering_file = Path(args.clustering)
    pca_file = Path(args.pca)
    embedding_dir = Path(args.embedding_dir)
    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print("Full Dataset Cluster Assignment")
    print(f"{'='*80}\n")

    # 1. Load subsample clustering
    labels, subsample_gene_ids, metadata = load_clustering(clustering_file)

    # 2. Load PCA embeddings and model
    pca_data = np.load(pca_file, allow_pickle=True)
    embeddings_pca_full = pca_data['embeddings_pca']
    pca_gene_ids_full = pca_data['gene_ids']
    n_components = embeddings_pca_full.shape[1]

    # Check if PCA model is in cache, otherwise fit from embeddings
    if 'pca_components' in pca_data and 'pca_mean' in pca_data:
        pca_components = pca_data['pca_components']
        pca_mean = pca_data['pca_mean']
        print(f"PCA model loaded from cache: {pca_components.shape[1]}D → {n_components}D")
    else:
        print("PCA model not in cache, fitting from embedding files...")
        pca_components, pca_mean = fit_pca_model(
            embedding_dir, n_genes_sample=len(pca_gene_ids_full),
            n_components=n_components
        )

    # 3. Align subsample clustering with PCA embeddings
    gene_to_idx = {gid: i for i, gid in enumerate(pca_gene_ids_full)}
    subsample_indices = [gene_to_idx[gid] for gid in subsample_gene_ids if gid in gene_to_idx]

    if len(subsample_indices) != len(subsample_gene_ids):
        print(f"Warning: Only {len(subsample_indices)}/{len(subsample_gene_ids)} genes found in PCA")

    embeddings_pca_subsample = embeddings_pca_full[subsample_indices]

    # 4. Compute cluster centroids
    cluster_centroids, cluster_ids = compute_cluster_centroids(
        embeddings_pca_subsample, labels, subsample_gene_ids
    )

    # 5. Find all embedding files
    embedding_files = sorted(embedding_dir.glob('*_embeddings.npz'))
    print(f"\nFound {len(embedding_files)} embedding files in {embedding_dir}")

    if len(embedding_files) == 0:
        raise ValueError(f"No embedding files found in {embedding_dir}")

    # 6. Process all embeddings and assign to clusters
    (all_gene_ids, all_genome_ids,
     all_cluster_assignments, all_distances) = process_embeddings_batch(
        embedding_files, pca_components, pca_mean,
        cluster_centroids, cluster_ids, args.batch_size
    )

    # 7. Compute cluster statistics if requested
    if args.compute_statistics:
        cluster_means, cluster_stds, cluster_sizes = compute_cluster_statistics(
            embedding_files, pca_components, pca_mean,
            all_cluster_assignments, all_gene_ids, cluster_ids
        )
    else:
        cluster_means = None
        cluster_stds = None
        cluster_sizes = None
        print("\nSkipping cluster statistics computation (use --compute-statistics to enable)")

    # 8. Save results
    print(f"\n{'='*80}")
    print(f"Saving results to {output_file}")
    print(f"{'='*80}\n")

    save_dict = {
        'gene_ids': all_gene_ids,
        'genome_ids': all_genome_ids,
        'cluster_assignments': all_cluster_assignments,
        'distances_to_centroid': all_distances,
        'cluster_centroids': cluster_centroids,
        'cluster_ids': cluster_ids,
        'n_genes': len(all_gene_ids),
        'n_clusters': len(cluster_ids),
        **metadata
    }

    if args.compute_statistics:
        save_dict.update({
            'cluster_means': cluster_means,
            'cluster_stds': cluster_stds,
            'cluster_sizes': cluster_sizes
        })

    np.savez_compressed(output_file, **save_dict)

    # Print summary
    print("Summary:")
    print(f"  Total genes assigned: {len(all_gene_ids):,}")
    print(f"  Unique genomes: {len(set(all_genome_ids)):,}")
    print(f"  Clusters: {len(cluster_ids):,}")
    print(f"  Mean distance to centroid: {all_distances.mean():.4f}")
    print(f"  Median distance to centroid: {np.median(all_distances):.4f}")

    if args.compute_statistics:
        print(f"\n  Cluster statistics computed:")
        print(f"    Mean size: {cluster_sizes.mean():.1f} genes")
        print(f"    Size range: {cluster_sizes.min():,} - {cluster_sizes.max():,}")

    print("\nDone!")


if __name__ == '__main__':
    main()
