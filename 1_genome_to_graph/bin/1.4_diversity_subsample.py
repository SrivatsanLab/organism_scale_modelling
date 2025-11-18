#!/usr/bin/env python
"""
Diversity-based subsampling for large datasets.

Instead of random sampling, select a diverse subset where no cluster dominates.
This ensures better representation of the full dataset in downstream analysis.

Strategies implemented:
1. Cluster-based sampling: Cluster first, then sample uniformly from each cluster
2. MaxMin sampling: Greedy farthest-point sampling
3. Hybrid: Cluster first, then farthest-point within each cluster
"""

import argparse
import numpy as np
from pathlib import Path
import time
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import euclidean_distances


def load_embeddings(embedding_file, use_pca=False):
    """Load embeddings from file."""
    print(f"Loading embeddings from {embedding_file}...")
    data = np.load(embedding_file, allow_pickle=True)

    if use_pca and 'embeddings_pca' in data:
        embeddings = data['embeddings_pca']
        print(f"  Using PCA embeddings")
    else:
        embeddings = data['embeddings']
        print(f"  Using full embeddings")

    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Shape: {embeddings.shape}")
    print(f"  Genes: {len(gene_ids):,}")

    return embeddings, gene_ids, genome_ids


def cluster_based_sampling(embeddings, gene_ids, genome_ids, n_samples, n_clusters=10000):
    """
    Cluster-based diversity sampling.

    1. Cluster data into k clusters
    2. Sample uniformly from each cluster (proportional to cluster size or uniform)

    This ensures all regions of the space are represented.
    """
    print(f"\nCluster-based sampling:")
    print(f"  Target samples: {n_samples:,}")
    print(f"  Clusters: {n_clusters:,}")

    # Cluster
    print(f"  Running MiniBatch K-means...")
    start = time.time()
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=42,
        batch_size=10000,
        max_iter=100,
        verbose=1
    )
    cluster_labels = kmeans.fit_predict(embeddings)
    print(f"  Clustering took {time.time() - start:.1f}s")

    # Count samples per cluster
    unique_clusters, cluster_counts = np.unique(cluster_labels, return_counts=True)
    print(f"  Cluster sizes: min={cluster_counts.min()}, max={cluster_counts.max()}, mean={cluster_counts.mean():.0f}")

    # Sample uniformly from each cluster
    samples_per_cluster = n_samples // n_clusters
    remainder = n_samples % n_clusters

    print(f"  Sampling {samples_per_cluster} per cluster (+ {remainder} from largest clusters)")

    selected_indices = []

    for i, cluster_id in enumerate(unique_clusters):
        cluster_mask = cluster_labels == cluster_id
        cluster_indices = np.where(cluster_mask)[0]

        # Larger clusters get the remainder samples
        n_to_sample = samples_per_cluster
        if i < remainder:
            n_to_sample += 1

        # Random sample from this cluster
        if len(cluster_indices) <= n_to_sample:
            # Take all if cluster is small
            selected_indices.extend(cluster_indices)
        else:
            sampled = np.random.choice(cluster_indices, size=n_to_sample, replace=False)
            selected_indices.extend(sampled)

    selected_indices = np.array(selected_indices)

    print(f"  Selected {len(selected_indices):,} samples")

    return selected_indices, cluster_labels


def maxmin_sampling(embeddings, gene_ids, genome_ids, n_samples, batch_size=10000):
    """
    MaxMin (farthest-point) sampling.

    Greedy algorithm:
    1. Start with a random point
    2. Iteratively add the point that's farthest from all selected points

    This is slow for large datasets, so we do it in batches.
    """
    print(f"\nMaxMin (farthest-point) sampling:")
    print(f"  Target samples: {n_samples:,}")
    print(f"  Note: This is slow, using batched approximation")

    selected_indices = []

    # Start with random point
    first_idx = np.random.randint(len(embeddings))
    selected_indices.append(first_idx)

    # Work on a subset for efficiency
    if len(embeddings) > 1000000:
        print(f"  Dataset too large, pre-sampling to 1M for MaxMin")
        candidate_indices = np.random.choice(len(embeddings), size=min(1000000, len(embeddings)), replace=False)
        candidate_embeddings = embeddings[candidate_indices]
    else:
        candidate_indices = np.arange(len(embeddings))
        candidate_embeddings = embeddings

    # Track minimum distances to selected set
    min_distances = np.full(len(candidate_embeddings), np.inf)

    for i in range(1, n_samples):
        if i % 100 == 0:
            print(f"  Selected {i}/{n_samples} samples...")

        # Get last selected point
        last_selected = embeddings[selected_indices[-1]]

        # Compute distances to last selected
        distances = euclidean_distances(
            candidate_embeddings,
            last_selected.reshape(1, -1)
        ).flatten()

        # Update minimum distances
        min_distances = np.minimum(min_distances, distances)

        # Select point with maximum minimum distance
        farthest_idx_in_candidates = np.argmax(min_distances)
        farthest_idx = candidate_indices[farthest_idx_in_candidates]

        selected_indices.append(farthest_idx)

        # Remove from candidates
        min_distances[farthest_idx_in_candidates] = -np.inf

    selected_indices = np.array(selected_indices)

    print(f"  Selected {len(selected_indices):,} samples")

    return selected_indices, None


def hybrid_sampling(embeddings, gene_ids, genome_ids, n_samples, n_clusters=1000):
    """
    Hybrid: Cluster first, then farthest-point within each cluster.

    This is faster than pure MaxMin while still being diverse.
    """
    print(f"\nHybrid sampling (cluster + farthest-point):")
    print(f"  Target samples: {n_samples:,}")
    print(f"  Clusters: {n_clusters:,}")

    # Cluster
    print(f"  Running MiniBatch K-means...")
    start = time.time()
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=42,
        batch_size=10000,
        max_iter=100,
        verbose=1
    )
    cluster_labels = kmeans.fit_predict(embeddings)
    print(f"  Clustering took {time.time() - start:.1f}s")

    # Samples per cluster
    samples_per_cluster = n_samples // n_clusters
    remainder = n_samples % n_clusters

    print(f"  Sampling ~{samples_per_cluster} per cluster using farthest-point")

    selected_indices = []

    for i, cluster_id in enumerate(range(n_clusters)):
        if i % 100 == 0:
            print(f"    Processing cluster {i}/{n_clusters}...")

        cluster_mask = cluster_labels == cluster_id
        cluster_indices = np.where(cluster_mask)[0]

        if len(cluster_indices) == 0:
            continue

        n_to_sample = samples_per_cluster
        if i < remainder:
            n_to_sample += 1

        if len(cluster_indices) <= n_to_sample:
            selected_indices.extend(cluster_indices)
            continue

        # Farthest-point sampling within cluster
        cluster_embeddings = embeddings[cluster_indices]

        # Start with cluster center
        center = kmeans.cluster_centers_[cluster_id]
        distances_to_center = euclidean_distances(
            cluster_embeddings,
            center.reshape(1, -1)
        ).flatten()

        # Pick farthest from center
        local_selected = [np.argmax(distances_to_center)]
        min_distances = distances_to_center.copy()

        # Iteratively add farthest points
        for _ in range(n_to_sample - 1):
            last_point = cluster_embeddings[local_selected[-1]]
            distances = euclidean_distances(
                cluster_embeddings,
                last_point.reshape(1, -1)
            ).flatten()

            min_distances = np.minimum(min_distances, distances)

            # Don't select already-selected points
            for idx in local_selected:
                min_distances[idx] = -np.inf

            next_idx = np.argmax(min_distances)
            local_selected.append(next_idx)

        # Convert local indices to global
        global_indices = cluster_indices[local_selected]
        selected_indices.extend(global_indices)

    selected_indices = np.array(selected_indices)

    print(f"  Selected {len(selected_indices):,} samples")

    return selected_indices, cluster_labels


def save_subsample(output_file, embeddings, gene_ids, genome_ids,
                   selected_indices, cluster_labels, method, n_clusters):
    """Save subsampled data."""
    print(f"\nSaving to {output_file}...")

    output_data = {
        'embeddings': embeddings[selected_indices],
        'gene_ids': gene_ids[selected_indices],
        'genome_ids': genome_ids[selected_indices],
        'selected_indices': selected_indices,
        'subsample_method': method,
        'n_original': len(embeddings),
        'n_subsampled': len(selected_indices)
    }

    if cluster_labels is not None:
        output_data['cluster_labels_full'] = cluster_labels
        output_data['cluster_labels_subsample'] = cluster_labels[selected_indices]
        output_data['n_clusters'] = n_clusters

    np.savez_compressed(output_file, **output_data)
    print(f"  Saved {len(selected_indices):,} samples")


def main():
    parser = argparse.ArgumentParser(description='Diversity-based subsampling')
    parser.add_argument('--input', required=True, help='Input embedding file')
    parser.add_argument('--output', required=True, help='Output subsampled file')
    parser.add_argument('--n-samples', type=int, required=True, help='Number of samples to select')
    parser.add_argument('--method', choices=['cluster', 'maxmin', 'hybrid'], default='hybrid',
                       help='Sampling method')
    parser.add_argument('--n-clusters', type=int, default=1000,
                       help='Number of clusters (for cluster/hybrid methods)')
    parser.add_argument('--use-pca', action='store_true', help='Use PCA embeddings if available')

    args = parser.parse_args()

    print("=" * 80)
    print("Diversity-Based Subsampling")
    print("=" * 80)
    print(f"Method: {args.method}")
    print(f"Target samples: {args.n_samples:,}")
    if args.method in ['cluster', 'hybrid']:
        print(f"Clusters: {args.n_clusters:,}")
    print()

    # Load embeddings
    embeddings, gene_ids, genome_ids = load_embeddings(args.input, args.use_pca)

    # Run sampling
    if args.method == 'cluster':
        selected_indices, cluster_labels = cluster_based_sampling(
            embeddings, gene_ids, genome_ids, args.n_samples, args.n_clusters
        )
    elif args.method == 'maxmin':
        selected_indices, cluster_labels = maxmin_sampling(
            embeddings, gene_ids, genome_ids, args.n_samples
        )
    elif args.method == 'hybrid':
        selected_indices, cluster_labels = hybrid_sampling(
            embeddings, gene_ids, genome_ids, args.n_samples, args.n_clusters
        )

    # Save
    save_subsample(
        args.output, embeddings, gene_ids, genome_ids,
        selected_indices, cluster_labels, args.method, args.n_clusters
    )

    print()
    print("=" * 80)
    print("✓ Diversity subsampling complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
