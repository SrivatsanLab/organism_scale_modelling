#!/usr/bin/env python
"""
GPU-accelerated dimensionality reduction and clustering using RAPIDS cuML.

Supports:
- GPU UMAP (faster than CPU)
- GPU HDBSCAN (much faster than CPU on large datasets)
- GPU K-means (very fast)
- GPU t-SNE (alternative to UMAP)

Uses RAPIDS cuML library which requires NVIDIA GPU.
"""

import argparse
import numpy as np
import time
from pathlib import Path
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="GPU-accelerated embeddings and clustering")

    parser.add_argument("--input", type=str, required=True,
                       help="Input embeddings (PCA cache .npz file)")
    parser.add_argument("--output", type=str, required=True,
                       help="Output file path (.npz)")

    # Embedding method
    parser.add_argument("--method", type=str, default="umap",
                       choices=["umap", "tsne"],
                       help="2D embedding method (GPU accelerated)")

    # UMAP parameters
    parser.add_argument("--n-neighbors", type=int, default=15,
                       help="UMAP n_neighbors")
    parser.add_argument("--min-dist", type=float, default=0.1,
                       help="UMAP min_dist")

    # t-SNE parameters
    parser.add_argument("--perplexity", type=float, default=30.0,
                       help="t-SNE perplexity")

    # Clustering (optional)
    parser.add_argument("--cluster", type=str, default=None,
                       choices=["hdbscan", "kmeans"],
                       help="Clustering method (optional)")
    parser.add_argument("--n-clusters", type=int, default=100,
                       help="Number of clusters (for kmeans)")
    parser.add_argument("--min-cluster-size", type=int, default=100,
                       help="Min cluster size (for HDBSCAN)")

    # Subsampling
    parser.add_argument("--subsample", type=int, default=None,
                       help="Subsample N genes (None = use all)")

    return parser.parse_args()


def load_embeddings(input_file, subsample=None):
    """Load PCA embeddings."""
    print(f"Loading embeddings from {input_file}...")
    data = np.load(input_file, allow_pickle=True)

    embeddings = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Loaded {len(gene_ids):,} genes")
    print(f"  Embedding shape: {embeddings.shape}")
    print(f"  Memory: {embeddings.nbytes / 1e9:.2f} GB")

    # Subsample if requested
    if subsample and subsample < len(gene_ids):
        print(f"\nSubsampling {subsample:,} genes...")
        np.random.seed(42)
        indices = np.random.choice(len(gene_ids), subsample, replace=False)
        indices = np.sort(indices)

        embeddings = embeddings[indices]
        gene_ids = gene_ids[indices]
        genome_ids = genome_ids[indices]

        print(f"  Subsampled shape: {embeddings.shape}")

    return embeddings, gene_ids, genome_ids


def compute_gpu_umap(embeddings, n_neighbors=15, min_dist=0.1):
    """Compute UMAP using GPU."""
    print("\n" + "="*70)
    print("GPU UMAP")
    print("="*70)

    try:
        from cuml.manifold import UMAP
        import cudf
    except ImportError:
        print("ERROR: cuML not available. Need RAPIDS cuML for GPU acceleration.")
        print("Run this script in a RAPIDS container or conda environment.")
        sys.exit(1)

    print(f"Parameters:")
    print(f"  n_neighbors: {n_neighbors}")
    print(f"  min_dist: {min_dist}")
    print(f"  metric: euclidean")

    # Convert to float32 for GPU
    embeddings = embeddings.astype(np.float32)

    start_time = time.time()

    # Create GPU UMAP model
    umap_model = UMAP(
        n_neighbors=n_neighbors,
        n_components=2,
        min_dist=min_dist,
        metric='euclidean',
        random_state=42,
        verbose=True
    )

    # Fit and transform
    print("\nFitting UMAP...")
    embedding_2d = umap_model.fit_transform(embeddings)

    # Convert back to numpy on CPU
    if hasattr(embedding_2d, 'to_numpy'):
        embedding_2d = embedding_2d.to_numpy()
    elif hasattr(embedding_2d, 'get'):
        embedding_2d = embedding_2d.get()

    elapsed = time.time() - start_time
    print(f"\nUMAP completed in {elapsed:.1f} seconds")
    print(f"Output shape: {embedding_2d.shape}")

    return embedding_2d


def compute_gpu_tsne(embeddings, perplexity=30.0):
    """Compute t-SNE using GPU."""
    print("\n" + "="*70)
    print("GPU t-SNE")
    print("="*70)

    try:
        from cuml.manifold import TSNE
    except ImportError:
        print("ERROR: cuML not available.")
        sys.exit(1)

    print(f"Parameters:")
    print(f"  perplexity: {perplexity}")

    embeddings = embeddings.astype(np.float32)

    start_time = time.time()

    tsne_model = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=42,
        verbose=True
    )

    print("\nFitting t-SNE...")
    embedding_2d = tsne_model.fit_transform(embeddings)

    if hasattr(embedding_2d, 'to_numpy'):
        embedding_2d = embedding_2d.to_numpy()
    elif hasattr(embedding_2d, 'get'):
        embedding_2d = embedding_2d.get()

    elapsed = time.time() - start_time
    print(f"\nt-SNE completed in {elapsed:.1f} seconds")

    return embedding_2d


def cluster_gpu_hdbscan(embeddings, min_cluster_size=100):
    """Cluster using GPU HDBSCAN."""
    print("\n" + "="*70)
    print("GPU HDBSCAN Clustering")
    print("="*70)

    try:
        from cuml.cluster import HDBSCAN
    except ImportError:
        print("ERROR: cuML not available.")
        sys.exit(1)

    print(f"Parameters:")
    print(f"  min_cluster_size: {min_cluster_size}")

    embeddings = embeddings.astype(np.float32)

    start_time = time.time()

    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        metric='euclidean',
        min_samples=None
    )

    print("\nFitting HDBSCAN...")
    labels = clusterer.fit_predict(embeddings)

    if hasattr(labels, 'to_numpy'):
        labels = labels.to_numpy()
    elif hasattr(labels, 'get'):
        labels = labels.get()

    elapsed = time.time() - start_time

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()

    print(f"\nHDBSCAN completed in {elapsed:.1f} seconds")
    print(f"  Clusters: {n_clusters}")
    print(f"  Noise points: {n_noise:,} ({100*n_noise/len(labels):.1f}%)")

    return labels


def cluster_gpu_kmeans(embeddings, n_clusters=100):
    """Cluster using GPU K-means."""
    print("\n" + "="*70)
    print("GPU K-means Clustering")
    print("="*70)

    try:
        from cuml.cluster import KMeans
    except ImportError:
        print("ERROR: cuML not available.")
        sys.exit(1)

    print(f"Parameters:")
    print(f"  n_clusters: {n_clusters}")

    embeddings = embeddings.astype(np.float32)

    start_time = time.time()

    kmeans = KMeans(
        n_clusters=n_clusters,
        random_state=42,
        max_iter=300
    )

    print("\nFitting K-means...")
    labels = kmeans.fit_predict(embeddings)

    if hasattr(labels, 'to_numpy'):
        labels = labels.to_numpy()
    elif hasattr(labels, 'get'):
        labels = labels.get()

    elapsed = time.time() - start_time

    print(f"\nK-means completed in {elapsed:.1f} seconds")
    print(f"  Clusters: {n_clusters}")

    return labels


def main():
    args = parse_args()

    print("="*70)
    print("GPU-Accelerated Embeddings and Clustering")
    print("="*70)
    print()
    print(f"Method: {args.method}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    if args.cluster:
        print(f"Clustering: {args.cluster}")
    print()

    # Load embeddings
    embeddings, gene_ids, genome_ids = load_embeddings(
        args.input,
        subsample=args.subsample
    )

    # Compute 2D embedding
    if args.method == "umap":
        embedding_2d = compute_gpu_umap(
            embeddings,
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist
        )
    elif args.method == "tsne":
        embedding_2d = compute_gpu_tsne(
            embeddings,
            perplexity=args.perplexity
        )

    # Optional clustering
    cluster_labels = None
    if args.cluster == "hdbscan":
        cluster_labels = cluster_gpu_hdbscan(
            embeddings,
            min_cluster_size=args.min_cluster_size
        )
    elif args.cluster == "kmeans":
        cluster_labels = cluster_gpu_kmeans(
            embeddings,
            n_clusters=args.n_clusters
        )

    # Save results
    print(f"\nSaving results to {args.output}...")
    output_data = {
        'embedding_2d': embedding_2d,
        'gene_ids': gene_ids,
        'genome_ids': genome_ids,
        'method': args.method
    }

    if cluster_labels is not None:
        output_data['cluster_labels'] = cluster_labels
        output_data['cluster_method'] = args.cluster

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **output_data)

    print("\n" + "="*70)
    print("Complete!")
    print("="*70)


if __name__ == "__main__":
    main()
