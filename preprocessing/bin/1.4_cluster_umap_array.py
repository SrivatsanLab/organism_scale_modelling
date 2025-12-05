#!/usr/bin/env python
"""
Clustering parameter sweep for UMAP embeddings.

Tests multiple clustering methods and parameters to find optimal clustering.

Usage:
    python cluster_umap_array.py --umap <umap.npz> --output <output.npz> --method <method> --params <params>
"""

import argparse
import numpy as np
import time
from pathlib import Path
import json


def parse_args():
    parser = argparse.ArgumentParser(description="Cluster UMAP embeddings")
    parser.add_argument("--umap", type=str, required=True, help="UMAP npz file")
    parser.add_argument("--output", type=str, required=True, help="Output npz file")
    parser.add_argument("--method", type=str, required=True,
                       choices=['hdbscan', 'kmeans', 'minibatch_kmeans', 'leiden', 'dbscan', 'gaussian_mixture'],
                       help="Clustering method")
    parser.add_argument("--params", type=str, required=True,
                       help="JSON string of method parameters")
    parser.add_argument("--use-pca", action='store_true',
                       help="Cluster on PCA space instead of UMAP")
    parser.add_argument("--pca-cache", type=str,
                       help="PCA cache file (if using --use-pca)")
    return parser.parse_args()


def load_umap(umap_file):
    """Load UMAP embeddings."""
    print(f"Loading UMAP from {umap_file}...")
    data = np.load(umap_file, allow_pickle=True)

    umap_embedding = data['umap_embedding']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Loaded {len(gene_ids):,} genes")
    print(f"  UMAP shape: {umap_embedding.shape}")

    return umap_embedding, gene_ids, genome_ids


def load_pca(pca_file):
    """Load PCA embeddings."""
    print(f"Loading PCA from {pca_file}...")
    data = np.load(pca_file, allow_pickle=True)

    pca_embedding = data['embeddings_pca']

    print(f"  PCA shape: {pca_embedding.shape}")

    return pca_embedding


def cluster_hdbscan(embeddings, params):
    """Cluster using HDBSCAN."""
    import hdbscan

    min_cluster_size = params.get('min_cluster_size', 100)
    min_samples = params.get('min_samples', None)
    cluster_selection_epsilon = params.get('cluster_selection_epsilon', 0.0)
    metric = params.get('metric', 'euclidean')

    print(f"Running HDBSCAN...")
    print(f"  min_cluster_size: {min_cluster_size}")
    print(f"  min_samples: {min_samples}")
    print(f"  cluster_selection_epsilon: {cluster_selection_epsilon}")
    print(f"  metric: {metric}")

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        metric=metric,
        core_dist_n_jobs=-1
    )

    labels = clusterer.fit_predict(embeddings)

    # Get additional metrics
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = np.sum(labels == -1)

    print(f"  Clusters found: {n_clusters}")
    print(f"  Noise points: {n_noise:,} ({100*n_noise/len(labels):.1f}%)")

    # Calculate cluster sizes
    unique, counts = np.unique(labels[labels != -1], return_counts=True)
    if len(unique) > 0:
        print(f"  Median cluster size: {np.median(counts):.0f}")
        print(f"  Mean cluster size: {np.mean(counts):.0f}")
        print(f"  Largest cluster: {np.max(counts):,}")
        print(f"  Smallest cluster: {np.min(counts):,}")

    return labels, {'n_clusters': n_clusters, 'n_noise': n_noise}


def cluster_kmeans(embeddings, params):
    """Cluster using K-means."""
    from sklearn.cluster import KMeans

    n_clusters = params.get('n_clusters', 10)

    print(f"Running K-means...")
    print(f"  n_clusters: {n_clusters}")

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(embeddings)

    # Calculate cluster sizes
    unique, counts = np.unique(labels, return_counts=True)
    print(f"  Median cluster size: {np.median(counts):.0f}")
    print(f"  Mean cluster size: {np.mean(counts):.0f}")
    print(f"  Inertia: {kmeans.inertia_:.2e}")

    return labels, {'n_clusters': n_clusters, 'inertia': float(kmeans.inertia_)}


def cluster_minibatch_kmeans(embeddings, params):
    """Cluster using MiniBatch K-means (faster than regular K-means)."""
    from sklearn.cluster import MiniBatchKMeans

    n_clusters = params.get('n_clusters', 10)
    batch_size = params.get('batch_size', 1024)

    print(f"Running MiniBatch K-means...")
    print(f"  n_clusters: {n_clusters}")
    print(f"  batch_size: {batch_size}")

    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=42,
        batch_size=batch_size,
        n_init=3
    )
    labels = kmeans.fit_predict(embeddings)

    # Calculate cluster sizes
    unique, counts = np.unique(labels, return_counts=True)
    print(f"  Median cluster size: {np.median(counts):.0f}")
    print(f"  Mean cluster size: {np.mean(counts):.0f}")
    print(f"  Inertia: {kmeans.inertia_:.2e}")

    return labels, {'n_clusters': n_clusters, 'inertia': float(kmeans.inertia_)}


def cluster_leiden(embeddings, params):
    """Cluster using Leiden algorithm on k-NN graph."""
    import scanpy as sc
    import anndata as ad

    resolution = params.get('resolution', 1.0)
    n_neighbors = params.get('n_neighbors', 15)

    print(f"Running Leiden clustering...")
    print(f"  resolution: {resolution}")
    print(f"  n_neighbors: {n_neighbors}")

    # Create AnnData object
    adata = ad.AnnData(X=embeddings)

    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X')

    # Run Leiden
    sc.tl.leiden(adata, resolution=resolution)

    labels = adata.obs['leiden'].astype(int).values

    n_clusters = len(set(labels))
    unique, counts = np.unique(labels, return_counts=True)

    print(f"  Clusters found: {n_clusters}")
    print(f"  Median cluster size: {np.median(counts):.0f}")
    print(f"  Mean cluster size: {np.mean(counts):.0f}")

    return labels, {'n_clusters': n_clusters}


def cluster_dbscan(embeddings, params):
    """Cluster using DBSCAN."""
    from sklearn.cluster import DBSCAN

    eps = params.get('eps', 0.5)
    min_samples = params.get('min_samples', 5)

    print(f"Running DBSCAN...")
    print(f"  eps: {eps}")
    print(f"  min_samples: {min_samples}")

    dbscan = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
    labels = dbscan.fit_predict(embeddings)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = np.sum(labels == -1)

    print(f"  Clusters found: {n_clusters}")
    print(f"  Noise points: {n_noise:,} ({100*n_noise/len(labels):.1f}%)")

    if n_clusters > 0:
        unique, counts = np.unique(labels[labels != -1], return_counts=True)
        print(f"  Median cluster size: {np.median(counts):.0f}")
        print(f"  Mean cluster size: {np.mean(counts):.0f}")

    return labels, {'n_clusters': n_clusters, 'n_noise': n_noise}


def cluster_gaussian_mixture(embeddings, params):
    """Cluster using Gaussian Mixture Model."""
    from sklearn.mixture import GaussianMixture

    n_components = params.get('n_components', 10)
    covariance_type = params.get('covariance_type', 'full')

    print(f"Running Gaussian Mixture...")
    print(f"  n_components: {n_components}")
    print(f"  covariance_type: {covariance_type}")

    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type=covariance_type,
        random_state=42
    )
    labels = gmm.fit_predict(embeddings)

    unique, counts = np.unique(labels, return_counts=True)
    print(f"  Median cluster size: {np.median(counts):.0f}")
    print(f"  Mean cluster size: {np.mean(counts):.0f}")
    print(f"  BIC: {gmm.bic(embeddings):.2e}")
    print(f"  AIC: {gmm.aic(embeddings):.2e}")

    return labels, {
        'n_clusters': n_components,
        'bic': float(gmm.bic(embeddings)),
        'aic': float(gmm.aic(embeddings))
    }


def main():
    args = parse_args()

    print("=" * 70)
    print("UMAP Clustering Parameter Sweep")
    print("=" * 70)
    print(f"Method: {args.method}")
    print(f"Parameters: {args.params}")
    print("")

    # Load UMAP
    umap_embedding, gene_ids, genome_ids = load_umap(args.umap)

    # Decide what to cluster on
    if args.use_pca:
        if not args.pca_cache:
            raise ValueError("--pca-cache required when using --use-pca")
        embeddings = load_pca(args.pca_cache)
    else:
        embeddings = umap_embedding

    # Parse parameters
    params = json.loads(args.params)

    # Run clustering
    start_time = time.time()

    if args.method == 'hdbscan':
        labels, metrics = cluster_hdbscan(embeddings, params)
    elif args.method == 'kmeans':
        labels, metrics = cluster_kmeans(embeddings, params)
    elif args.method == 'minibatch_kmeans':
        labels, metrics = cluster_minibatch_kmeans(embeddings, params)
    elif args.method == 'leiden':
        labels, metrics = cluster_leiden(embeddings, params)
    elif args.method == 'dbscan':
        labels, metrics = cluster_dbscan(embeddings, params)
    elif args.method == 'gaussian_mixture':
        labels, metrics = cluster_gaussian_mixture(embeddings, params)
    else:
        raise ValueError(f"Unknown method: {args.method}")

    elapsed = time.time() - start_time
    print(f"\nClustering complete in {elapsed:.1f} seconds")

    # Save results
    print(f"\nSaving results to {args.output}...")
    np.savez_compressed(
        args.output,
        cluster_labels=labels,
        gene_ids=gene_ids,
        genome_ids=genome_ids,
        method=args.method,
        params=json.dumps(params),
        metrics=json.dumps(metrics),
        elapsed_seconds=elapsed,
        use_pca=args.use_pca
    )

    print("=" * 70)
    print("Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
