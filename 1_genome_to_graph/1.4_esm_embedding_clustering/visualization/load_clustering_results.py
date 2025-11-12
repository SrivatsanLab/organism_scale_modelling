#!/usr/bin/env python
"""
Helper script to load and explore clustering results.

Usage in notebook:
    from scripts.embeddings.load_clustering_results import load_cluster_labels, get_clustering_summary

    # Load specific clustering
    labels = load_cluster_labels('hdbscan_umap_minclust500')

    # Get summary of all results
    summary = get_clustering_summary()
    print(summary)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import json


def load_cluster_labels(suffix, results_dir='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering'):
    """
    Load cluster labels from a specific clustering result.

    Args:
        suffix: Clustering result suffix (e.g., 'hdbscan_umap_minclust500')
        results_dir: Directory containing clustering results

    Returns:
        cluster_labels: numpy array of cluster labels
    """
    file_path = Path(results_dir) / f"clusters_{suffix}.npz"

    if not file_path.exists():
        raise FileNotFoundError(f"Clustering result not found: {file_path}")

    data = np.load(file_path, allow_pickle=True)

    return data['cluster_labels']


def load_cluster_data(suffix, results_dir='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering'):
    """
    Load full clustering data including metadata.

    Args:
        suffix: Clustering result suffix
        results_dir: Directory containing clustering results

    Returns:
        dict with cluster_labels, gene_ids, genome_ids, method, params, metrics
    """
    file_path = Path(results_dir) / f"clusters_{suffix}.npz"

    if not file_path.exists():
        raise FileNotFoundError(f"Clustering result not found: {file_path}")

    data = np.load(file_path, allow_pickle=True)

    return {
        'cluster_labels': data['cluster_labels'],
        'gene_ids': data['gene_ids'],
        'genome_ids': data['genome_ids'],
        'method': str(data['method']),
        'params': json.loads(str(data['params'])),
        'metrics': json.loads(str(data['metrics'])),
        'elapsed_seconds': float(data['elapsed_seconds']),
        'use_pca': bool(data['use_pca'])
    }


def get_clustering_summary(results_dir='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering'):
    """
    Get summary statistics for all clustering results.

    Returns:
        DataFrame with columns: suffix, method, n_clusters, n_noise, params, metrics
    """
    results_dir = Path(results_dir)

    if not results_dir.exists():
        print(f"Results directory not found: {results_dir}")
        return pd.DataFrame()

    cluster_files = sorted(results_dir.glob("clusters_*.npz"))

    if not cluster_files:
        print(f"No clustering results found in {results_dir}")
        return pd.DataFrame()

    summaries = []

    for file_path in cluster_files:
        suffix = file_path.stem.replace('clusters_', '')

        try:
            data = load_cluster_data(suffix, results_dir)

            labels = data['cluster_labels']
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise = np.sum(labels == -1) if -1 in labels else 0
            n_genes = len(labels)

            # Get cluster size statistics
            if n_clusters > 0:
                unique, counts = np.unique(labels[labels != -1], return_counts=True)
                median_size = np.median(counts)
                mean_size = np.mean(counts)
                min_size = np.min(counts)
                max_size = np.max(counts)
            else:
                median_size = mean_size = min_size = max_size = 0

            summaries.append({
                'suffix': suffix,
                'method': data['method'],
                'use_pca': data['use_pca'],
                'n_clusters': n_clusters,
                'n_noise': n_noise,
                'pct_noise': 100 * n_noise / n_genes,
                'median_cluster_size': median_size,
                'mean_cluster_size': mean_size,
                'min_cluster_size': min_size,
                'max_cluster_size': max_size,
                'params': str(data['params']),
                'elapsed_sec': data['elapsed_seconds']
            })

        except Exception as e:
            print(f"Error loading {suffix}: {e}")

    df = pd.DataFrame(summaries)

    if len(df) > 0:
        # Sort by method and n_clusters
        df = df.sort_values(['method', 'n_clusters'])

    return df


def plot_clustering_comparison(umap_embedding, cluster_labels_dict, figsize=(20, 12)):
    """
    Plot multiple clustering results side-by-side.

    Args:
        umap_embedding: (N, 2) array of UMAP coordinates
        cluster_labels_dict: dict mapping names to cluster label arrays
        figsize: Figure size

    Example:
        plot_clustering_comparison(
            umap_embedding,
            {
                'HDBSCAN min=100': labels1,
                'HDBSCAN min=500': labels2,
                'K-means k=50': labels3
            }
        )
    """
    import matplotlib.pyplot as plt
    import numpy as np

    n_plots = len(cluster_labels_dict)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, (name, labels) in enumerate(cluster_labels_dict.items()):
        ax = axes[idx]

        # Handle noise points
        mask_noise = labels == -1
        mask_clustered = ~mask_noise

        # Plot noise points first (gray, small)
        if np.any(mask_noise):
            ax.scatter(
                umap_embedding[mask_noise, 0],
                umap_embedding[mask_noise, 1],
                c='lightgray',
                s=0.1,
                alpha=0.5,
                rasterized=True,
                label='Noise'
            )

        # Plot clustered points
        if np.any(mask_clustered):
            ax.scatter(
                umap_embedding[mask_clustered, 0],
                umap_embedding[mask_clustered, 1],
                c=labels[mask_clustered],
                s=0.5,
                alpha=0.5,
                cmap='tab20',
                rasterized=True
            )

        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = np.sum(mask_noise)

        ax.set_title(f'{name}\n{n_clusters} clusters, {n_noise:,} noise ({100*n_noise/len(labels):.1f}%)',
                    fontsize=10)
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.grid(True, alpha=0.2)

    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Example usage
    print("Clustering Results Summary")
    print("=" * 80)

    summary = get_clustering_summary()

    if len(summary) > 0:
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', 50)

        print(summary.to_string(index=False))
        print(f"\nTotal clustering results: {len(summary)}")
    else:
        print("No clustering results found.")
        print("Run the clustering sweep first:")
        print("  sbatch scripts/embeddings/submit_clustering_sweep.sh")
