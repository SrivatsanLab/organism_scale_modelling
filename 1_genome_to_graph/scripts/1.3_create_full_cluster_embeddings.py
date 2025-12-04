#!/usr/bin/env python3
"""
Create full cluster embeddings (1152D) for all 388k high-confidence clusters.

This script uses the same direct approach as the toy graph but keeps full dimensionality.

Input:
  - MMseqs cluster statistics CSV (388k clusters)
  - MMseqs cluster assignments TSV
  - ESM embeddings (7,664 genome files)

Output:
  - graph_outputs/protein_graph/cluster_embeddings.npz (full 1152D)
  - graph_outputs/protein_graph/cluster_index.csv (cluster metadata)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from collections import defaultdict
from tqdm import tqdm


def load_used_clusters(cluster_stats_file):
    """
    Load the set of cluster IDs from the statistics file.

    These are the 388k clusters that passed filtering criteria.
    """
    print(f"\nLoading cluster statistics from {cluster_stats_file}")
    df = pd.read_csv(cluster_stats_file)

    used_clusters = set(df['cluster_representative'].values)
    print(f"Found {len(used_clusters):,} clusters in statistics file")

    return used_clusters


def load_cluster_assignments(cluster_file, used_clusters):
    """
    Load protein->cluster assignments for the used clusters only.

    Returns:
        dict: protein_id -> cluster_id
        dict: cluster_id -> size
    """
    print(f"\nLoading cluster assignments from {cluster_file}")

    protein_to_cluster = {}
    cluster_sizes = defaultdict(int)

    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Reading clusters"):
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_rep, protein_id = parts

                # Only keep if this cluster is in the used set
                if cluster_rep in used_clusters:
                    protein_to_cluster[protein_id] = cluster_rep
                    cluster_sizes[cluster_rep] += 1

    print(f"Loaded {len(protein_to_cluster):,} proteins in {len(cluster_sizes):,} clusters")

    return protein_to_cluster, cluster_sizes


def load_embeddings_for_clusters(embedding_dir, protein_to_cluster):
    """
    Load ESM embeddings only for proteins in the used clusters.

    Returns:
        dict: cluster_id -> list of embeddings (1152D)
    """
    embedding_dir = Path(embedding_dir)
    embedding_files = sorted(embedding_dir.glob("*.npz"))

    print(f"\nFound {len(embedding_files)} embedding files")

    cluster_embeddings = defaultdict(list)
    proteins_found = 0
    proteins_skipped = 0

    for emb_file in tqdm(embedding_files, desc="Loading embeddings"):
        # Extract genome ID from filename
        genome_id = emb_file.stem.split('_esmc_embeddings')[0]

        data = np.load(emb_file, allow_pickle=True)
        embeddings = data['embeddings']  # Shape: (n_proteins, 1152)
        seq_ids = data['seq_ids']

        for protein_id, embedding in zip(seq_ids, embeddings):
            # Prepend genome ID to match cluster file format
            protein_id_full = f"{genome_id}_{protein_id}"

            if protein_id_full in protein_to_cluster:
                cluster_id = protein_to_cluster[protein_id_full]
                cluster_embeddings[cluster_id].append(embedding)
                proteins_found += 1
            else:
                proteins_skipped += 1

    print(f"\nProteins with embeddings: {proteins_found:,}")
    print(f"Proteins skipped (not in used clusters): {proteins_skipped:,}")
    print(f"Clusters with embeddings: {len(cluster_embeddings):,}")

    return cluster_embeddings


def compute_cluster_means(cluster_embeddings):
    """
    Compute mean embedding for each cluster.

    Returns:
        np.ndarray: (n_clusters, 1152) mean embeddings
        list: cluster IDs (in same order as matrix rows)
    """
    print("\nComputing mean embeddings per cluster...")

    cluster_ids = sorted(cluster_embeddings.keys())
    n_clusters = len(cluster_ids)
    embedding_dim = len(cluster_embeddings[cluster_ids[0]][0])

    print(f"Computing means for {n_clusters:,} clusters")
    print(f"Embedding dimension: {embedding_dim}")

    mean_embeddings = np.zeros((n_clusters, embedding_dim), dtype=np.float32)

    for i, cluster_id in enumerate(tqdm(cluster_ids, desc="Computing means")):
        embeddings = np.array(cluster_embeddings[cluster_id])
        mean_embeddings[i] = embeddings.mean(axis=0)

    return mean_embeddings, cluster_ids


def save_outputs(mean_embeddings, cluster_ids, cluster_sizes, output_dir):
    """
    Save cluster embeddings and metadata.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nSaving outputs to {output_dir}")

    # Save embeddings
    embeddings_file = output_dir / "cluster_embeddings.npz"
    np.savez_compressed(
        embeddings_file,
        embeddings=mean_embeddings,
        cluster_ids=np.array(cluster_ids, dtype=object)
    )
    print(f"Saved embeddings: {embeddings_file}")
    print(f"  Shape: {mean_embeddings.shape}")
    print(f"  Size: {embeddings_file.stat().st_size / 1e6:.1f} MB")

    # Save cluster index
    cluster_df = pd.DataFrame({
        'cluster_id': cluster_ids,
        'cluster_size': [cluster_sizes[c] for c in cluster_ids]
    })

    index_file = output_dir / "cluster_index.csv"
    cluster_df.to_csv(index_file, index=False)
    print(f"Saved cluster index: {index_file}")

    # Print statistics
    print("\nCluster Statistics:")
    print(f"  Total clusters: {len(cluster_ids):,}")
    print(f"  Embedding dimension: {mean_embeddings.shape[1]}")
    print(f"  Mean cluster size: {cluster_df['cluster_size'].mean():.1f}")
    print(f"  Median cluster size: {cluster_df['cluster_size'].median():.0f}")
    print(f"  Min cluster size: {cluster_df['cluster_size'].min()}")
    print(f"  Max cluster size: {cluster_df['cluster_size'].max()}")


def main():
    parser = argparse.ArgumentParser(
        description='Create full cluster embeddings (1152D) using direct approach'
    )
    parser.add_argument(
        '--cluster-stats',
        type=str,
        default='1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv',
        help='CSV file with used cluster IDs (388k clusters)'
    )
    parser.add_argument(
        '--clusters',
        type=str,
        default='1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv',
        help='MMseqs2 cluster TSV file'
    )
    parser.add_argument(
        '--embeddings',
        type=str,
        default='data/esm_embeddings',
        help='Directory with ESM embedding .npz files'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph',
        help='Output directory'
    )

    args = parser.parse_args()

    print("="*70)
    print("Creating Full Cluster-Level ESM Embeddings (1152D)")
    print("="*70)

    # Load the set of clusters to use (388k)
    used_clusters = load_used_clusters(args.cluster_stats)

    # Load cluster assignments (filtered to used clusters only)
    protein_to_cluster, cluster_sizes = load_cluster_assignments(
        args.clusters,
        used_clusters
    )

    # Load embeddings for proteins in used clusters
    cluster_embeddings = load_embeddings_for_clusters(
        args.embeddings,
        protein_to_cluster
    )

    # Compute mean embeddings
    mean_embeddings, cluster_ids = compute_cluster_means(cluster_embeddings)

    # Save outputs
    save_outputs(mean_embeddings, cluster_ids, cluster_sizes, args.output)

    print("\n" + "="*70)
    print("DONE!")
    print("="*70)


if __name__ == "__main__":
    main()
