#!/usr/bin/env python3
"""
Create mean ESM embeddings for each MMseqs2 protein cluster.

This script:
1. Loads MMseqs2 cluster assignments (~30M proteins → ~400k clusters)
2. Loads ESM embeddings from all genomes (7,664 files, 1152D embeddings)
3. Computes mean embedding per cluster
4. Saves cluster × embedding matrix in compressed format

Input:
  - MMseqs2 clusters: intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv
  - ESM embeddings: data/esm_embeddings/*.npz (7,664 files)

Output:
  - graph_outputs/protein_graph/cluster_embeddings.npz (cluster × 1152D matrix)
  - graph_outputs/protein_graph/cluster_index.csv (cluster IDs and metadata)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
import argparse


def load_cluster_assignments(cluster_file):
    """
    Load MMseqs2 cluster assignments.

    Returns:
        dict: protein_id -> cluster_representative
    """
    print(f"Loading cluster assignments from {cluster_file}...")

    protein_to_cluster = {}
    cluster_sizes = defaultdict(int)

    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Reading clusters"):
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_rep, protein_id = parts
                protein_to_cluster[protein_id] = cluster_rep
                cluster_sizes[cluster_rep] += 1

    n_proteins = len(protein_to_cluster)
    n_clusters = len(cluster_sizes)

    print(f"Loaded {n_proteins:,} proteins in {n_clusters:,} clusters")
    print(f"Mean cluster size: {n_proteins/n_clusters:.1f}")
    print(f"Median cluster size: {np.median(list(cluster_sizes.values())):.0f}")
    print(f"Max cluster size: {max(cluster_sizes.values()):,}")

    return protein_to_cluster, cluster_sizes


def load_all_embeddings(embedding_dir, protein_to_cluster):
    """
    Load all ESM embeddings and organize by cluster.

    Returns:
        dict: cluster_id -> list of embeddings (1152D)
    """
    embedding_dir = Path(embedding_dir)
    embedding_files = sorted(embedding_dir.glob("*.npz"))

    print(f"\nFound {len(embedding_files)} embedding files")

    cluster_embeddings = defaultdict(list)
    proteins_found = 0
    proteins_missing = 0

    for emb_file in tqdm(embedding_files, desc="Loading embeddings"):
        data = np.load(emb_file)
        embeddings = data['embeddings']  # Shape: (n_proteins, 1152)
        seq_ids = data['seq_ids']

        for protein_id, embedding in zip(seq_ids, embeddings):
            protein_id_str = str(protein_id)

            if protein_id_str in protein_to_cluster:
                cluster_id = protein_to_cluster[protein_id_str]
                cluster_embeddings[cluster_id].append(embedding)
                proteins_found += 1
            else:
                proteins_missing += 1

    print(f"\nProteins with embeddings: {proteins_found:,}")
    print(f"Proteins missing embeddings: {proteins_missing:,}")
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
    """Save cluster embeddings and metadata."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save embeddings matrix
    embedding_file = output_dir / 'cluster_embeddings.npz'
    print(f"\nSaving embeddings to {embedding_file}...")
    np.savez_compressed(
        embedding_file,
        embeddings=mean_embeddings,
        cluster_ids=cluster_ids
    )

    # Create cluster index with metadata
    cluster_data = []
    for cluster_id in cluster_ids:
        # Parse genome ID from cluster representative
        genome_id = '_'.join(cluster_id.split('_')[:2])  # e.g., GCF_000006985.1

        cluster_data.append({
            'cluster_id': cluster_id,
            'genome_id': genome_id,
            'cluster_size': cluster_sizes[cluster_id]
        })

    cluster_df = pd.DataFrame(cluster_data)

    # Save cluster index
    index_file = output_dir / 'cluster_index.csv'
    print(f"Saving cluster index to {index_file}...")
    cluster_df.to_csv(index_file, index=False)

    # Print summary stats
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total clusters: {len(cluster_ids):,}")
    print(f"Embedding dimension: {mean_embeddings.shape[1]}")
    print(f"Matrix shape: {mean_embeddings.shape}")
    print(f"Matrix size: {embedding_file.stat().st_size / 1024**2:.1f} MB")
    print(f"\nCluster size distribution:")
    print(f"  Mean: {cluster_df['cluster_size'].mean():.1f}")
    print(f"  Median: {cluster_df['cluster_size'].median():.0f}")
    print(f"  Min: {cluster_df['cluster_size'].min()}")
    print(f"  Max: {cluster_df['cluster_size'].max():,}")
    print("="*70)


def main():
    parser = argparse.ArgumentParser(
        description='Create mean ESM embeddings for MMseqs2 protein clusters'
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
        help='Directory containing ESM embedding NPZ files'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph',
        help='Output directory for cluster embeddings'
    )

    args = parser.parse_args()

    print("="*70)
    print("Creating Cluster-Level ESM Embeddings")
    print("="*70)

    # Load cluster assignments
    protein_to_cluster, cluster_sizes = load_cluster_assignments(args.clusters)

    # Load all embeddings organized by cluster
    cluster_embeddings = load_all_embeddings(args.embeddings, protein_to_cluster)

    # Compute mean embeddings
    mean_embeddings, cluster_ids = compute_cluster_means(cluster_embeddings)

    # Save outputs
    save_outputs(mean_embeddings, cluster_ids, cluster_sizes, args.output_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
