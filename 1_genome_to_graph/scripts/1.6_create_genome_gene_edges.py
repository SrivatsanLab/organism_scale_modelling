#!/usr/bin/env python3
"""
Create genome-gene edge matrix (presence/absence of genes in genomes).

This script creates a binary matrix indicating which protein clusters (genes)
are present in which genomes. This forms the bipartite edges between genome
nodes and gene nodes in the heterogeneous graph.

Input:
  - MMseqs2 cluster assignments (TSV)
  - Cluster statistics CSV (388k filtered clusters)
  - Genome metadata

Output:
  - graph_outputs/protein_graph/genome_gene_edges.csv - genome-gene edges
  - graph_outputs/protein_graph/genome_gene_matrix.npz - sparse binary matrix
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from collections import defaultdict
from tqdm import tqdm


def load_cluster_list(cluster_stats_file):
    """
    Load the set of cluster IDs to use.

    Handles both cluster_statistics.csv (cluster_representative column)
    and cluster_index.csv (cluster_id column) formats.

    Returns:
        set: Cluster IDs
    """
    print(f"Loading cluster list from {cluster_stats_file}")
    df = pd.read_csv(cluster_stats_file)

    # Handle both formats
    if 'cluster_representative' in df.columns:
        cluster_ids = set(df['cluster_representative'].values)
    elif 'cluster_id' in df.columns:
        cluster_ids = set(df['cluster_id'].values)
    else:
        raise ValueError(f"Cluster file must have 'cluster_representative' or 'cluster_id' column")

    print(f"Found {len(cluster_ids):,} clusters")

    return cluster_ids


def load_genome_list(genome_metadata_file):
    """
    Load the list of genomes.

    Returns:
        list: Genome IDs
    """
    print(f"\nLoading genome list from {genome_metadata_file}")
    df = pd.read_csv(genome_metadata_file)

    genome_ids = df['genome_id'].tolist()
    print(f"Found {len(genome_ids):,} genomes")

    return genome_ids


def build_genome_gene_matrix(cluster_file, used_clusters, genome_list):
    """
    Build genome-gene presence/absence matrix from cluster assignments.

    Args:
        cluster_file: Path to MMseqs2 clusters TSV
        used_clusters: Set of cluster IDs to include
        genome_list: List of genome IDs

    Returns:
        dict: genome_id -> set of cluster_ids present in that genome
        set: All cluster IDs that appear in the data
    """
    print(f"\nBuilding genome-gene matrix from {cluster_file}")

    # Map genome -> set of clusters present
    genome_to_clusters = defaultdict(set)

    # Track which clusters actually appear
    clusters_seen = set()

    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Reading clusters"):
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue

            cluster_rep, protein_id = parts

            # Only process clusters in our filtered set
            if cluster_rep not in used_clusters:
                continue

            # Extract genome ID from protein ID
            # Format: GCF_XXXXXX.X_contig_gene
            genome_id = '_'.join(protein_id.split('_')[:2])  # GCF_XXXXXX.X

            genome_to_clusters[genome_id].add(cluster_rep)
            clusters_seen.add(cluster_rep)

    print(f"\nGenomes with clusters: {len(genome_to_clusters):,}")
    print(f"Clusters found in data: {len(clusters_seen):,}")

    return genome_to_clusters, clusters_seen


def create_sparse_matrix(genome_to_clusters, genome_list, cluster_list):
    """
    Create sparse binary matrix of genome-gene presence/absence.

    Args:
        genome_to_clusters: Dict mapping genome_id -> set of cluster_ids
        genome_list: Ordered list of genome IDs
        cluster_list: Ordered list of cluster IDs

    Returns:
        csr_matrix: Sparse binary matrix (n_genomes x n_clusters)
        list: Genome IDs (row indices)
        list: Cluster IDs (column indices)
    """
    print("\nCreating sparse matrix...")

    # Create index mappings
    genome_to_idx = {g: i for i, g in enumerate(genome_list)}
    cluster_to_idx = {c: i for i, c in enumerate(cluster_list)}

    n_genomes = len(genome_list)
    n_clusters = len(cluster_list)

    print(f"Matrix dimensions: {n_genomes:,} genomes x {n_clusters:,} clusters")

    # Build sparse matrix
    row_indices = []
    col_indices = []

    for genome_id, clusters in tqdm(genome_to_clusters.items(), desc="Building matrix"):
        if genome_id not in genome_to_idx:
            continue

        genome_idx = genome_to_idx[genome_id]

        for cluster_id in clusters:
            if cluster_id not in cluster_to_idx:
                continue

            cluster_idx = cluster_to_idx[cluster_id]
            row_indices.append(genome_idx)
            col_indices.append(cluster_idx)

    # Create sparse matrix with binary values
    data = np.ones(len(row_indices), dtype=np.int8)
    matrix = csr_matrix(
        (data, (row_indices, col_indices)),
        shape=(n_genomes, n_clusters),
        dtype=np.int8
    )

    print(f"Total edges: {len(row_indices):,}")
    print(f"Matrix sparsity: {100 * (1 - matrix.nnz / (n_genomes * n_clusters)):.2f}%")

    return matrix, genome_list, cluster_list


def create_edge_list(genome_to_clusters, genome_list, cluster_list):
    """
    Create edge list format (genome_id, cluster_id pairs).

    Returns:
        DataFrame: columns [genome_id, cluster_id]
    """
    print("\nCreating edge list...")

    edges = []

    # Filter to only genomes and clusters in our lists
    genome_set = set(genome_list)
    cluster_set = set(cluster_list)

    for genome_id, clusters in tqdm(genome_to_clusters.items(), desc="Building edges"):
        if genome_id not in genome_set:
            continue

        for cluster_id in clusters:
            if cluster_id not in cluster_set:
                continue

            edges.append({
                'genome_id': genome_id,
                'cluster_id': cluster_id
            })

    df = pd.DataFrame(edges)
    print(f"Created {len(df):,} genome-gene edges")

    return df


def save_outputs(matrix, genome_list, cluster_list, edge_df, output_dir):
    """
    Save genome-gene matrix and edge list.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nSaving outputs to {output_dir}")

    # Save sparse matrix
    matrix_file = output_dir / "genome_gene_matrix.npz"
    np.savez_compressed(
        matrix_file,
        data=matrix.data,
        indices=matrix.indices,
        indptr=matrix.indptr,
        shape=matrix.shape,
        genome_ids=np.array(genome_list, dtype=object),
        cluster_ids=np.array(cluster_list, dtype=object)
    )
    print(f"Saved sparse matrix: {matrix_file}")
    print(f"  Size: {matrix_file.stat().st_size / 1e6:.1f} MB")

    # Save edge list
    edges_file = output_dir / "genome_gene_edges.csv"
    edge_df.to_csv(edges_file, index=False)
    print(f"Saved edge list: {edges_file}")
    print(f"  Size: {edges_file.stat().st_size / 1e6:.1f} MB")

    # Print statistics
    print("\nGenome-Gene Matrix Statistics:")
    print(f"  Genomes: {len(genome_list):,}")
    print(f"  Clusters: {len(cluster_list):,}")
    print(f"  Edges: {matrix.nnz:,}")
    print(f"  Sparsity: {100 * (1 - matrix.nnz / (len(genome_list) * len(cluster_list))):.2f}%")
    print(f"  Mean clusters per genome: {matrix.sum(axis=1).mean():.1f}")
    print(f"  Mean genomes per cluster: {matrix.sum(axis=0).mean():.1f}")


def main():
    parser = argparse.ArgumentParser(
        description='Create genome-gene edge matrix'
    )
    parser.add_argument(
        '--clusters',
        type=str,
        default='1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv',
        help='MMseqs2 cluster TSV file'
    )
    parser.add_argument(
        '--cluster-stats',
        type=str,
        default='1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv',
        help='Cluster statistics CSV (filtered 388k clusters)'
    )
    parser.add_argument(
        '--genome-list',
        type=str,
        default='data/genome_metadata.csv',
        help='Genome metadata CSV'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph',
        help='Output directory'
    )

    args = parser.parse_args()

    print("="*70)
    print("Creating Genome-Gene Edge Matrix")
    print("="*70)

    # Load cluster and genome lists
    used_clusters = load_cluster_list(args.cluster_stats)
    genome_list = load_genome_list(args.genome_list)

    # Build genome-gene associations
    genome_to_clusters, clusters_seen = build_genome_gene_matrix(
        args.clusters,
        used_clusters,
        genome_list
    )

    # Create ordered cluster list (only clusters that appear in data)
    cluster_list = sorted(clusters_seen)

    # Create sparse matrix
    matrix, genome_ids, cluster_ids = create_sparse_matrix(
        genome_to_clusters,
        genome_list,
        cluster_list
    )

    # Create edge list
    edge_df = create_edge_list(
        genome_to_clusters,
        genome_list,
        cluster_list
    )

    # Save outputs
    save_outputs(matrix, genome_ids, cluster_ids, edge_df, args.output)

    print("\n" + "="*70)
    print("DONE!")
    print("="*70)


if __name__ == '__main__':
    main()
