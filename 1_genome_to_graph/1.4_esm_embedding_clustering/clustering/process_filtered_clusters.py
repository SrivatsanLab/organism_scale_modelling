#!/usr/bin/env python
"""
Process filtered MMseqs2 clusters (70% identity, 10+ members).

1. Extract protein lists for filtered clusters
2. Get ESM embeddings for these proteins
3. Compute cluster statistics (mean, SD per dimension)
4. Sample 10 proteins per cluster for UMAP
"""

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import argparse


def filter_clusters(cluster_summary_path, cluster_assignments_path, min_size=10):
    """Filter clusters to those with minimum size."""
    print(f"Loading cluster data from {cluster_summary_path}")

    # Load cluster summary
    summary_df = pd.read_csv(cluster_summary_path)

    print(f"  Total clusters: {len(summary_df):,}")
    print(f"  Total proteins: {summary_df['size'].sum():,}")

    # Filter to minimum size
    filtered_summary = summary_df[summary_df['size'] >= min_size].copy()

    print(f"\nFiltered to clusters with {min_size}+ members:")
    print(f"  Clusters remaining: {len(filtered_summary):,}")
    print(f"  Proteins remaining: {filtered_summary['size'].sum():,}")
    print(f"  Mean cluster size: {filtered_summary['size'].mean():.1f}")
    print(f"  Median cluster size: {filtered_summary['size'].median():.0f}")

    # Get list of valid cluster representatives
    valid_clusters = set(filtered_summary['representative'].values)

    # Load cluster assignments and filter
    print(f"\nLoading cluster assignments from {cluster_assignments_path}")
    assignments_df = pd.read_csv(cluster_assignments_path)

    print(f"  Total assignments: {len(assignments_df):,}")

    # Filter to valid clusters
    filtered_assignments = assignments_df[
        assignments_df['representative'].isin(valid_clusters)
    ].copy()

    print(f"  Assignments after filtering: {len(filtered_assignments):,}")

    return filtered_summary, filtered_assignments, valid_clusters


def get_esm_embeddings_for_proteins(gene_ids, pca_cache_path):
    """Load ESM embeddings (original 1152D) for specified proteins."""
    print(f"\nLoading ESM embeddings from {pca_cache_path}")

    data = np.load(pca_cache_path, allow_pickle=True)

    # Check what's in the cache
    print(f"  Cache contains: {list(data.keys())}")

    # The PCA cache has reduced dimensions, we need original embeddings
    # Let's check if we have access to original embeddings
    if 'embeddings_pca' in data:
        print("  WARNING: Only PCA embeddings (50D) available in cache")
        print("  Using PCA embeddings instead of original 1152D")
        all_embeddings = data['embeddings_pca']
        embedding_dim = all_embeddings.shape[1]
    else:
        raise ValueError("No embeddings found in cache")

    all_gene_ids = data['gene_ids']

    print(f"  Total genes in cache: {len(all_gene_ids):,}")
    print(f"  Embedding dimensions: {embedding_dim}")

    # Create mapping
    gene_to_idx = {gene_id: idx for idx, gene_id in enumerate(all_gene_ids)}

    # Find indices for our genes
    gene_ids_set = set(gene_ids)
    found_genes = []
    found_indices = []

    for gene_id in tqdm(gene_ids, desc="Finding gene embeddings"):
        if gene_id in gene_to_idx:
            found_genes.append(gene_id)
            found_indices.append(gene_to_idx[gene_id])

    print(f"\n  Found embeddings for {len(found_genes):,} / {len(gene_ids):,} genes")
    print(f"  Coverage: {len(found_genes) / len(gene_ids) * 100:.1f}%")

    # Extract embeddings
    embeddings = all_embeddings[found_indices]

    return np.array(found_genes), embeddings, embedding_dim


def compute_cluster_statistics(filtered_summary, filtered_assignments,
                               gene_ids, embeddings):
    """Compute mean and SD for each cluster."""
    print("\nComputing cluster statistics...")

    # Create gene -> embedding mapping
    gene_to_embedding = {gene_id: emb for gene_id, emb in zip(gene_ids, embeddings)}

    embedding_dim = embeddings.shape[1]

    cluster_stats = []

    for _, row in tqdm(filtered_summary.iterrows(),
                       total=len(filtered_summary),
                       desc="Processing clusters"):
        cluster_rep = row['representative']

        # Get all genes in this cluster
        cluster_genes = filtered_assignments[
            filtered_assignments['representative'] == cluster_rep
        ]['member'].values

        # Get embeddings for these genes
        cluster_embeddings = []
        for gene in cluster_genes:
            if gene in gene_to_embedding:
                cluster_embeddings.append(gene_to_embedding[gene])

        if len(cluster_embeddings) == 0:
            continue

        cluster_embeddings = np.array(cluster_embeddings)

        # Compute statistics
        mean_embedding = cluster_embeddings.mean(axis=0)
        std_embedding = cluster_embeddings.std(axis=0)

        cluster_stats.append({
            'cluster_rep': cluster_rep,
            'n_proteins': len(cluster_embeddings),
            'mean_embedding': mean_embedding,
            'std_embedding': std_embedding
        })

    print(f"  Computed statistics for {len(cluster_stats):,} clusters")

    return cluster_stats


def sample_proteins_for_umap(filtered_summary, filtered_assignments,
                             gene_ids, embeddings, n_per_cluster=10, seed=42):
    """Sample up to n proteins per cluster for UMAP visualization."""
    print(f"\nSampling {n_per_cluster} proteins per cluster...")

    np.random.seed(seed)

    # Create gene -> embedding mapping
    gene_to_embedding = {gene_id: emb for gene_id, emb in zip(gene_ids, embeddings)}

    sampled_data = []

    for _, row in tqdm(filtered_summary.iterrows(),
                       total=len(filtered_summary),
                       desc="Sampling clusters"):
        cluster_rep = row['representative']

        # Get all genes in this cluster
        cluster_genes = filtered_assignments[
            filtered_assignments['representative'] == cluster_rep
        ]['member'].values

        # Filter to genes with embeddings
        cluster_genes_with_emb = [g for g in cluster_genes if g in gene_to_embedding]

        if len(cluster_genes_with_emb) == 0:
            continue

        # Sample (with replacement if needed)
        if len(cluster_genes_with_emb) < n_per_cluster:
            sampled_genes = cluster_genes_with_emb
        else:
            sampled_genes = np.random.choice(cluster_genes_with_emb,
                                            size=n_per_cluster,
                                            replace=False)

        # Store sampled data
        for gene in sampled_genes:
            sampled_data.append({
                'gene_id': gene,
                'cluster_rep': cluster_rep,
                'embedding': gene_to_embedding[gene]
            })

    print(f"  Sampled {len(sampled_data):,} proteins total")
    print(f"  From {len(filtered_summary):,} clusters")
    print(f"  Average: {len(sampled_data) / len(filtered_summary):.1f} proteins/cluster")

    return sampled_data


def main():
    parser = argparse.ArgumentParser(
        description='Process filtered MMseqs2 clusters for UMAP analysis'
    )
    parser.add_argument('--mmseqs-dir', type=str, default='data/mmseqs_seqid_0p7',
                       help='MMseqs2 clustering directory')
    parser.add_argument('--pca-cache', type=str, default='results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz',
                       help='PCA cache with embeddings')
    parser.add_argument('--min-size', type=int, default=10,
                       help='Minimum cluster size')
    parser.add_argument('--n-sample', type=int, default=10,
                       help='Number of proteins to sample per cluster')
    parser.add_argument('--output-dir', type=str, default='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7',
                       help='Output directory')

    args = parser.parse_args()

    print("=" * 80)
    print("PROCESSING FILTERED MMseqs2 CLUSTERS")
    print("=" * 80)
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Filter clusters
    cluster_summary = Path(args.mmseqs_dir) / 'cluster_summary.csv'
    cluster_assignments = Path(args.mmseqs_dir) / 'cluster_summary_assignments.csv'

    filtered_summary, filtered_assignments, valid_clusters = filter_clusters(
        cluster_summary, cluster_assignments, min_size=args.min_size
    )

    # Save filtered cluster info
    filtered_summary.to_csv(output_dir / 'filtered_cluster_summary.csv', index=False)
    filtered_assignments.to_csv(output_dir / 'filtered_cluster_assignments.csv', index=False)

    print(f"\nSaved filtered cluster data to {output_dir}")

    # Step 2: Get ESM embeddings for all proteins in filtered clusters
    all_gene_ids = filtered_assignments['member'].unique()
    print(f"\nTotal unique genes in filtered clusters: {len(all_gene_ids):,}")

    gene_ids, embeddings, embedding_dim = get_esm_embeddings_for_proteins(
        all_gene_ids, args.pca_cache
    )

    # Step 3: Compute cluster statistics
    cluster_stats = compute_cluster_statistics(
        filtered_summary, filtered_assignments, gene_ids, embeddings
    )

    # Save cluster statistics
    print(f"\nSaving cluster statistics to {output_dir / 'cluster_statistics.npz'}")
    np.savez(
        output_dir / 'cluster_statistics.npz',
        cluster_reps=[s['cluster_rep'] for s in cluster_stats],
        n_proteins=[s['n_proteins'] for s in cluster_stats],
        mean_embeddings=np.array([s['mean_embedding'] for s in cluster_stats]),
        std_embeddings=np.array([s['std_embedding'] for s in cluster_stats]),
        embedding_dim=embedding_dim
    )

    # Step 4: Sample proteins for UMAP
    sampled_data = sample_proteins_for_umap(
        filtered_summary, filtered_assignments, gene_ids, embeddings,
        n_per_cluster=args.n_sample
    )

    # Save sampled data
    print(f"\nSaving sampled data to {output_dir / 'umap_sample.npz'}")
    np.savez(
        output_dir / 'umap_sample.npz',
        gene_ids=[s['gene_id'] for s in sampled_data],
        cluster_reps=[s['cluster_rep'] for s in sampled_data],
        embeddings=np.array([s['embedding'] for s in sampled_data])
    )

    print()
    print("=" * 80)
    print("PROCESSING COMPLETE")
    print("=" * 80)
    print()
    print(f"Output files in {output_dir}:")
    print(f"  1. filtered_cluster_summary.csv - Cluster metadata")
    print(f"  2. filtered_cluster_assignments.csv - Gene -> cluster mapping")
    print(f"  3. cluster_statistics.npz - Mean/SD per cluster")
    print(f"  4. umap_sample.npz - Sampled proteins for UMAP")
    print()


if __name__ == '__main__':
    main()
