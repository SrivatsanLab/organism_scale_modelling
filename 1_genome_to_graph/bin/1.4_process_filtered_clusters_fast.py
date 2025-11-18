#!/usr/bin/env python
"""
Fast version: Process filtered MMseqs2 clusters using only genes with embeddings.

Strategy:
1. Load the 1M genes with embeddings from PCA cache
2. Filter 70% clusters to those with 10+ members
3. For each cluster, keep only members that have embeddings
4. Refilter to clusters with at least 10 members WITH embeddings
5. Compute stats and sample for UMAP
"""

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm


def main():
    print("=" * 80)
    print("FAST PROCESSING OF FILTERED MMseqs2 CLUSTERS")
    print("=" * 80)
    print()

    # Paths
    mmseqs_dir = Path('data/mmseqs_seqid_0p7')
    pca_cache = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz')
    output_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load PCA cache
    print(f"Loading ESM embeddings from {pca_cache}")
    data = np.load(pca_cache, allow_pickle=True)
    cache_gene_ids_short = data['gene_ids']
    cache_genome_ids = data['genome_ids']
    cache_embeddings = data['embeddings_pca']

    print(f"  Genes in cache: {len(cache_gene_ids_short):,}")
    print(f"  Embedding dimensions: {cache_embeddings.shape[1]}")

    # Reconstruct full gene IDs (genome_contig_gene format to match MMseqs)
    cache_gene_ids_full = [f"{genome}_{gene}" for genome, gene in zip(cache_genome_ids, cache_gene_ids_short)]

    print(f"  Sample reconstructed IDs:")
    for i in range(3):
        print(f"    {cache_gene_ids_full[i]}")

    # Create gene -> embedding mapping
    gene_to_idx = {gene_id: idx for idx, gene_id in enumerate(cache_gene_ids_full)}
    cache_genes_set = set(cache_gene_ids_full)

    # Load MMseqs2 cluster assignments
    print(f"\nLoading MMseqs2 cluster assignments...")
    assignments_df = pd.read_csv(mmseqs_dir / 'cluster_summary_assignments.csv')

    print(f"  Total assignments: {len(assignments_df):,}")

    # Filter to genes with embeddings
    assignments_df = assignments_df[assignments_df['member'].isin(cache_genes_set)].copy()

    print(f"  After filtering to genes with embeddings: {len(assignments_df):,}")
    print(f"  Coverage: {len(assignments_df) / 30098843 * 100:.1f}% of all proteins")

    # Count members per cluster
    cluster_sizes = assignments_df.groupby('representative').size().reset_index(name='size')

    # Filter to clusters with 10+ members
    min_size = 10
    cluster_sizes_filtered = cluster_sizes[cluster_sizes['size'] >= min_size]

    print(f"\nClusters with {min_size}+ members (with embeddings):")
    print(f"  Total clusters: {len(cluster_sizes_filtered):,}")
    print(f"  Total proteins: {cluster_sizes_filtered['size'].sum():,}")
    print(f"  Mean cluster size: {cluster_sizes_filtered['size'].mean():.1f}")
    print(f"  Median cluster size: {cluster_sizes_filtered['size'].median():.0f}")
    print(f"  Max cluster size: {cluster_sizes_filtered['size'].max():,}")

    # Save filtered cluster summary
    cluster_sizes_filtered.to_csv(output_dir / 'filtered_cluster_summary.csv', index=False)

    valid_clusters = set(cluster_sizes_filtered['representative'])
    assignments_filtered = assignments_df[assignments_df['representative'].isin(valid_clusters)]
    assignments_filtered.to_csv(output_dir / 'filtered_cluster_assignments.csv', index=False)

    print(f"\nSaved filtered cluster data to {output_dir}")

    # Compute cluster statistics
    print(f"\nComputing cluster statistics...")

    cluster_stats = []

    for cluster_rep in tqdm(cluster_sizes_filtered['representative'], desc="Processing clusters"):
        # Get members
        members = assignments_filtered[
            assignments_filtered['representative'] == cluster_rep
        ]['member'].values

        # Get embeddings
        indices = [gene_to_idx[gene] for gene in members]
        cluster_embeddings = cache_embeddings[indices]

        # Compute stats
        mean_emb = cluster_embeddings.mean(axis=0)
        std_emb = cluster_embeddings.std(axis=0)

        cluster_stats.append({
            'cluster_rep': cluster_rep,
            'n_proteins': len(members),
            'mean_emb': mean_emb,
            'std_emb': std_emb
        })

    print(f"  Computed statistics for {len(cluster_stats):,} clusters")

    # Save cluster statistics
    print(f"\nSaving cluster statistics...")
    np.savez(
        output_dir / 'cluster_statistics.npz',
        cluster_reps=[s['cluster_rep'] for s in cluster_stats],
        n_proteins=[s['n_proteins'] for s in cluster_stats],
        mean_embeddings=np.array([s['mean_emb'] for s in cluster_stats]),
        std_embeddings=np.array([s['std_emb'] for s in cluster_stats]),
        embedding_dim=cache_embeddings.shape[1]
    )

    # Sample proteins for UMAP
    print(f"\nSampling 10 proteins per cluster for UMAP...")

    np.random.seed(42)
    sampled_data = []

    for cluster_rep in tqdm(cluster_sizes_filtered['representative'], desc="Sampling"):
        members = assignments_filtered[
            assignments_filtered['representative'] == cluster_rep
        ]['member'].values

        # Sample up to 10
        if len(members) <= 10:
            sampled = members
        else:
            sampled = np.random.choice(members, size=10, replace=False)

        # Get embeddings
        for gene in sampled:
            sampled_data.append({
                'gene_id': gene,
                'cluster_rep': cluster_rep,
                'embedding': cache_embeddings[gene_to_idx[gene]]
            })

    print(f"  Sampled {len(sampled_data):,} proteins")

    # Save sampled data
    print(f"\nSaving sampled data...")
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
    print(f"  1. filtered_cluster_summary.csv - {len(cluster_sizes_filtered):,} clusters")
    print(f"  2. filtered_cluster_assignments.csv - {len(assignments_filtered):,} assignments")
    print(f"  3. cluster_statistics.npz - Mean/SD for {len(cluster_stats):,} clusters")
    print(f"  4. umap_sample.npz - {len(sampled_data):,} proteins for UMAP")
    print()


if __name__ == '__main__':
    main()
