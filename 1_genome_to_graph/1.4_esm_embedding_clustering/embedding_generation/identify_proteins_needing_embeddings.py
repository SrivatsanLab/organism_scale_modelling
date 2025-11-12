#!/usr/bin/env python
"""
Identify proteins in filtered clusters that need ESM embeddings generated.

For 70% seq ID clusters with 10+ members, find all proteins that don't
already have embeddings in the PCA cache.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm


def main():
    print("=" * 80)
    print("IDENTIFY PROTEINS NEEDING EMBEDDINGS")
    print("=" * 80)
    print()

    # Load MMseqs2 cluster summary
    print("Loading MMseqs2 cluster data...")
    summary = pd.read_csv('data/mmseqs_seqid_0p7/cluster_summary.csv')
    print(f"  Total clusters: {len(summary):,}")

    # Filter to clusters with 10+ members
    large_clusters = summary[summary['size'] >= 10]
    print(f"  Clusters with 10+ members: {len(large_clusters):,}")
    print(f"  Total proteins: {large_clusters['size'].sum():,}")

    # Load cluster assignments
    print("\nLoading cluster assignments...")
    assignments = pd.read_csv('data/mmseqs_seqid_0p7/cluster_summary_assignments.csv')
    print(f"  Total assignments: {len(assignments):,}")

    # Filter to large clusters
    valid_reps = set(large_clusters['representative'])
    large_cluster_assignments = assignments[assignments['representative'].isin(valid_reps)]
    print(f"  Proteins in large clusters: {len(large_cluster_assignments):,}")

    # Load existing embeddings
    print("\nLoading existing embeddings from PCA cache...")
    cache = np.load('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz', allow_pickle=True)
    cache_gene_ids_short = cache['gene_ids']
    cache_genome_ids = cache['genome_ids']

    # Reconstruct full gene IDs
    cache_gene_ids_full = [f"{genome}_{gene}" for genome, gene in zip(cache_genome_ids, cache_gene_ids_short)]
    cache_set = set(cache_gene_ids_full)

    print(f"  Genes with embeddings: {len(cache_set):,}")

    # Find proteins needing embeddings
    print("\nIdentifying proteins needing embeddings...")
    all_proteins = large_cluster_assignments['member'].unique()
    proteins_needing_embeddings = [p for p in tqdm(all_proteins) if p not in cache_set]

    print(f"\n  Proteins WITH embeddings: {len(all_proteins) - len(proteins_needing_embeddings):,}")
    print(f"  Proteins NEEDING embeddings: {len(proteins_needing_embeddings):,}")
    print(f"  Coverage: {(len(all_proteins) - len(proteins_needing_embeddings)) / len(all_proteins) * 100:.2f}%")

    # Save list of proteins needing embeddings
    output_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7')
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / 'proteins_needing_embeddings.txt'
    print(f"\nSaving list to {output_file}...")

    with open(output_file, 'w') as f:
        for protein in proteins_needing_embeddings:
            f.write(f"{protein}\n")

    print(f"  Saved {len(proteins_needing_embeddings):,} protein IDs")

    # Also save summary statistics
    summary_file = output_dir / 'embedding_coverage_summary.txt'
    print(f"\nSaving summary to {summary_file}...")

    with open(summary_file, 'w') as f:
        f.write("EMBEDDING COVERAGE SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"70% sequence identity clustering:\n")
        f.write(f"  Total clusters: {len(summary):,}\n")
        f.write(f"  Clusters with 10+ members: {len(large_clusters):,}\n")
        f.write(f"  Total proteins in large clusters: {len(all_proteins):,}\n\n")
        f.write(f"Embedding coverage:\n")
        f.write(f"  Proteins WITH embeddings: {len(all_proteins) - len(proteins_needing_embeddings):,}\n")
        f.write(f"  Proteins NEEDING embeddings: {len(proteins_needing_embeddings):,}\n")
        f.write(f"  Coverage: {(len(all_proteins) - len(proteins_needing_embeddings)) / len(all_proteins) * 100:.2f}%\n\n")
        f.write(f"Output files:\n")
        f.write(f"  Protein list: {output_file}\n")

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print(f"Next steps:")
    print(f"  1. Extract sequences for {len(proteins_needing_embeddings):,} proteins")
    print(f"  2. Generate ESM embeddings in batches")
    print(f"  3. Merge with existing embeddings")
    print()


if __name__ == '__main__':
    main()
