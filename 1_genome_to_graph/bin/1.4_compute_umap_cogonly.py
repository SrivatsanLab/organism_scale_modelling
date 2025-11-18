#!/usr/bin/env python
"""
Compute UMAP embeddings for COG-annotated genes only.

This creates UMAPs that match the gene set used in COG-filtered clusterings,
providing better visualizations since both clustering and UMAP use the same genes.
"""

import argparse
import numpy as np
from pathlib import Path
import time
import sys

# Add project root to path
sys.path.insert(0, '/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling')

from scripts.embeddings.plot_all_clusterings import load_cog_annotations_for_genes


def load_pca_and_filter_to_cog(pca_file):
    """Load PCA embeddings and filter to COG-annotated genes."""
    print(f"Loading PCA from {pca_file}...")
    data = np.load(pca_file, allow_pickle=True)

    pca_embedding = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Total genes: {len(gene_ids):,}")
    print()

    # Load COG annotations
    cog_lookup = load_cog_annotations_for_genes(gene_ids, genome_ids)
    print()

    # Filter to annotated genes
    print("Filtering to COG-annotated genes...")
    annotated_mask = np.array([gid in cog_lookup for gid in gene_ids])

    filtered_pca = pca_embedding[annotated_mask]
    filtered_gene_ids = gene_ids[annotated_mask]
    filtered_genome_ids = genome_ids[annotated_mask]

    print(f"  Kept: {len(filtered_gene_ids):,} genes ({100*len(filtered_gene_ids)/len(gene_ids):.1f}%)")
    print()

    return filtered_pca, filtered_gene_ids, filtered_genome_ids, cog_lookup


def compute_umap(embeddings, n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42):
    """Compute UMAP embedding."""
    import umap

    print(f"Computing UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})...")
    print(f"  Input shape: {embeddings.shape}")

    start = time.time()

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        n_components=2,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
        verbose=True
    )

    embedding_2d = reducer.fit_transform(embeddings)

    elapsed = time.time() - start
    print(f"  UMAP completed in {elapsed:.1f}s")
    print(f"  Output shape: {embedding_2d.shape}")

    return embedding_2d


def main():
    parser = argparse.ArgumentParser(description='Compute UMAP for COG-annotated genes')
    parser.add_argument('--pca-cache', required=True, help='PCA cache file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--n-neighbors', type=int, default=15, help='UMAP n_neighbors')
    parser.add_argument('--min-dist', type=float, default=0.1, help='UMAP min_dist')

    args = parser.parse_args()

    print("=" * 80)
    print("UMAP for COG-Annotated Genes Only")
    print("=" * 80)
    print(f"N neighbors: {args.n_neighbors}")
    print(f"Min dist: {args.min_dist}")
    print()

    # Load PCA and filter to COG-annotated genes
    pca_embedding, gene_ids, genome_ids, cog_lookup = load_pca_and_filter_to_cog(args.pca_cache)

    # Compute UMAP
    umap_embedding = compute_umap(
        pca_embedding,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist
    )
    print()

    # Save results
    print(f"Saving to {args.output}...")
    np.savez_compressed(
        args.output,
        umap_embedding=umap_embedding,
        gene_ids=gene_ids,
        genome_ids=genome_ids,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        cog_only=True
    )

    print()
    print("=" * 80)
    print("✓ UMAP complete!")
    print(f"Genes: {len(gene_ids):,}")
    print(f"Output: {args.output}")
    print("=" * 80)


if __name__ == '__main__':
    main()
