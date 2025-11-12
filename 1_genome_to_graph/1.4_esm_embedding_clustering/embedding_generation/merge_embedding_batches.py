#!/usr/bin/env python
"""
Merge batched embeddings into a single file and combine with existing PCA cache.

Combines:
1. Existing embeddings from PCA cache (1M proteins, 50D)
2. Newly generated embeddings (11.8M proteins, original dimensions)

Then applies PCA to new embeddings to match 50D, and creates final combined cache.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from sklearn.decomposition import IncrementalPCA


def load_existing_embeddings():
    """Load existing embeddings from PCA cache."""
    print("Loading existing embeddings from PCA cache...")
    cache = np.load('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz', allow_pickle=True)

    gene_ids_short = cache['gene_ids']
    genome_ids = cache['genome_ids']
    embeddings_pca = cache['embeddings_pca']

    # Reconstruct full gene IDs
    gene_ids_full = [f"{genome}_{gene}" for genome, gene in zip(genome_ids, gene_ids_short)]

    print(f"  Loaded {len(gene_ids_full):,} proteins")
    print(f"  Embedding dimensions: {embeddings_pca.shape[1]}")

    return gene_ids_full, embeddings_pca


def load_batch_embeddings(batch_dir):
    """Load all batch embeddings."""
    print(f"\nLoading batch embeddings from {batch_dir}...")

    batch_files = sorted(Path(batch_dir).glob('embeddings_batch_*.npz'))
    print(f"  Found {len(batch_files)} batch files")

    if len(batch_files) == 0:
        raise ValueError(f"No batch files found in {batch_dir}")

    all_gene_ids = []
    all_embeddings = []

    for batch_file in tqdm(batch_files, desc="Loading batches"):
        data = np.load(batch_file, allow_pickle=True)
        all_gene_ids.extend(data['gene_ids'])
        all_embeddings.append(data['embeddings'])

    all_embeddings = np.vstack(all_embeddings)

    print(f"  Loaded {len(all_gene_ids):,} proteins")
    print(f"  Embedding dimensions: {all_embeddings.shape[1]}")

    return all_gene_ids, all_embeddings


def apply_pca_to_new_embeddings(embeddings, n_components=50):
    """Apply PCA to new embeddings to reduce to 50D."""
    print(f"\nApplying PCA to reduce from {embeddings.shape[1]}D to {n_components}D...")

    # Use Incremental PCA for large datasets
    pca = IncrementalPCA(n_components=n_components)

    batch_size = 10000
    for i in tqdm(range(0, len(embeddings), batch_size), desc="Fitting PCA"):
        batch = embeddings[i:i + batch_size]
        pca.partial_fit(batch)

    # Transform all embeddings
    embeddings_pca = pca.transform(embeddings)

    print(f"  Explained variance ratio: {pca.explained_variance_ratio_.sum():.4f}")
    print(f"  Output shape: {embeddings_pca.shape}")

    return embeddings_pca, pca


def main():
    print("=" * 80)
    print("MERGE EMBEDDING BATCHES")
    print("=" * 80)
    print()

    # Load existing embeddings
    existing_gene_ids, existing_embeddings = load_existing_embeddings()

    # Load batch embeddings
    batch_dir = Path('/fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches')
    new_gene_ids, new_embeddings = load_batch_embeddings(batch_dir)

    # Apply PCA to new embeddings
    new_embeddings_pca, pca = apply_pca_to_new_embeddings(new_embeddings, n_components=50)

    # Combine with existing embeddings
    print("\nCombining embeddings...")
    combined_gene_ids = existing_gene_ids + new_gene_ids
    combined_embeddings = np.vstack([existing_embeddings, new_embeddings_pca])

    print(f"  Total proteins: {len(combined_gene_ids):,}")
    print(f"  Combined shape: {combined_embeddings.shape}")

    # Parse gene IDs back into genome and gene components
    print("\nParsing gene IDs...")
    genome_ids = []
    gene_ids_short = []

    for gene_id in tqdm(combined_gene_ids, desc="Parsing"):
        parts = gene_id.split('_', 1)
        if len(parts) == 2:
            genome_ids.append(parts[0])
            gene_ids_short.append(parts[1])
        else:
            # Fallback
            genome_ids.append('')
            gene_ids_short.append(gene_id)

    # Save combined cache
    output_file = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache_full.npz')
    print(f"\nSaving combined embeddings to {output_file}...")

    np.savez_compressed(
        output_file,
        gene_ids=gene_ids_short,
        genome_ids=genome_ids,
        embeddings_pca=combined_embeddings,
        explained_variance_ratio=pca.explained_variance_ratio_
    )

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print(f"Combined embeddings saved to: {output_file}")
    print(f"  Total proteins: {len(combined_gene_ids):,}")
    print(f"  Dimensions: {combined_embeddings.shape[1]}")
    print(f"  File size: {output_file.stat().st_size / 1e9:.2f} GB")
    print()


if __name__ == '__main__':
    main()
