#!/usr/bin/env python
"""
Merge batched embeddings into a single PCA cache file.

Processes all 11.8M newly generated embeddings:
1. Loads all batch files
2. Applies PCA to reduce from 1152D to 50D
3. Saves final PCA cache for downstream analysis
"""

import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from sklearn.decomposition import IncrementalPCA


def load_existing_embeddings():
    """Deprecated - no longer merging with old cache."""
    # Old testing cache has been moved to testing/ directory
    # This function is kept for reference but not used
    pass


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

    # Load batch embeddings
    batch_dir = Path('/fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches')
    gene_ids, embeddings = load_batch_embeddings(batch_dir)

    # Apply PCA to reduce dimensions
    embeddings_pca, pca = apply_pca_to_new_embeddings(embeddings, n_components=50)

    # Parse gene IDs into genome and gene components
    print("\nParsing gene IDs...")
    genome_ids = []
    gene_ids_short = []

    for gene_id in tqdm(gene_ids, desc="Parsing"):
        parts = gene_id.split('_', 1)
        if len(parts) == 2:
            genome_ids.append(parts[0])
            gene_ids_short.append(parts[1])
        else:
            # Fallback
            genome_ids.append('')
            gene_ids_short.append(gene_id)

    # Save PCA cache
    output_file = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz')
    print(f"\nSaving PCA cache to {output_file}...")

    np.savez_compressed(
        output_file,
        gene_ids=np.array(gene_ids_short),
        genome_ids=np.array(genome_ids),
        embeddings_pca=embeddings_pca,
        explained_variance_ratio=pca.explained_variance_ratio_
    )

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print(f"PCA cache saved to: {output_file}")
    print(f"  Total proteins: {len(gene_ids):,}")
    print(f"  Dimensions: {embeddings_pca.shape[1]}")
    print(f"  File size: {output_file.stat().st_size / 1e9:.2f} GB")
    print()


if __name__ == '__main__':
    main()
