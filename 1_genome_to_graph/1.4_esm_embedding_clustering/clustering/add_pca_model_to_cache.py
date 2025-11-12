#!/usr/bin/env python
"""
Add PCA model parameters to existing PCA cache.

The existing cache has embeddings_pca but not the PCA model itself.
We need to recompute PCA to extract the components and mean.

Usage:
    python add_pca_model_to_cache.py \
        --input results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache_with_model.npz \
        --n-components 50
"""

import numpy as np
from sklearn.decomposition import PCA
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Add PCA model to cache')
    parser.add_argument('--input', type=str, required=True,
                        help='Input PCA cache (without model)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output PCA cache (with model)')
    parser.add_argument('--n-components', type=int, default=50,
                        help='Number of PCA components')
    parser.add_argument('--embedding-dir', type=str,
                        default='data/refseq_esm_embeddings',
                        help='Directory with original embeddings (if needed to recompute)')

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    print(f"Loading PCA cache from {input_file}...")
    data = np.load(input_file, allow_pickle=True)

    # Check what we have
    print(f"Keys in cache: {list(data.keys())}")

    if 'pca_components' in data and 'pca_mean' in data:
        print("PCA model already present!")
        return

    # Need to recompute PCA from embeddings
    embeddings_pca = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"\nRecomputing PCA model from {len(gene_ids):,} genes...")
    print(f"  Embedding dimension: {embeddings_pca.shape[1]}")

    # We need original embeddings to compute PCA properly
    # Load a sample to get original dimension
    embedding_dir = Path(args.embedding_dir)
    emb_files = list(embedding_dir.glob('*_embeddings.npz'))

    if len(emb_files) == 0:
        raise ValueError(f"No embeddings found in {embedding_dir}")

    # Get original dimension
    sample_data = np.load(emb_files[0], allow_pickle=True)
    orig_dim = sample_data['embeddings'].shape[1]
    print(f"  Original dimension: {orig_dim}")

    # Collect original embeddings for genes in cache
    print("  Loading original embeddings for PCA model...")
    gene_id_set = set(gene_ids)
    embeddings_orig = []
    collected_gene_ids = []

    from tqdm import tqdm
    for emb_file in tqdm(emb_files, desc="Loading embeddings"):
        try:
            data_file = np.load(emb_file, allow_pickle=True)
            file_embeddings = data_file['embeddings']
            file_gene_ids = data_file['gene_ids']

            for i, gid in enumerate(file_gene_ids):
                if gid in gene_id_set:
                    embeddings_orig.append(file_embeddings[i])
                    collected_gene_ids.append(gid)
        except Exception as e:
            print(f"\nError loading {emb_file.name}: {e}")
            continue

    embeddings_orig = np.array(embeddings_orig)
    print(f"\n  Collected {len(embeddings_orig):,} embeddings")

    # Reorder to match original gene_ids order
    gene_id_to_idx = {gid: i for i, gid in enumerate(collected_gene_ids)}
    reorder_indices = [gene_id_to_idx[gid] for gid in gene_ids if gid in gene_id_to_idx]
    embeddings_orig = embeddings_orig[reorder_indices]

    print(f"  Reordered to {len(embeddings_orig):,} embeddings")

    # Fit PCA
    print(f"\n  Fitting PCA with {args.n_components} components...")
    pca = PCA(n_components=args.n_components, random_state=42)
    embeddings_pca_new = pca.fit_transform(embeddings_orig)

    print(f"  Explained variance: {pca.explained_variance_ratio_.sum():.4f}")
    print(f"  PCA components shape: {pca.components_.shape}")
    print(f"  PCA mean shape: {pca.mean_.shape}")

    # Verify PCA is consistent
    diff = np.abs(embeddings_pca_new - embeddings_pca).max()
    print(f"  Max difference from cached PCA: {diff:.6f}")

    if diff > 1e-3:
        print("  WARNING: PCA differs from cache! This may indicate different data.")

    # Save updated cache
    print(f"\nSaving updated cache to {output_file}...")
    np.savez_compressed(
        output_file,
        embeddings_pca=embeddings_pca,
        gene_ids=gene_ids,
        genome_ids=genome_ids,
        pca_components=pca.components_,
        pca_mean=pca.mean_,
        explained_variance_ratio=pca.explained_variance_ratio_,
        n_components=args.n_components
    )

    print("Done!")


if __name__ == '__main__':
    main()
