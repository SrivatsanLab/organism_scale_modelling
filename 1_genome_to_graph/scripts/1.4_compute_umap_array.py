#!/usr/bin/env python
"""
Compute UMAP embeddings for gene data with specified parameters.

This script is designed to be run as part of a SLURM array job,
with each job processing a different n_neighbors value.

Usage:
    python compute_umap_array.py --n-neighbors 15 --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap_n15.npz
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import time
from tqdm import tqdm
import gc

# Import ML libraries
from sklearn.decomposition import PCA
import umap


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute UMAP embeddings for gene data"
    )

    parser.add_argument(
        "--embeddings-dir",
        type=str,
        default="data/esm_embeddings",
        help="Directory containing embedding files",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file for UMAP results (.npz)",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        required=True,
        help="UMAP n_neighbors parameter",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of PCA components (default: 50)",
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.1,
        help="UMAP min_dist parameter (default: 0.1)",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=None,
        help="Subsample N genes (default: use all)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Use GPU acceleration (requires cuML)",
    )
    parser.add_argument(
        "--load-cached-pca",
        type=str,
        default=None,
        help="Load pre-computed PCA embeddings from file",
    )
    parser.add_argument(
        "--save-pca",
        type=str,
        default=None,
        help="Save PCA embeddings to file",
    )

    return parser.parse_args()


def load_all_embeddings(embeddings_dir, subsample=None, seed=42):
    """
    Load all gene embeddings from directory.

    Returns:
        embeddings: np.ndarray of shape (n_genes, embedding_dim)
        gene_ids: np.ndarray of gene IDs
        genome_ids: np.ndarray of genome IDs
    """
    embeddings_dir = Path(embeddings_dir)
    embedding_files = sorted(embeddings_dir.glob("*_esmc_embeddings_esmc_600m.npz"))

    print(f"Loading embeddings from {len(embedding_files)} genomes...")

    all_embeddings = []
    all_gene_ids = []
    all_genome_ids = []

    for emb_file in tqdm(embedding_files, desc="Loading"):
        genome_id = emb_file.stem.split("_esmc_embeddings")[0]

        try:
            data = np.load(emb_file, allow_pickle=True)
            embeddings = data['embeddings']
            seq_ids = data['seq_ids']

            all_embeddings.append(embeddings)
            all_gene_ids.extend(seq_ids)
            all_genome_ids.extend([genome_id] * len(seq_ids))

        except Exception as e:
            print(f"Error loading {genome_id}: {e}", file=sys.stderr)
            continue

    # Concatenate
    print("Concatenating embeddings...")
    all_embeddings = np.vstack(all_embeddings)
    all_gene_ids = np.array(all_gene_ids)
    all_genome_ids = np.array(all_genome_ids)

    print(f"Loaded {len(all_embeddings):,} genes from {len(embedding_files)} genomes")
    print(f"Embedding shape: {all_embeddings.shape}")
    print(f"Memory: {all_embeddings.nbytes / 1e9:.2f} GB")

    # Subsample if requested
    if subsample is not None and subsample < len(all_embeddings):
        print(f"\nSubsampling {subsample:,} genes...")
        np.random.seed(seed)
        indices = np.random.choice(len(all_embeddings), subsample, replace=False)

        all_embeddings = all_embeddings[indices]
        all_gene_ids = all_gene_ids[indices]
        all_genome_ids = all_genome_ids[indices]

        print(f"Subset shape: {all_embeddings.shape}")

    gc.collect()

    return all_embeddings, all_gene_ids, all_genome_ids


def compute_pca(embeddings, n_components=50):
    """
    Compute PCA on embeddings.

    Returns:
        embeddings_pca: PCA-transformed embeddings
        pca: Fitted PCA object
    """
    print(f"\nRunning PCA: {embeddings.shape[1]}D → {n_components}D")

    pca = PCA(n_components=n_components)
    embeddings_pca = pca.fit_transform(embeddings)

    variance_explained = np.sum(pca.explained_variance_ratio_)
    print(f"Variance explained: {variance_explained:.2%}")
    print(f"PCA shape: {embeddings_pca.shape}")

    return embeddings_pca, pca


def compute_umap(embeddings, n_neighbors, min_dist=0.1, seed=42, use_gpu=False):
    """
    Compute UMAP embedding.

    Returns:
        umap_embedding: 2D UMAP coordinates
        umap_reducer: Fitted UMAP object
    """
    print(f"\nRunning UMAP...")
    print(f"  n_neighbors: {n_neighbors}")
    print(f"  min_dist: {min_dist}")
    print(f"  metric: cosine")
    print(f"  GPU: {use_gpu}")

    if use_gpu:
        try:
            from cuml import UMAP as cumlUMAP
            import cupy as cp

            # Test GPU
            cp.cuda.Device(0).compute_capability

            print("  Using cuML GPU UMAP")

            # Convert to GPU array
            embeddings_gpu = cp.asarray(embeddings)

            # Run UMAP
            umap_reducer = cumlUMAP(
                n_neighbors=n_neighbors,
                n_components=2,
                min_dist=min_dist,
                metric='cosine',
                random_state=seed,
                verbose=True
            )

            umap_embedding_gpu = umap_reducer.fit_transform(embeddings_gpu)

            # Convert back to CPU
            umap_embedding = cp.asnumpy(umap_embedding_gpu)

        except Exception as e:
            print(f"  GPU UMAP failed: {e}", file=sys.stderr)
            print("  Falling back to CPU UMAP")
            use_gpu = False

    if not use_gpu:
        print("  Using CPU UMAP")

        umap_reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            n_components=2,
            min_dist=min_dist,
            metric='cosine',
            random_state=seed,
            verbose=True
        )

        umap_embedding = umap_reducer.fit_transform(embeddings)

    print(f"UMAP complete! Shape: {umap_embedding.shape}")

    return umap_embedding, umap_reducer


def main():
    args = parse_args()

    print("=" * 70)
    print("UMAP Computation")
    print("=" * 70)
    print(f"n_neighbors: {args.n_neighbors}")
    print(f"n_pcs: {args.n_pcs}")
    print(f"min_dist: {args.min_dist}")
    print(f"Output: {args.output}")
    print("=" * 70)

    start_time = time.time()

    # Step 1: Load or compute PCA embeddings
    if args.load_cached_pca:
        print(f"\nLoading cached PCA from {args.load_cached_pca}")
        cache = np.load(args.load_cached_pca, allow_pickle=True)
        embeddings_pca = cache['embeddings_pca']
        gene_ids = cache['gene_ids']
        genome_ids = cache['genome_ids']
        print(f"Loaded PCA embeddings: {embeddings_pca.shape}")
    else:
        # Load raw embeddings
        embeddings, gene_ids, genome_ids = load_all_embeddings(
            args.embeddings_dir,
            subsample=args.subsample,
            seed=args.seed
        )

        # Compute PCA
        embeddings_pca, pca = compute_pca(embeddings, n_components=args.n_pcs)

        # Save PCA if requested
        if args.save_pca:
            print(f"\nSaving PCA to {args.save_pca}")
            np.savez_compressed(
                args.save_pca,
                embeddings_pca=embeddings_pca,
                gene_ids=gene_ids,
                genome_ids=genome_ids,
                explained_variance_ratio=pca.explained_variance_ratio_
            )

        # Free memory
        del embeddings
        gc.collect()

    # Step 2: Compute UMAP
    umap_embedding, umap_reducer = compute_umap(
        embeddings_pca,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        seed=args.seed,
        use_gpu=args.use_gpu
    )

    # Step 3: Save results
    print(f"\nSaving results to {args.output}")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        args.output,
        umap_embedding=umap_embedding,
        gene_ids=gene_ids,
        genome_ids=genome_ids,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        n_pcs=args.n_pcs,
        seed=args.seed
    )

    elapsed = time.time() - start_time
    print(f"\n✓ Complete in {elapsed/60:.1f} minutes")
    print(f"Output saved to: {args.output}")


if __name__ == "__main__":
    main()
