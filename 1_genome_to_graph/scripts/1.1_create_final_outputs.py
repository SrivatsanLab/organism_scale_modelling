#!/usr/bin/env python3
"""
Create final outputs for Component 1.1: K-mer profiling

Outputs:
1. Compressed 6-mer matrix (genome × k-mer) with genome and k-mer indices
2. K-nearest neighbor graph adjacency matrix optimized for genus-level structure
"""

import numpy as np
import pandas as pd
from scipy.sparse import load_npz, save_npz, csr_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
import argparse
from pathlib import Path
import re


def extract_taxonomy_from_ids(genome_ids):
    """Extract genus and species from genome IDs."""
    taxonomy = {}
    for genome_id in genome_ids:
        parts = genome_id.split('_')
        if len(parts) >= 4:
            genus = parts[2]
            species = parts[3] if len(parts) > 3 else 'unknown'
        else:
            genus = 'unknown'
            species = 'unknown'

        taxonomy[genome_id] = {
            'genus': genus,
            'species': species,
            'full_name': f"{genus}_{species}"
        }

    return taxonomy


def evaluate_graph_quality(adjacency, taxonomy_df, level='genus'):
    """
    Evaluate how well a graph captures taxonomic structure.

    Returns genus purity: fraction of neighbors sharing the same genus.
    """
    n_samples = adjacency.shape[0]
    purities = []

    for i in range(n_samples):
        # For CSR matrix, nonzero() returns (row_indices, col_indices)
        # We want column indices (the actual neighbors)
        neighbors = adjacency[i].nonzero()[1]
        if len(neighbors) == 0:
            continue

        my_taxon = taxonomy_df.iloc[i][level]
        if my_taxon == 'unknown':
            continue

        neighbor_taxa = taxonomy_df.iloc[neighbors][level].values
        purity = (neighbor_taxa == my_taxon).mean()
        purities.append(purity)

    return np.mean(purities) if purities else 0.0


def create_knn_graph(X, n_neighbors, metric='cosine'):
    """
    Create k-NN graph using specified distance metric.

    Returns binary adjacency matrix (unweighted, undirected).
    """
    print(f"Building k-NN graph with k={n_neighbors}, metric={metric}...")

    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric, n_jobs=-1)
    nbrs.fit(X)

    # Get k-NN graph (includes self-connections)
    adjacency = nbrs.kneighbors_graph(X, mode='connectivity')

    # Make symmetric (undirected)
    adjacency = adjacency + adjacency.T
    adjacency = (adjacency > 0).astype(int)

    # Remove self-loops
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()

    return adjacency


def find_optimal_k(X, taxonomy_df, k_range, metric='cosine'):
    """
    Find optimal k that maximizes genus-level taxonomic purity.
    """
    print("\nSearching for optimal k...")
    results = []

    for k in k_range:
        adjacency = create_knn_graph(X, k, metric=metric)
        purity = evaluate_graph_quality(adjacency, taxonomy_df, level='genus')

        print(f"  k={k:3d}: genus purity = {purity:.3f}")
        results.append({'k': k, 'genus_purity': purity})

    results_df = pd.DataFrame(results)
    optimal_k = results_df.loc[results_df['genus_purity'].idxmax(), 'k']

    print(f"\nOptimal k = {optimal_k} (genus purity = {results_df['genus_purity'].max():.3f})")

    return int(optimal_k), results_df


def main():
    parser = argparse.ArgumentParser(description='Create final outputs for Component 1.1')
    parser.add_argument('--output-dir', type=str,
                        default='results/1_genome_to_graph/1.1_kmer_profiling/final',
                        help='Output directory')
    parser.add_argument('--n-pca', type=int, default=50,
                        help='Number of PCA components for graph construction')
    parser.add_argument('--metric', type=str, default='cosine',
                        choices=['cosine', 'euclidean', 'correlation'],
                        help='Distance metric for k-NN graph')
    parser.add_argument('--k-min', type=int, default=10,
                        help='Minimum k to test')
    parser.add_argument('--k-max', type=int, default=50,
                        help='Maximum k to test')
    parser.add_argument('--k-step', type=int, default=5,
                        help='Step size for k search')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("Creating final outputs for Component 1.1: K-mer profiling")
    print("="*70)

    # ========================================================================
    # Step 1: Load 6-mer matrix
    # ========================================================================
    print("\n[1] Loading 6-mer matrix...")
    matrix_path = 'results/1_genome_to_graph/1.1_kmer_profiling/6mer/6mer_matrix.npz'
    meta_path = 'results/1_genome_to_graph/1.1_kmer_profiling/6mer/6mer_matrix.meta.npz'

    matrix = load_npz(matrix_path)
    meta = np.load(meta_path, allow_pickle=True)
    genome_ids = meta['genome_ids']
    kmers = meta['kmer_list']

    print(f"Matrix shape: {matrix.shape}")
    print(f"Matrix format: {type(matrix).__name__}")
    print(f"Matrix dtype: {matrix.dtype}")
    print(f"Sparsity: {1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1]):.4f}")
    print(f"Matrix size on disk: {Path(matrix_path).stat().st_size / 1024**2:.1f} MB")

    # ========================================================================
    # Step 2: Compress matrix maximally
    # ========================================================================
    print("\n[2] Compressing matrix...")

    # Convert to CSR format (optimal for row slicing)
    # Scipy sparse matrices only support float32 as the smallest float type
    matrix_csr = matrix.tocsr().astype(np.float32)

    print("\nMatrix already in optimal format (CSR, float32)")
    matrix_compressed = matrix_csr

    # Save compressed matrix
    output_matrix_path = output_dir / '6mer_matrix_compressed.npz'
    save_npz(output_matrix_path, matrix_compressed)
    final_size = output_matrix_path.stat().st_size / 1024**2
    print(f"\nSaved compressed matrix: {output_matrix_path}")
    print(f"Final size: {final_size:.1f} MB")

    # ========================================================================
    # Step 3: Save genome and k-mer indices
    # ========================================================================
    print("\n[3] Saving genome and k-mer indices...")

    # Genome index (with taxonomy)
    taxonomy = extract_taxonomy_from_ids(genome_ids)
    genome_df = pd.DataFrame([
        {
            'genome_index': i,
            'genome_id': gid,
            **taxonomy[gid]
        }
        for i, gid in enumerate(genome_ids)
    ])

    genome_index_path = output_dir / 'genome_index.csv'
    genome_df.to_csv(genome_index_path, index=False)
    print(f"Saved genome index: {genome_index_path}")
    print(f"  {len(genome_df)} genomes")
    print(f"  {genome_df['genus'].nunique()} unique genera")

    # K-mer index
    kmer_df = pd.DataFrame({
        'kmer_index': range(len(kmers)),
        'kmer': kmers
    })

    kmer_index_path = output_dir / 'kmer_index.csv'
    kmer_df.to_csv(kmer_index_path, index=False)
    print(f"Saved k-mer index: {kmer_index_path}")
    print(f"  {len(kmer_df)} canonical 6-mers")

    # ========================================================================
    # Step 4: Determine best representation for graph construction
    # ========================================================================
    print("\n[4] Determining optimal representation for graph construction...")

    # Convert to dense for PCA
    X_full = matrix_compressed.toarray().astype(np.float32)

    # Run PCA
    print(f"\nRunning PCA (n_components={args.n_pca})...")
    pca = PCA(n_components=args.n_pca, random_state=42)
    X_pca = pca.fit_transform(X_full)

    explained_var = pca.explained_variance_ratio_.sum()
    print(f"Explained variance: {explained_var:.3f}")
    print(f"PC1 variance: {pca.explained_variance_ratio_[0]:.3f}")

    # Test which representation captures genus structure better
    print("\nComparing representations for genus-level structure...")

    # Sample k for quick comparison
    k_test = 20

    print(f"\nTesting full 6-mer space (2080D)...")
    adj_full = create_knn_graph(X_full, k_test, metric=args.metric)
    purity_full = evaluate_graph_quality(adj_full, genome_df, level='genus')
    print(f"  Genus purity: {purity_full:.3f}")

    print(f"\nTesting PCA space ({args.n_pca}D)...")
    adj_pca = create_knn_graph(X_pca, k_test, metric=args.metric)
    purity_pca = evaluate_graph_quality(adj_pca, genome_df, level='genus')
    print(f"  Genus purity: {purity_pca:.3f}")

    # Choose best representation
    if purity_pca >= purity_full:
        print(f"\nUsing PCA space ({args.n_pca}D) for graph construction")
        X_graph = X_pca
    else:
        print(f"\nUsing full 6-mer space (2080D) for graph construction")
        X_graph = X_full

    # ========================================================================
    # Step 5: Find optimal k for genus-level structure
    # ========================================================================
    print("\n[5] Finding optimal k for genus-level neighborhood graph...")

    k_range = range(args.k_min, args.k_max + 1, args.k_step)
    optimal_k, k_results = find_optimal_k(X_graph, genome_df, k_range, metric=args.metric)

    # Save k optimization results
    k_results_path = output_dir / 'k_optimization.csv'
    k_results.to_csv(k_results_path, index=False)
    print(f"\nSaved k optimization results: {k_results_path}")

    # ========================================================================
    # Step 6: Create final adjacency matrix
    # ========================================================================
    print("\n[6] Creating final k-NN adjacency matrix...")

    adjacency = create_knn_graph(X_graph, optimal_k, metric=args.metric)

    print(f"\nAdjacency matrix properties:")
    print(f"  Shape: {adjacency.shape}")
    print(f"  Format: {type(adjacency).__name__}")
    print(f"  Edges: {adjacency.nnz // 2}")  # Divide by 2 because symmetric
    print(f"  Avg degree: {adjacency.nnz / adjacency.shape[0]:.1f}")
    print(f"  Sparsity: {1 - adjacency.nnz / (adjacency.shape[0] * adjacency.shape[1]):.6f}")

    # Evaluate final graph
    final_purity = evaluate_graph_quality(adjacency, genome_df, level='genus')
    print(f"\nFinal genus purity: {final_purity:.3f}")

    # Save adjacency matrix (as CSR sparse matrix)
    adjacency_path = output_dir / 'adjacency_matrix.npz'
    save_npz(adjacency_path, adjacency)
    print(f"\nSaved adjacency matrix: {adjacency_path}")
    print(f"Size: {adjacency_path.stat().st_size / 1024**2:.1f} MB")

    # ========================================================================
    # Step 7: Save metadata
    # ========================================================================
    print("\n[7] Saving metadata...")

    metadata = {
        'n_genomes': matrix.shape[0],
        'n_kmers': matrix.shape[1],
        'kmer_size': 6,
        'matrix_dtype': 'float32',
        'matrix_format': 'CSR',
        'graph_representation': 'PCA' if purity_pca >= purity_full else 'full',
        'graph_dimensionality': X_graph.shape[1],
        'distance_metric': args.metric,
        'n_neighbors': int(optimal_k),
        'n_edges': int(adjacency.nnz // 2),
        'genus_purity': float(final_purity),
        'pca_explained_variance': float(explained_var) if purity_pca >= purity_full else None,
    }

    metadata_path = output_dir / 'metadata.txt'
    with open(metadata_path, 'w') as f:
        f.write("Component 1.1 Final Outputs - Metadata\n")
        f.write("="*60 + "\n\n")
        for key, value in metadata.items():
            f.write(f"{key}: {value}\n")

    print(f"Saved metadata: {metadata_path}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*70)
    print("SUMMARY - Final outputs created:")
    print("="*70)
    print(f"\n1. Compressed 6-mer matrix:")
    print(f"   {output_matrix_path}")
    print(f"   Shape: {matrix_compressed.shape}, Size: {final_size:.1f} MB")
    print(f"\n2. Genome index:")
    print(f"   {genome_index_path}")
    print(f"   {len(genome_df)} genomes, {genome_df['genus'].nunique()} genera")
    print(f"\n3. K-mer index:")
    print(f"   {kmer_index_path}")
    print(f"   {len(kmer_df)} canonical 6-mers")
    print(f"\n4. K-NN adjacency matrix:")
    print(f"   {adjacency_path}")
    print(f"   k={optimal_k}, metric={args.metric}, {adjacency.nnz // 2} edges")
    print(f"   Genus purity: {final_purity:.3f}")
    print(f"\n5. Metadata:")
    print(f"   {metadata_path}")
    print("\n" + "="*70)


if __name__ == '__main__':
    main()
