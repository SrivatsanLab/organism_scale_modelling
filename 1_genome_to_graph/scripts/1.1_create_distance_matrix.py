#!/usr/bin/env python3
"""
Create pairwise cosine distance matrix for all genomes.

This is the full genome × genome distance matrix that was used to construct
the k-NN adjacency graph.
"""

import numpy as np
from scipy.sparse import load_npz, save_npz
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import cosine_distances
from pathlib import Path
import argparse


def main():
    parser = argparse.ArgumentParser(description='Create pairwise cosine distance matrix')
    parser.add_argument('--output-dir', type=str,
                        default='results/1_genome_to_graph/1.1_kmer_profiling/final',
                        help='Output directory')
    parser.add_argument('--format', type=str, default='npz', choices=['npz', 'npy'],
                        help='Output format: npz (compressed) or npy (uncompressed)')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    print("="*70)
    print("Creating pairwise cosine distance matrix")
    print("="*70)

    # Load 6-mer matrix
    print("\n[1] Loading 6-mer matrix...")
    matrix_path = 'results/1_genome_to_graph/1.1_kmer_profiling/6mer/6mer_matrix.npz'
    matrix = load_npz(matrix_path)

    print(f"Matrix shape: {matrix.shape}")
    print(f"Converting to dense array for distance computation...")
    X = matrix.toarray()

    # Compute pairwise cosine distances
    print("\n[2] Computing pairwise cosine distances...")
    print(f"This will create a {matrix.shape[0]} × {matrix.shape[0]} distance matrix")
    print("(This may take a few minutes...)")

    # Use sklearn's cosine_distances which is efficient for dense arrays
    distance_matrix = cosine_distances(X)

    print(f"Distance matrix shape: {distance_matrix.shape}")
    print(f"Distance matrix dtype: {distance_matrix.dtype}")

    # Verify it's symmetric and has zeros on diagonal
    print("\n[3] Verifying distance matrix properties...")
    diagonal_check = np.allclose(np.diag(distance_matrix), 0, atol=1e-6)
    symmetry_check = np.allclose(distance_matrix, distance_matrix.T, atol=1e-6)

    print(f"  Diagonal is zero: {diagonal_check}")
    print(f"  Matrix is symmetric: {symmetry_check}")
    print(f"  Min distance: {distance_matrix.min():.6f}")
    print(f"  Max distance: {distance_matrix.max():.6f}")
    print(f"  Mean distance: {distance_matrix.mean():.6f}")
    print(f"  Median distance: {np.median(distance_matrix):.6f}")

    # Save distance matrix
    print("\n[4] Saving distance matrix...")

    if args.format == 'npz':
        # Save as compressed numpy archive
        output_path = output_dir / 'distance_matrix.npz'
        np.savez_compressed(output_path, distance_matrix=distance_matrix)
    else:
        # Save as uncompressed numpy array
        output_path = output_dir / 'distance_matrix.npy'
        np.save(output_path, distance_matrix)

    file_size_mb = output_path.stat().st_size / 1024**2
    print(f"Saved to: {output_path}")
    print(f"File size: {file_size_mb:.1f} MB")

    # Memory size comparison
    memory_size_mb = distance_matrix.nbytes / 1024**2
    compression_ratio = memory_size_mb / file_size_mb if args.format == 'npz' else 1.0
    print(f"Memory size: {memory_size_mb:.1f} MB")
    if args.format == 'npz':
        print(f"Compression ratio: {compression_ratio:.2f}x")

    # Update metadata
    print("\n[5] Updating metadata...")
    metadata_path = output_dir / 'metadata.txt'

    with open(metadata_path, 'a') as f:
        f.write("\n# Distance Matrix\n")
        f.write(f"distance_matrix_shape: {distance_matrix.shape}\n")
        f.write(f"distance_matrix_dtype: {distance_matrix.dtype}\n")
        f.write(f"distance_metric: cosine\n")
        f.write(f"distance_min: {distance_matrix.min():.6f}\n")
        f.write(f"distance_max: {distance_matrix.max():.6f}\n")
        f.write(f"distance_mean: {distance_matrix.mean():.6f}\n")
        f.write(f"distance_matrix_size_mb: {file_size_mb:.1f}\n")

    print("Updated metadata file")

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nDistance matrix: {output_path}")
    print(f"Shape: {distance_matrix.shape}")
    print(f"Size: {file_size_mb:.1f} MB")
    print(f"\nLoad with: np.load('{output_path}')['distance_matrix']")
    print("="*70)


if __name__ == '__main__':
    main()
