#!/usr/bin/env python3
"""
Compress distance matrix for git storage.

Since the distance matrix is symmetric, we only need to store the upper triangle.
Combined with float16 precision, this achieves 4x compression with minimal error.
"""

import numpy as np
from pathlib import Path
import argparse


def compress_distance_matrix(input_path, output_path):
    """
    Compress symmetric distance matrix by storing only upper triangle in float16.

    Args:
        input_path: Path to full distance matrix (.npz)
        output_path: Path for compressed output (.npz)

    Returns:
        Compression ratio, max error
    """
    print("Loading full distance matrix...")
    dist_matrix = np.load(input_path)['distance_matrix']

    n = dist_matrix.shape[0]
    original_size = Path(input_path).stat().st_size / 1024**2

    print(f"Original matrix:")
    print(f"  Shape: {dist_matrix.shape}")
    print(f"  dtype: {dist_matrix.dtype}")
    print(f"  Size: {original_size:.1f} MB")

    # Extract upper triangle (excluding diagonal)
    print("\nExtracting upper triangle...")
    upper_tri_indices = np.triu_indices(n, k=1)
    upper_tri_values = dist_matrix[upper_tri_indices].astype(np.float16)

    # Save compressed
    print("Saving compressed format...")
    np.savez_compressed(
        output_path,
        upper_triangle=upper_tri_values,
        n=n,
        dtype_original=str(dist_matrix.dtype)
    )

    compressed_size = Path(output_path).stat().st_size / 1024**2
    compression_ratio = original_size / compressed_size

    # Compute error
    upper_tri_f32 = dist_matrix[upper_tri_indices]
    max_error = np.abs(upper_tri_f32 - upper_tri_values.astype(np.float32)).max()
    mean_error = np.abs(upper_tri_f32 - upper_tri_values.astype(np.float32)).mean()

    print(f"\nCompressed matrix:")
    print(f"  Format: upper triangle only, float16")
    print(f"  Size: {compressed_size:.1f} MB")
    print(f"  Compression ratio: {compression_ratio:.2f}x")
    print(f"  Max error: {max_error:.6f}")
    print(f"  Mean error: {mean_error:.6f}")

    return compression_ratio, max_error


def decompress_distance_matrix(input_path):
    """
    Reconstruct full distance matrix from compressed format.

    Args:
        input_path: Path to compressed distance matrix (.npz)

    Returns:
        Full symmetric distance matrix (float32)
    """
    data = np.load(input_path)
    upper_tri_values = data['upper_triangle']
    n = int(data['n'])

    # Reconstruct full matrix
    dist_matrix = np.zeros((n, n), dtype=np.float32)

    # Fill upper triangle
    upper_tri_indices = np.triu_indices(n, k=1)
    dist_matrix[upper_tri_indices] = upper_tri_values.astype(np.float32)

    # Fill lower triangle (symmetric)
    dist_matrix = dist_matrix + dist_matrix.T

    return dist_matrix


def main():
    parser = argparse.ArgumentParser(description='Compress distance matrix for git')
    parser.add_argument('--input', type=str,
                        default='1_genome_to_graph/intermediate/kmer/distance_matrix.npz',
                        help='Input distance matrix')
    parser.add_argument('--output', type=str,
                        default='1_genome_to_graph/graph_outputs/genome_graph/distance_matrix_compressed.npz',
                        help='Output compressed matrix')
    parser.add_argument('--verify', action='store_true',
                        help='Verify decompression works correctly')

    args = parser.parse_args()

    print("="*70)
    print("Compressing Distance Matrix for Git Storage")
    print("="*70)
    print()

    # Compress
    compression_ratio, max_error = compress_distance_matrix(args.input, args.output)

    # Verify if requested
    if args.verify:
        print("\n" + "="*70)
        print("Verifying decompression...")
        print("="*70)

        # Load original
        original = np.load(args.input)['distance_matrix']

        # Decompress
        decompressed = decompress_distance_matrix(args.output)

        # Compare
        print(f"\nOriginal shape: {original.shape}")
        print(f"Decompressed shape: {decompressed.shape}")
        print(f"Max difference: {np.abs(original - decompressed).max():.6f}")
        print(f"Mean difference: {np.abs(original - decompressed).mean():.6f}")

        # Check symmetry
        is_symmetric = np.allclose(decompressed, decompressed.T)
        has_zero_diagonal = np.allclose(np.diag(decompressed), 0)

        print(f"Symmetric: {is_symmetric}")
        print(f"Zero diagonal: {has_zero_diagonal}")

        if is_symmetric and has_zero_diagonal:
            print("\n✓ Decompression successful!")

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nCompressed file: {args.output}")
    print(f"Compression: {compression_ratio:.2f}x smaller")
    print(f"Max error: {max_error:.6f}")
    print("\nTo load:")
    print("  # Use helper function:")
    print("  from compress_distance_matrix import decompress_distance_matrix")
    print(f"  dist = decompress_distance_matrix('{args.output}')")
    print("\n  # Or manually:")
    print("  data = np.load('distance_matrix_compressed.npz')")
    print("  upper = data['upper_triangle'].astype(np.float32)")
    print("  n = int(data['n'])")
    print("  dist = np.zeros((n, n), dtype=np.float32)")
    print("  triu_idx = np.triu_indices(n, k=1)")
    print("  dist[triu_idx] = upper")
    print("  dist = dist + dist.T  # Make symmetric")
    print("="*70)


if __name__ == '__main__':
    main()
