#!/home/dmullane/micromamba/envs/esm3_env/bin/python
"""
create_kmer_matrix.py

Convert individual 5-mer profiles into a single matrix for downstream analysis.

Usage:
    python create_kmer_matrix.py --input-dir results/kmer_profiles/5mer \
                                   --output results/kmer_profiles/5mer_matrix.npz

Output formats:
    - .npz: Sparse matrix (recommended for large datasets)
    - .csv: Full matrix (use --format csv, only for small subsets)
"""

import argparse
import numpy as np
from pathlib import Path
from scipy.sparse import lil_matrix, save_npz
import pandas as pd
from collections import defaultdict
from tqdm import tqdm


def parse_kmer_file(filepath):
    """Parse jellyfish dump file and return k-mer counts as dict."""
    kmer_counts = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                kmer, count = line.strip().split()
                kmer_counts[kmer] = int(count)
    return kmer_counts


def normalize_counts(counts, method='frequency'):
    """
    Normalize k-mer counts.

    Args:
        counts: dict of k-mer -> count
        method: 'frequency', 'presence', or 'none'

    Returns:
        Normalized counts dict
    """
    if method == 'none':
        return counts

    total = sum(counts.values())

    if method == 'frequency':
        # Convert to frequencies (0-1)
        return {k: v/total for k, v in counts.items()}
    elif method == 'presence':
        # Binary: present (1) or absent (0)
        return {k: 1 for k in counts.keys()}
    else:
        raise ValueError(f"Unknown normalization method: {method}")


def create_kmer_matrix(input_dir, output_file, format='npz', normalize='frequency',
                       genome_pattern='*_5mer.txt'):
    """
    Create k-mer frequency matrix from individual genome profiles.

    Args:
        input_dir: Directory containing *_5mer.txt files
        output_file: Output file path
        format: 'npz' (sparse) or 'csv' (dense)
        normalize: 'frequency', 'presence', or 'none'
        genome_pattern: File pattern for k-mer files
    """
    input_path = Path(input_dir)
    kmer_files = sorted(input_path.glob(genome_pattern))

    if not kmer_files:
        raise ValueError(f"No k-mer files found in {input_dir} with pattern {genome_pattern}")

    print(f"Found {len(kmer_files)} k-mer profile files")

    # Step 1: Collect all unique k-mers across all genomes
    print("Collecting unique k-mers...")
    all_kmers = set()
    for kmer_file in tqdm(kmer_files):
        counts = parse_kmer_file(kmer_file)
        all_kmers.update(counts.keys())

    # Create k-mer index
    kmer_to_idx = {kmer: idx for idx, kmer in enumerate(sorted(all_kmers))}
    n_kmers = len(kmer_to_idx)
    n_genomes = len(kmer_files)

    print(f"Total unique k-mers: {n_kmers}")
    print(f"Total genomes: {n_genomes}")
    print(f"Normalization: {normalize}")

    # Step 2: Build matrix
    print("Building k-mer matrix...")

    if format == 'npz':
        # Use sparse matrix for efficiency
        matrix = lil_matrix((n_genomes, n_kmers), dtype=np.float32)
    else:
        # Dense matrix
        matrix = np.zeros((n_genomes, n_kmers), dtype=np.float32)

    genome_ids = []

    for genome_idx, kmer_file in enumerate(tqdm(kmer_files)):
        # Extract genome ID
        genome_id = kmer_file.stem.replace('_5mer', '')
        genome_ids.append(genome_id)

        # Parse k-mer counts
        counts = parse_kmer_file(kmer_file)

        # Normalize if requested
        counts = normalize_counts(counts, method=normalize)

        # Fill matrix
        for kmer, count in counts.items():
            kmer_idx = kmer_to_idx[kmer]
            matrix[genome_idx, kmer_idx] = count

    # Step 3: Save output
    print(f"Saving to {output_file}...")
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == 'npz':
        # Convert to CSR for efficient storage
        matrix_csr = matrix.tocsr()
        save_npz(output_path, matrix_csr)

        # Save metadata separately
        metadata = {
            'genome_ids': genome_ids,
            'kmer_list': sorted(kmer_to_idx.keys())
        }
        np.savez(output_path.with_suffix('.meta.npz'), **metadata)

        print(f"Sparse matrix saved: {output_path}")
        print(f"Metadata saved: {output_path.with_suffix('.meta.npz')}")
        print(f"Matrix shape: {matrix_csr.shape}")
        print(f"Non-zero elements: {matrix_csr.nnz} ({100*matrix_csr.nnz/matrix_csr.size:.2f}%)")

    else:  # CSV
        # Create DataFrame
        df = pd.DataFrame(
            matrix,
            index=genome_ids,
            columns=sorted(kmer_to_idx.keys())
        )
        df.to_csv(output_path)
        print(f"CSV matrix saved: {output_path}")
        print(f"Matrix shape: {df.shape}")

    print("Done!")

    return matrix, genome_ids, kmer_to_idx


def main():
    parser = argparse.ArgumentParser(
        description='Create k-mer frequency matrix from individual profiles',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--input-dir', required=True,
                        help='Directory containing *_5mer.txt files')
    parser.add_argument('--output', required=True,
                        help='Output file path (.npz or .csv)')
    parser.add_argument('--format', choices=['npz', 'csv'], default='npz',
                        help='Output format (default: npz)')
    parser.add_argument('--normalize', choices=['frequency', 'presence', 'none'],
                        default='frequency',
                        help='Normalization method (default: frequency)')
    parser.add_argument('--pattern', default='*_5mer.txt',
                        help='File pattern for k-mer files (default: *_5mer.txt)')

    args = parser.parse_args()

    create_kmer_matrix(
        args.input_dir,
        args.output,
        format=args.format,
        normalize=args.normalize,
        genome_pattern=args.pattern
    )


if __name__ == '__main__':
    main()
