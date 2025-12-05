#!/usr/bin/env python
"""
Utility script to inspect generated ESM-C embeddings.

Usage:
    python inspect_embeddings.py embeddings/esmc_embeddings.h5
    python inspect_embeddings.py embeddings/ --format npz
"""

import argparse
from pathlib import Path
import h5py
import numpy as np
import pandas as pd


def inspect_hdf5(file_path: Path):
    """Inspect HDF5 embeddings file."""
    print(f"Inspecting HDF5 file: {file_path}")
    print("=" * 70)

    with h5py.File(file_path, "r") as f:
        genome_ids = list(f.keys())
        print(f"\nTotal genomes: {len(genome_ids)}")

        # Collect summary statistics
        stats = []
        for genome_id in genome_ids:
            grp = f[genome_id]
            embeddings = grp["embeddings"]
            gene_ids = grp["gene_ids"]
            lengths = grp["lengths"][:]

            stats.append({
                "genome_id": genome_id,
                "num_genes": len(gene_ids),
                "embedding_dim": embeddings.shape[1],
                "avg_length": lengths.mean(),
                "min_length": lengths.min(),
                "max_length": lengths.max(),
            })

        df = pd.DataFrame(stats)

        print("\nSummary Statistics:")
        print("-" * 70)
        print(f"Total genes: {df['num_genes'].sum():,}")
        print(f"Embedding dimension: {df['embedding_dim'].iloc[0]}")
        print(f"Average genes per genome: {df['num_genes'].mean():.1f}")
        print(f"Average sequence length: {df['avg_length'].mean():.1f}")

        print("\nGenomes with most genes:")
        print(df.nlargest(5, "num_genes")[["genome_id", "num_genes"]])

        print("\nGenomes with longest average sequences:")
        print(df.nlargest(5, "avg_length")[["genome_id", "avg_length"]])

        # Show example
        if genome_ids:
            example_id = genome_ids[0]
            print(f"\nExample genome: {example_id}")
            print("-" * 70)
            grp = f[example_id]

            print(f"Number of genes: {len(grp['gene_ids'])}")
            print(f"Embedding shape: {grp['embeddings'].shape}")
            print(f"First 3 gene IDs:")
            for i, gene_id in enumerate(grp["gene_ids"][:3]):
                print(f"  {i+1}. {gene_id.decode('utf-8')}")

            # Show embedding statistics for first gene
            first_embedding = grp["embeddings"][0]
            print(f"\nFirst gene embedding statistics:")
            print(f"  Mean: {first_embedding.mean():.6f}")
            print(f"  Std:  {first_embedding.std():.6f}")
            print(f"  Min:  {first_embedding.min():.6f}")
            print(f"  Max:  {first_embedding.max():.6f}")


def inspect_npz_directory(directory: Path):
    """Inspect directory of NPZ embedding files."""
    print(f"Inspecting NPZ directory: {directory}")
    print("=" * 70)

    npz_files = list(directory.glob("*_embeddings.npz"))
    print(f"\nTotal genome files: {len(npz_files)}")

    if not npz_files:
        print("No NPZ files found!")
        return

    # Collect summary statistics
    stats = []
    for npz_file in npz_files:
        genome_id = npz_file.stem.replace("_embeddings", "")

        try:
            data = np.load(npz_file, allow_pickle=True)
            embeddings = data["embeddings"]
            gene_ids = data["gene_ids"]
            lengths = data["lengths"]

            stats.append({
                "genome_id": genome_id,
                "num_genes": len(gene_ids),
                "embedding_dim": embeddings.shape[1],
                "avg_length": lengths.mean(),
                "min_length": lengths.min(),
                "max_length": lengths.max(),
            })
        except Exception as e:
            print(f"Error reading {npz_file.name}: {e}")

    df = pd.DataFrame(stats)

    print("\nSummary Statistics:")
    print("-" * 70)
    print(f"Total genes: {df['num_genes'].sum():,}")
    print(f"Embedding dimension: {df['embedding_dim'].iloc[0]}")
    print(f"Average genes per genome: {df['num_genes'].mean():.1f}")
    print(f"Average sequence length: {df['avg_length'].mean():.1f}")

    print("\nGenomes with most genes:")
    print(df.nlargest(5, "num_genes")[["genome_id", "num_genes"]])

    # Show example
    if npz_files:
        example_file = npz_files[0]
        print(f"\nExample file: {example_file.name}")
        print("-" * 70)

        data = np.load(example_file, allow_pickle=True)
        print(f"Number of genes: {len(data['gene_ids'])}")
        print(f"Embedding shape: {data['embeddings'].shape}")
        print(f"First 3 gene IDs:")
        for i, gene_id in enumerate(data["gene_ids"][:3]):
            print(f"  {i+1}. {gene_id}")


def main():
    parser = argparse.ArgumentParser(
        description="Inspect ESM-C embedding files"
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to HDF5 file or directory with NPZ files",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["hdf5", "npz", "auto"],
        default="auto",
        help="Format of embeddings (auto-detected by default)",
    )

    args = parser.parse_args()
    path = Path(args.path)

    if not path.exists():
        print(f"Error: Path does not exist: {path}")
        return

    # Auto-detect format
    if args.format == "auto":
        if path.is_file() and path.suffix == ".h5":
            args.format = "hdf5"
        elif path.is_dir():
            args.format = "npz"
        else:
            print("Error: Could not auto-detect format. Specify --format")
            return

    # Inspect based on format
    if args.format == "hdf5":
        if not path.is_file():
            print("Error: HDF5 format requires a file path")
            return
        inspect_hdf5(path)
    else:  # npz
        if not path.is_dir():
            print("Error: NPZ format requires a directory path")
            return
        inspect_npz_directory(path)


if __name__ == "__main__":
    main()
