#!/usr/bin/env python
"""
Cluster protein sequences using MMseqs2 to identify homologous groups.

This script runs MMseqs2 to cluster proteins by sequence similarity,
then analyzes the resulting clusters to enable balanced sampling.

Usage:
    python cluster_proteins_mmseqs.py --input proteins.faa --output clusters.tsv
"""

import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import shutil


def run_mmseqs_clustering(input_fasta, output_dir, min_seq_id=0.5, coverage=0.8,
                          cluster_mode=0, threads=8):
    """
    Run MMseqs2 clustering on protein sequences.

    Args:
        input_fasta: Input protein FASTA file
        output_dir: Output directory for clustering results
        min_seq_id: Minimum sequence identity (0-1)
        coverage: Minimum coverage (0-1)
        cluster_mode: Clustering mode (0=greedy, 1=connected component, 2=greedy incremental)
        threads: Number of threads to use

    Returns:
        Path to cluster TSV file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create temporary directory for MMseqs2
    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(exist_ok=True)

    # Output files
    db_file = output_dir / 'seq_db'
    cluster_db = output_dir / 'cluster_db'
    cluster_tsv = output_dir / 'clusters.tsv'

    print(f"Running MMseqs2 clustering:")
    print(f"  Input: {input_fasta}")
    print(f"  Min sequence identity: {min_seq_id}")
    print(f"  Min coverage: {coverage}")
    print(f"  Cluster mode: {cluster_mode}")
    print(f"  Threads: {threads}")
    print()

    try:
        # Load MMseqs2 module
        print("[1/4] Loading MMseqs2 module...")
        # Note: Module loading needs to be done in the shell, not here

        # Create sequence database
        print("[2/4] Creating sequence database...")
        cmd = f"mmseqs createdb {input_fasta} {db_file}"
        subprocess.run(cmd, shell=True, check=True)

        # Run clustering
        print("[3/4] Running clustering (this may take a while)...")
        cmd = (f"mmseqs cluster {db_file} {cluster_db} {tmp_dir} "
               f"--min-seq-id {min_seq_id} -c {coverage} "
               f"--cluster-mode {cluster_mode} --threads {threads}")
        subprocess.run(cmd, shell=True, check=True)

        # Convert to TSV
        print("[4/4] Converting results to TSV...")
        cmd = f"mmseqs createtsv {db_file} {db_file} {cluster_db} {cluster_tsv}"
        subprocess.run(cmd, shell=True, check=True)

        # Clean up temporary files
        print("\nCleaning up temporary files...")
        shutil.rmtree(tmp_dir)
        for f in [db_file, f"{db_file}.index", f"{db_file}.dbtype",
                  cluster_db, f"{cluster_db}.index", f"{cluster_db}.dbtype"]:
            if Path(f).exists():
                Path(f).unlink()

        print(f"\nClustering complete! Results saved to {cluster_tsv}")
        return cluster_tsv

    except subprocess.CalledProcessError as e:
        print(f"\nError running MMseqs2: {e}")
        return None


def analyze_clusters(cluster_tsv, output_summary):
    """
    Analyze MMseqs2 clustering results.

    Args:
        cluster_tsv: Path to MMseqs2 cluster TSV file
        output_summary: Path to output summary CSV

    Returns:
        DataFrame with cluster statistics
    """
    print(f"\nAnalyzing clusters from {cluster_tsv}")

    # Load cluster assignments
    # TSV format: representative_id, member_id
    clusters = pd.read_csv(cluster_tsv, sep='\t', header=None,
                           names=['representative', 'member'])

    print(f"  Total proteins: {len(clusters)}")
    print(f"  Total clusters: {clusters['representative'].nunique()}")

    # Count cluster sizes
    cluster_sizes = clusters.groupby('representative').size().reset_index(name='size')
    cluster_sizes = cluster_sizes.sort_values('size', ascending=False)

    print(f"\nCluster size distribution:")
    print(f"  Mean: {cluster_sizes['size'].mean():.1f}")
    print(f"  Median: {cluster_sizes['size'].median():.1f}")
    print(f"  Max: {cluster_sizes['size'].max()}")
    print(f"  Min: {cluster_sizes['size'].min()}")
    print()

    # Print size histogram
    print("Cluster size histogram:")
    bins = [1, 2, 5, 10, 20, 50, 100, 500, np.inf]
    labels = ['1', '2-4', '5-9', '10-19', '20-49', '50-99', '100-499', '500+']
    cluster_sizes['size_bin'] = pd.cut(cluster_sizes['size'], bins=bins, labels=labels, right=False)
    size_hist = cluster_sizes['size_bin'].value_counts().sort_index()

    for size_range, count in size_hist.items():
        pct = count / len(cluster_sizes) * 100
        print(f"  {size_range:>10}: {count:>8,} clusters ({pct:>5.1f}%)")

    # Save cluster assignments and summary
    clusters.to_csv(output_summary.replace('.csv', '_assignments.csv'), index=False)
    cluster_sizes.to_csv(output_summary, index=False)

    print(f"\nSaved cluster assignments to {output_summary.replace('.csv', '_assignments.csv')}")
    print(f"Saved cluster sizes to {output_summary}")

    return cluster_sizes


def main():
    parser = argparse.ArgumentParser(description='Cluster proteins using MMseqs2')
    parser.add_argument('--input', type=str, required=True,
                        help='Input protein FASTA file')
    parser.add_argument('--output-dir', type=str, default='data/mmseqs_clusters',
                        help='Output directory (default: data/mmseqs_clusters)')
    parser.add_argument('--min-seq-id', type=float, default=0.5,
                        help='Minimum sequence identity (default: 0.5)')
    parser.add_argument('--coverage', type=float, default=0.8,
                        help='Minimum coverage (default: 0.8)')
    parser.add_argument('--cluster-mode', type=int, default=0, choices=[0, 1, 2],
                        help='Clustering mode: 0=greedy, 1=connected component, 2=greedy incremental (default: 0)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads (default: 8)')

    args = parser.parse_args()

    # Run clustering
    cluster_tsv = run_mmseqs_clustering(
        args.input,
        args.output_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        cluster_mode=args.cluster_mode,
        threads=args.threads
    )

    if cluster_tsv is None:
        print("\nClustering failed!")
        return 1

    # Analyze results
    output_summary = Path(args.output_dir) / 'cluster_summary.csv'
    analyze_clusters(cluster_tsv, str(output_summary))

    print("\nDone!")
    return 0


if __name__ == '__main__':
    exit(main())
