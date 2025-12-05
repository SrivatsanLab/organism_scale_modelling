#!/usr/bin/env python
"""
Analyze MMseqs2 clustering results.

Usage:
    python analyze_mmseqs_clusters.py --clusters clusters.tsv --output summary.csv
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Analyze MMseqs2 clustering results')
    parser.add_argument('--clusters', type=str, required=True,
                        help='MMseqs2 clusters TSV file')
    parser.add_argument('--output', type=str, required=True,
                        help='Output summary CSV file')

    args = parser.parse_args()

    print(f"Analyzing clusters from {args.clusters}")

    # Load cluster assignments
    # TSV format: representative_id, member_id
    clusters = pd.read_csv(args.clusters, sep='\t', header=None,
                           names=['representative', 'member'])

    print(f"  Total proteins: {len(clusters):,}")
    print(f"  Total clusters: {clusters['representative'].nunique():,}")

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

    # Save results
    output_path = Path(args.output)
    assignments_path = output_path.parent / f"{output_path.stem}_assignments.csv"

    clusters.to_csv(assignments_path, index=False)
    cluster_sizes.to_csv(output_path, index=False)

    print(f"\nSaved cluster assignments to {assignments_path}")
    print(f"Saved cluster summary to {output_path}")


if __name__ == '__main__':
    main()
