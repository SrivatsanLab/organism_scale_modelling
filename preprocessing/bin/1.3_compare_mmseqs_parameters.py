#!/usr/bin/env python
"""
Compare MMseqs2 clustering results across different sequence identity thresholds.

Analyzes cluster size distributions to determine optimal parameters for balanced sampling.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def analyze_clustering(cluster_summary_path, seq_id):
    """Analyze a single clustering result."""
    df = pd.read_csv(cluster_summary_path)

    total_clusters = len(df)
    total_proteins = df['size'].sum()

    # Cluster size statistics
    mean_size = df['size'].mean()
    median_size = df['size'].median()
    max_size = df['size'].max()

    # Size distribution
    singletons = (df['size'] == 1).sum()
    small = ((df['size'] >= 2) & (df['size'] < 10)).sum()
    medium = ((df['size'] >= 10) & (df['size'] < 100)).sum()
    large = (df['size'] >= 100).sum()

    pct_singletons = singletons / total_clusters * 100
    pct_small = small / total_clusters * 100
    pct_medium = medium / total_clusters * 100
    pct_large = large / total_clusters * 100

    # Proteins in each category
    proteins_singleton = df[df['size'] == 1]['size'].sum()
    proteins_small = df[(df['size'] >= 2) & (df['size'] < 10)]['size'].sum()
    proteins_medium = df[(df['size'] >= 10) & (df['size'] < 100)]['size'].sum()
    proteins_large = df[df['size'] >= 100]['size'].sum()

    pct_proteins_singleton = proteins_singleton / total_proteins * 100
    pct_proteins_small = proteins_small / total_proteins * 100
    pct_proteins_medium = proteins_medium / total_proteins * 100
    pct_proteins_large = proteins_large / total_proteins * 100

    return {
        'seq_id': seq_id,
        'total_clusters': total_clusters,
        'total_proteins': total_proteins,
        'mean_size': mean_size,
        'median_size': median_size,
        'max_size': max_size,
        'n_singletons': singletons,
        'n_small_2-9': small,
        'n_medium_10-99': medium,
        'n_large_100+': large,
        'pct_singletons': pct_singletons,
        'pct_small': pct_small,
        'pct_medium': pct_medium,
        'pct_large': pct_large,
        'proteins_in_singletons': proteins_singleton,
        'proteins_in_small': proteins_small,
        'proteins_in_medium': proteins_medium,
        'proteins_in_large': proteins_large,
        'pct_proteins_singleton': pct_proteins_singleton,
        'pct_proteins_small': pct_proteins_small,
        'pct_proteins_medium': pct_proteins_medium,
        'pct_proteins_large': pct_proteins_large,
    }


def main():
    print("=" * 80)
    print("MMseqs2 PARAMETER COMPARISON")
    print("=" * 80)
    print()

    # Define parameter runs
    runs = [
        ('data/mmseqs_seqid_0p3/cluster_summary.csv', 0.3),
        ('data/mmseqs_seqid_0p4/cluster_summary.csv', 0.4),
        ('data/mmseqs_seqid_0p5/cluster_summary.csv', 0.5),
        ('data/mmseqs_full_dataset/cluster_summary.csv', 0.5),  # Same as 0.5
        ('data/mmseqs_seqid_0p6/cluster_summary.csv', 0.6),
        ('data/mmseqs_seqid_0p7/cluster_summary.csv', 0.7),
    ]

    # Analyze all runs
    results = []
    for path, seq_id in runs:
        if Path(path).exists():
            result = analyze_clustering(path, seq_id)
            results.append(result)

    # Remove duplicate 0.5
    results_df = pd.DataFrame(results)
    results_df = results_df.drop_duplicates(subset=['seq_id']).sort_values('seq_id')

    print("CLUSTER COUNT BY SIZE CATEGORY")
    print("-" * 80)
    header1 = "{:<10} {:<12} {:<12} {:<12} {:<12} {:<12}".format(
        "Seq ID", "Total", "Singletons", "Small (2-9)", "Med (10-99)", "Large (100+)")
    header2 = "{:<10} {:<12} {:<12} {:<12} {:<12} {:<12}".format(
        "", "Clusters", "(% clust)", "(% clust)", "(% clust)", "(% clust)")
    print(header1)
    print(header2)
    print("-" * 80)

    for _, row in results_df.iterrows():
        print(f"{row['seq_id']:<10.1f} {row['total_clusters']:<12,} "
              f"{row['n_singletons']:>7,} ({row['pct_singletons']:>4.1f}%) "
              f"{row['n_small_2-9']:>7,} ({row['pct_small']:>4.1f}%) "
              f"{row['n_medium_10-99']:>7,} ({row['pct_medium']:>4.1f}%) "
              f"{row['n_large_100+']:>7,} ({row['pct_large']:>4.1f}%)")

    print()
    print("PROTEIN DISTRIBUTION BY CLUSTER SIZE")
    print("-" * 80)
    header3 = "{:<10} {:<12} {:<15} {:<15} {:<15} {:<15}".format(
        "Seq ID", "Total", "In Singleton", "In Small", "In Medium", "In Large")
    header4 = "{:<10} {:<12} {:<15} {:<15} {:<15} {:<15}".format(
        "", "Proteins", "(% proteins)", "(% proteins)", "(% proteins)", "(% proteins)")
    print(header3)
    print(header4)
    print("-" * 80)

    for _, row in results_df.iterrows():
        print(f"{row['seq_id']:<10.1f} {row['total_proteins']:<12,} "
              f"{row['proteins_in_singletons']:>9,} ({row['pct_proteins_singleton']:>4.1f}%) "
              f"{row['proteins_in_small']:>9,} ({row['pct_proteins_small']:>4.1f}%) "
              f"{row['proteins_in_medium']:>9,} ({row['pct_proteins_medium']:>4.1f}%) "
              f"{row['proteins_in_large']:>9,} ({row['pct_proteins_large']:>4.1f}%)")

    print()
    print("CLUSTER SIZE STATISTICS")
    print("-" * 80)
    header5 = "{:<10} {:<10} {:<10} {:<10}".format("Seq ID", "Mean", "Median", "Max")
    print(header5)
    print("-" * 80)

    for _, row in results_df.iterrows():
        print(f"{row['seq_id']:<10.1f} {row['mean_size']:<10.2f} {row['median_size']:<10.0f} {row['max_size']:<10,}")

    print()
    print("=" * 80)
    print("RECOMMENDATIONS FOR BALANCED SAMPLING")
    print("=" * 80)
    print()

    # Find best parameter
    # We want to minimize singletons while maintaining reasonable cluster homogeneity
    # Prioritize having proteins in medium/large clusters

    results_df['score'] = (
        results_df['pct_proteins_medium'] +
        results_df['pct_proteins_large'] * 2 -  # Weight large clusters more
        results_df['pct_proteins_singleton'] * 0.5  # Penalize singletons
    )

    best_idx = results_df['score'].idxmax()
    best_row = results_df.loc[best_idx]

    print(f"Best parameter: Sequence identity = {best_row['seq_id']:.1f}")
    print(f"  Total clusters: {best_row['total_clusters']:,}")
    print(f"  Mean cluster size: {best_row['mean_size']:.1f}")
    print(f"  Median cluster size: {best_row['median_size']:.0f}")
    print(f"  Singletons: {best_row['pct_singletons']:.1f}% of clusters ({best_row['pct_proteins_singleton']:.1f}% of proteins)")
    print(f"  Medium clusters (10-99): {best_row['n_medium_10-99']:,} clusters ({best_row['pct_proteins_medium']:.1f}% of proteins)")
    print(f"  Large clusters (100+): {best_row['n_large_100+']:,} clusters ({best_row['pct_proteins_large']:.1f}% of proteins)")
    print()

    print("Why this is optimal:")
    print(f"  - {100 - best_row['pct_proteins_singleton']:.1f}% of proteins are in multi-member clusters")
    print(f"  - {best_row['n_medium_10-99'] + best_row['n_large_100+']:,} clusters have 10+ members (good for sampling)")
    print(f"  - Balance between cluster granularity and sample diversity")
    print()

    # Save comparison
    output_file = 'results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/mmseqs_parameter_comparison.csv'
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_file, index=False)
    print(f"Detailed comparison saved to: {output_file}")
    print()


if __name__ == '__main__':
    main()
