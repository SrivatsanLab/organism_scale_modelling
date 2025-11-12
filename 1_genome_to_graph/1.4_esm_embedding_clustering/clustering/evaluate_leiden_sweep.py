#!/usr/bin/env python
"""
Evaluate Leiden clustering sweep results.

Loads all clustering results and compares them based on:
- Number of clusters
- COG annotation rate
- Mean COG homogeneity
- Weighted COG homogeneity
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys

def load_clustering_results():
    """Load all Leiden sweep clustering results with distribution metrics."""
    results_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering')

    results = []

    for npz_file in sorted(results_dir.glob('clusters_leiden_*.npz')):
        try:
            data = np.load(npz_file, allow_pickle=True)

            result = {
                'file': npz_file.name,
                'resolution': float(data.get('resolution', 0)),
                'n_neighbors': int(data.get('n_neighbors', 0)),
                'cog_only': bool(data.get('cog_only', False)),
                'n_clusters': int(data.get('n_clusters', 0)),
                'n_genes': len(data['gene_ids'])
            }

            # Extract evaluation metrics if available
            if 'evaluation' in data:
                eval_data = data['evaluation'].item()
                result['mean_homogeneity'] = eval_data['mean_homogeneity']
                result['weighted_homogeneity'] = eval_data['weighted_homogeneity']
                result['annotation_rate'] = eval_data['annotation_rate']

                # Calculate distribution metrics
                cluster_metrics = eval_data['cluster_metrics']
                homogeneities = [m['homogeneity'] for m in cluster_metrics]
                sizes = [m['size'] for m in cluster_metrics]
                total_genes = sum(sizes)

                # Compute counts and percentages at different thresholds
                for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
                    n_clusters = sum(1 for h in homogeneities if h >= thresh)
                    n_genes = sum(s for h, s in zip(homogeneities, sizes) if h >= thresh)

                    result[f'n_clusters_hom>={thresh}'] = n_clusters
                    result[f'pct_clusters_hom>={thresh}'] = 100 * n_clusters / len(cluster_metrics)
                    result[f'n_genes_hom>={thresh}'] = n_genes
                    result[f'pct_genes_hom>={thresh}'] = 100 * n_genes / total_genes

                # Add median homogeneity
                result['median_homogeneity'] = float(np.median(homogeneities))

            else:
                result['mean_homogeneity'] = None
                result['weighted_homogeneity'] = None
                result['annotation_rate'] = None
                result['median_homogeneity'] = None

                # Set distribution metrics to None
                for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
                    result[f'n_clusters_hom>={thresh}'] = None
                    result[f'pct_clusters_hom>={thresh}'] = None
                    result[f'n_genes_hom>={thresh}'] = None
                    result[f'pct_genes_hom>={thresh}'] = None

            results.append(result)

        except Exception as e:
            print(f"Warning: Could not load {npz_file.name}: {e}", file=sys.stderr)
            continue

    return pd.DataFrame(results)


def print_summary(df):
    """Print summary statistics."""
    print("=" * 100)
    print("Leiden Clustering Sweep Summary")
    print("=" * 100)
    print()

    print(f"Total clusterings: {len(df)}")
    print()

    # Split by version
    df_all = df[df['cog_only'] == False]
    df_cog = df[df['cog_only'] == True]

    print(f"All genes: {len(df_all)} clusterings")
    print(f"COG-only: {len(df_cog)} clusterings")
    print()

    # Evaluate completeness
    n_with_eval = df['mean_homogeneity'].notna().sum()
    print(f"With evaluation: {n_with_eval}/{len(df)}")
    print()


def print_top_results(df, n=10):
    """Print top results by weighted homogeneity."""
    print("=" * 100)
    print(f"Top {n} by Weighted COG Homogeneity")
    print("=" * 100)
    print()

    # Filter to those with evaluation
    df_eval = df[df['weighted_homogeneity'].notna()].copy()

    if len(df_eval) == 0:
        print("No evaluated clusterings found.")
        return

    # Sort by weighted homogeneity
    df_eval = df_eval.sort_values('weighted_homogeneity', ascending=False)

    # Format for display
    for i, (idx, row) in enumerate(df_eval.head(n).iterrows()):
        print(f"{i+1}. {row['file']}")
        print(f"   Resolution: {row['resolution']:.1f}, N neighbors: {row['n_neighbors']}")
        print(f"   Version: {'COG-only' if row['cog_only'] else 'All genes'}")
        print(f"   Clusters: {row['n_clusters']:,}, Genes: {row['n_genes']:,}")
        print(f"   Weighted homogeneity: {row['weighted_homogeneity']:.4f}")
        print(f"   Mean homogeneity: {row['mean_homogeneity']:.4f}")
        print(f"   Annotation rate: {100*row['annotation_rate']:.1f}%")
        print()


def plot_parameter_effects(df):
    """Analyze effects of different parameters."""
    print("=" * 100)
    print("Parameter Effects on COG Homogeneity")
    print("=" * 100)
    print()

    df_eval = df[df['weighted_homogeneity'].notna()].copy()

    if len(df_eval) == 0:
        print("No evaluated clusterings found.")
        return

    # Group by parameters
    print("By Resolution (averaged over n_neighbors):")
    res_summary = df_eval.groupby(['resolution', 'cog_only']).agg({
        'weighted_homogeneity': 'mean',
        'n_clusters': 'mean'
    }).round(4)
    print(res_summary)
    print()

    print("By N Neighbors (averaged over resolution):")
    nn_summary = df_eval.groupby(['n_neighbors', 'cog_only']).agg({
        'weighted_homogeneity': 'mean',
        'n_clusters': 'mean'
    }).round(4)
    print(nn_summary)
    print()

    print("By Version:")
    version_summary = df_eval.groupby('cog_only').agg({
        'weighted_homogeneity': ['mean', 'std'],
        'mean_homogeneity': ['mean', 'std'],
        'n_clusters': ['mean', 'min', 'max']
    }).round(4)
    print(version_summary)
    print()


def print_distribution_analysis(df, top_n=5):
    """Print homogeneity distribution analysis for top configs."""
    print("=" * 100)
    print(f"Homogeneity Distribution Analysis (Top {top_n} Configs)")
    print("=" * 100)
    print()

    # Filter to evaluated results and sort by weighted homogeneity
    df_eval = df[df['weighted_homogeneity'].notna()].copy()
    df_eval = df_eval.sort_values('weighted_homogeneity', ascending=False)

    for i, (idx, row) in enumerate(df_eval.head(top_n).iterrows()):
        print(f"{i+1}. {row['file']}")
        print(f"   Resolution: {row['resolution']:.1f}, N neighbors: {row['n_neighbors']}, "
              f"Version: {'COG-only' if row['cog_only'] else 'All genes'}")
        print(f"   Total clusters: {row['n_clusters']:,}, Total genes: {row['n_genes']:,}")
        print()

        # Print distribution
        print(f"   {'Threshold':<12} {'# Clusters':<12} {'% Clusters':<12} {'# Genes':<12} {'% Genes':<12}")
        print(f"   {'-' * 60}")

        for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
            n_clust = int(row[f'n_clusters_hom>={thresh}'])
            pct_clust = row[f'pct_clusters_hom>={thresh}']
            n_genes = int(row[f'n_genes_hom>={thresh}'])
            pct_genes = row[f'pct_genes_hom>={thresh}']

            print(f"   ≥{thresh:<11.2f} {n_clust:<12,} {pct_clust:<12.1f} {n_genes:<12,} {pct_genes:<12.1f}")

        print()
        print(f"   Mean homogeneity: {row['mean_homogeneity']:.3f}")
        print(f"   Median homogeneity: {row['median_homogeneity']:.3f}")
        print(f"   Weighted homogeneity: {row['weighted_homogeneity']:.3f}")
        print()


def export_results(df):
    """Export results to CSV."""
    output_file = 'results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/leiden_sweep_summary.csv'
    df.to_csv(output_file, index=False)
    print(f"Exported results to {output_file}")
    print(f"  Total rows: {len(df)}")
    print(f"  Columns: {len(df.columns)}")
    print(f"    Base metrics: file, resolution, n_neighbors, cog_only, n_clusters, n_genes")
    print(f"    Homogeneity: mean, median, weighted, annotation_rate")
    print(f"    Distribution: n_clusters_hom>=X, pct_clusters_hom>=X, n_genes_hom>=X, pct_genes_hom>=X")
    print(f"                  (for X = 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)")


def main():
    print("Loading Leiden clustering results...\n")
    df = load_clustering_results()

    if len(df) == 0:
        print("No clustering results found!")
        return

    print_summary(df)
    print()

    print_top_results(df, n=15)
    print()

    plot_parameter_effects(df)
    print()

    print_distribution_analysis(df, top_n=5)
    print()

    export_results(df)


if __name__ == '__main__':
    main()
