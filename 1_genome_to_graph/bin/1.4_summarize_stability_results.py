#!/usr/bin/env python
"""
Summarize stability evaluation results across all configurations.

Creates a comparison table showing stability metrics across resolutions.

Usage:
    python summarize_stability_results.py
"""

import numpy as np
import pandas as pd
from pathlib import Path


def main():
    stability_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability')

    # Find all stability result files
    stability_files = sorted(stability_dir.glob('stability_eff_*.npz'))

    if len(stability_files) == 0:
        print(f"No stability files found in {stability_dir}")
        print("Run the stability evaluation first:")
        print("  sbatch scripts/embeddings/submit_stability_efficient.sh")
        return

    print(f"Found {len(stability_files)} stability result files\n")

    # Collect results
    results = []

    for stability_file in stability_files:
        data = np.load(stability_file, allow_pickle=True)

        result = {
            'file': stability_file.name,
            'resolution': float(data['resolution']),
            'n_neighbors': int(data['n_neighbors']),
            'cog_only': bool(data['cog_only']),
            'n_subsamples': int(data['n_subsamples']),
            'subsample_size': int(data['subsample_size']),
            'n_total_genes': int(data['n_total_genes']),
            # ARI metrics
            'ari_mean': float(data['ari_mean']),
            'ari_std': float(data['ari_std']),
            'ari_min': float(data['ari_min']),
            'ari_max': float(data['ari_max']),
            # Gene stability
            'gene_stability_mean': float(data['gene_stability_mean']),
            'gene_stability_median': float(data['gene_stability_median']),
            'pct_genes_stable_90': float(data['pct_genes_stable_90']),
            'pct_genes_stable_80': float(data['pct_genes_stable_80']),
            'pct_genes_stable_70': float(data['pct_genes_stable_70']),
            'n_genes_tracked': int(data['n_genes_tracked']),
            # Cluster count stability
            'n_clusters_mean': float(data['n_clusters_mean']),
            'n_clusters_std': float(data['n_clusters_std']),
            'n_clusters_min': int(data['n_clusters_min']),
            'n_clusters_max': int(data['n_clusters_max']),
        }

        results.append(result)

    # Create DataFrame
    df = pd.DataFrame(results)

    # Sort by resolution and dataset
    df = df.sort_values(['cog_only', 'resolution'], ascending=[False, True])

    # Save full results
    output_file = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability_summary.csv')
    df.to_csv(output_file, index=False)
    print(f"Saved full results to {output_file}\n")

    # Print summary table
    print("=" * 100)
    print("STABILITY EVALUATION SUMMARY")
    print("=" * 100)
    print()

    # COG-only results
    print("COG-ONLY DATASET (798K genes)")
    print("-" * 100)
    cog_df = df[df['cog_only']]
    print(cog_df[['resolution', 'ari_mean', 'ari_std', 'gene_stability_mean',
                   'pct_genes_stable_90', 'n_clusters_mean', 'n_clusters_std']].to_string(index=False))
    print()

    # All genes results
    print("ALL GENES DATASET (1M genes)")
    print("-" * 100)
    all_df = df[~df['cog_only']]
    print(all_df[['resolution', 'ari_mean', 'ari_std', 'gene_stability_mean',
                   'pct_genes_stable_90', 'n_clusters_mean', 'n_clusters_std']].to_string(index=False))
    print()

    # Interpretation
    print("=" * 100)
    print("INTERPRETATION GUIDE")
    print("=" * 100)
    print()
    print("ARI (Adjusted Rand Index between subsamples):")
    print("  > 0.8: Excellent stability - clusters highly reproducible")
    print("  0.6-0.8: Good stability - clusters generally reproducible")
    print("  0.4-0.6: Moderate stability - some variation across samples")
    print("  0.2-0.4: Poor stability - substantial variation")
    print("  < 0.2: Very poor stability - essentially random")
    print()
    print("Gene Stability (fraction assigned to same cluster):")
    print("  % genes with stability > 0.9: Very stable assignments")
    print("  Mean stability: Average consistency across all genes")
    print()
    print("Cluster Count Stability:")
    print("  Low std: Consistent number of clusters across samples")
    print("  High std: Variable cluster counts (less stable)")
    print()

    # Find best configuration
    print("=" * 100)
    print("RECOMMENDATIONS")
    print("=" * 100)
    print()

    # Best by ARI
    best_ari_idx = df['ari_mean'].idxmax()
    best_ari = df.loc[best_ari_idx]
    print(f"Highest ARI: res={best_ari['resolution']:.0f}, "
          f"{'COG-only' if best_ari['cog_only'] else 'all genes'}, "
          f"ARI={best_ari['ari_mean']:.3f}")

    # Best by gene stability
    best_gene_idx = df['gene_stability_mean'].idxmax()
    best_gene = df.loc[best_gene_idx]
    print(f"Highest gene stability: res={best_gene['resolution']:.0f}, "
          f"{'COG-only' if best_gene['cog_only'] else 'all genes'}, "
          f"stability={best_gene['gene_stability_mean']:.3f}")

    # Best by % highly stable genes
    best_pct_idx = df['pct_genes_stable_90'].idxmax()
    best_pct = df.loc[best_pct_idx]
    print(f"Most highly stable genes: res={best_pct['resolution']:.0f}, "
          f"{'COG-only' if best_pct['cog_only'] else 'all genes'}, "
          f"{best_pct['pct_genes_stable_90']:.1f}% genes with stability > 0.9")

    print()

    # Resolution sweet spot (COG-only)
    print("Resolution trade-off (COG-only):")
    cog_sorted = cog_df.sort_values('resolution')
    for _, row in cog_sorted.iterrows():
        status = "✅" if row['ari_mean'] > 0.6 else "⚠️" if row['ari_mean'] > 0.4 else "❌"
        print(f"  {status} res={row['resolution']:>4.0f}: "
              f"ARI={row['ari_mean']:.3f}, "
              f"{row['n_clusters_mean']:>6.0f} clusters, "
              f"{row['pct_genes_stable_90']:>4.1f}% highly stable genes")

    print()
    print("Recommended configuration:")

    # Find highest resolution with good stability
    good_stability = cog_df[cog_df['ari_mean'] > 0.6]
    if len(good_stability) > 0:
        best = good_stability.sort_values('resolution', ascending=False).iloc[0]
        print(f"  Use res={best['resolution']:.0f}, nn={best['n_neighbors']}, COG-only")
        print(f"  This gives ~{best['n_clusters_mean']:.0f} clusters with good stability (ARI={best['ari_mean']:.3f})")
    else:
        moderate = cog_df[cog_df['ari_mean'] > 0.4]
        if len(moderate) > 0:
            best = moderate.sort_values('resolution', ascending=False).iloc[0]
            print(f"  Use res={best['resolution']:.0f}, nn={best['n_neighbors']}, COG-only")
            print(f"  This gives ~{best['n_clusters_mean']:.0f} clusters with moderate stability (ARI={best['ari_mean']:.3f})")
            print(f"  Note: No configurations have high stability (ARI > 0.6)")
        else:
            print("  WARNING: All configurations have low stability (ARI < 0.4)")
            print("  Consider using lower resolution or different approach")

    print()


if __name__ == '__main__':
    main()
