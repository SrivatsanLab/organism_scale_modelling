#!/usr/bin/env python
"""
Compare different subsampling strategies.

Evaluates:
- Coverage of the full embedding space
- Cluster representation
- Distance distribution
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial.distance import pdist


def load_subsample(subsample_file):
    """Load subsampled data."""
    data = np.load(subsample_file, allow_pickle=True)
    return data


def compute_coverage_metrics(full_embeddings, subsample_indices):
    """
    Compute how well the subsample covers the full space.

    Metrics:
    - Mean distance from each full point to nearest subsample point
    - Max distance from each full point to nearest subsample point
    - Cluster representation (if available)
    """
    print("Computing coverage metrics...")

    subsample_embeddings = full_embeddings[subsample_indices]

    # Sample 100k points from full set for efficiency
    if len(full_embeddings) > 100000:
        test_indices = np.random.choice(len(full_embeddings), size=100000, replace=False)
        test_embeddings = full_embeddings[test_indices]
    else:
        test_embeddings = full_embeddings

    # Compute distances to nearest subsample point (in batches)
    print("  Computing nearest-neighbor distances...")
    batch_size = 10000
    min_distances = []

    for i in range(0, len(test_embeddings), batch_size):
        batch = test_embeddings[i:i+batch_size]
        distances = euclidean_distances(batch, subsample_embeddings)
        min_dist_batch = distances.min(axis=1)
        min_distances.extend(min_dist_batch)

    min_distances = np.array(min_distances)

    metrics = {
        'mean_distance_to_subsample': min_distances.mean(),
        'max_distance_to_subsample': min_distances.max(),
        'median_distance_to_subsample': np.median(min_distances),
        'p95_distance_to_subsample': np.percentile(min_distances, 95)
    }

    print(f"  Mean distance to nearest subsample: {metrics['mean_distance_to_subsample']:.3f}")
    print(f"  Max distance to nearest subsample: {metrics['max_distance_to_subsample']:.3f}")
    print(f"  95th percentile: {metrics['p95_distance_to_subsample']:.3f}")

    return metrics


def compute_diversity_metrics(embeddings):
    """
    Compute diversity within the subsample.

    Metrics:
    - Mean pairwise distance
    - Minimum pairwise distance
    - Distance distribution
    """
    print("Computing diversity metrics...")

    # Sample pairs for efficiency
    if len(embeddings) > 10000:
        sample_indices = np.random.choice(len(embeddings), size=10000, replace=False)
        sample_embeddings = embeddings[sample_indices]
    else:
        sample_embeddings = embeddings

    pairwise_distances = pdist(sample_embeddings, metric='euclidean')

    metrics = {
        'mean_pairwise_distance': pairwise_distances.mean(),
        'min_pairwise_distance': pairwise_distances.min(),
        'median_pairwise_distance': np.median(pairwise_distances),
        'std_pairwise_distance': pairwise_distances.std()
    }

    print(f"  Mean pairwise distance: {metrics['mean_pairwise_distance']:.3f}")
    print(f"  Min pairwise distance: {metrics['min_pairwise_distance']:.3f}")

    return metrics, pairwise_distances


def plot_comparison(subsamples_info, output_dir):
    """Create comparison plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Coverage plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    methods = list(subsamples_info.keys())
    mean_dists = [subsamples_info[m]['coverage']['mean_distance_to_subsample'] for m in methods]
    max_dists = [subsamples_info[m]['coverage']['max_distance_to_subsample'] for m in methods]

    ax = axes[0]
    x = np.arange(len(methods))
    ax.bar(x, mean_dists)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45)
    ax.set_ylabel('Mean Distance to Nearest Subsample Point')
    ax.set_title('Coverage: Lower is Better')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.bar(x, max_dists)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45)
    ax.set_ylabel('Max Distance to Nearest Subsample Point')
    ax.set_title('Worst-Case Coverage: Lower is Better')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'coverage_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Diversity plot
    fig, ax = plt.subplots(figsize=(8, 5))

    for method in methods:
        distances = subsamples_info[method]['pairwise_distances']
        ax.hist(distances, bins=50, alpha=0.5, label=method, density=True)

    ax.set_xlabel('Pairwise Distance')
    ax.set_ylabel('Density')
    ax.set_title('Pairwise Distance Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'diversity_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved comparison plots to {output_dir}/")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Compare subsampling strategies')
    parser.add_argument('--full-data', required=True, help='Full dataset file')
    parser.add_argument('--subsamples', nargs='+', required=True, help='Subsample files to compare')
    parser.add_argument('--output-dir', default='results/subsampling/comparison', help='Output directory')

    args = parser.parse_args()

    print("=" * 80)
    print("Comparing Subsampling Strategies")
    print("=" * 80)
    print()

    # Load full data
    print("Loading full dataset...")
    full_data = np.load(args.full_data, allow_pickle=True)
    if 'embeddings_pca' in full_data:
        full_embeddings = full_data['embeddings_pca']
    else:
        full_embeddings = full_data['embeddings']
    print(f"  Full dataset: {full_embeddings.shape[0]:,} samples")
    print()

    # Analyze each subsample
    subsamples_info = {}

    for subsample_file in args.subsamples:
        method_name = Path(subsample_file).stem
        print(f"\n{'='*80}")
        print(f"Analyzing: {method_name}")
        print('='*80)

        subsample_data = load_subsample(subsample_file)
        subsample_indices = subsample_data['selected_indices']
        subsample_embeddings = subsample_data['embeddings']

        print(f"  Subsample size: {len(subsample_indices):,}")

        # Coverage
        coverage_metrics = compute_coverage_metrics(full_embeddings, subsample_indices)

        # Diversity
        diversity_metrics, pairwise_distances = compute_diversity_metrics(subsample_embeddings)

        subsamples_info[method_name] = {
            'coverage': coverage_metrics,
            'diversity': diversity_metrics,
            'pairwise_distances': pairwise_distances,
            'n_samples': len(subsample_indices)
        }

    # Create comparison plots
    print(f"\n{'='*80}")
    print("Creating comparison plots...")
    print('='*80)
    plot_comparison(subsamples_info, args.output_dir)

    # Print summary
    print(f"\n{'='*80}")
    print("Summary")
    print('='*80)
    print()

    for method, info in subsamples_info.items():
        print(f"{method}:")
        print(f"  Coverage (mean dist to nearest): {info['coverage']['mean_distance_to_subsample']:.3f}")
        print(f"  Diversity (mean pairwise dist): {info['diversity']['mean_pairwise_distance']:.3f}")
        print()


if __name__ == '__main__':
    main()
