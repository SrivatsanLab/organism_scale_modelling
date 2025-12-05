#!/usr/bin/env python
"""
Analyze correlation between ESM-C and gLM2 pairwise distances.

This script:
1. Loads distance matrices from both models
2. Computes correlation (Pearson, Spearman)
3. Creates visualization plots
4. Analyzes genome-specific patterns

For Issue #2: Compare protein embeddings with/without genomic context
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import squareform


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze correlation between embedding distance matrices"
    )
    parser.add_argument(
        "--embedding-dir",
        type=str,
        default="data/embeddings/protein_only",
        help="Directory containing embeddings and distance matrices",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/embedding_analysis",
        help="Output directory for plots and results",
    )
    return parser.parse_args()


def load_data(embedding_dir: Path):
    """Load embeddings, distance matrices, and metadata."""
    print(f"Loading data from {embedding_dir}...")

    # Load distance matrices
    esmc_dist = np.load(embedding_dir / "esmc_distance_matrix.npy")
    glm2_dist = np.load(embedding_dir / "glm2_protein_only_distance_matrix.npy")

    # Load metadata
    metadata = pd.read_csv(embedding_dir / "protein_metadata.csv")

    print(f"  ESM-C distance matrix: {esmc_dist.shape}")
    print(f"  gLM2 distance matrix: {glm2_dist.shape}")
    print(f"  Metadata: {len(metadata)} proteins")

    return esmc_dist, glm2_dist, metadata


def compute_correlations(esmc_dist, glm2_dist):
    """Compute correlation between distance matrices."""
    print("\n" + "="*80)
    print("Computing Correlations")
    print("="*80)

    # Get upper triangle (excluding diagonal) to avoid redundancy
    n = esmc_dist.shape[0]
    triu_indices = np.triu_indices(n, k=1)

    esmc_distances = esmc_dist[triu_indices]
    glm2_distances = glm2_dist[triu_indices]

    print(f"\nNumber of pairwise distances: {len(esmc_distances):,}")

    # Pearson correlation
    pearson_r, pearson_p = pearsonr(esmc_distances, glm2_distances)
    print(f"\nPearson correlation:")
    print(f"  r = {pearson_r:.4f}")
    print(f"  p-value = {pearson_p:.4e}")

    # Spearman correlation
    spearman_r, spearman_p = spearmanr(esmc_distances, glm2_distances)
    print(f"\nSpearman correlation:")
    print(f"  ρ = {spearman_r:.4f}")
    print(f"  p-value = {spearman_p:.4e}")

    # R-squared
    r_squared = pearson_r ** 2
    print(f"\nR² = {r_squared:.4f}")

    return {
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'r_squared': r_squared,
        'n_pairs': len(esmc_distances),
    }


def plot_distance_correlation(esmc_dist, glm2_dist, output_dir: Path, stats: dict):
    """Create scatter plot of distance correlation."""
    print("\nCreating distance correlation plot...")

    # Get upper triangle
    n = esmc_dist.shape[0]
    triu_indices = np.triu_indices(n, k=1)
    esmc_distances = esmc_dist[triu_indices]
    glm2_distances = glm2_dist[triu_indices]

    # Subsample for plotting (if too many points)
    max_points = 50000
    if len(esmc_distances) > max_points:
        idx = np.random.choice(len(esmc_distances), max_points, replace=False)
        esmc_plot = esmc_distances[idx]
        glm2_plot = glm2_distances[idx]
        print(f"  Subsampling {max_points:,} points for visualization")
    else:
        esmc_plot = esmc_distances
        glm2_plot = glm2_distances

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 10))

    # Hexbin plot for density
    hexbin = ax.hexbin(esmc_plot, glm2_plot, gridsize=50, cmap='Blues',
                       mincnt=1, alpha=0.8)

    # Add diagonal line (perfect correlation)
    lim_max = max(esmc_plot.max(), glm2_plot.max())
    ax.plot([0, lim_max], [0, lim_max], 'r--', alpha=0.5, label='y=x')

    # Labels and title
    ax.set_xlabel('ESM-C Cosine Distance', fontsize=14)
    ax.set_ylabel('gLM2 Cosine Distance', fontsize=14)
    ax.set_title(
        f'ESM-C vs gLM2 Pairwise Distance Correlation\n'
        f'Pearson r = {stats["pearson_r"]:.4f}, R² = {stats["r_squared"]:.4f}',
        fontsize=16
    )

    # Add colorbar
    cbar = plt.colorbar(hexbin, ax=ax)
    cbar.set_label('Count', fontsize=12)

    # Add stats text box
    textstr = '\n'.join([
        f'Pearson r = {stats["pearson_r"]:.4f}',
        f'Spearman ρ = {stats["spearman_r"]:.4f}',
        f'R² = {stats["r_squared"]:.4f}',
        f'n = {stats["n_pairs"]:,} pairs',
    ])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / "distance_correlation_scatter.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def analyze_genome_specificity(esmc_dist, glm2_dist, metadata, output_dir: Path):
    """Analyze within-genome vs between-genome distances."""
    print("\n" + "="*80)
    print("Genome-Specific Analysis")
    print("="*80)

    # Get genome IDs
    genome_ids = metadata['genome_id'].values
    unique_genomes = np.unique(genome_ids)

    print(f"\nGenomes: {list(unique_genomes)}")
    print(f"Proteins per genome:")
    for genome in unique_genomes:
        count = (genome_ids == genome).sum()
        print(f"  {genome}: {count}")

    # Calculate within-genome and between-genome distances
    n = len(genome_ids)

    within_esmc = []
    within_glm2 = []
    between_esmc = []
    between_glm2 = []

    for i in range(n):
        for j in range(i+1, n):
            if genome_ids[i] == genome_ids[j]:
                # Within genome
                within_esmc.append(esmc_dist[i, j])
                within_glm2.append(glm2_dist[i, j])
            else:
                # Between genomes
                between_esmc.append(esmc_dist[i, j])
                between_glm2.append(glm2_dist[i, j])

    within_esmc = np.array(within_esmc)
    within_glm2 = np.array(within_glm2)
    between_esmc = np.array(between_esmc)
    between_glm2 = np.array(between_glm2)

    print(f"\nWithin-genome pairs: {len(within_esmc):,}")
    print(f"Between-genome pairs: {len(between_esmc):,}")

    # Statistics
    print(f"\nESM-C distances:")
    print(f"  Within-genome:  {within_esmc.mean():.4f} ± {within_esmc.std():.4f}")
    print(f"  Between-genome: {between_esmc.mean():.4f} ± {between_esmc.std():.4f}")

    print(f"\ngLM2 distances:")
    print(f"  Within-genome:  {within_glm2.mean():.4f} ± {within_glm2.std():.4f}")
    print(f"  Between-genome: {between_glm2.mean():.4f} ± {between_glm2.std():.4f}")

    # Create violin plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ESM-C
    ax = axes[0]
    data = pd.DataFrame({
        'Distance': np.concatenate([within_esmc, between_esmc]),
        'Type': ['Within-genome']*len(within_esmc) + ['Between-genome']*len(between_esmc)
    })
    sns.violinplot(data=data, x='Type', y='Distance', ax=ax)
    ax.set_title('ESM-C Distances', fontsize=14)
    ax.set_ylabel('Cosine Distance', fontsize=12)

    # gLM2
    ax = axes[1]
    data = pd.DataFrame({
        'Distance': np.concatenate([within_glm2, between_glm2]),
        'Type': ['Within-genome']*len(within_glm2) + ['Between-genome']*len(between_glm2)
    })
    sns.violinplot(data=data, x='Type', y='Distance', ax=ax)
    ax.set_title('gLM2 Distances', fontsize=14)
    ax.set_ylabel('Cosine Distance', fontsize=12)

    plt.tight_layout()
    output_file = output_dir / "genome_specificity_violin.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    plt.close()


def save_summary(stats: dict, output_dir: Path):
    """Save summary statistics to file."""
    output_file = output_dir / "correlation_summary.txt"

    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("ESM-C vs gLM2 Distance Correlation Analysis\n")
        f.write("="*80 + "\n\n")

        f.write(f"Number of pairwise distances: {stats['n_pairs']:,}\n\n")

        f.write("Pearson Correlation:\n")
        f.write(f"  r = {stats['pearson_r']:.6f}\n")
        f.write(f"  p-value = {stats['pearson_p']:.4e}\n\n")

        f.write("Spearman Correlation:\n")
        f.write(f"  ρ = {stats['spearman_r']:.6f}\n")
        f.write(f"  p-value = {stats['spearman_p']:.4e}\n\n")

        f.write(f"R² = {stats['r_squared']:.6f}\n")

    print(f"\nSaved summary: {output_file}")


def main():
    args = parse_args()

    print("="*80)
    print("ESM-C vs gLM2 Distance Correlation Analysis")
    print("="*80)
    print(f"\nEmbedding directory: {args.embedding_dir}")
    print(f"Output directory: {args.output_dir}")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load data
    esmc_dist, glm2_dist, metadata = load_data(Path(args.embedding_dir))

    # Compute correlations
    stats = compute_correlations(esmc_dist, glm2_dist)

    # Create plots
    plot_distance_correlation(esmc_dist, glm2_dist, output_dir, stats)

    # Genome-specific analysis
    analyze_genome_specificity(esmc_dist, glm2_dist, metadata, output_dir)

    # Save summary
    save_summary(stats, output_dir)

    print("\n" + "="*80)
    print("Analysis Complete!")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
