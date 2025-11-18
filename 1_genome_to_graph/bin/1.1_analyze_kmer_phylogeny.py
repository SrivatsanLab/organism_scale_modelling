#!/home/dmullane/micromamba/envs/esm3_env/bin/python
"""
analyze_kmer_phylogeny.py

Analyze k-mer profiles to assess phylogenetic signal:
1. Run PCA on k-mer matrices
2. Generate UMAP visualizations
3. Extract taxonomic information from genome IDs
4. Compute clustering metrics by taxonomy

Usage:
    python analyze_kmer_phylogeny.py --kmer-size 5
    python analyze_kmer_phylogeny.py --kmer-size 6
    python analyze_kmer_phylogeny.py --both  # Analyze both 5-mer and 6-mer
"""

import argparse
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict, Counter
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score


def extract_taxonomy_from_ids(genome_ids):
    """
    Extract taxonomic information from RefSeq genome IDs.

    Format: GCF_XXXXXX.Y_Genus_species_strain

    Returns dict with genus, species for each genome.
    """
    taxonomy = {}

    for gid in genome_ids:
        # Remove GCF prefix
        name_parts = gid.split('_', 2)
        if len(name_parts) >= 3:
            # Extract organism name
            organism = name_parts[2]

            # Split into components
            parts = organism.split('_')

            # First part is usually genus, second is species
            genus = parts[0] if len(parts) > 0 else 'Unknown'
            species = parts[1] if len(parts) > 1 else 'sp'

            # Clean up special characters
            genus = re.sub(r'[^a-zA-Z]', '', genus)
            species = re.sub(r'[^a-zA-Z]', '', species)

            taxonomy[gid] = {
                'genus': genus,
                'species': species,
                'full_name': f"{genus}_{species}"
            }
        else:
            taxonomy[gid] = {
                'genus': 'Unknown',
                'species': 'sp',
                'full_name': 'Unknown_sp'
            }

    return taxonomy


def compute_phylogenetic_metrics(embeddings, taxonomy_df, level='genus'):
    """
    Compute metrics quantifying phylogenetic clustering.

    Args:
        embeddings: 2D array of coordinates (N x D)
        taxonomy_df: DataFrame with taxonomy info
        level: 'genus' or 'species'

    Returns:
        dict of metrics
    """
    labels = taxonomy_df[level].values

    # Filter out rare taxa (< 3 genomes)
    label_counts = Counter(labels)
    common_labels = {l for l, c in label_counts.items() if c >= 3}

    mask = np.array([l in common_labels for l in labels])
    filtered_embeddings = embeddings[mask]
    filtered_labels = labels[mask]

    if len(np.unique(filtered_labels)) < 2:
        return {
            'n_taxa': len(label_counts),
            'n_common_taxa': len(common_labels),
            'n_genomes_analyzed': 0,
            'silhouette': np.nan,
            'taxonomic_purity_mean': np.nan,
            'taxonomic_purity_median': np.nan
        }

    # Silhouette score (how well-separated are taxonomic groups?)
    try:
        sil_score = silhouette_score(filtered_embeddings, filtered_labels,
                                     metric='euclidean', sample_size=min(5000, len(filtered_labels)))
    except:
        sil_score = np.nan

    # Taxonomic purity using k-NN (for each genome, what fraction of nearest neighbors are same taxon?)
    from sklearn.neighbors import NearestNeighbors

    k = min(20, len(filtered_embeddings) - 1)  # Adjust k if not enough samples
    nbrs = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nbrs.fit(filtered_embeddings)
    distances, indices = nbrs.kneighbors(filtered_embeddings)

    purities = []
    for i, label in enumerate(filtered_labels):
        # Get labels of k nearest neighbors (excluding self)
        neighbor_labels = filtered_labels[indices[i, 1:]]
        # Fraction that match this genome's taxon
        purity = np.mean(neighbor_labels == label)
        purities.append(purity)

    return {
        'n_taxa': len(label_counts),
        'n_common_taxa': len(common_labels),
        'n_genomes_analyzed': len(filtered_labels),
        'silhouette': sil_score,
        'taxonomic_purity_mean': np.mean(purities),
        'taxonomic_purity_median': np.median(purities),
        'taxonomic_purity_q25': np.percentile(purities, 25),
        'taxonomic_purity_q75': np.percentile(purities, 75)
    }


def analyze_kmer_matrix(kmer_size, output_dir, n_pca_components=50, umap_neighbors=15):
    """
    Complete analysis pipeline for one k-mer size.
    """
    print(f"\n{'='*60}")
    print(f"Analyzing {kmer_size}-mer profiles")
    print(f"{'='*60}\n")

    # Load matrix
    if kmer_size == 5:
        matrix_path = 'results/1_genome_to_graph/1.1_kmer_profiling/5mer/5mer_matrix.npz'
        meta_path = 'results/1_genome_to_graph/1.1_kmer_profiling/5mer/5mer_matrix.meta.npz'
    else:
        matrix_path = 'results/1_genome_to_graph/1.1_kmer_profiling/6mer/6mer_matrix.npz'
        meta_path = 'results/1_genome_to_graph/1.1_kmer_profiling/6mer/6mer_matrix.meta.npz'

    print(f"Loading matrix from: {matrix_path}")
    matrix = load_npz(matrix_path)
    meta = np.load(meta_path, allow_pickle=True)
    genome_ids = meta['genome_ids']

    print(f"Matrix shape: {matrix.shape}")
    print(f"K-mer dimensionality: {matrix.shape[1]}")

    # Convert to dense
    print("Converting to dense array...")
    X = matrix.toarray()

    # Extract taxonomy
    print("Extracting taxonomy from genome IDs...")
    taxonomy = extract_taxonomy_from_ids(genome_ids)
    taxonomy_df = pd.DataFrame([
        {'genome_id': gid, **taxonomy[gid]}
        for gid in genome_ids
    ])

    print(f"Unique genera: {taxonomy_df['genus'].nunique()}")
    print(f"Unique species: {taxonomy_df['full_name'].nunique()}")

    # Save taxonomy
    taxonomy_path = output_dir / f'{kmer_size}mer_taxonomy.csv'
    taxonomy_df.to_csv(taxonomy_path, index=False)
    print(f"Saved taxonomy to: {taxonomy_path}")

    # PCA
    print(f"\nRunning PCA (n_components={n_pca_components})...")
    pca = PCA(n_components=n_pca_components, random_state=42)
    X_pca = pca.fit_transform(X)

    print(f"Explained variance: {pca.explained_variance_ratio_.sum():.3f}")
    print(f"PC1 variance: {pca.explained_variance_ratio_[0]:.3f}")
    print(f"PC2 variance: {pca.explained_variance_ratio_[1]:.3f}")

    # Save PCA results
    pca_path = output_dir / f'{kmer_size}mer_pca.npz'
    np.savez(pca_path,
             pca_coords=X_pca,
             explained_variance_ratio=pca.explained_variance_ratio_,
             genome_ids=genome_ids)
    print(f"Saved PCA to: {pca_path}")

    # UMAP
    print(f"\nRunning UMAP (n_neighbors={umap_neighbors})...")
    reducer = umap.UMAP(n_neighbors=umap_neighbors, min_dist=0.1,
                       n_components=2, random_state=42, metric='euclidean')
    X_umap = reducer.fit_transform(X_pca)

    # Save UMAP results
    umap_path = output_dir / f'{kmer_size}mer_umap.npz'
    np.savez(umap_path,
             umap_coords=X_umap,
             genome_ids=genome_ids)
    print(f"Saved UMAP to: {umap_path}")

    # Compute phylogenetic metrics
    print("\nComputing phylogenetic clustering metrics...")

    metrics_genus = compute_phylogenetic_metrics(X_umap, taxonomy_df, level='genus')
    metrics_species = compute_phylogenetic_metrics(X_umap, taxonomy_df, level='full_name')

    print(f"\nGenus-level metrics:")
    print(f"  Total genera: {metrics_genus['n_taxa']}")
    print(f"  Common genera (≥3 genomes): {metrics_genus['n_common_taxa']}")
    print(f"  Silhouette score: {metrics_genus['silhouette']:.3f}")
    print(f"  Taxonomic purity (mean): {metrics_genus['taxonomic_purity_mean']:.3f}")
    print(f"  Taxonomic purity (median): {metrics_genus['taxonomic_purity_median']:.3f}")

    print(f"\nSpecies-level metrics:")
    print(f"  Total species: {metrics_species['n_taxa']}")
    print(f"  Common species (≥3 genomes): {metrics_species['n_common_taxa']}")
    print(f"  Silhouette score: {metrics_species['silhouette']:.3f}")
    print(f"  Taxonomic purity (mean): {metrics_species['taxonomic_purity_mean']:.3f}")
    print(f"  Taxonomic purity (median): {metrics_species['taxonomic_purity_median']:.3f}")

    # Save metrics
    metrics_df = pd.DataFrame({
        'kmer_size': [kmer_size, kmer_size],
        'taxonomic_level': ['genus', 'species'],
        **{k: [metrics_genus[k], metrics_species[k]]
           for k in metrics_genus.keys()}
    })
    metrics_path = output_dir / f'{kmer_size}mer_metrics.csv'
    metrics_df.to_csv(metrics_path, index=False)
    print(f"\nSaved metrics to: {metrics_path}")

    # Create visualizations
    print("\nGenerating visualizations...")
    create_umap_plots(X_umap, taxonomy_df, kmer_size, output_dir)

    return {
        'X_pca': X_pca,
        'X_umap': X_umap,
        'taxonomy_df': taxonomy_df,
        'metrics_genus': metrics_genus,
        'metrics_species': metrics_species
    }


def create_umap_plots(X_umap, taxonomy_df, kmer_size, output_dir):
    """
    Create UMAP visualization plots colored by taxonomy.
    """
    # Get top 20 most common genera for coloring
    genus_counts = taxonomy_df['genus'].value_counts()
    top_genera = genus_counts.head(20).index.tolist()

    # Create color map
    taxonomy_df['genus_plot'] = taxonomy_df['genus'].apply(
        lambda x: x if x in top_genera else 'Other'
    )

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # Plot 1: Colored by genus
    ax = axes[0]

    # Plot "Other" first in gray
    other_mask = taxonomy_df['genus_plot'] == 'Other'
    ax.scatter(X_umap[other_mask, 0], X_umap[other_mask, 1],
               c='lightgray', s=1, alpha=0.3, label='Other')

    # Plot top genera with distinct colors
    colors = plt.cm.tab20(np.linspace(0, 1, 20))
    for i, genus in enumerate(top_genera):
        mask = taxonomy_df['genus'] == genus
        count = mask.sum()
        ax.scatter(X_umap[mask, 0], X_umap[mask, 1],
                  c=[colors[i]], s=5, alpha=0.6,
                  label=f'{genus} (n={count})')

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'{kmer_size}-mer UMAP - Colored by Genus (Top 20)', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, markerscale=2)

    # Plot 2: Density plot
    ax = axes[1]
    ax.hexbin(X_umap[:, 0], X_umap[:, 1], gridsize=50, cmap='viridis', mincnt=1)
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'{kmer_size}-mer UMAP - Density', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plot_path = output_dir / f'{kmer_size}mer_umap_taxonomy.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved plot to: {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze k-mer profiles for phylogenetic signal',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--kmer-size', type=int, choices=[5, 6],
                       help='K-mer size to analyze (5 or 6)')
    parser.add_argument('--both', action='store_true',
                       help='Analyze both 5-mer and 6-mer')
    parser.add_argument('--output-dir', type=str,
                       default='results/1_genome_to_graph/1.1_kmer_profiling/analysis',
                       help='Output directory for results')
    parser.add_argument('--n-pca', type=int, default=50,
                       help='Number of PCA components (default: 50)')
    parser.add_argument('--umap-neighbors', type=int, default=15,
                       help='UMAP n_neighbors parameter (default: 15)')

    args = parser.parse_args()

    if not args.both and args.kmer_size is None:
        parser.error("Must specify either --kmer-size or --both")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {output_dir}")

    # Analyze requested k-mer sizes
    kmer_sizes = [args.kmer_size] if args.kmer_size else [5, 6]

    all_metrics = []

    for kmer_size in kmer_sizes:
        result = analyze_kmer_matrix(
            kmer_size,
            output_dir,
            n_pca_components=args.n_pca,
            umap_neighbors=args.umap_neighbors
        )

        all_metrics.append({
            'kmer_size': kmer_size,
            'genus_silhouette': result['metrics_genus']['silhouette'],
            'genus_purity_mean': result['metrics_genus']['taxonomic_purity_mean'],
            'species_silhouette': result['metrics_species']['silhouette'],
            'species_purity_mean': result['metrics_species']['taxonomic_purity_mean']
        })

    # Compare k-mer sizes if both analyzed
    if len(kmer_sizes) == 2:
        print(f"\n{'='*60}")
        print("Comparison: 5-mer vs 6-mer")
        print(f"{'='*60}\n")

        comparison_df = pd.DataFrame(all_metrics)
        print(comparison_df.to_string(index=False))

        comparison_path = output_dir / 'kmer_comparison.csv'
        comparison_df.to_csv(comparison_path, index=False)
        print(f"\nSaved comparison to: {comparison_path}")

    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
