#!/usr/bin/env python
"""
Plot all UMAP x Clustering combinations.

For each UMAP embedding (n_neighbors = 15, 25, 50, 100, 200):
    For each clustering result:
        Create a figure with 2 subplots:
            - Left: UMAP colored by cluster labels
            - Right: UMAP colored by COG categories

Saves plots to: results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/umap_n{N}/{clustering_method}.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
sys.path.insert(0, '/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling')

from scripts.embeddings.load_clustering_results import get_clustering_summary, load_cluster_labels

# COG color scheme
COG_COLORS = {
    'J': '#e41a1c',  # Translation
    'K': '#377eb8',  # Transcription
    'L': '#4daf4a',  # Replication
    'D': '#984ea3',  # Cell cycle
    'O': '#ff7f00',  # Post-translational modification
    'M': '#ffff33',  # Cell wall/membrane
    'N': '#a65628',  # Cell motility
    'P': '#f781bf',  # Inorganic ion transport
    'T': '#999999',  # Signal transduction
    'C': '#fc8d62',  # Energy production
    'G': '#8da0cb',  # Carbohydrate transport
    'E': '#e78ac3',  # Amino acid transport
    'F': '#a6d854',  # Nucleotide transport
    'H': '#ffd92f',  # Coenzyme transport
    'I': '#e5c494',  # Lipid transport
    'Q': '#b3b3b3',  # Secondary metabolites
    'R': '#bebada',  # General function prediction
    'S': '#fb8072',  # Function unknown
    'U': '#80b1d3',  # Intracellular trafficking
    'V': '#fdb462',  # Defense mechanisms
    'W': '#b3de69',  # Extracellular structures
    'Y': '#fccde5',  # Nuclear structure
    'Z': '#d9d9d9',  # Cytoskeleton
    'X': '#bc80bd',  # Mobilome
    'No COG': '#d3d3d3'  # No annotation
}


def load_cog_annotations_for_genes(gene_ids, genome_ids):
    """Load COG annotations only for specified genes (much faster!)."""
    import pandas as pd

    print("Loading COG annotations for UMAP subset...")

    # Load COG database files
    cog_db_dir = Path('/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/cog')
    cog_def_file = cog_db_dir / 'cog-20.def.tab'
    cog_csv_file = cog_db_dir / 'cog-20.cog.csv'

    # Load protein → COG mapping
    print(f"  Loading protein→COG mapping from {cog_csv_file.name}...")
    prot_to_cog = {}
    with open(cog_csv_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            fields = line.strip().split(',')
            if len(fields) >= 7:
                prot_id = fields[2]  # Protein ID
                cog_id = fields[6]   # COG ID
                if prot_id not in prot_to_cog:
                    prot_to_cog[prot_id] = cog_id
    print(f"  Loaded {len(prot_to_cog):,} protein→COG mappings")

    # Load COG → category mapping
    print(f"  Loading COG→category mapping from {cog_def_file.name}...")
    cog_def = pd.read_csv(cog_def_file, sep='\t', header=None, encoding='latin-1',
                         names=['cog_id', 'category', 'description', 'gene', 'pathway', 'extra1', 'extra2'])
    cog_to_category = dict(zip(cog_def['cog_id'], cog_def['category']))

    # Get unique genomes in this dataset
    unique_genomes = set(genome_ids)
    print(f"  Need annotations for {len(gene_ids):,} genes from {len(unique_genomes)} genomes")

    # Create gene_id set for fast lookup
    gene_id_set = set(gene_ids)

    # Load DIAMOND hits only for needed genomes
    diamond_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation')
    cog_lookup = {}

    for genome_id in unique_genomes:
        hit_file = diamond_dir / f'{genome_id}_cog_diamond' / f'{genome_id}_cog_hits.tsv'

        if not hit_file.exists():
            continue

        try:
            # Read DIAMOND hits
            hits_df = pd.read_csv(hit_file, sep='\t',
                                 names=['gene_id', 'protein_id', 'pident', 'length', 'mismatch',
                                       'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                       'evalue', 'bitscore', 'qcovhsp', 'scovhsp'])

            # Only process genes that are in our UMAP subset
            hits_df = hits_df[hits_df['gene_id'].isin(gene_id_set)]

            # For each gene, take best hit (already sorted by bitscore in DIAMOND output)
            for gene_id, group in hits_df.groupby('gene_id'):
                best_hit = group.iloc[0]  # First hit is best
                protein_id = best_hit['protein_id']

                if pd.notna(protein_id):
                    # Convert underscore to dot for matching (WP_123_1 → WP_123.1)
                    if '_' in protein_id:
                        parts = protein_id.rsplit('_', 1)
                        protein_id_dot = f"{parts[0]}.{parts[1]}"
                    else:
                        protein_id_dot = protein_id

                    # Map protein → COG → category
                    if protein_id_dot in prot_to_cog:
                        cog_id = prot_to_cog[protein_id_dot]
                        if cog_id in cog_to_category:
                            categories = cog_to_category[cog_id]
                            # Take first category (some COGs have multiple)
                            if len(categories) > 0 and categories[0] != '-':
                                cog_lookup[gene_id] = categories[0]

        except Exception as e:
            print(f"  Warning: Could not load {genome_id}: {e}")
            continue

    n_annotated = len(cog_lookup)
    pct_annotated = 100 * n_annotated / len(gene_ids)
    print(f"  Loaded {n_annotated:,} COG annotations ({pct_annotated:.1f}%)")

    return cog_lookup


def load_umap_embedding(n_neighbors=15, subsample=1000000, path=None):
    """Load UMAP embedding for specific n_neighbors."""
    if path == None:
        umap_file = f'results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n{n_neighbors}_subsample{subsample}.npz'
    else:
        umap_file = path

    print(f"Loading UMAP n={n_neighbors}...")
    data = np.load(umap_file, allow_pickle=True)

    umap_coords = data['umap_embedding']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Shape: {umap_coords.shape}")

    return umap_coords, gene_ids, genome_ids


def map_cogs_to_genes(gene_ids, cog_lookup):
    """Map COG categories to gene IDs."""
    cog_categories = np.array([cog_lookup.get(gid, 'No COG') for gid in gene_ids])

    n_annotated = (cog_categories != 'No COG').sum()
    pct_annotated = 100 * n_annotated / len(cog_categories)

    print(f"  COG annotation rate: {pct_annotated:.1f}% ({n_annotated:,}/{len(cog_categories):,})")

    return cog_categories


def plot_clustering_with_cog(umap_coords, cluster_labels, cog_categories,
                              clustering_name, n_neighbors, output_file):
    """
    Create a figure with 2 subplots: clusters and COG annotations.
    """
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))

    # === LEFT PLOT: Clusters ===
    ax = axes[0]

    # Plot noise points first (gray, small)
    mask_noise = cluster_labels == -1
    mask_clustered = ~mask_noise

    if np.any(mask_noise):
        ax.scatter(
            umap_coords[mask_noise, 0],
            umap_coords[mask_noise, 1],
            c='lightgray',
            s=0.1,
            alpha=0.3,
            rasterized=True,
            label='Noise'
        )

    # Plot clustered points
    if np.any(mask_clustered):
        scatter = ax.scatter(
            umap_coords[mask_clustered, 0],
            umap_coords[mask_clustered, 1],
            c=cluster_labels[mask_clustered],
            s=0.5,
            alpha=0.5,
            cmap='tab20',
            rasterized=True
        )

    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = np.sum(mask_noise)
    pct_noise = 100 * n_noise / len(cluster_labels)

    ax.set_title(
        f'{clustering_name}\n{n_clusters} clusters, {n_noise:,} noise ({pct_noise:.1f}%)',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel(f'UMAP 1 (n_neighbors={n_neighbors})', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.grid(True, alpha=0.2)

    # === RIGHT PLOT: COG Categories ===
    ax = axes[1]

    # Plot each COG category
    for cog, color in COG_COLORS.items():
        mask = cog_categories == cog
        if mask.sum() > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=color,
                s=0.5,
                alpha=0.3,
                label=cog,
                rasterized=True
            )

    n_annotated = (cog_categories != 'No COG').sum()
    pct_annotated = 100 * n_annotated / len(cog_categories)

    ax.set_title(
        f'COG Categories\n{pct_annotated:.1f}% annotated ({n_annotated:,}/{len(cog_categories):,})',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel(f'UMAP 1 (n_neighbors={n_neighbors})', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        markerscale=10,
        fontsize=9,
        framealpha=0.9
    )
    ax.grid(True, alpha=0.2)

    plt.tight_layout()

    # Save figure
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"  Saved: {output_file}")


def main():
    print("=" * 80)
    print("Plotting All UMAP × Clustering Combinations")
    print("=" * 80)
    print()

    # Configuration
    n_neighbors_values = [15, 25, 50, 100, 200]
    subsample = 1000000
    results_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering')
    plots_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/plots')

    # Load one UMAP file to get the gene list (all UMAPs use same genes)
    print("Loading gene list from UMAP n=15...")
    sample_umap_coords, sample_gene_ids, sample_genome_ids = load_umap_embedding(15, subsample)
    print()

    # Load COG annotations only for these genes
    cog_lookup = load_cog_annotations_for_genes(sample_gene_ids, sample_genome_ids)
    print()

    # Get available clustering results
    summary = get_clustering_summary()

    if len(summary) == 0:
        print("No clustering results found yet. Waiting for jobs to complete...")
        return

    print(f"Found {len(summary)} clustering results")
    print()

    # Filter to PCA-based clustering only
    summary = summary[summary['use_pca'] == True]
    print(f"Using {len(summary)} PCA-based clustering results")
    print()

    # Loop through UMAP embeddings
    for n_neighbors in n_neighbors_values:
        print(f"\n{'='*80}")
        print(f"Processing UMAP n_neighbors={n_neighbors}")
        print(f"{'='*80}\n")

        try:
            # Load UMAP
            umap_coords, gene_ids, genome_ids = load_umap_embedding(n_neighbors, subsample)

            # Map COG categories to these genes
            cog_categories = map_cogs_to_genes(gene_ids, cog_lookup)

            # Output directory for this UMAP
            umap_plots_dir = plots_dir / f'umap_n{n_neighbors}'

            # Loop through clustering results
            for idx, row in summary.iterrows():
                suffix = row['suffix']
                method = row['method']

                print(f"\n  Processing {suffix}...")

                try:
                    # Load cluster labels
                    cluster_labels = load_cluster_labels(suffix)

                    # Check if lengths match
                    if len(cluster_labels) != len(gene_ids):
                        print(f"    WARNING: Length mismatch! Clusters: {len(cluster_labels)}, Genes: {len(gene_ids)}")
                        print(f"    Skipping {suffix}")
                        continue

                    # Create plot
                    output_file = umap_plots_dir / f'{suffix}.png'

                    plot_clustering_with_cog(
                        umap_coords,
                        cluster_labels,
                        cog_categories,
                        clustering_name=suffix.replace('_', ' ').title(),
                        n_neighbors=n_neighbors,
                        output_file=output_file
                    )

                except Exception as e:
                    print(f"    ERROR: {e}")
                    continue

            print(f"\n  Completed UMAP n={n_neighbors}")
            print(f"  Plots saved to: {umap_plots_dir}")

        except Exception as e:
            print(f"  ERROR loading UMAP n={n_neighbors}: {e}")
            continue

    print("\n" + "=" * 80)
    print("All plots generated!")
    print(f"Output directory: {plots_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
