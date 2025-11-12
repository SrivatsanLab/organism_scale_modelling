#!/usr/bin/env python
"""
Plot individual COG categories on UMAP embeddings.

For each UMAP embedding (n_neighbors = 15, 25, 50, 100, 200):
    For each COG category (J, K, L, C, E, etc.):
        Create a UMAP plot highlighting genes from that category

Output: results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cog_categories/umap_n{N}_cog_{category}.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# COG categories and their names
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'D': 'Cell cycle control, cell division',
    'O': 'Post-translational modification, protein turnover',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'P': 'Inorganic ion transport and metabolism',
    'T': 'Signal transduction mechanisms',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'U': 'Intracellular trafficking, secretion',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton',
    'X': 'Mobilome: prophages, transposons',
}

# Colors for foreground/background
CATEGORY_COLOR = '#e41a1c'  # Red for highlighted category
BACKGROUND_COLOR = '#d3d3d3'  # Light gray for other genes


def load_cog_annotations_for_genes(gene_ids, genome_ids):
    """Load COG annotations only for specified genes."""
    import pandas as pd

    print("Loading COG annotations...")

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

    # Get unique genomes
    unique_genomes = set(genome_ids)
    print(f"  Loading from {len(unique_genomes)} genomes...")

    gene_id_set = set(gene_ids)
    diamond_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation')
    cog_lookup = {}

    for genome_id in unique_genomes:
        hit_file = diamond_dir / f'{genome_id}_cog_diamond' / f'{genome_id}_cog_hits.tsv'

        if not hit_file.exists():
            continue

        try:
            hits_df = pd.read_csv(hit_file, sep='\t',
                                 names=['gene_id', 'protein_id', 'pident', 'length', 'mismatch',
                                       'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                       'evalue', 'bitscore', 'qcovhsp', 'scovhsp'])

            hits_df = hits_df[hits_df['gene_id'].isin(gene_id_set)]

            for gene_id, group in hits_df.groupby('gene_id'):
                best_hit = group.iloc[0]
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
            continue

    n_annotated = len(cog_lookup)
    pct_annotated = 100 * n_annotated / len(gene_ids)
    print(f"  Annotated: {n_annotated:,} / {len(gene_ids):,} ({pct_annotated:.1f}%)")

    return cog_lookup


def load_umap_embedding(n_neighbors, subsample=1000000):
    """Load UMAP embedding."""
    umap_file = f'results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n{n_neighbors}_subsample{subsample}.npz'

    print(f"Loading UMAP n={n_neighbors}...")
    data = np.load(umap_file, allow_pickle=True)

    umap_coords = data['umap_embedding']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  Shape: {umap_coords.shape}")

    return umap_coords, gene_ids, genome_ids


def plot_cog_category(umap_coords, gene_ids, cog_lookup, category, n_neighbors, output_file):
    """
    Plot UMAP highlighting a specific COG category.
    """
    # Map genes to COG categories
    gene_cogs = np.array([cog_lookup.get(gid, 'None') for gid in gene_ids])

    # Create masks
    mask_category = gene_cogs == category
    mask_other = ~mask_category

    n_category = mask_category.sum()
    pct_category = 100 * n_category / len(gene_ids)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot background genes (light gray, small)
    ax.scatter(
        umap_coords[mask_other, 0],
        umap_coords[mask_other, 1],
        c=BACKGROUND_COLOR,
        s=0.1,
        alpha=0.2,
        rasterized=True,
        label='Other genes'
    )

    # Plot category genes (red, larger, on top)
    if n_category > 0:
        ax.scatter(
            umap_coords[mask_category, 0],
            umap_coords[mask_category, 1],
            c=CATEGORY_COLOR,
            s=1.0,
            alpha=0.6,
            rasterized=True,
            label=f'COG {category}'
        )

    # Title and labels
    category_name = COG_CATEGORIES.get(category, 'Unknown')
    ax.set_title(
        f'COG Category {category}: {category_name}\n'
        f'{n_category:,} genes ({pct_category:.1f}%) | UMAP n_neighbors={n_neighbors}',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel(f'UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.legend(markerscale=10, fontsize=10)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()

    # Save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"  Saved: {output_file.name} ({n_category:,} genes)")


def main():
    print("=" * 80)
    print("Plotting Individual COG Categories on UMAP")
    print("=" * 80)
    print()

    # Configuration
    n_neighbors_values = [15, 25, 50, 100, 200]
    subsample = 1000000
    plots_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cog_categories')

    # Load gene list from one UMAP
    print("Loading gene list from UMAP n=15...")
    sample_umap_coords, sample_gene_ids, sample_genome_ids = load_umap_embedding(15, subsample)
    print()

    # Load COG annotations
    cog_lookup = load_cog_annotations_for_genes(sample_gene_ids, sample_genome_ids)
    print()

    # Count genes per category
    print("COG category distribution:")
    gene_cogs = np.array([cog_lookup.get(gid, 'None') for gid in sample_gene_ids])
    for cat in sorted(COG_CATEGORIES.keys()):
        n_genes = (gene_cogs == cat).sum()
        pct = 100 * n_genes / len(sample_gene_ids)
        if n_genes > 0:
            print(f"  {cat}: {n_genes:,} genes ({pct:.1f}%) - {COG_CATEGORIES[cat]}")
    n_unannotated = (gene_cogs == 'None').sum()
    print(f"  Unannotated: {n_unannotated:,} genes ({100*n_unannotated/len(sample_gene_ids):.1f}%)")
    print()

    # Plot each UMAP x COG category combination
    for n_neighbors in n_neighbors_values:
        print(f"\n{'='*80}")
        print(f"Processing UMAP n_neighbors={n_neighbors}")
        print(f"{'='*80}\n")

        # Load UMAP
        umap_coords, gene_ids, genome_ids = load_umap_embedding(n_neighbors, subsample)

        # Plot each COG category
        for category in sorted(COG_CATEGORIES.keys()):
            # Check if this category has genes
            gene_cogs = np.array([cog_lookup.get(gid, 'None') for gid in gene_ids])
            if (gene_cogs == category).sum() == 0:
                continue  # Skip categories with no genes

            output_file = plots_dir / f'umap_n{n_neighbors}_cog_{category}.png'

            plot_cog_category(
                umap_coords,
                gene_ids,
                cog_lookup,
                category,
                n_neighbors,
                output_file
            )

    print("\n" + "=" * 80)
    print("All COG category plots generated!")
    print(f"Output directory: {plots_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
