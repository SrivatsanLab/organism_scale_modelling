#!/usr/bin/env python
"""
Visualize UMAP embeddings colored by COG functional categories.

This script creates publication-quality plots of UMAP embeddings with
genes colored by their COG functional annotations.

Usage:
    # Option 1: Load UMAP from .npz and merge with COG metadata
    python plot_cog_umap.py \
        --umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample100000.npz \
        --cog-metadata results/1_genome_to_graph/1.4_esm_embedding_clustering/metadata/cog_annotations \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_cog_visualization.png

    # Option 2: Use pre-merged CSV
    python plot_cog_umap.py \
        --annotated-csv results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_with_cogs.csv \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_cog_visualization.png
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import sys

# Add script directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from load_cog_metadata import load_all_cog_metadata, merge_with_cog

# COG category color palette (distinct colors for visualization)
COG_COLORS = {
    'J': '#e41a1c',  # Translation - red
    'K': '#377eb8',  # Transcription - blue
    'L': '#4daf4a',  # Replication - green
    'D': '#984ea3',  # Cell cycle - purple
    'V': '#ff7f00',  # Defense - orange
    'T': '#ffff33',  # Signal transduction - yellow
    'M': '#a65628',  # Cell wall - brown
    'N': '#f781bf',  # Cell motility - pink
    'U': '#999999',  # Intracellular trafficking - gray
    'O': '#66c2a5',  # Posttranslational modification - teal
    'C': '#fc8d62',  # Energy production - coral
    'G': '#8da0cb',  # Carbohydrate metabolism - lavender
    'E': '#e78ac3',  # Amino acid metabolism - magenta
    'F': '#a6d854',  # Nucleotide metabolism - lime
    'H': '#ffd92f',  # Coenzyme metabolism - gold
    'I': '#e5c494',  # Lipid metabolism - tan
    'P': '#b3b3b3',  # Inorganic ion transport - light gray
    'Q': '#8dd3c7',  # Secondary metabolites - cyan
    'R': '#bebada',  # General function - light purple
    'S': '#fb8072',  # Function unknown - salmon
    'No COG': '#d3d3d3',  # No annotation - light gray
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize UMAP with COG annotations"
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--annotated-csv",
        type=str,
        help="Pre-merged CSV file with UMAP + COG annotations",
    )
    input_group.add_argument(
        "--umap",
        type=str,
        help="UMAP .npz file (requires --cog-metadata)",
    )

    parser.add_argument(
        "--cog-metadata",
        type=str,
        help="Directory with COG metadata files (required with --umap)",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output plot file (.png or .pdf)",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=None,
        help="Subsample to this many points for plotting (default: plot all)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Plot DPI (default: 300)",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=0.5,
        help="Point size for scatter plot (default: 0.5)",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.3,
        help="Point transparency (default: 0.3)",
    )

    args = parser.parse_args()

    # Validate: --umap requires --cog-metadata
    if args.umap and not args.cog_metadata:
        parser.error("--umap requires --cog-metadata")

    return args


def load_umap_from_npz(npz_file):
    """Load UMAP coordinates from .npz file."""
    print(f"Loading UMAP from {npz_file}...")

    data = np.load(npz_file, allow_pickle=True)

    df = pd.DataFrame({
        'gene_id': data['gene_ids'],
        'genome_id': data['genome_ids'],
        'umap_1': data['umap_embedding'][:, 0],
        'umap_2': data['umap_embedding'][:, 1]
    })

    print(f"  Loaded {len(df):,} genes")
    print(f"  From {df['genome_id'].nunique()} genomes")

    return df


def load_and_merge_data(umap_file=None, cog_metadata_dir=None,
                        annotated_csv=None, max_points=None):
    """
    Load UMAP and COG data, either from separate files or pre-merged CSV.

    Args:
        umap_file: Path to UMAP .npz file
        cog_metadata_dir: Directory with COG metadata files
        annotated_csv: Pre-merged CSV file
        max_points: Subsample to this many points

    Returns:
        DataFrame with umap_1, umap_2, gene_id, cog_category, category_name
    """

    if annotated_csv:
        # Load pre-merged CSV
        print(f"Loading pre-merged data from {annotated_csv}...")
        df = pd.read_csv(annotated_csv)
        print(f"  Loaded {len(df):,} genes")

    else:
        # Load UMAP and merge with COG metadata
        umap_df = load_umap_from_npz(umap_file)

        print(f"\nLoading COG metadata from {cog_metadata_dir}...")
        cog_df = load_all_cog_metadata(cog_metadata_dir)

        print("\nMerging UMAP with COG annotations...")
        df = merge_with_cog(umap_df, cog_df, on='gene_id', how='left')

    # Subsample if requested
    if max_points and len(df) > max_points:
        print(f"\nSubsampling to {max_points:,} points for plotting...")
        df = df.sample(n=max_points, random_state=42)

    return df


def plot_umap_by_cog(df, output_file, point_size=0.5, alpha=0.3, dpi=300):
    """
    Create comprehensive UMAP visualization with COG categories.

    Creates a figure with:
    1. Main UMAP colored by COG category
    2. Bar chart of category distributions
    """

    print(f"\nCreating visualization...")

    # Count categories
    category_counts = df['cog_category'].value_counts()

    # Create figure with main plot + side panels
    fig = plt.figure(figsize=(18, 8))

    # Main UMAP plot
    ax_main = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)

    # Plot each COG category separately for better color control
    categories_to_plot = [cat for cat in category_counts.index if cat != 'No COG']

    # Plot unannotated genes first (background)
    if 'No COG' in category_counts.index:
        no_cog_df = df[df['cog_category'] == 'No COG']
        ax_main.scatter(
            no_cog_df['umap_1'],
            no_cog_df['umap_2'],
            c=COG_COLORS.get('No COG', '#d3d3d3'),
            s=point_size,
            alpha=alpha * 0.5,  # More transparent for background
            rasterized=True,
            label='No COG'
        )

    # Plot annotated genes by category
    for category in sorted(categories_to_plot):
        cat_df = df[df['cog_category'] == category]
        cat_name = cat_df['category_name'].iloc[0] if len(cat_df) > 0 else 'Unknown'

        ax_main.scatter(
            cat_df['umap_1'],
            cat_df['umap_2'],
            c=COG_COLORS.get(category, '#999999'),
            s=point_size,
            alpha=alpha,
            rasterized=True,
            label=f'{category}: {cat_name[:30]}'
        )

    ax_main.set_xlabel('UMAP 1', fontsize=12)
    ax_main.set_ylabel('UMAP 2', fontsize=12)
    ax_main.set_title(f'Gene Embeddings Colored by COG Category\n{len(df):,} genes',
                     fontsize=14, fontweight='bold')
    ax_main.grid(True, alpha=0.2)

    # Legend outside plot
    ax_main.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                  fontsize=8, frameon=True, ncol=1)

    # Bar chart of category distribution
    ax_bar = plt.subplot2grid((2, 3), (0, 2))

    top_categories = category_counts.head(15)
    colors = [COG_COLORS.get(cat, '#999999') for cat in top_categories.index]

    ax_bar.barh(range(len(top_categories)), top_categories.values, color=colors)
    ax_bar.set_yticks(range(len(top_categories)))
    ax_bar.set_yticklabels([f'{cat}' for cat in top_categories.index], fontsize=8)
    ax_bar.set_xlabel('Gene Count', fontsize=10)
    ax_bar.set_title('Top COG Categories', fontsize=11, fontweight='bold')
    ax_bar.invert_yaxis()

    # Add count labels
    for i, count in enumerate(top_categories.values):
        pct = 100 * count / len(df)
        ax_bar.text(count, i, f' {count:,} ({pct:.1f}%)',
                   va='center', fontsize=7)

    # Summary statistics
    ax_stats = plt.subplot2grid((2, 3), (1, 2))
    ax_stats.axis('off')

    n_annotated = (df['cog_category'] != 'No COG').sum()
    n_total = len(df)
    pct_annotated = 100 * n_annotated / n_total
    n_genomes = df['genome_id'].nunique() if 'genome_id' in df.columns else 'N/A'
    n_categories = len(category_counts[category_counts.index != 'No COG'])

    stats_text = f"""
Summary Statistics
{'='*30}

Total genes:        {n_total:,}
Annotated:          {n_annotated:,} ({pct_annotated:.1f}%)
Unannotated:        {n_total - n_annotated:,}

Genomes:            {n_genomes}
COG categories:     {n_categories}

Top 3 Categories:
"""

    for i, (cat, count) in enumerate(category_counts.head(3).items(), 1):
        if cat in df.columns or 'category_name' in df.columns:
            cat_name = df[df['cog_category'] == cat]['category_name'].iloc[0] if len(df[df['cog_category'] == cat]) > 0 else 'Unknown'
        else:
            cat_name = cat
        pct = 100 * count / n_total
        stats_text += f"\n{i}. {cat}: {count:,} ({pct:.1f}%)\n   {cat_name[:35]}"

    ax_stats.text(0.05, 0.95, stats_text,
                 transform=ax_stats.transAxes,
                 fontsize=9,
                 verticalalignment='top',
                 fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()

    # Save
    print(f"\nSaving to {output_file}...")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved visualization")

    plt.close()


def plot_individual_categories(df, output_dir, point_size=0.5, alpha=0.5, dpi=150):
    """
    Create individual plots for each major COG category.

    Useful for detailed inspection of specific functional classes.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nCreating individual category plots...")

    category_counts = df['cog_category'].value_counts()

    # Plot top categories (excluding 'No COG')
    top_categories = [cat for cat in category_counts.head(10).index if cat != 'No COG']

    for category in top_categories:
        cat_df = df[df['cog_category'] == category]
        other_df = df[df['cog_category'] != category]

        cat_name = cat_df['category_name'].iloc[0] if len(cat_df) > 0 else category

        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot other genes in background
        ax.scatter(
            other_df['umap_1'],
            other_df['umap_2'],
            c='lightgray',
            s=point_size * 0.5,
            alpha=alpha * 0.3,
            rasterized=True
        )

        # Highlight this category
        ax.scatter(
            cat_df['umap_1'],
            cat_df['umap_2'],
            c=COG_COLORS.get(category, '#999999'),
            s=point_size * 2,
            alpha=alpha,
            rasterized=True
        )

        ax.set_xlabel('UMAP 1', fontsize=12)
        ax.set_ylabel('UMAP 2', fontsize=12)

        count = len(cat_df)
        pct = 100 * count / len(df)
        ax.set_title(f'COG {category}: {cat_name}\n{count:,} genes ({pct:.1f}%)',
                    fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.2)

        output_file = output_dir / f'umap_cog_{category}.png'
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Saved {output_file.name}")


def main():
    args = parse_args()

    print("=" * 70)
    print("UMAP COG Visualization")
    print("=" * 70)

    # Load and merge data
    df = load_and_merge_data(
        umap_file=args.umap,
        cog_metadata_dir=args.cog_metadata,
        annotated_csv=args.annotated_csv,
        max_points=args.max_points
    )

    # Create main visualization
    plot_umap_by_cog(
        df,
        args.output,
        point_size=args.point_size,
        alpha=args.alpha,
        dpi=args.dpi
    )

    # Create individual category plots
    output_dir = Path(args.output).parent / "cog_individual"
    plot_individual_categories(
        df,
        output_dir,
        point_size=args.point_size,
        alpha=args.alpha,
        dpi=150
    )

    print("\n" + "=" * 70)
    print("Complete!")
    print("=" * 70)
    print(f"\nMain plot: {args.output}")
    print(f"Individual plots: {output_dir}/")


if __name__ == "__main__":
    main()
