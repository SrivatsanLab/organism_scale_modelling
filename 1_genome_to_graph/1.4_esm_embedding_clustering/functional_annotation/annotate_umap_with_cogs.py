#!/usr/bin/env python
"""
Annotate UMAP embeddings with COG functional categories.

This script merges COG annotations from eggNOG-mapper with UMAP coordinates
to enable functional visualization of gene clusters.

Usage:
    python annotate_umap_with_cogs.py \
        --umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample100000.npz \
        --annotations results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_with_cogs.csv
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict

# COG category definitions
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion',
    'O': 'Posttranslational modification, protein turnover',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'X': 'Mobilome (prophages, transposons)',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate UMAP embeddings with COG categories"
    )
    parser.add_argument(
        "--umap",
        type=str,
        required=True,
        help="UMAP .npz file from compute_umap_array.py",
    )
    parser.add_argument(
        "--annotations",
        type=str,
        required=True,
        help="Directory with eggNOG-mapper results",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV file with UMAP + COG annotations",
    )
    return parser.parse_args()


def load_umap_data(umap_file):
    """Load UMAP embeddings and metadata."""
    print(f"Loading UMAP from {umap_file}...")

    data = np.load(umap_file, allow_pickle=True)

    df = pd.DataFrame({
        'gene_id': data['gene_ids'],
        'genome_id': data['genome_ids'],
        'umap_1': data['umap_embedding'][:, 0],
        'umap_2': data['umap_embedding'][:, 1]
    })

    print(f"  Loaded {len(df):,} genes")
    print(f"  From {df['genome_id'].nunique()} genomes")

    return df


def load_cog_annotations(annotations_dir):
    """
    Load COG annotations from all eggNOG-mapper result files.

    Returns:
        dict: {gene_id: {cog_category, cog_id, description}}
    """
    annotations_dir = Path(annotations_dir)
    annotation_files = list(annotations_dir.glob("*_eggnog/*.emapper.annotations"))

    print(f"\nLoading COG annotations from {len(annotation_files)} genomes...")

    cog_dict = {}
    total_annotated = 0

    for annot_file in annotation_files:
        genome_id = annot_file.parent.name.replace('_eggnog', '')

        try:
            with open(annot_file, 'r') as f:
                for line in f:
                    # Skip headers and comments
                    if line.startswith('#') or line.strip() == '':
                        continue

                    fields = line.strip().split('\t')

                    if len(fields) < 7:
                        continue

                    # Parse fields
                    query_id = fields[0]  # Gene ID
                    cog_category = fields[6]  # COG category letters
                    description = fields[7] if len(fields) > 7 else ""

                    # Skip unannotated
                    if cog_category == '-' or not cog_category:
                        continue

                    # Get primary category (first letter)
                    primary_cog = cog_category[0] if cog_category else None

                    if primary_cog and primary_cog in COG_CATEGORIES:
                        cog_dict[query_id] = {
                            'cog_category': primary_cog,
                            'cog_categories_all': cog_category,
                            'description': description,
                            'category_name': COG_CATEGORIES[primary_cog]
                        }
                        total_annotated += 1

        except Exception as e:
            print(f"  Error reading {annot_file.name}: {e}")
            continue

    print(f"  Loaded {total_annotated:,} COG annotations")

    return cog_dict


def merge_umap_with_cogs(umap_df, cog_dict):
    """Merge UMAP coordinates with COG annotations."""
    print("\nMerging UMAP with COG annotations...")

    # Add COG annotations
    umap_df['cog_category'] = umap_df['gene_id'].map(
        lambda x: cog_dict.get(x, {}).get('cog_category', 'No COG')
    )
    umap_df['category_name'] = umap_df['gene_id'].map(
        lambda x: cog_dict.get(x, {}).get('category_name', 'No annotation')
    )
    umap_df['description'] = umap_df['gene_id'].map(
        lambda x: cog_dict.get(x, {}).get('description', '')
    )

    # Count annotations
    n_annotated = (umap_df['cog_category'] != 'No COG').sum()
    n_total = len(umap_df)

    print(f"  Annotated: {n_annotated:,} / {n_total:,} ({100*n_annotated/n_total:.1f}%)")

    # Show category distribution
    print("\nCOG category distribution:")
    category_counts = umap_df['cog_category'].value_counts()
    for cat, count in category_counts.head(15).items():
        if cat == 'No COG':
            cat_name = 'No annotation'
        else:
            cat_name = COG_CATEGORIES.get(cat, 'Unknown')
        pct = 100 * count / n_total
        print(f"  {cat}: {cat_name[:50]:50s} {count:7,} ({pct:5.1f}%)")

    return umap_df


def main():
    args = parse_args()

    print("=" * 70)
    print("Annotate UMAP with COG Categories")
    print("=" * 70)

    # Load UMAP data
    umap_df = load_umap_data(args.umap)

    # Load COG annotations
    cog_dict = load_cog_annotations(args.annotations)

    # Merge
    annotated_df = merge_umap_with_cogs(umap_df, cog_dict)

    # Save
    print(f"\nSaving to {args.output}...")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    annotated_df.to_csv(args.output, index=False)

    print(f"✓ Saved {len(annotated_df):,} annotated genes")
    print("\nOutput columns:")
    for col in annotated_df.columns:
        print(f"  - {col}")

    print("\n" + "=" * 70)
    print("Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
