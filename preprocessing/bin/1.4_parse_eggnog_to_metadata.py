#!/usr/bin/env python
"""
Parse eggNOG-mapper results into per-genome COG metadata files.

Creates reusable metadata files that can be merged with any analysis
(UMAP, clustering, differential analysis, etc.)

Usage:
    python parse_eggnog_to_metadata.py \
        --annotations results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/metadata/cog_annotations

Output:
    results/1_genome_to_graph/1.4_esm_embedding_clustering/metadata/cog_annotations/
        ├── GCF_000005845_cog.tsv
        ├── GCF_000006945_cog.tsv
        └── ... (one file per genome)

Each TSV has columns: gene_id, cog_category, cog_categories_all, cog_description, category_name
"""

import argparse
from pathlib import Path
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
        description="Parse eggNOG results to per-genome COG metadata"
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
        help="Output directory for COG metadata files",
    )
    return parser.parse_args()


def parse_eggnog_file(annot_file):
    """
    Parse a single eggNOG-mapper annotation file.

    Returns:
        DataFrame with columns: gene_id, cog_category, cog_categories_all,
                                cog_description, category_name
    """
    annotations = []

    with open(annot_file, 'r') as f:
        for line in f:
            # Skip headers and comments
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')

            if len(fields) < 7:
                continue

            # Parse fields from eggNOG-mapper output
            gene_id = fields[0]              # Query gene ID
            cog_category = fields[6]         # COG category letters (e.g., "J" or "KL")
            description = fields[7] if len(fields) > 7 else ""

            # Handle unannotated genes
            if cog_category == '-' or not cog_category:
                annotations.append({
                    'gene_id': gene_id,
                    'cog_category': None,
                    'cog_categories_all': None,
                    'cog_description': description,
                    'category_name': None
                })
            else:
                # Get primary category (first letter)
                primary_cog = cog_category[0]

                # Get category name
                category_name = COG_CATEGORIES.get(primary_cog, 'Unknown')

                annotations.append({
                    'gene_id': gene_id,
                    'cog_category': primary_cog,
                    'cog_categories_all': cog_category,  # Keep all categories (e.g., "KL")
                    'cog_description': description,
                    'category_name': category_name
                })

    return pd.DataFrame(annotations)


def process_all_genomes(annotations_dir, output_dir):
    """Process all eggNOG-mapper results into per-genome metadata files."""

    annotations_dir = Path(annotations_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all annotation files
    annotation_files = list(annotations_dir.glob("*_eggnog/*.emapper.annotations"))

    print(f"\nProcessing {len(annotation_files)} genomes...")
    print("=" * 70)

    total_genes = 0
    total_annotated = 0
    processed = 0

    for annot_file in sorted(annotation_files):
        genome_id = annot_file.parent.name.replace('_eggnog', '')

        try:
            # Parse eggNOG results
            df = parse_eggnog_file(annot_file)

            if len(df) == 0:
                print(f"⚠ {genome_id}: No genes found")
                continue

            # Count annotations
            n_genes = len(df)
            n_annotated = df['cog_category'].notna().sum()
            pct = 100 * n_annotated / n_genes if n_genes > 0 else 0

            # Save metadata
            output_file = output_dir / f"{genome_id}_cog.tsv"
            df.to_csv(output_file, sep='\t', index=False)

            # Update totals
            total_genes += n_genes
            total_annotated += n_annotated
            processed += 1

            print(f"✓ {genome_id:25s} {n_genes:5,} genes, {n_annotated:5,} annotated ({pct:5.1f}%)")

        except Exception as e:
            print(f"✗ {genome_id}: Error - {e}")
            continue

    print("=" * 70)
    print(f"\nProcessed {processed} genomes")
    print(f"Total genes: {total_genes:,}")
    print(f"Total annotated: {total_annotated:,} ({100*total_annotated/total_genes:.1f}%)")
    print(f"\nMetadata files saved to: {output_dir}/")

    return processed


def create_summary(output_dir):
    """Create a summary file with annotation statistics per genome."""

    output_dir = Path(output_dir)
    metadata_files = list(output_dir.glob("*_cog.tsv"))

    if not metadata_files:
        return

    print("\nCreating summary statistics...")

    summary = []

    for metadata_file in sorted(metadata_files):
        genome_id = metadata_file.name.replace('_cog.tsv', '')

        df = pd.read_csv(metadata_file, sep='\t')

        n_genes = len(df)
        n_annotated = df['cog_category'].notna().sum()

        # Count categories
        category_counts = df['cog_category'].value_counts()

        summary.append({
            'genome_id': genome_id,
            'total_genes': n_genes,
            'annotated_genes': n_annotated,
            'annotation_rate': n_annotated / n_genes if n_genes > 0 else 0,
            'n_categories': df['cog_category'].nunique() - 1,  # Exclude NA
            'most_common_category': category_counts.index[0] if len(category_counts) > 0 else None,
        })

    summary_df = pd.DataFrame(summary)
    summary_file = output_dir / "annotation_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)

    print(f"✓ Saved summary to: {summary_file}")

    # Print overall stats
    print("\nOverall Statistics:")
    print(f"  Mean annotation rate: {summary_df['annotation_rate'].mean():.1%}")
    print(f"  Median annotation rate: {summary_df['annotation_rate'].median():.1%}")
    print(f"  Min annotation rate: {summary_df['annotation_rate'].min():.1%}")
    print(f"  Max annotation rate: {summary_df['annotation_rate'].max():.1%}")


def main():
    args = parse_args()

    print("=" * 70)
    print("Parse eggNOG Results to COG Metadata")
    print("=" * 70)

    # Process all genomes
    n_processed = process_all_genomes(args.annotations, args.output)

    if n_processed > 0:
        # Create summary
        create_summary(args.output)

    print("\n" + "=" * 70)
    print("Complete!")
    print("=" * 70)
    print(f"\nOutput: {args.output}/")
    print("\nUsage:")
    print("  Load metadata for a genome:")
    print(f"    pd.read_csv('{args.output}/GCF_000005845_cog.tsv', sep='\\t')")
    print("\n  Merge with any analysis:")
    print("    umap_df.merge(cog_df, on='gene_id', how='left')")


if __name__ == "__main__":
    main()
