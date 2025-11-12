#!/usr/bin/env python3
"""
Parse COG Categories from eggNOG-mapper Results

Version: 1.0.0
Author: ShotgunDomestication Project
Description: Extract and summarize COG functional categories from eggNOG-mapper annotations

Usage: python parse_cog_categories.py <results_dir> -o <output.csv>
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import Counter
import re

SCRIPT_VERSION = "1.0.0"

# COG category definitions
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
}


def parse_emapper_annotations(annot_file):
    """
    Parse eggNOG-mapper annotation file and extract COG categories

    Returns:
        dict: Counts of proteins per COG category
        int: Total proteins
        int: Annotated proteins
    """
    cog_counts = Counter()
    total_proteins = 0
    annotated_proteins = 0

    with open(annot_file, 'r') as f:
        for line in f:
            # Skip header and comment lines
            if line.startswith('#') or line.strip() == '':
                continue

            total_proteins += 1

            fields = line.strip().split('\t')

            # COG category is column 6 (0-indexed)
            if len(fields) > 6:
                cog_str = fields[6]

                # Check if has annotation (not '-' or empty)
                if cog_str and cog_str != '-':
                    annotated_proteins += 1

                    # COG categories can be multiple letters (e.g., "COG")
                    # Count each letter separately
                    for cog_letter in cog_str:
                        if cog_letter in COG_CATEGORIES:
                            cog_counts[cog_letter] += 1

    return cog_counts, total_proteins, annotated_proteins


def process_all_samples(results_dir, output_file):
    """
    Process all eggNOG-mapper results in a directory

    Args:
        results_dir: Directory containing sample_eggnog subdirectories
        output_file: Path to output CSV file
    """
    results_path = Path(results_dir)

    # Find all annotation files
    annot_files = list(results_path.glob('*_eggnog/*.emapper.annotations'))

    if not annot_files:
        print(f"Error: No annotation files found in {results_dir}")
        sys.exit(1)

    print(f"Found {len(annot_files)} annotation files")
    print()

    # Initialize results list
    results = []

    # Process each sample
    for i, annot_file in enumerate(annot_files, 1):
        sample_name = annot_file.parent.name.replace('_eggnog', '')

        print(f"[{i}/{len(annot_files)}] Processing {sample_name}...", end=' ')

        try:
            cog_counts, total_proteins, annotated_proteins = parse_emapper_annotations(annot_file)

            # Calculate annotation rate
            annotation_rate = (annotated_proteins / total_proteins * 100) if total_proteins > 0 else 0

            # Create result dictionary
            result = {
                'sample_name': sample_name,
                'total_proteins': total_proteins,
                'annotated_proteins': annotated_proteins,
                'annotation_rate': annotation_rate,
            }

            # Add COG category counts
            for cog in sorted(COG_CATEGORIES.keys()):
                result[f'COG_{cog}'] = cog_counts.get(cog, 0)

            # Add unannotated count
            result['unannotated'] = total_proteins - annotated_proteins

            results.append(result)

            print(f"✓ ({total_proteins} proteins, {annotation_rate:.1f}% annotated)")

        except Exception as e:
            print(f"✗ Error: {e}")
            continue

    print()

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Sort by sample name
    df = df.sort_values('sample_name')

    # Write to CSV
    df.to_csv(output_file, index=False)

    print(f"Results written to: {output_file}")
    print()

    # Print summary statistics
    print("=" * 70)
    print("Summary Statistics")
    print("=" * 70)
    print()
    print(f"Total samples: {len(df)}")
    print(f"Total proteins: {df['total_proteins'].sum():,}")
    print(f"Annotated proteins: {df['annotated_proteins'].sum():,}")
    print(f"Mean annotation rate: {df['annotation_rate'].mean():.1f}%")
    print()

    print("Mean COG category distribution:")
    cog_cols = [c for c in df.columns if c.startswith('COG_')]
    for cog_col in cog_cols:
        cog_letter = cog_col.split('_')[1]
        mean_count = df[cog_col].mean()
        mean_pct = (df[cog_col].sum() / df['total_proteins'].sum() * 100)
        print(f"  {cog_letter} ({COG_CATEGORIES[cog_letter][:50]:50s}): {mean_count:6.1f} ({mean_pct:5.1f}%)")
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Parse COG categories from eggNOG-mapper results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Example usage:
  python {sys.argv[0]} results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation -o cog_summary.csv

Version: {SCRIPT_VERSION}
        """
    )

    parser.add_argument(
        'results_dir',
        help='Directory containing eggNOG-mapper results (sample_eggnog subdirectories)'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file path'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {SCRIPT_VERSION}'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.results_dir):
        print(f"Error: Results directory not found: {args.results_dir}")
        sys.exit(1)

    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("=" * 70)
    print("Parse COG Categories from eggNOG-mapper")
    print(f"Version: {SCRIPT_VERSION}")
    print("=" * 70)
    print()

    # Process samples
    process_all_samples(args.results_dir, args.output)

    print("Done!")


if __name__ == '__main__':
    main()
