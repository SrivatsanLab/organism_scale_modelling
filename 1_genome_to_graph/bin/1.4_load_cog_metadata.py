#!/usr/bin/env python
"""
Utility functions to load and merge COG metadata with any analysis.

This module provides helper functions to work with per-genome COG metadata files.

Example Usage:

    from load_cog_metadata import load_all_cog_metadata, merge_with_cog

    # Load all COG metadata
    cog_df = load_all_cog_metadata('results/1_genome_to_graph/1.4_esm_embedding_clustering/metadata/cog_annotations')

    # Merge with UMAP or any other analysis
    umap_df = pd.DataFrame({
        'gene_id': ['GCF_000005845_00001', 'GCF_000006945_00523'],
        'umap_1': [1.2, 3.4],
        'umap_2': [5.6, 7.8]
    })

    annotated_df = merge_with_cog(umap_df, cog_df)
"""

from pathlib import Path
import pandas as pd
from typing import Union, List, Optional


def load_genome_cog_metadata(metadata_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load COG metadata for a single genome.

    Args:
        metadata_file: Path to *_cog.tsv file

    Returns:
        DataFrame with columns: gene_id, cog_category, cog_categories_all,
                                cog_description, category_name
    """
    return pd.read_csv(metadata_file, sep='\t')


def load_all_cog_metadata(metadata_dir: Union[str, Path]) -> pd.DataFrame:
    """
    Load and concatenate all COG metadata files.

    Args:
        metadata_dir: Directory containing *_cog.tsv files

    Returns:
        DataFrame with all genes from all genomes
    """
    metadata_dir = Path(metadata_dir)
    metadata_files = sorted(metadata_dir.glob("*_cog.tsv"))

    if not metadata_files:
        raise FileNotFoundError(f"No COG metadata files found in {metadata_dir}")

    print(f"Loading COG metadata from {len(metadata_files)} genomes...")

    dfs = []
    for metadata_file in metadata_files:
        df = pd.read_csv(metadata_file, sep='\t')
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)

    print(f"  Loaded {len(combined_df):,} genes")
    print(f"  Annotated: {combined_df['cog_category'].notna().sum():,}")

    return combined_df


def load_cog_for_genomes(metadata_dir: Union[str, Path],
                         genome_ids: List[str]) -> pd.DataFrame:
    """
    Load COG metadata for specific genomes only.

    Args:
        metadata_dir: Directory containing *_cog.tsv files
        genome_ids: List of genome IDs to load

    Returns:
        DataFrame with genes from specified genomes
    """
    metadata_dir = Path(metadata_dir)

    dfs = []
    for genome_id in genome_ids:
        metadata_file = metadata_dir / f"{genome_id}_cog.tsv"

        if not metadata_file.exists():
            print(f"Warning: {genome_id} metadata not found")
            continue

        df = pd.read_csv(metadata_file, sep='\t')
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError(f"No metadata found for genomes: {genome_ids}")

    return pd.concat(dfs, ignore_index=True)


def merge_with_cog(data_df: pd.DataFrame,
                   cog_df: pd.DataFrame,
                   on: str = 'gene_id',
                   how: str = 'left') -> pd.DataFrame:
    """
    Merge any analysis DataFrame with COG annotations.

    Args:
        data_df: DataFrame with gene_id column (e.g., UMAP, clusters, DE results)
        cog_df: COG metadata DataFrame (from load_all_cog_metadata)
        on: Column to merge on (default: 'gene_id')
        how: Merge type (default: 'left' to keep all genes from data_df)

    Returns:
        Merged DataFrame with COG annotations
    """
    merged_df = data_df.merge(cog_df, on=on, how=how)

    # Fill missing annotations
    merged_df['cog_category'] = merged_df['cog_category'].fillna('No COG')
    merged_df['category_name'] = merged_df['category_name'].fillna('No annotation')

    n_annotated = (merged_df['cog_category'] != 'No COG').sum()
    n_total = len(merged_df)

    print(f"\nMerge summary:")
    print(f"  Total genes: {n_total:,}")
    print(f"  Annotated: {n_annotated:,} ({100*n_annotated/n_total:.1f}%)")

    return merged_df


def filter_by_cog_category(data_df: pd.DataFrame,
                           categories: Union[str, List[str]]) -> pd.DataFrame:
    """
    Filter DataFrame to genes with specific COG categories.

    Args:
        data_df: DataFrame with cog_category column
        categories: Single category ('J') or list of categories (['J', 'K', 'L'])

    Returns:
        Filtered DataFrame
    """
    if isinstance(categories, str):
        categories = [categories]

    filtered_df = data_df[data_df['cog_category'].isin(categories)]

    print(f"Filtered to {len(filtered_df):,} genes in categories: {', '.join(categories)}")

    return filtered_df


def get_cog_summary(data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get summary statistics of COG categories in a DataFrame.

    Args:
        data_df: DataFrame with cog_category column

    Returns:
        DataFrame with columns: cog_category, count, percentage, category_name
    """
    category_counts = data_df['cog_category'].value_counts()

    summary = []
    for category, count in category_counts.items():
        # Get category name
        if category == 'No COG':
            cat_name = 'No annotation'
        else:
            cat_name = data_df[data_df['cog_category'] == category]['category_name'].iloc[0]

        summary.append({
            'cog_category': category,
            'count': count,
            'percentage': 100 * count / len(data_df),
            'category_name': cat_name
        })

    summary_df = pd.DataFrame(summary)
    summary_df = summary_df.sort_values('count', ascending=False)

    return summary_df


# Example usage
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python load_cog_metadata.py <metadata_dir>")
        print("\nExample:")
        print("  python load_cog_metadata.py results/1_genome_to_graph/1.4_esm_embedding_clustering/metadata/cog_annotations")
        sys.exit(1)

    metadata_dir = sys.argv[1]

    print("=" * 70)
    print("COG Metadata Loader - Example Usage")
    print("=" * 70)

    # Load all metadata
    cog_df = load_all_cog_metadata(metadata_dir)

    # Show summary
    print("\nCOG Category Distribution:")
    summary = get_cog_summary(cog_df)
    print(summary.to_string(index=False))

    # Example: Filter to translation genes
    print("\n" + "=" * 70)
    print("Example: Filter to Translation genes (category J)")
    print("=" * 70)
    translation_genes = filter_by_cog_category(cog_df, 'J')
    print(f"\nFound {len(translation_genes):,} translation genes")
    print(translation_genes.head(10).to_string(index=False))

    print("\n" + "=" * 70)
    print("Example complete!")
    print("=" * 70)
