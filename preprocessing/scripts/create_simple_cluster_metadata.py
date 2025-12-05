#!/usr/bin/env python3
"""
Create Simple Cluster Metadata (fast version)

Just parse cluster IDs and add basic information without loading the full clusters.tsv.
This is much faster than the full version.

Usage:
    python create_simple_cluster_metadata.py \
        --cluster-index graph_outputs/protein_graph/cluster_index.csv \
        --output graph_outputs/protein_graph/metadata/cluster_metadata_simple.csv
"""

import argparse
import pandas as pd
from pathlib import Path

def parse_cluster_id(cluster_id):
    """
    Parse cluster ID into components.

    Format: {genome_accession}_{contig_id}_{gene_number}
    Example: GCF_000006985.1_NC_002932.3_1196

    Returns:
        dict with genome, contig, gene_num
    """
    parts = cluster_id.split('_')

    # Handle case where contig ID has underscores
    # Genome is GCF_XXXXXX.X (first 2 parts)
    # Gene number is last part
    # Everything in between is contig

    genome = '_'.join(parts[:2])  # GCF_000006985.1
    gene_num = parts[-1]           # 1196
    contig = '_'.join(parts[2:-1]) # NC_002932.3

    return {
        'genome': genome,
        'contig': contig,
        'gene_num': gene_num
    }


def create_simple_metadata(cluster_index_file, output_file):
    """Create simple metadata by just parsing cluster IDs."""

    print("=" * 70)
    print("Creating Simple Cluster Metadata")
    print("=" * 70)

    # Load cluster index
    print(f"\nLoading cluster index: {cluster_index_file}")
    cluster_index = pd.read_csv(cluster_index_file)
    print(f"  Found {len(cluster_index):,} clusters")

    # Parse cluster IDs
    print("\nParsing cluster IDs...")
    metadata = []

    for _, row in cluster_index.iterrows():
        cluster_id = row['cluster_id']
        cluster_size = row['cluster_size']

        # Parse ID
        parsed = parse_cluster_id(cluster_id)

        metadata.append({
            'cluster_id': cluster_id,
            'cluster_size': cluster_size,
            'representative_genome': parsed['genome'],
            'representative_contig': parsed['contig'],
            'representative_gene': parsed['gene_num'],
        })

    # Create DataFrame
    metadata_df = pd.DataFrame(metadata)

    # Save
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_df.to_csv(output_file, index=False)

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Total clusters: {len(metadata_df):,}")
    print(f"Unique genomes: {metadata_df['representative_genome'].nunique():,}")
    print(f"Unique contigs: {metadata_df['representative_contig'].nunique():,}")

    print(f"\n✓ Metadata saved to: {output_file}")
    print(f"  File size: {output_path.stat().st_size / 1024:.2f} KB")

    print("\nUsage:")
    print("  metadata_df = pd.read_csv('{}')".format(output_file))
    print("  cluster_data.merge(metadata_df, on='cluster_id')")

    return metadata_df


def main():
    parser = argparse.ArgumentParser(
        description='Create simple cluster metadata (fast version)'
    )

    parser.add_argument(
        '--cluster-index',
        required=True,
        help='Path to cluster_index.csv'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output path for metadata CSV'
    )

    args = parser.parse_args()

    create_simple_metadata(args.cluster_index, args.output)


if __name__ == "__main__":
    main()
