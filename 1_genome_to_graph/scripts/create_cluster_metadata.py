#!/usr/bin/env python3
"""
Create Cluster Metadata CSV

Generates comprehensive metadata for each gene cluster including:
- Top COG functional category (when diamond/eggnog annotations available)
- COG category distribution within cluster
- Cluster size
- Representative gene information
- Gene names from RefSeq (when available)

This creates a reusable metadata file that can be merged with any cluster-level analysis.

Usage:
    python create_cluster_metadata.py \
        --cluster-index graph_outputs/protein_graph/cluster_index.csv \
        --clusters intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
        --output graph_outputs/protein_graph/metadata/cluster_metadata.csv

    # If you have COG annotations from diamond:
    python create_cluster_metadata.py \
        --cluster-index graph_outputs/protein_graph/cluster_index.csv \
        --clusters intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
        --cog-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation \
        --output graph_outputs/protein_graph/metadata/cluster_metadata.csv

Output CSV columns:
    - cluster_id: Representative gene (e.g., GCF_000006985.1_NC_002932.3_1196)
    - cluster_size: Number of genes in cluster
    - top_cog_category: Most common COG category in cluster
    - top_cog_name: Name of top COG category
    - cog_counts: JSON dict of all COG categories and counts
    - annotation_rate: Fraction of cluster members with COG annotation
    - representative_genome: Genome accession of representative
    - representative_contig: Contig ID of representative
    - representative_gene: Gene number of representative
"""

import argparse
import pandas as pd
from pathlib import Path
from collections import Counter
import json
from tqdm import tqdm

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


def load_gene_to_cog_mapping(cog_dir):
    """
    Load COG annotations from diamond results.

    Returns:
        dict: {gene_id: cog_category}
    """
    cog_dir = Path(cog_dir)
    gene_to_cog = {}

    if not cog_dir.exists():
        print(f"Warning: COG directory not found: {cog_dir}")
        return gene_to_cog

    # Load COG database definitions
    cog_db = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/cog")

    prot_to_cog = {}
    cog_definitions = {}

    if cog_db.exists():
        # Load protein-to-COG mapping
        cog_csv = cog_db / "cog-20.cog.csv"
        if cog_csv.exists():
            print("Loading COG protein mapping...")
            with open(cog_csv, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    fields = line.strip().split(',')
                    if len(fields) >= 7:
                        prot_id = fields[2]
                        cog_id = fields[6]
                        prot_to_cog[prot_id] = cog_id

        # Load COG definitions
        cog_def = cog_db / "cog-20.def.tab"
        if cog_def.exists():
            print("Loading COG definitions...")
            with open(cog_def, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        cog_id = fields[0]
                        cog_category = fields[1]
                        cog_definitions[cog_id] = cog_category

    # Process diamond results for each genome
    diamond_dirs = list(cog_dir.glob("*_cog_diamond"))

    if not diamond_dirs:
        print(f"Warning: No diamond results found in {cog_dir}")
        return gene_to_cog

    print(f"Loading COG annotations from {len(diamond_dirs)} genomes...")

    for diamond_dir in tqdm(diamond_dirs):
        genome_id = diamond_dir.name.replace('_cog_diamond', '')
        hits_file = diamond_dir / f"{genome_id}_cog_hits.tsv"

        if not hits_file.exists():
            continue

        with open(hits_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue

                query_id = fields[0]   # Gene ID (e.g., NC_002932.3_1)
                subject_id = fields[1] # COG protein ID

                # Convert gene ID format: NC_002932.3_1 -> GCF_000006985.1_NC_002932.3_1
                # We need to prepend genome_id
                full_gene_id = f"{genome_id}_{query_id}"

                # Convert underscore to dot in subject ID for matching
                if '_' in subject_id:
                    parts = subject_id.rsplit('_', 1)
                    subject_id_dot = f"{parts[0]}.{parts[1]}"
                else:
                    subject_id_dot = subject_id

                # Map to COG ID and then to category
                if subject_id_dot in prot_to_cog:
                    cog_id = prot_to_cog[subject_id_dot]
                    if cog_id in cog_definitions:
                        categories = cog_definitions[cog_id]
                        # Use first category as primary
                        if categories and categories[0] != '-':
                            gene_to_cog[full_gene_id] = categories[0]

    print(f"  Loaded COG annotations for {len(gene_to_cog):,} genes")

    return gene_to_cog


def create_cluster_metadata(cluster_index_file, clusters_file, gene_to_cog=None, output_file=None):
    """
    Create comprehensive cluster metadata.

    Args:
        cluster_index_file: cluster_index.csv with cluster_id and cluster_size
        clusters_file: MMseqs2 clusters.tsv with representative -> member mappings
        gene_to_cog: Optional dict mapping gene IDs to COG categories
        output_file: Path to save metadata CSV
    """

    print("=" * 70)
    print("Creating Cluster Metadata")
    print("=" * 70)

    # Load cluster index
    print(f"\nLoading cluster index from: {cluster_index_file}")
    cluster_index = pd.read_csv(cluster_index_file)
    print(f"  Found {len(cluster_index):,} clusters")

    # Load cluster memberships
    print(f"\nLoading cluster memberships from: {clusters_file}")
    clusters_df = pd.read_csv(clusters_file, sep='\t', header=None,
                              names=['representative', 'member'])
    print(f"  Found {len(clusters_df):,} cluster assignments")

    # Group members by representative
    print("\nGrouping cluster members...")
    cluster_members = clusters_df.groupby('representative')['member'].apply(list).to_dict()

    # Create metadata for each cluster
    print("\nGenerating metadata for each cluster...")
    metadata = []

    for _, row in tqdm(cluster_index.iterrows(), total=len(cluster_index)):
        cluster_id = row['cluster_id']
        cluster_size = row['cluster_size']

        # Parse cluster ID
        parsed = parse_cluster_id(cluster_id)

        # Get cluster members
        members = cluster_members.get(cluster_id, [cluster_id])

        # Calculate COG statistics if available
        top_cog = None
        top_cog_name = None
        cog_counts = {}
        annotation_rate = 0.0

        if gene_to_cog:
            # Count COG categories in cluster
            cog_counter = Counter()
            annotated_count = 0

            for member in members:
                if member in gene_to_cog:
                    cog = gene_to_cog[member]
                    cog_counter[cog] += 1
                    annotated_count += 1

            # Get top COG
            if cog_counter:
                top_cog = cog_counter.most_common(1)[0][0]
                top_cog_name = COG_CATEGORIES.get(top_cog, 'Unknown')
                cog_counts = dict(cog_counter)

            # Calculate annotation rate
            annotation_rate = annotated_count / len(members) if members else 0.0

        metadata.append({
            'cluster_id': cluster_id,
            'cluster_size': cluster_size,
            'top_cog_category': top_cog,
            'top_cog_name': top_cog_name,
            'cog_counts': json.dumps(cog_counts) if cog_counts else None,
            'annotation_rate': annotation_rate,
            'representative_genome': parsed['genome'],
            'representative_contig': parsed['contig'],
            'representative_gene': parsed['gene_num'],
        })

    # Create DataFrame
    metadata_df = pd.DataFrame(metadata)

    # Print summary
    print("\n" + "=" * 70)
    print("Metadata Summary")
    print("=" * 70)
    print(f"\nTotal clusters: {len(metadata_df):,}")

    if gene_to_cog:
        annotated_clusters = metadata_df['top_cog_category'].notna().sum()
        print(f"Clusters with COG annotation: {annotated_clusters:,} ({100*annotated_clusters/len(metadata_df):.1f}%)")
        print(f"Mean annotation rate per cluster: {metadata_df['annotation_rate'].mean():.1%}")

        print("\nTop 10 COG categories:")
        top_cogs = metadata_df['top_cog_category'].value_counts().head(10)
        for cog, count in top_cogs.items():
            cog_name = COG_CATEGORIES.get(cog, 'Unknown')
            print(f"  {cog} ({cog_name[:50]:50s}): {count:6,} clusters ({100*count/len(metadata_df):5.1f}%)")

    # Save to file
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_df.to_csv(output_file, index=False)
        print(f"\n✓ Metadata saved to: {output_file}")
        print(f"  File size: {output_path.stat().st_size / 1024 / 1024:.2f} MB")

    return metadata_df


def main():
    parser = argparse.ArgumentParser(
        description='Create cluster metadata CSV with COG annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--cluster-index',
        required=True,
        help='Path to cluster_index.csv'
    )

    parser.add_argument(
        '--clusters',
        required=True,
        help='Path to MMseqs2 clusters.tsv'
    )

    parser.add_argument(
        '--cog-dir',
        help='Optional directory with diamond COG results (*_cog_diamond/ subdirs)'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output path for cluster_metadata.csv'
    )

    args = parser.parse_args()

    # Load COG annotations if provided
    gene_to_cog = None
    if args.cog_dir:
        gene_to_cog = load_gene_to_cog_mapping(args.cog_dir)
        if not gene_to_cog:
            print("\nWarning: No COG annotations loaded. Metadata will not include COG info.")
            print("Run diamond COG annotation first, or omit --cog-dir to skip COG metadata.")

    # Create metadata
    metadata_df = create_cluster_metadata(
        cluster_index_file=args.cluster_index,
        clusters_file=args.clusters,
        gene_to_cog=gene_to_cog,
        output_file=args.output
    )

    print("\n" + "=" * 70)
    print("Complete!")
    print("=" * 70)
    print("\nUsage:")
    print("  Load metadata:")
    print(f"    metadata_df = pd.read_csv('{args.output}')")
    print("\n  Merge with cluster embeddings:")
    print("    cluster_embeddings.merge(metadata_df, on='cluster_id', how='left')")
    print("\n  Filter to specific COG category:")
    print("    translation_clusters = metadata_df[metadata_df['top_cog_category'] == 'J']")


if __name__ == "__main__":
    main()
