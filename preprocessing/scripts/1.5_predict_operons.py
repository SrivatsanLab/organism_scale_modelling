#!/usr/bin/env python3
"""
Predict operons from bacterial genomes using gene positions and distances.

This script predicts operons based on:
1. Gene proximity (intergenic distance < threshold, typically 50-150bp)
2. Same strand orientation
3. Prodigal gene predictions (coordinates from GFF/GBK files)

Operons define protein-protein adjacency in the graph - proteins in the same
operon are connected, reflecting their coordinated expression and function.

Input:
  - Genome sequences (FASTA format)
  - Gene predictions from Prodigal (regenerate if needed)

Output:
  - graph_outputs/protein_graph/operon_edges.csv (protein pairs in same operon)
  - intermediate/protein/operons/operon_predictions.tsv (per-genome operon assignments)
"""

import argparse
from pathlib import Path
import subprocess
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import re


def run_prodigal(genome_fasta, output_dir, genome_id):
    """
    Run Prodigal to predict genes and get coordinates.

    Returns:
        Path to GFF file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gff_file = output_dir / f"{genome_id}_genes.gff"
    faa_file = output_dir / f"{genome_id}_genes.faa"

    # Skip if already exists
    if gff_file.exists() and faa_file.exists():
        return gff_file

    cmd = [
        'prodigal',
        '-i', str(genome_fasta),
        '-f', 'gff',
        '-o', str(gff_file),
        '-a', str(faa_file),
        '-q'  # Quiet mode
    ]

    subprocess.run(cmd, check=True, capture_output=True)

    return gff_file


def parse_gff(gff_file, genome_id):
    """
    Parse Prodigal GFF file to extract gene coordinates.

    Returns:
        DataFrame with columns: protein_id, contig, start, end, strand
    """
    genes = []

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Extract ID from attributes (9th column)
            # Format: ID=1_1;partial=00;start_type=ATG;...
            attrs = parts[8]
            id_match = re.search(r'ID=([^;]+)', attrs)
            if id_match:
                gene_id = id_match.group(1)
                # Create protein ID in same format as other components
                # Format: GCF_XXXXXX.X_contig_genenum
                protein_id = f"{genome_id}_{contig}_{gene_id.split('_')[-1]}"

                genes.append({
                    'protein_id': protein_id,
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_num': int(gene_id.split('_')[-1])
                })

    return pd.DataFrame(genes)


def predict_operons(genes_df, max_intergenic_distance=150):
    """
    Predict operons based on gene proximity and strand.

    Args:
        genes_df: DataFrame with gene coordinates
        max_intergenic_distance: Maximum distance (bp) between genes in an operon

    Returns:
        DataFrame with operon assignments
    """
    # Sort by contig and position
    genes_df = genes_df.sort_values(['contig', 'start']).reset_index(drop=True)

    # Assign operon IDs
    operon_id = 1
    operon_assignments = []

    for contig in genes_df['contig'].unique():
        contig_genes = genes_df[genes_df['contig'] == contig].copy()

        current_operon = None
        prev_gene = None

        for idx, gene in contig_genes.iterrows():
            if prev_gene is None:
                # First gene on contig
                current_operon = operon_id
                operon_id += 1
            else:
                # Check if gene should be in same operon as previous
                same_strand = gene['strand'] == prev_gene['strand']
                distance = gene['start'] - prev_gene['end']

                if same_strand and 0 <= distance <= max_intergenic_distance:
                    # Same operon
                    pass
                else:
                    # New operon
                    current_operon = operon_id
                    operon_id += 1

            operon_assignments.append({
                'protein_id': gene['protein_id'],
                'contig': gene['contig'],
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'operon_id': current_operon,
                'intergenic_dist': gene['start'] - prev_gene['end'] if prev_gene is not None else None
            })

            prev_gene = gene

    return pd.DataFrame(operon_assignments)


def create_operon_edges(operon_df):
    """
    Create protein-protein edges for proteins in the same operon.

    Returns:
        DataFrame with columns: protein_id_1, protein_id_2, operon_id
    """
    edges = []

    for operon_id, group in operon_df.groupby('operon_id'):
        proteins = group['protein_id'].tolist()

        # Create edges between adjacent genes in operon
        for i in range(len(proteins) - 1):
            edges.append({
                'protein_id_1': proteins[i],
                'protein_id_2': proteins[i + 1],
                'operon_id': operon_id,
                'edge_type': 'operon_adjacency'
            })

    return pd.DataFrame(edges)


def process_genome(genome_fasta, genome_id, prodigal_output_dir, max_intergenic_distance):
    """
    Process a single genome to predict operons.

    Returns:
        tuple: (operon_df, edges_df)
    """
    # Run Prodigal
    gff_file = run_prodigal(genome_fasta, prodigal_output_dir, genome_id)

    # Parse genes
    genes_df = parse_gff(gff_file, genome_id)

    # Predict operons
    operon_df = predict_operons(genes_df, max_intergenic_distance)
    operon_df['genome_id'] = genome_id

    # Create edges
    edges_df = create_operon_edges(operon_df)

    return operon_df, edges_df


def main():
    parser = argparse.ArgumentParser(
        description='Predict operons from bacterial genomes'
    )
    parser.add_argument(
        '--genomes-dir',
        type=str,
        default='data/refseq_genomes',
        help='Directory containing genome FASTA files'
    )
    parser.add_argument(
        '--genome-list',
        type=str,
        default='data/genome_metadata.csv',
        help='CSV file with genome IDs (column: genome_id)'
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        default=150,
        help='Maximum intergenic distance (bp) for operon prediction'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph',
        help='Output directory for operon edges'
    )
    parser.add_argument(
        '--intermediate-dir',
        type=str,
        default='1_genome_to_graph/intermediate/protein/operons',
        help='Directory for intermediate operon predictions'
    )
    parser.add_argument(
        '--prodigal-dir',
        type=str,
        default='1_genome_to_graph/intermediate/protein/gene_predictions',
        help='Directory for Prodigal outputs (GFF files)'
    )

    args = parser.parse_args()

    print("="*70)
    print("Operon Prediction")
    print("="*70)
    print(f"Max intergenic distance: {args.max_distance} bp")

    # Load genome list
    genome_df = pd.read_csv(args.genome_list)
    genome_ids = genome_df['genome_id'].tolist()

    print(f"\nProcessing {len(genome_ids):,} genomes...")

    # Process each genome
    all_operons = []
    all_edges = []

    genomes_dir = Path(args.genomes_dir)

    for genome_id in tqdm(genome_ids, desc="Predicting operons"):
        # Try both .fasta and .fna extensions
        genome_fasta = genomes_dir / f"{genome_id}.fasta"
        if not genome_fasta.exists():
            genome_fasta = genomes_dir / f"{genome_id}.fna"

        if not genome_fasta.exists():
            print(f"Warning: Genome file not found: {genome_id}")
            continue

        try:
            operon_df, edges_df = process_genome(
                genome_fasta,
                genome_id,
                args.prodigal_dir,
                args.max_distance
            )

            all_operons.append(operon_df)
            all_edges.append(edges_df)

        except Exception as e:
            print(f"Error processing {genome_id}: {e}")
            continue

    # Combine results
    print("\nCombining results...")
    operons_combined = pd.concat(all_operons, ignore_index=True)
    edges_combined = pd.concat(all_edges, ignore_index=True)

    # Save outputs
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    intermediate_dir = Path(args.intermediate_dir)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Save edges (for graph construction)
    edges_file = output_dir / 'operon_edges.csv'
    print(f"\nSaving operon edges to {edges_file}...")
    edges_combined.to_csv(edges_file, index=False)

    # Save full operon predictions (for analysis)
    operon_file = intermediate_dir / 'operon_predictions.tsv'
    print(f"Saving operon predictions to {operon_file}...")
    operons_combined.to_csv(operon_file, sep='\t', index=False)

    # Print summary stats
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total genomes processed: {len(all_operons):,}")
    print(f"Total genes: {len(operons_combined):,}")
    print(f"Total operons predicted: {operons_combined['operon_id'].nunique():,}")
    print(f"Mean operon size: {operons_combined.groupby('operon_id').size().mean():.1f} genes")
    print(f"Median operon size: {operons_combined.groupby('operon_id').size().median():.0f} genes")
    print(f"\nOperon edges created: {len(edges_combined):,}")
    print(f"Mean intergenic distance: {operons_combined['intergenic_dist'].dropna().mean():.1f} bp")
    print("="*70)

    print("\nDone!")


if __name__ == '__main__':
    main()
