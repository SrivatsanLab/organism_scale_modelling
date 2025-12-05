#!/usr/bin/env python3
"""
Create gene-gene (cluster-cluster) adjacency matrix based on operon predictions.

This script:
1. Predicts operons from bacterial genomes using gene proximity
2. Maps genes to protein clusters
3. Creates cluster-cluster edges when genes from two clusters are in operons together
4. Computes edge weights based on co-occurrence frequency

Input:
  - Genome sequences (FASTA format)
  - MMseqs2 cluster assignments
  - Cluster statistics (filtered clusters to use)

Output:
  - graph_outputs/protein_graph/gene_gene_edges.csv - cluster pair edges with weights
  - graph_outputs/protein_graph/gene_gene_matrix.npz - sparse adjacency matrix
  - intermediate/protein/operons/operon_predictions.tsv - per-genome operon assignments
"""

import argparse
from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
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
            attrs = parts[8]
            id_match = re.search(r'ID=([^;]+)', attrs)
            if id_match:
                gene_id = id_match.group(1)
                # Create protein ID matching cluster file format
                # Format: GCF_XXXXXX.X_contig_genenum
                protein_id = f"{genome_id}_{contig}_{gene_id.split('_')[-1]}"

                genes.append({
                    'protein_id': protein_id,
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand
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
                'operon_id': current_operon
            })

            prev_gene = gene

    return pd.DataFrame(operon_assignments)


def load_protein_clusters(cluster_file, used_clusters):
    """
    Load protein -> cluster mapping for used clusters only.

    Returns:
        dict: protein_id -> cluster_id
    """
    print(f"\nLoading cluster assignments from {cluster_file}")

    protein_to_cluster = {}

    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Reading clusters"):
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue

            cluster_rep, protein_id = parts

            # Only keep if this cluster is in the used set
            if cluster_rep in used_clusters:
                protein_to_cluster[protein_id] = cluster_rep

    print(f"Loaded {len(protein_to_cluster):,} proteins in used clusters")

    return protein_to_cluster


def create_cluster_edges_from_operons(operon_df, protein_to_cluster, min_cooccurrence=1):
    """
    Create cluster-cluster edges from operon predictions.

    Two clusters are connected if genes from those clusters appear in operons together
    across genomes. Edge weight = number of times they co-occur in operons.

    Args:
        operon_df: DataFrame with operon assignments (protein_id, operon_id, genome_id)
        protein_to_cluster: Dict mapping protein_id -> cluster_id
        min_cooccurrence: Minimum co-occurrence count to create an edge

    Returns:
        DataFrame: cluster_id_1, cluster_id_2, weight
    """
    print("\nCreating cluster-cluster edges from operons...")

    # Map proteins to clusters
    operon_df['cluster_id'] = operon_df['protein_id'].map(protein_to_cluster)

    # Remove proteins not in clusters
    operon_df = operon_df.dropna(subset=['cluster_id'])

    print(f"Proteins mapped to clusters: {len(operon_df):,}")

    # Count cluster co-occurrences in operons
    edge_counts = defaultdict(int)

    for (genome_id, operon_id), group in tqdm(
        operon_df.groupby(['genome_id', 'operon_id']),
        desc="Processing operons"
    ):
        clusters = group['cluster_id'].unique()

        # Skip single-cluster operons
        if len(clusters) < 2:
            continue

        # Create edges between all cluster pairs in this operon
        for i, c1 in enumerate(clusters):
            for c2 in clusters[i+1:]:
                # Create canonical edge (sorted)
                edge = tuple(sorted([c1, c2]))
                edge_counts[edge] += 1

    # Filter by minimum co-occurrence
    edges = [
        {'cluster_id_1': edge[0], 'cluster_id_2': edge[1], 'weight': count}
        for edge, count in edge_counts.items()
        if count >= min_cooccurrence
    ]

    edge_df = pd.DataFrame(edges)

    print(f"\nCluster-cluster edges created: {len(edge_df):,}")
    print(f"Unique clusters with edges: {len(set(edge_df['cluster_id_1']) | set(edge_df['cluster_id_2'])):,}")

    return edge_df


def create_adjacency_matrix(edge_df, cluster_list):
    """
    Create sparse adjacency matrix from edge list.

    Args:
        edge_df: DataFrame with cluster_id_1, cluster_id_2, weight
        cluster_list: Ordered list of all cluster IDs

    Returns:
        csr_matrix: Symmetric sparse matrix (n_clusters x n_clusters)
    """
    print("\nCreating adjacency matrix...")

    cluster_to_idx = {c: i for i, c in enumerate(cluster_list)}
    n_clusters = len(cluster_list)

    print(f"Matrix dimensions: {n_clusters:,} x {n_clusters:,}")

    # Build sparse matrix (symmetric)
    row_indices = []
    col_indices = []
    weights = []

    for _, row in tqdm(edge_df.iterrows(), total=len(edge_df), desc="Building matrix"):
        c1 = row['cluster_id_1']
        c2 = row['cluster_id_2']
        weight = row['weight']

        if c1 not in cluster_to_idx or c2 not in cluster_to_idx:
            continue

        idx1 = cluster_to_idx[c1]
        idx2 = cluster_to_idx[c2]

        # Add both directions (symmetric)
        row_indices.extend([idx1, idx2])
        col_indices.extend([idx2, idx1])
        weights.extend([weight, weight])

    # Create sparse matrix
    matrix = csr_matrix(
        (weights, (row_indices, col_indices)),
        shape=(n_clusters, n_clusters),
        dtype=np.float32
    )

    print(f"Total edges (undirected): {len(edge_df):,}")
    print(f"Matrix sparsity: {100 * (1 - matrix.nnz / (n_clusters * n_clusters)):.4f}%")

    return matrix


def process_genome(genome_fasta, genome_id, prodigal_output_dir, max_intergenic_distance):
    """
    Process a single genome to predict operons.

    Returns:
        DataFrame: operon assignments
    """
    # Run Prodigal
    gff_file = run_prodigal(genome_fasta, prodigal_output_dir, genome_id)

    # Parse genes
    genes_df = parse_gff(gff_file, genome_id)

    if len(genes_df) == 0:
        return pd.DataFrame()

    # Predict operons
    operon_df = predict_operons(genes_df, max_intergenic_distance)
    operon_df['genome_id'] = genome_id

    return operon_df


def main():
    parser = argparse.ArgumentParser(
        description='Create gene-gene (cluster-cluster) edges from operon predictions'
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
        '--clusters',
        type=str,
        default='1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv',
        help='MMseqs2 cluster TSV file'
    )
    parser.add_argument(
        '--cluster-stats',
        type=str,
        default='1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv',
        help='Cluster statistics CSV (filtered clusters to use)'
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        default=150,
        help='Maximum intergenic distance (bp) for operon prediction'
    )
    parser.add_argument(
        '--min-cooccurrence',
        type=int,
        default=1,
        help='Minimum co-occurrence count for cluster-cluster edge'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph',
        help='Output directory for edges'
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
    print("Gene-Gene Edge Prediction (Operon-Based)")
    print("="*70)
    print(f"Max intergenic distance: {args.max_distance} bp")
    print(f"Min co-occurrence: {args.min_cooccurrence}")

    # Load cluster list
    print(f"\nLoading cluster list from {args.cluster_stats}")
    cluster_df = pd.read_csv(args.cluster_stats)

    # Handle both cluster_statistics.csv (cluster_representative column)
    # and cluster_index.csv (cluster_id column) formats
    if 'cluster_representative' in cluster_df.columns:
        used_clusters = set(cluster_df['cluster_representative'].values)
    elif 'cluster_id' in cluster_df.columns:
        used_clusters = set(cluster_df['cluster_id'].values)
    else:
        raise ValueError(f"Cluster file must have 'cluster_representative' or 'cluster_id' column")

    print(f"Using {len(used_clusters):,} clusters")

    # Load protein -> cluster mapping
    protein_to_cluster = load_protein_clusters(args.clusters, used_clusters)

    # Load genome list
    genome_df = pd.read_csv(args.genome_list)
    genome_ids = genome_df['genome_id'].tolist()

    print(f"\nProcessing {len(genome_ids):,} genomes...")

    # Process each genome
    all_operons = []
    genomes_dir = Path(args.genomes_dir)

    for genome_id in tqdm(genome_ids, desc="Predicting operons"):
        # Find genome file (files have full name like GCF_000006985.1_Organism_name.fasta)
        genome_files = list(genomes_dir.glob(f"{genome_id}*.fasta"))
        if not genome_files:
            genome_files = list(genomes_dir.glob(f"{genome_id}*.fna"))

        if not genome_files:
            continue

        genome_fasta = genome_files[0]  # Take first match

        try:
            operon_df = process_genome(
                genome_fasta,
                genome_id,
                args.prodigal_dir,
                args.max_distance
            )

            if len(operon_df) > 0:
                all_operons.append(operon_df)

        except Exception as e:
            print(f"Error processing {genome_id}: {e}")
            continue

    # Combine results
    print("\nCombining operon predictions...")

    if len(all_operons) == 0:
        print("ERROR: No operons were predicted from any genome!")
        print("This likely means no genome files were found or all predictions failed.")
        return

    operons_combined = pd.concat(all_operons, ignore_index=True)

    print(f"Total genes: {len(operons_combined):,}")
    print(f"Total operons: {operons_combined['operon_id'].nunique():,}")

    # Create cluster-cluster edges
    edge_df = create_cluster_edges_from_operons(
        operons_combined,
        protein_to_cluster,
        args.min_cooccurrence
    )

    # Get ordered cluster list
    cluster_list = sorted(used_clusters)

    # Create adjacency matrix
    adj_matrix = create_adjacency_matrix(edge_df, cluster_list)

    # Save outputs
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    intermediate_dir = Path(args.intermediate_dir)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Save edges (for graph construction)
    edges_file = output_dir / 'gene_gene_edges.csv'
    print(f"\nSaving gene-gene edges to {edges_file}...")
    edge_df.to_csv(edges_file, index=False)

    # Save adjacency matrix
    matrix_file = output_dir / 'gene_gene_matrix.npz'
    print(f"Saving adjacency matrix to {matrix_file}...")
    np.savez_compressed(
        matrix_file,
        data=adj_matrix.data,
        indices=adj_matrix.indices,
        indptr=adj_matrix.indptr,
        shape=adj_matrix.shape,
        cluster_ids=np.array(cluster_list, dtype=object)
    )

    # Save full operon predictions (for analysis)
    operon_file = intermediate_dir / 'operon_predictions.tsv'
    print(f"Saving operon predictions to {operon_file}...")
    operons_combined.to_csv(operon_file, sep='\t', index=False)

    # Print summary stats
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Genomes processed: {len(all_operons):,}")
    print(f"Total genes: {len(operons_combined):,}")
    print(f"Total operons: {operons_combined['operon_id'].nunique():,}")
    print(f"Mean operon size: {operons_combined.groupby('operon_id').size().mean():.1f} genes")
    print(f"\nCluster-cluster edges: {len(edge_df):,}")
    print(f"Mean edge weight: {edge_df['weight'].mean():.1f}")
    print(f"Median edge weight: {edge_df['weight'].median():.0f}")
    print(f"Max edge weight: {edge_df['weight'].max()}")
    print("="*70)

    print("\nDone!")


if __name__ == '__main__':
    main()
