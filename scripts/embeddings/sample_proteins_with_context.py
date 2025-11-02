#!/usr/bin/env python
"""
Sample proteins with genomic context for embedding comparison.

This script:
1. Randomly selects 2 genomes from the GTDB database
2. Samples 1000 proteins from each genome
3. Extracts 5' and 3' flanking DNA sequences for each protein
4. Creates datasets for comparing ESM-C vs gLM2 embeddings

For Issue #2: Compare protein embeddings with/without genomic context
For Issue #3: Test generalization to molecular interactions
"""

import argparse
import random
from pathlib import Path
from typing import List, Dict, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Sample proteins with genomic context"
    )
    parser.add_argument(
        "--genome-dir",
        type=str,
        default="data/refseq_genomes",
        help="Directory containing genome .fna files",
    )
    parser.add_argument(
        "--annotation-dir",
        type=str,
        default="data/refseq_gene_annotations",
        help="Directory containing gene annotation files (.gff)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/protein_samples",
        help="Output directory for sampled proteins",
    )
    parser.add_argument(
        "--n-genomes",
        type=int,
        default=2,
        help="Number of genomes to sample",
    )
    parser.add_argument(
        "--n-proteins-per-genome",
        type=int,
        default=1000,
        help="Number of proteins to sample per genome",
    )
    parser.add_argument(
        "--flanking-length",
        type=int,
        default=500,
        help="Length of flanking DNA sequence (bp) on each side",
    )
    parser.add_argument(
        "--min-protein-length",
        type=int,
        default=50,
        help="Minimum protein length (aa)",
    )
    parser.add_argument(
        "--max-protein-length",
        type=int,
        default=1000,
        help="Maximum protein length (aa)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )
    return parser.parse_args()


def get_genome_files(genome_dir: Path) -> List[Path]:
    """Get list of all genome files."""
    genome_dir = Path(genome_dir)

    # Try .fna first
    genome_files = sorted(genome_dir.glob("*.fna"))
    if len(genome_files) == 0:
        # Try .fasta
        genome_files = sorted(genome_dir.glob("*.fasta"))

    return genome_files


def extract_genome_id(filename: str) -> str:
    """Extract genome ID from filename (GCF_XXXXXXXXX.X format)."""
    name = Path(filename).stem
    # Extract just the GCF ID (e.g., GCF_024496245.1)
    if name.startswith("GCF_") or name.startswith("GCA_"):
        # Split on underscore and take first two parts (GCF_XXXXXXXXX.X)
        parts = name.split("_")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
    # Fallback: remove common suffixes
    name = name.replace("_genomic", "").replace("_complete", "")
    return name


def load_gff_annotations(gff_file: Path) -> List[Dict]:
    """
    Load gene annotations from GFF file.

    Returns:
        List of dicts with gene information
    """
    annotations = []

    if not gff_file.exists():
        print(f"  Warning: GFF file not found: {gff_file}")
        return annotations

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'CDS':
                continue

            seqid = fields[0]
            start = int(fields[3]) - 1  # GFF is 1-indexed, convert to 0-indexed
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key] = value

            gene_id = attributes.get('ID', f"{seqid}_{start}_{end}")

            annotations.append({
                'gene_id': gene_id,
                'seqid': seqid,
                'start': start,
                'end': end,
                'strand': strand,
                'length': end - start,
            })

    return annotations


def extract_protein_with_context(
    genome_seq: str,
    annotation: Dict,
    flanking_length: int,
) -> Dict:
    """
    Extract protein sequence with flanking genomic context.

    Returns:
        Dict with protein sequence, 5' and 3' flanking DNA
    """
    start = annotation['start']
    end = annotation['end']
    strand = annotation['strand']

    # Get CDS sequence
    cds_seq = genome_seq[start:end]

    # Get 5' and 3' flanking regions
    upstream_start = max(0, start - flanking_length)
    downstream_end = min(len(genome_seq), end + flanking_length)

    flanking_5prime = genome_seq[upstream_start:start]
    flanking_3prime = genome_seq[end:downstream_end]

    # Handle reverse strand
    if strand == '-':
        cds_seq = str(Seq(cds_seq).reverse_complement())
        # For reverse strand, flip 5' and 3'
        flanking_5prime, flanking_3prime = (
            str(Seq(flanking_3prime).reverse_complement()),
            str(Seq(flanking_5prime).reverse_complement())
        )

    # Translate CDS to protein
    try:
        protein_seq = str(Seq(cds_seq).translate(to_stop=True))
        protein_seq = protein_seq.replace('*', '')  # Remove stop codons
    except Exception as e:
        protein_seq = None

    return {
        'gene_id': annotation['gene_id'],
        'protein_seq': protein_seq,
        'cds_seq': cds_seq,
        'flanking_5prime': flanking_5prime,
        'flanking_3prime': flanking_3prime,
        'strand': strand,
        'start': start,
        'end': end,
        'protein_length': len(protein_seq) if protein_seq else 0,
    }


def sample_proteins_from_genome(
    genome_file: Path,
    annotation_file: Path,
    n_proteins: int,
    flanking_length: int,
    min_length: int,
    max_length: int,
    seed: int,
) -> List[Dict]:
    """
    Sample proteins from a genome with genomic context.

    Returns:
        List of protein data dicts
    """
    genome_id = extract_genome_id(genome_file.name)
    print(f"\nProcessing genome: {genome_id}")

    # Load genome sequence
    print(f"  Loading genome sequence...")
    genome_record = next(SeqIO.parse(genome_file, "fasta"))
    genome_seq = str(genome_record.seq)
    print(f"  Genome length: {len(genome_seq):,} bp")

    # Load annotations
    print(f"  Loading annotations...")
    annotations = load_gff_annotations(annotation_file)
    print(f"  Total genes: {len(annotations)}")

    if len(annotations) == 0:
        print(f"  ERROR: No annotations found!")
        return []

    # Filter by protein length
    print(f"  Extracting proteins...")
    all_proteins = []
    for annot in annotations:
        protein_data = extract_protein_with_context(
            genome_seq, annot, flanking_length
        )

        if protein_data['protein_seq'] is None:
            continue

        prot_len = protein_data['protein_length']
        if min_length <= prot_len <= max_length:
            protein_data['genome_id'] = genome_id
            all_proteins.append(protein_data)

    print(f"  Valid proteins ({min_length}-{max_length} aa): {len(all_proteins)}")

    # Sample proteins
    if len(all_proteins) > n_proteins:
        random.seed(seed)
        sampled = random.sample(all_proteins, n_proteins)
        print(f"  Sampled: {n_proteins} proteins")
    else:
        sampled = all_proteins
        print(f"  Using all {len(sampled)} proteins (fewer than requested)")

    return sampled


def save_protein_data(proteins: List[Dict], output_dir: Path):
    """Save protein data to files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    print(f"\nSaving data to {output_dir}...")

    # Create DataFrame
    df = pd.DataFrame(proteins)

    # Save metadata CSV
    metadata_cols = ['genome_id', 'gene_id', 'strand', 'start', 'end', 'protein_length']
    metadata_df = df[metadata_cols]
    metadata_file = output_dir / "protein_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    print(f"  Saved metadata: {metadata_file}")

    # Save protein sequences (FASTA)
    protein_fasta = output_dir / "proteins.faa"
    with open(protein_fasta, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row['genome_id']}|{row['gene_id']}\n")
            f.write(f"{row['protein_seq']}\n")
    print(f"  Saved proteins: {protein_fasta}")

    # Save CDS sequences (FASTA)
    cds_fasta = output_dir / "cds_sequences.fna"
    with open(cds_fasta, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row['genome_id']}|{row['gene_id']}\n")
            f.write(f"{row['cds_seq']}\n")
    print(f"  Saved CDS: {cds_fasta}")

    # Save 5' flanking sequences
    flanking_5_fasta = output_dir / "flanking_5prime.fna"
    with open(flanking_5_fasta, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row['genome_id']}|{row['gene_id']}\n")
            f.write(f"{row['flanking_5prime']}\n")
    print(f"  Saved 5' flanking: {flanking_5_fasta}")

    # Save 3' flanking sequences
    flanking_3_fasta = output_dir / "flanking_3prime.fna"
    with open(flanking_3_fasta, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row['genome_id']}|{row['gene_id']}\n")
            f.write(f"{row['flanking_3prime']}\n")
    print(f"  Saved 3' flanking: {flanking_3_fasta}")

    # Save combined genomic context (5' + CDS + 3')
    genomic_fasta = output_dir / "genomic_context.fna"
    with open(genomic_fasta, 'w') as f:
        for _, row in df.iterrows():
            combined = row['flanking_5prime'] + row['cds_seq'] + row['flanking_3prime']
            f.write(f">{row['genome_id']}|{row['gene_id']}\n")
            f.write(f"{combined}\n")
    print(f"  Saved genomic context: {genomic_fasta}")

    # Summary statistics
    print(f"\n  Summary:")
    print(f"    Total proteins: {len(df)}")
    print(f"    Genomes: {df['genome_id'].nunique()}")
    print(f"    Avg protein length: {df['protein_length'].mean():.1f} aa")
    print(f"    Avg 5' flank length: {df['flanking_5prime'].str.len().mean():.1f} bp")
    print(f"    Avg 3' flank length: {df['flanking_3prime'].str.len().mean():.1f} bp")


def main():
    args = parse_args()

    print("="*80)
    print("Sample Proteins with Genomic Context")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Genomes to sample: {args.n_genomes}")
    print(f"  Proteins per genome: {args.n_proteins_per_genome}")
    print(f"  Flanking DNA length: {args.flanking_length} bp")
    print(f"  Protein length range: {args.min_protein_length}-{args.max_protein_length} aa")
    print(f"  Random seed: {args.seed}")

    # Get list of genome files
    genome_dir = Path(args.genome_dir)
    annotation_dir = Path(args.annotation_dir)

    genome_files = get_genome_files(genome_dir)
    print(f"\nFound {len(genome_files)} genome files in {genome_dir}")

    if len(genome_files) == 0:
        print("ERROR: No genome files found!")
        return

    # Randomly select genomes
    random.seed(args.seed)
    selected_genomes = random.sample(genome_files, min(args.n_genomes, len(genome_files)))

    print(f"\nSelected {len(selected_genomes)} genomes:")
    for gf in selected_genomes:
        print(f"  - {extract_genome_id(gf.name)}")

    # Sample proteins from each genome
    all_proteins = []

    for genome_file in selected_genomes:
        genome_id = extract_genome_id(genome_file.name)

        # Find annotation file
        annotation_patterns = [
            annotation_dir / f"{genome_id}_prodigal_genes.gff",
            annotation_dir / f"{genome_id}_prodigal.gff",
            annotation_dir / f"{genome_id}.gff",
            annotation_dir / f"{genome_id}_genes.gff",
        ]

        annotation_file = None
        for pattern in annotation_patterns:
            if pattern.exists():
                annotation_file = pattern
                break

        if annotation_file is None:
            print(f"  Warning: No annotation file found for {genome_id}")
            continue

        # Sample proteins
        proteins = sample_proteins_from_genome(
            genome_file,
            annotation_file,
            args.n_proteins_per_genome,
            args.flanking_length,
            args.min_protein_length,
            args.max_protein_length,
            args.seed,
        )

        all_proteins.extend(proteins)

    # Save all data
    if len(all_proteins) > 0:
        save_protein_data(all_proteins, Path(args.output_dir))

        print("\n" + "="*80)
        print("Complete!")
        print("="*80)
        print(f"\nTotal proteins sampled: {len(all_proteins)}")
        print(f"Output directory: {args.output_dir}")

        print("\nNext steps:")
        print("1. Generate ESM-C embeddings (protein-only):")
        print(f"   python scripts/embeddings/compare_glm2_esmc.py --protein-dir {args.output_dir}")
        print("\n2. Generate gLM2 embeddings (with genomic context):")
        print(f"   python scripts/embeddings/get_glm2_embeddings.py --mode genes")
        print("\n3. Compare embeddings with/without context (Issue #2)")
        print("4. Test interaction predictions (Issue #3)")
    else:
        print("\nERROR: No proteins sampled!")


if __name__ == "__main__":
    main()
