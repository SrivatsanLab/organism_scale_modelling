#!/usr/bin/env python
"""
Extract sequences for proteins that need ESM embeddings.

Reads the list of protein IDs and extracts their sequences from the
concatenated FASTA file.
"""

import argparse
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO


def load_protein_list(protein_list_file):
    """Load list of proteins needing embeddings."""
    print(f"Loading protein list from {protein_list_file}...")
    proteins = set()
    with open(protein_list_file) as f:
        for line in f:
            proteins.add(line.strip())
    print(f"  Loaded {len(proteins):,} protein IDs")
    return proteins


def extract_sequences(fasta_file, protein_ids, output_file):
    """Extract sequences for specified proteins."""
    print(f"\nExtracting sequences from {fasta_file}...")
    print(f"  This may take a while...")

    found_count = 0
    missing_count = 0

    with open(output_file, 'w') as out_f:
        for record in tqdm(SeqIO.parse(fasta_file, 'fasta'), desc="Scanning FASTA"):
            if record.id in protein_ids:
                SeqIO.write(record, out_f, 'fasta')
                found_count += 1

                if found_count % 100000 == 0:
                    print(f"    Found {found_count:,} sequences so far...")

    missing_count = len(protein_ids) - found_count

    print(f"\n  Found: {found_count:,} sequences")
    print(f"  Missing: {missing_count:,} sequences")
    print(f"  Coverage: {found_count / len(protein_ids) * 100:.2f}%")

    return found_count, missing_count


def main():
    parser = argparse.ArgumentParser(
        description='Extract sequences for proteins needing embeddings'
    )
    parser.add_argument('--protein-list', type=str,
                       default='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/filtered_0p7/proteins_needing_embeddings.txt',
                       help='File with list of protein IDs')
    parser.add_argument('--fasta', type=str,
                       default='data/all_proteins.faa',
                       help='Input FASTA file with all proteins')
    parser.add_argument('--output', type=str,
                       default='data/proteins_for_embedding.faa',
                       help='Output FASTA file')

    args = parser.parse_args()

    print("=" * 80)
    print("EXTRACT SEQUENCES FOR EMBEDDING")
    print("=" * 80)
    print()

    # Load protein list
    protein_ids = load_protein_list(args.protein_list)

    # Extract sequences
    found, missing = extract_sequences(args.fasta, protein_ids, args.output)

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print(f"Output: {args.output}")
    print(f"  Sequences extracted: {found:,}")
    print(f"  File size: {Path(args.output).stat().st_size / 1e9:.2f} GB")
    print()


if __name__ == '__main__':
    main()
