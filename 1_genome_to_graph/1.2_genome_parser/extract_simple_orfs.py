#!/usr/bin/env python
"""
Simple ORF (Open Reading Frame) finder and translator.

This is a lightweight alternative to Prodigal for basic ORF extraction.
Use this if you just need simple translation of all possible ORFs.

For better gene prediction, use predict_genes.py with Prodigal/pyrodigal.

Usage:
    python extract_simple_orfs.py genome.fasta --output proteins.faa
"""

import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract and translate ORFs from genome sequences"
    )
    parser.add_argument(
        "genome",
        type=str,
        help="Input genome FASTA file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="orfs.faa",
        help="Output protein FASTA file",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=30,
        help="Minimum ORF length in amino acids (default: 30)",
    )
    parser.add_argument(
        "--start-codons",
        type=str,
        default="ATG,GTG,TTG",
        help="Comma-separated start codons (default: ATG,GTG,TTG)",
    )
    parser.add_argument(
        "--table",
        type=int,
        default=11,
        help="Genetic code table (default: 11 for bacteria)",
    )
    return parser.parse_args()


def find_orfs(
    sequence: str,
    min_length: int = 30,
    start_codons: list = ["ATG", "GTG", "TTG"],
    table: int = 11,
):
    """
    Find all ORFs in a DNA sequence (both strands).

    Args:
        sequence: DNA sequence string
        min_length: Minimum ORF length in amino acids
        start_codons: List of valid start codons
        table: Genetic code table number

    Yields:
        Tuples of (start, end, strand, protein_sequence)
    """
    sequence = sequence.upper()

    # Forward strand
    for frame in range(3):
        seq_obj = Seq(sequence[frame:])
        for start_codon in start_codons:
            # Find all occurrences of start codon
            for match in re.finditer(start_codon, str(seq_obj)):
                start_pos = match.start()
                subseq = seq_obj[start_pos:]

                # Translate until stop codon
                try:
                    protein = subseq.translate(table=table, to_stop=True)
                    if len(protein) >= min_length:
                        actual_start = frame + start_pos
                        actual_end = actual_start + len(protein) * 3
                        yield (actual_start, actual_end, 1, str(protein))
                except:
                    continue

    # Reverse strand
    rev_seq_obj = Seq(sequence).reverse_complement()
    for frame in range(3):
        seq_obj = Seq(str(rev_seq_obj)[frame:])
        for start_codon in start_codons:
            for match in re.finditer(start_codon, str(seq_obj)):
                start_pos = match.start()
                subseq = seq_obj[start_pos:]

                try:
                    protein = subseq.translate(table=table, to_stop=True)
                    if len(protein) >= min_length:
                        # Convert to original sequence coordinates
                        actual_end = len(sequence) - (frame + start_pos)
                        actual_start = actual_end - len(protein) * 3
                        yield (actual_start, actual_end, -1, str(protein))
                except:
                    continue


def main():
    args = parse_args()

    # Parse start codons
    start_codons = [s.strip().upper() for s in args.start_codons.split(",")]

    print(f"Reading genome: {args.genome}")
    print(f"Minimum ORF length: {args.min_length} aa")
    print(f"Start codons: {', '.join(start_codons)}")
    print(f"Genetic code table: {args.table}")

    # Process genome
    output_file = Path(args.output)
    orf_count = 0

    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(args.genome, "fasta"):
            print(f"\nProcessing: {record.id} ({len(record.seq)} bp)")

            orfs = list(find_orfs(
                str(record.seq),
                min_length=args.min_length,
                start_codons=start_codons,
                table=args.table,
            ))

            print(f"  Found {len(orfs)} ORFs")

            for i, (start, end, strand, protein) in enumerate(orfs, start=1):
                orf_id = f"{record.id}_ORF_{orf_count + 1}"
                strand_str = "+" if strand == 1 else "-"

                out_f.write(
                    f">{orf_id} {record.description} | "
                    f"{start}..{end} ({strand_str}) | {len(protein)}aa\n"
                )
                out_f.write(f"{protein}\n")

                orf_count += 1

    print(f"\n{'=' * 60}")
    print(f"Total ORFs found: {orf_count}")
    print(f"Output saved to: {output_file}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
