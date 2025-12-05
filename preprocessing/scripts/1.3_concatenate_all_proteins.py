#!/usr/bin/env python
"""
Concatenate all protein FASTA files from 7,664 genomes into a single file.

This creates a master protein FASTA file for MMseqs2 clustering.
The gene IDs are modified to include genome information for traceability.
"""

import argparse
from pathlib import Path
from tqdm import tqdm
import gzip


def process_fasta_file(fasta_file, output_handle, genome_id):
    """
    Process a single FASTA file and write to output with modified headers.

    Args:
        fasta_file: Input FASTA file path
        output_handle: Output file handle
        genome_id: Genome ID (e.g., GCF_000006985.1)
    """
    sequence_count = 0

    with open(fasta_file, 'r') as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # Write previous sequence if exists
                if current_header is not None:
                    seq_str = ''.join(current_seq)
                    output_handle.write(f"{current_header}\n{seq_str}\n")
                    sequence_count += 1

                # Parse new header
                # Original format: >NC_002932.3_1 # 2 # 1126 # -1 # ID=1_1;...
                # Extract contig_gene ID
                parts = line[1:].split()
                contig_gene = parts[0]  # e.g., NC_002932.3_1

                # Create new header with genome_contig_gene format
                # This matches the format in the PCA cache
                new_id = f"{genome_id}_{contig_gene}"
                current_header = f">{new_id}"
                current_seq = []
            else:
                # Accumulate sequence
                current_seq.append(line)

        # Write last sequence
        if current_header is not None:
            seq_str = ''.join(current_seq)
            output_handle.write(f"{current_header}\n{seq_str}\n")
            sequence_count += 1

    return sequence_count


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate all protein FASTA files from RefSeq genomes'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='/fh/working/srivatsan_s/ShotgunDomestication/data/refseq_genomes/gene_annotations',
        help='Directory containing protein FASTA files'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='data/all_proteins.faa',
        help='Output concatenated FASTA file'
    )
    parser.add_argument(
        '--pattern',
        type=str,
        default='*_prodigal_proteins.faa',
        help='File pattern to match (default: *_prodigal_proteins.faa)'
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_file = Path(args.output)

    print("=" * 80)
    print("CONCATENATING PROTEIN SEQUENCES")
    print("=" * 80)
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print()

    # Find all protein FASTA files
    fasta_files = sorted(input_dir.glob(args.pattern))
    print(f"Found {len(fasta_files):,} protein FASTA files")
    print()

    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Process all files
    total_proteins = 0

    with open(output_file, 'w') as out_handle:
        for fasta_file in tqdm(fasta_files, desc="Processing genomes"):
            # Extract genome ID from filename
            # e.g., GCF_000006985.1_prodigal_proteins.faa -> GCF_000006985.1
            genome_id = fasta_file.stem.replace('_prodigal_proteins', '')

            # Process this genome's proteins
            n_proteins = process_fasta_file(fasta_file, out_handle, genome_id)
            total_proteins += n_proteins

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Genomes processed: {len(fasta_files):,}")
    print(f"Total proteins: {total_proteins:,}")
    print(f"Average proteins per genome: {total_proteins / len(fasta_files):.0f}")
    print(f"Output file: {output_file}")
    print(f"File size: {output_file.stat().st_size / 1e9:.2f} GB")
    print()


if __name__ == '__main__':
    main()
