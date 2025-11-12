#!/usr/bin/env python
"""
Predict genes from bacterial genome sequences and extract protein sequences.

This script takes bacterial genome FASTA files and:
1. Predicts coding sequences (CDS) using Prodigal
2. Extracts the DNA sequences of predicted genes
3. Translates them to amino acid sequences
4. Saves outputs in FASTA and GFF format

Usage:
    python predict_genes.py genome.fasta --output-dir genes/
    python predict_genes.py --genome-dir genomes/ --output-dir genes/
"""

import argparse
import subprocess
from pathlib import Path
from typing import Optional, List, Dict
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Predict genes and extract proteins from bacterial genomes"
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--genome",
        type=str,
        help="Single genome FASTA file to process",
    )
    input_group.add_argument(
        "--genome-dir",
        type=str,
        help="Directory containing genome FASTA files",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="predicted_genes",
        help="Output directory for gene predictions",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="pyrodigal",
        choices=["pyrodigal", "prodigal"],
        help="Gene prediction method (default: pyrodigal)",
    )
    parser.add_argument(
        "--closed",
        action="store_true",
        help="Treat sequences as closed ends (complete genomes)",
    )
    parser.add_argument(
        "--min-gene-length",
        type=int,
        default=90,
        help="Minimum gene length in nucleotides (default: 90)",
    )
    parser.add_argument(
        "--num-genomes",
        type=int,
        default=None,
        help="Limit number of genomes to process (for testing)",
    )

    return parser.parse_args()


def check_dependencies(method: str) -> bool:
    """Check if required dependencies are available."""
    if method == "prodigal":
        try:
            result = subprocess.run(
                ["prodigal", "-v"],
                capture_output=True,
                text=True
            )
            print(f"Using Prodigal: {result.stderr.split()[1]}")
            return True
        except FileNotFoundError:
            print("Error: Prodigal not found in PATH")
            print("Install with: conda install -c bioconda prodigal")
            return False

    elif method == "pyrodigal":
        try:
            import pyrodigal
            print(f"Using pyrodigal version: {pyrodigal.__version__}")
            return True
        except ImportError:
            print("Error: pyrodigal not installed")
            print("Install with: pip install pyrodigal")
            return False

    return False


def run_prodigal_external(
    genome_file: Path,
    output_prefix: Path,
    closed: bool = False,
) -> Dict[str, Path]:
    """
    Run Prodigal using external command line tool.

    Returns:
        Dictionary with paths to output files
    """
    # Prodigal command
    cmd = [
        "prodigal",
        "-i", str(genome_file),
        "-a", str(output_prefix) + "_proteins.faa",
        "-d", str(output_prefix) + "_genes.fna",
        "-f", "gff",
        "-o", str(output_prefix) + "_genes.gff",
        "-s", str(output_prefix) + "_stats.txt",
    ]

    if closed:
        cmd.append("-c")

    # Run Prodigal
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Prodigal: {e}")
        print(f"stderr: {e.stderr.decode()}")
        raise

    return {
        "proteins": output_prefix.with_name(output_prefix.name + "_proteins.faa"),
        "genes": output_prefix.with_name(output_prefix.name + "_genes.fna"),
        "gff": output_prefix.with_name(output_prefix.name + "_genes.gff"),
        "stats": output_prefix.with_name(output_prefix.name + "_stats.txt"),
    }


def run_pyrodigal(
    genome_file: Path,
    output_prefix: Path,
    closed: bool = False,
    min_gene_length: int = 90,
) -> Dict[str, Path]:
    """
    Run Prodigal using pyrodigal Python library.

    Returns:
        Dictionary with paths to output files
    """
    import pyrodigal

    # Initialize gene finder
    orf_finder = pyrodigal.GeneFinder(
        closed=closed,
        min_gene=min_gene_length,
    )

    # Read genome
    genome_records = list(SeqIO.parse(genome_file, "fasta"))

    # Output files
    proteins_file = output_prefix.with_name(output_prefix.name + "_proteins.faa")
    genes_file = output_prefix.with_name(output_prefix.name + "_genes.fna")
    gff_file = output_prefix.with_name(output_prefix.name + "_genes.gff")
    stats_file = output_prefix.with_name(output_prefix.name + "_stats.txt")

    all_proteins = []
    all_genes = []
    all_gff_lines = ["##gff-version 3"]
    stats = []

    gene_counter = 0

    for record in genome_records:
        sequence = str(record.seq)

        # Train or use pre-trained model for the sequence
        genes = orf_finder.find_genes(sequence)

        for i, gene in enumerate(genes, start=1):
            gene_counter += 1

            # Create gene ID
            gene_id = f"{record.id}_{gene_counter}"

            # Get gene sequence
            gene_seq = gene.sequence()

            # Get protein sequence
            protein_seq = gene.translate()

            # Add to collections
            all_genes.append(f">{gene_id} # {gene.begin} # {gene.end} # {gene.strand} # partial={'1' if gene.partial_begin else '0'}{'1' if gene.partial_end else '0'}\n{gene_seq}\n")
            all_proteins.append(f">{gene_id} # {gene.begin} # {gene.end} # {gene.strand} # partial={'1' if gene.partial_begin else '0'}{'1' if gene.partial_end else '0'}\n{protein_seq}\n")

            # GFF line
            gff_line = f"{record.id}\tProdigal\tCDS\t{gene.begin}\t{gene.end}\t{gene.score:.2f}\t{'+' if gene.strand == 1 else '-'}\t0\tID={gene_id}"
            all_gff_lines.append(gff_line)

            # Stats
            stats.append({
                "gene_id": gene_id,
                "contig": record.id,
                "start": gene.begin,
                "end": gene.end,
                "strand": gene.strand,
                "length": len(gene_seq),
                "score": gene.score,
            })

    # Write outputs
    with open(proteins_file, "w") as f:
        f.writelines(all_proteins)

    with open(genes_file, "w") as f:
        f.writelines(all_genes)

    with open(gff_file, "w") as f:
        f.write("\n".join(all_gff_lines) + "\n")

    # Write stats
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(stats_file, sep="\t", index=False)

    print(f"  Found {gene_counter} genes")
    print(f"  Average gene length: {stats_df['length'].mean():.1f} bp")

    return {
        "proteins": proteins_file,
        "genes": genes_file,
        "gff": gff_file,
        "stats": stats_file,
    }


def predict_genes_from_genome(
    genome_file: Path,
    output_dir: Path,
    method: str = "pyrodigal",
    closed: bool = False,
    min_gene_length: int = 90,
) -> Optional[Dict[str, Path]]:
    """
    Predict genes from a genome file.

    Args:
        genome_file: Path to genome FASTA file
        output_dir: Directory to save outputs
        method: Gene prediction method ('pyrodigal' or 'prodigal')
        closed: Treat sequences as closed ends
        min_gene_length: Minimum gene length

    Returns:
        Dictionary with paths to output files
    """
    # Create genome ID from filename
    genome_id = genome_file.stem

    # Output prefix
    output_prefix = output_dir / genome_id

    print(f"Processing: {genome_id}")

    try:
        if method == "prodigal":
            results = run_prodigal_external(
                genome_file, output_prefix, closed
            )
        else:  # pyrodigal
            results = run_pyrodigal(
                genome_file, output_prefix, closed, min_gene_length
            )

        return results

    except Exception as e:
        print(f"  Error processing {genome_id}: {e}")
        return None


def get_genome_files(path: Path) -> List[Path]:
    """Get list of genome FASTA files."""
    if path.is_file():
        return [path]
    elif path.is_dir():
        # Look for FASTA files
        fasta_files = []
        for ext in ["*.fasta", "*.fa", "*.fna"]:
            fasta_files.extend(path.glob(ext))
        return sorted(fasta_files)
    else:
        return []


def main():
    args = parse_args()

    print("=" * 70)
    print("Bacterial Gene Prediction Pipeline")
    print("=" * 70)

    # Check dependencies
    if not check_dependencies(args.method):
        sys.exit(1)

    # Get genome files
    if args.genome:
        genome_files = [Path(args.genome)]
    else:
        genome_files = get_genome_files(Path(args.genome_dir))

    if not genome_files:
        print("Error: No genome files found")
        sys.exit(1)

    print(f"\nFound {len(genome_files)} genome file(s)")

    if args.num_genomes:
        genome_files = genome_files[:args.num_genomes]
        print(f"Processing first {args.num_genomes} genomes")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Process each genome
    results = []
    for genome_file in genome_files:
        result = predict_genes_from_genome(
            genome_file,
            output_dir,
            method=args.method,
            closed=args.closed,
            min_gene_length=args.min_gene_length,
        )
        if result:
            results.append({
                "genome": genome_file.stem,
                **{k: str(v) for k, v in result.items()}
            })

    # Save summary
    if results:
        summary_df = pd.DataFrame(results)
        summary_file = output_dir / "prediction_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"\n{'=' * 70}")
        print(f"Successfully processed {len(results)} genome(s)")
        print(f"Output directory: {output_dir}")
        print(f"Summary saved to: {summary_file}")
        print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
