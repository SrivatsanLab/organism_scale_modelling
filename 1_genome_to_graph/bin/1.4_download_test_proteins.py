#!/usr/bin/env python
"""
Download protein sequences for a test subset of genomes from NCBI.

This script downloads protein FASTA files for a small set of genomes
to enable MMseqs2 clustering tests.

Usage:
    python download_test_proteins.py --n-genomes 10 --output test_proteins.faa
"""

import argparse
import numpy as np
import time
from pathlib import Path
from Bio import Entrez
from io import StringIO


# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = "your_email@example.com"


def download_genome_proteins(genome_accession, retries=3, delay=2):
    """
    Download protein sequences for a genome from NCBI using Entrez.

    Args:
        genome_accession: GCF accession (e.g., 'GCF_000006985.1')
        retries: Number of retry attempts
        delay: Delay between retries in seconds

    Returns:
        List of (header, sequence) tuples, or None if failed
    """
    for attempt in range(retries):
        try:
            print(f"  Searching NCBI for {genome_accession}... (attempt {attempt + 1}/{retries})")

            # Search for the assembly
            search_handle = Entrez.esearch(db="assembly", term=genome_accession, retmax=1)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            if not search_results['IdList']:
                print(f"    Assembly not found: {genome_accession}")
                return None

            assembly_id = search_results['IdList'][0]

            # Get assembly summary
            summary_handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(summary_handle, validate=False)
            summary_handle.close()

            # Get FTP path from assembly summary
            try:
                ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
                if not ftp_path:
                    ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']

                # Convert FTP to HTTPS
                ftp_path = ftp_path.replace('ftp://', 'https://')

                # Construct protein FASTA URL
                assembly_name = ftp_path.split('/')[-1]
                protein_url = f"{ftp_path}/{assembly_name}_protein.faa.gz"

                print(f"    Downloading from: {protein_url}")

                # Download protein FASTA
                import urllib.request
                import gzip
                import io

                response = urllib.request.urlopen(protein_url, timeout=60)

                # Decompress gzip data
                with gzip.GzipFile(fileobj=io.BytesIO(response.read())) as f:
                    content = f.read().decode('utf-8')

                # Parse FASTA
                proteins = []
                current_header = None
                current_seq = []

                for line in content.split('\n'):
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        # Save previous protein
                        if current_header is not None:
                            proteins.append((current_header, ''.join(current_seq)))

                        # Start new protein
                        current_header = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line)

                # Save last protein
                if current_header is not None:
                    proteins.append((current_header, ''.join(current_seq)))

                print(f"    Downloaded {len(proteins)} proteins")
                return proteins

            except (KeyError, IndexError) as e:
                print(f"    Error extracting FTP path: {e}")
                return None

        except Exception as e:
            print(f"    Error: {e}")
            if attempt < retries - 1:
                print(f"    Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print(f"    Failed to download {genome_accession}")
                return None

    return None


def main():
    parser = argparse.ArgumentParser(description='Download test protein sequences from NCBI')
    parser.add_argument('--n-genomes', type=int, default=10,
                        help='Number of genomes to download (default: 10)')
    parser.add_argument('--output', type=str, default='data/test_proteins.faa',
                        help='Output FASTA file')
    parser.add_argument('--pca-cache', type=str, default='results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz',
                        help='PCA cache with gene IDs (default: results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz)')

    args = parser.parse_args()

    # Load PCA cache to get genome accessions
    print(f"Loading PCA cache from {args.pca_cache}")
    data = np.load(args.pca_cache, allow_pickle=True)

    genome_ids = data['genome_ids']
    gene_ids = data['gene_ids']

    # Get unique genomes
    unique_genomes = np.unique(genome_ids)
    print(f"Found {len(unique_genomes)} unique genomes in cache")

    # Select subset of genomes
    selected_genomes = unique_genomes[:args.n_genomes]
    print(f"\nDownloading proteins for {len(selected_genomes)} genomes:")
    print(selected_genomes)
    print()

    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Download proteins for each genome
    all_proteins = []
    downloaded_genomes = []

    for i, genome_acc in enumerate(selected_genomes):
        print(f"[{i+1}/{len(selected_genomes)}] Processing {genome_acc}")

        proteins = download_genome_proteins(genome_acc)

        if proteins is not None:
            all_proteins.extend(proteins)
            downloaded_genomes.append(genome_acc)

        # Be nice to NCBI servers
        time.sleep(1)

    print(f"\nSuccessfully downloaded {len(downloaded_genomes)}/{len(selected_genomes)} genomes")
    print(f"Total proteins: {len(all_proteins)}")

    # Write to FASTA
    print(f"\nWriting proteins to {output_path}")
    with open(output_path, 'w') as f:
        for header, seq in all_proteins:
            f.write(f">{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")

    print(f"\nDone! Wrote {len(all_proteins)} proteins to {output_path}")

    # Print statistics
    print("\nStatistics:")
    print(f"  Total proteins: {len(all_proteins)}")
    if len(downloaded_genomes) > 0:
        print(f"  Average per genome: {len(all_proteins) / len(downloaded_genomes):.1f}")

    # Count how many genes from PCA cache are in downloaded genomes
    genes_in_downloaded = 0
    for genome_acc in downloaded_genomes:
        genes_in_downloaded += np.sum(genome_ids == genome_acc)

    print(f"  Genes in PCA cache from downloaded genomes: {genes_in_downloaded:,}")
    print(f"  Coverage: {genes_in_downloaded / len(gene_ids) * 100:.2f}%")


if __name__ == '__main__':
    main()
