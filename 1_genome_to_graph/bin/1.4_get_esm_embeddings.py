#!/usr/bin/env python
"""
Generate ESM-C embeddings for all genes in RefSeq genomes.

This script processes protein sequences from Prodigal gene annotations
and generates embeddings using the ESM-C 600M model.
"""

import os
import sys
from pathlib import Path
import argparse
from typing import List, Dict, Tuple
import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm
from Bio import SeqIO
import torch

# ESM imports
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate ESM-C embeddings for RefSeq genes"
    )
    parser.add_argument(
        "--gene-dir",
        type=str,
        default="data/refseq_gene_annotations",
        help="Directory containing .faa protein files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="embeddings",
        help="Output directory for embeddings",
    )
    parser.add_argument(
        "--model-name",
        type=str,
        default="esmc_600m",
        help="ESM-C model to use (default: esmc_600m)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for processing sequences",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=1024,
        help="Maximum sequence length to process",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to use (cuda/cpu)",
    )
    parser.add_argument(
        "--num-genomes",
        type=int,
        default=None,
        help="Number of genomes to process (for testing)",
    )
    parser.add_argument(
        "--save-format",
        type=str,
        default="hdf5",
        choices=["hdf5", "npz"],
        help="Format to save embeddings",
    )
    return parser.parse_args()


def get_protein_files(gene_dir: str) -> List[Path]:
    """Get list of all protein .faa files."""
    gene_dir = Path(gene_dir)
    protein_files = sorted(gene_dir.glob("*_prodigal_proteins.faa"))
    print(f"Found {len(protein_files)} protein files")
    return protein_files


def extract_genome_id(filename: str) -> str:
    """Extract genome ID from filename (e.g., GCF_000006985.1)."""
    return filename.split("_prodigal_proteins")[0]


def load_sequences(faa_file: Path, max_length: int = 1024) -> List[Dict]:
    """
    Load protein sequences from a FASTA file.

    Returns:
        List of dicts with 'gene_id', 'sequence', and 'length'
    """
    sequences = []
    for record in SeqIO.parse(faa_file, "fasta"):
        seq_str = str(record.seq)
        # Filter out stop codons and very long sequences
        seq_str = seq_str.replace("*", "")

        if len(seq_str) > 0 and len(seq_str) <= max_length:
            sequences.append({
                "gene_id": record.id,
                "sequence": seq_str,
                "length": len(seq_str),
                "description": record.description,
            })
    return sequences


def get_embeddings_batch(
    model: ESMC,
    sequences: List[str],
    device: str,
) -> np.ndarray:
    """
    Get embeddings for a batch of sequences using ESM-C.

    Args:
        model: ESM-C model
        sequences: List of protein sequences
        device: Device to use

    Returns:
        Array of embeddings (batch_size, embedding_dim)
    """
    # Convert sequences to ESMProtein objects
    proteins = [ESMProtein(sequence=seq) for seq in sequences]

    # Get embeddings - ESM-C returns per-residue embeddings
    # We'll use mean pooling to get sequence-level embeddings
    embeddings = []

    with torch.no_grad():
        for protein in proteins:
            # Get protein embedding
            embedding_result = model.encode(protein)

            # Mean pool over sequence length to get fixed-size embedding
            # embedding_result should be shape (seq_len, embed_dim)
            embedding = embedding_result.mean(axis=0)
            embeddings.append(embedding.cpu().numpy())

    return np.array(embeddings)


def process_genome(
    faa_file: Path,
    model: ESMC,
    args,
) -> Tuple[np.ndarray, List[str], List[int]]:
    """
    Process all genes in a genome and return embeddings.

    Returns:
        embeddings: Array of shape (n_genes, embedding_dim)
        gene_ids: List of gene IDs
        lengths: List of sequence lengths
    """
    # Load sequences
    sequences = load_sequences(faa_file, max_length=args.max_length)

    if len(sequences) == 0:
        return None, None, None

    # Process in batches
    all_embeddings = []
    gene_ids = []
    lengths = []

    for i in range(0, len(sequences), args.batch_size):
        batch = sequences[i:i + args.batch_size]
        batch_seqs = [s["sequence"] for s in batch]

        # Get embeddings
        embeddings = get_embeddings_batch(model, batch_seqs, args.device)
        all_embeddings.append(embeddings)

        gene_ids.extend([s["gene_id"] for s in batch])
        lengths.extend([s["length"] for s in batch])

    # Concatenate all batches
    all_embeddings = np.vstack(all_embeddings)

    return all_embeddings, gene_ids, lengths


def save_embeddings_hdf5(
    output_file: Path,
    genome_id: str,
    embeddings: np.ndarray,
    gene_ids: List[str],
    lengths: List[int],
):
    """Save embeddings to HDF5 file."""
    with h5py.File(output_file, "a") as f:
        # Create group for this genome
        grp = f.create_group(genome_id)

        # Save embeddings
        grp.create_dataset("embeddings", data=embeddings, compression="gzip")

        # Save gene IDs as strings
        dt = h5py.string_dtype(encoding='utf-8')
        grp.create_dataset("gene_ids", data=gene_ids, dtype=dt)

        # Save sequence lengths
        grp.create_dataset("lengths", data=np.array(lengths))


def save_embeddings_npz(
    output_dir: Path,
    genome_id: str,
    embeddings: np.ndarray,
    gene_ids: List[str],
    lengths: List[int],
):
    """Save embeddings to individual NPZ files per genome."""
    output_file = output_dir / f"{genome_id}_embeddings.npz"
    np.savez_compressed(
        output_file,
        embeddings=embeddings,
        gene_ids=gene_ids,
        lengths=lengths,
    )


def main():
    args = parse_args()

    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    print(f"Using device: {args.device}")
    print(f"Loading ESM-C model: {args.model_name}")

    # Load ESM-C model
    model = ESMC.from_pretrained(args.model_name).to(args.device)
    model.eval()

    print(f"Model loaded successfully")

    # Get list of protein files
    protein_files = get_protein_files(args.gene_dir)

    if args.num_genomes is not None:
        protein_files = protein_files[:args.num_genomes]
        print(f"Processing first {args.num_genomes} genomes")

    # Setup output file
    if args.save_format == "hdf5":
        output_file = output_dir / "esmc_embeddings.h5"
        print(f"Saving embeddings to: {output_file}")
    else:
        print(f"Saving embeddings to: {output_dir}")

    # Process each genome
    metadata = []

    for faa_file in tqdm(protein_files, desc="Processing genomes"):
        genome_id = extract_genome_id(faa_file.name)

        try:
            # Get embeddings for this genome
            embeddings, gene_ids, lengths = process_genome(
                faa_file, model, args
            )

            if embeddings is None:
                print(f"  No valid sequences for {genome_id}")
                continue

            # Save embeddings
            if args.save_format == "hdf5":
                save_embeddings_hdf5(
                    output_file, genome_id, embeddings, gene_ids, lengths
                )
            else:
                save_embeddings_npz(
                    output_dir, genome_id, embeddings, gene_ids, lengths
                )

            # Track metadata
            metadata.append({
                "genome_id": genome_id,
                "num_genes": len(gene_ids),
                "embedding_dim": embeddings.shape[1],
                "avg_length": np.mean(lengths),
            })

        except Exception as e:
            print(f"  Error processing {genome_id}: {e}")
            continue

    # Save metadata
    metadata_df = pd.DataFrame(metadata)
    metadata_file = output_dir / "embedding_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)

    print(f"\nProcessing complete!")
    print(f"Processed {len(metadata)} genomes")
    print(f"Total genes: {metadata_df['num_genes'].sum()}")
    print(f"Metadata saved to: {metadata_file}")


if __name__ == "__main__":
    main()
