#!/usr/bin/env python
"""
Generate ESM embeddings for a batch of proteins.

This script is designed to be run as part of a SLURM array job,
where each array task processes a subset of proteins.
"""

import argparse
import torch
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from esm.models.esmc import ESMC
from esm.sdk.api import ESM3InferenceClient, ESMProtein


def load_sequences_for_batch(fasta_file, batch_idx, batch_size):
    """Load sequences for a specific batch."""
    print(f"Loading sequences for batch {batch_idx}...")

    start_idx = batch_idx * batch_size
    end_idx = start_idx + batch_size

    sequences = []
    current_idx = 0

    for record in SeqIO.parse(fasta_file, 'fasta'):
        if current_idx >= start_idx and current_idx < end_idx:
            sequences.append((record.id, str(record.seq)))
        current_idx += 1

        if current_idx >= end_idx:
            break

    print(f"  Loaded {len(sequences):,} sequences (indices {start_idx:,}-{end_idx-1:,})")
    return sequences


def generate_embeddings_esmc(sequences, device='cuda'):
    """Generate embeddings using ESM-C model."""
    print(f"Generating embeddings using ESM-C...")
    print(f"  Device: {device}")

    # Convert device string to torch.device object
    device = torch.device(device)

    # Load model
    print("  Loading model...")
    client = ESMC.from_pretrained("esmc_600m", device=device)
    # Convert model to float32 for GPU compatibility (BFloat16 not supported on older GPUs)
    client = client.float()

    embeddings = []
    gene_ids = []

    for gene_id, sequence in tqdm(sequences, desc="Generating embeddings"):
        try:
            # Create protein and encode
            protein = ESMProtein(sequence=sequence)
            protein_tensor = client.encode(protein)

            # Generate embedding (mean pool across sequence length)
            with torch.no_grad():
                output = client.forward(protein_tensor.sequence.unsqueeze(0))
                embedding = output.embeddings.mean(dim=1).squeeze(0).cpu().numpy()

            embeddings.append(embedding)
            gene_ids.append(gene_id)

        except Exception as e:
            print(f"  WARNING: Failed to process {gene_id}: {e}")
            continue

    print(f"  Generated {len(embeddings):,} embeddings")

    return gene_ids, np.array(embeddings)


def main():
    parser = argparse.ArgumentParser(
        description='Generate ESM embeddings for a batch of proteins'
    )
    parser.add_argument('--fasta', type=str, required=True,
                       help='Input FASTA file')
    parser.add_argument('--batch-idx', type=int, required=True,
                       help='Batch index (0-based)')
    parser.add_argument('--batch-size', type=int, default=10000,
                       help='Number of proteins per batch')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    parser.add_argument('--device', type=str, default='cuda',
                       help='Device (cuda or cpu)')

    args = parser.parse_args()

    print("=" * 80)
    print(f"BATCH EMBEDDING GENERATION - Batch {args.batch_idx}")
    print("=" * 80)
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load sequences for this batch
    sequences = load_sequences_for_batch(args.fasta, args.batch_idx, args.batch_size)

    if len(sequences) == 0:
        print("  No sequences in this batch, exiting...")
        return

    # Generate embeddings
    gene_ids, embeddings = generate_embeddings_esmc(sequences, device=args.device)

    # Save embeddings
    output_file = output_dir / f'embeddings_batch_{args.batch_idx:05d}.npz'
    print(f"\nSaving embeddings to {output_file}...")

    np.savez_compressed(
        output_file,
        gene_ids=gene_ids,
        embeddings=embeddings,
        batch_idx=args.batch_idx
    )

    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print(f"Batch {args.batch_idx}:")
    print(f"  Sequences processed: {len(gene_ids):,}")
    if len(embeddings.shape) > 1:
        print(f"  Embedding dimensions: {embeddings.shape[1]}")
    else:
        print(f"  Embedding shape: {embeddings.shape}")
    print(f"  Output file: {output_file}")
    print()


if __name__ == '__main__':
    main()
