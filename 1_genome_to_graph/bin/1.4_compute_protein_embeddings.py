#!/usr/bin/env python
"""
Compute protein-only embeddings and pairwise distances for ESM-C and gLM2.

This script:
1. Loads protein sequences from sampled data
2. Generates embeddings using ESM-C (protein-only)
3. Generates embeddings using gLM2 (protein-only, no genomic context)
4. Computes pairwise distance matrices for both
5. Saves embeddings and distance matrices

For Issue #2: Compare protein embeddings with/without genomic context
"""

import argparse
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from typing import List, Dict, Tuple
from Bio import SeqIO
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute protein-only embeddings for ESM-C and gLM2"
    )
    parser.add_argument(
        "--protein-dir",
        type=str,
        default="data/protein_samples",
        help="Directory containing protein samples",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/embeddings",
        help="Output directory for embeddings",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for embedding computation",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to use (cuda/cpu)",
    )
    return parser.parse_args()


def load_proteins(protein_file: Path) -> List[Dict]:
    """Load protein sequences from FASTA file."""
    proteins = []
    for record in SeqIO.parse(protein_file, "fasta"):
        genome_id, gene_id = record.id.split("|")
        proteins.append({
            "id": record.id,
            "genome_id": genome_id,
            "gene_id": gene_id,
            "sequence": str(record.seq),
            "length": len(record.seq),
        })
    return proteins


def compute_esmc_embeddings(
    proteins: List[Dict],
    batch_size: int = 32,
    device: str = "cuda",
) -> np.ndarray:
    """
    Compute ESM-C embeddings for proteins.

    Returns:
        np.ndarray of shape (n_proteins, embedding_dim)
    """
    print("\n" + "="*80)
    print("Computing ESM-C Embeddings (Protein-Only)")
    print("="*80)

    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    print("Loading ESM-C model...")
    # Convert device string to torch.device
    device_obj = torch.device(device)
    model = ESMC.from_pretrained("esmc_300m", device=device_obj)

    print(f"Model loaded on {device}")
    print(f"Expected embedding dimension: 960")

    embeddings = []

    print(f"\nProcessing {len(proteins)} proteins in batches of {batch_size}...")

    for i in tqdm(range(0, len(proteins), batch_size)):
        batch = proteins[i:i + batch_size]

        batch_embeddings = []
        for prot in batch:
            # Create ESMProtein object
            protein = ESMProtein(sequence=prot["sequence"])

            # Generate embedding using ESM-C API
            # Step 1: Encode protein to get protein_tensor
            protein_tensor = model.encode(protein)

            # Step 2: Get embeddings via logits with return_embeddings=True
            logits_output = model.logits(
                protein_tensor,
                LogitsConfig(sequence=True, return_embeddings=True)
            )

            # Step 3: Extract embeddings (shape: (1, seq_len, embedding_dim))
            embedding = logits_output.embeddings.detach().cpu().numpy()

            # Mean pool across sequence length
            embedding = embedding.mean(axis=1)[0]  # Shape: (embedding_dim,)
            batch_embeddings.append(embedding)

        embeddings.extend(batch_embeddings)

    embeddings = np.array(embeddings)
    print(f"\nESM-C embeddings shape: {embeddings.shape}")
    print(f"  Expected: ({len(proteins)}, 960)")

    return embeddings


def compute_glm2_embeddings(
    proteins: List[Dict],
    batch_size: int = 32,
    device: str = "cuda",
) -> np.ndarray:
    """
    Compute gLM2 embeddings for proteins (protein-only, no genomic context).

    Returns:
        np.ndarray of shape (n_proteins, embedding_dim)
    """
    print("\n" + "="*80)
    print("Computing gLM2 Embeddings (Protein-Only)")
    print("="*80)

    from transformers import AutoModel, AutoTokenizer

    print("Loading gLM2 model...")
    model = AutoModel.from_pretrained(
        "tattabio/gLM2_650M",
        torch_dtype=torch.bfloat16,
        trust_remote_code=True,
    )
    model = model.to(device)
    model.eval()

    tokenizer = AutoTokenizer.from_pretrained(
        "tattabio/gLM2_650M",
        trust_remote_code=True,
    )

    print(f"Model loaded on {device}")
    # gLM2 config may not have hidden_size attribute, get it from first forward pass or use default
    hidden_dim = getattr(model.config, 'hidden_size', 1280)
    print(f"Expected embedding dimension: {hidden_dim}")

    embeddings = []

    print(f"\nProcessing {len(proteins)} proteins in batches of {batch_size}...")

    for i in tqdm(range(0, len(proteins), batch_size)):
        batch = proteins[i:i + batch_size]

        # Format sequences for gLM2 (protein-only, use <+> strand indicator)
        batch_seqs = [f"<+>{prot['sequence']}" for prot in batch]

        # Tokenize
        inputs = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=1024,
        )
        # Remove token_type_ids if present (gLM2 doesn't use them)
        if 'token_type_ids' in inputs:
            del inputs['token_type_ids']
        inputs = {k: v.to(device) for k, v in inputs.items()}

        # Generate embeddings
        with torch.no_grad():
            outputs = model(**inputs)
            # Get mean pooled embeddings
            # outputs.last_hidden_state shape: (batch_size, seq_len, hidden_size)
            batch_embeddings = outputs.last_hidden_state.mean(dim=1)
            batch_embeddings = batch_embeddings.cpu().float().numpy()

        embeddings.append(batch_embeddings)

    embeddings = np.vstack(embeddings)
    print(f"\ngLM2 embeddings shape: {embeddings.shape}")
    print(f"  Expected: ({len(proteins)}, {hidden_dim})")

    return embeddings


def compute_distance_matrix(
    embeddings: np.ndarray,
    metric: str = "cosine",
) -> np.ndarray:
    """
    Compute pairwise distance matrix.

    Args:
        embeddings: Array of shape (n_samples, embedding_dim)
        metric: Distance metric (cosine, euclidean, correlation)

    Returns:
        Distance matrix of shape (n_samples, n_samples)
    """
    print(f"\nComputing {metric} distance matrix...")
    distances = pdist(embeddings, metric=metric)
    distance_matrix = squareform(distances)
    print(f"  Distance matrix shape: {distance_matrix.shape}")
    return distance_matrix


def save_results(
    proteins: List[Dict],
    esmc_embeddings: np.ndarray,
    glm2_embeddings: np.ndarray,
    esmc_distances: np.ndarray,
    glm2_distances: np.ndarray,
    output_dir: Path,
):
    """Save embeddings and distance matrices."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    print("\n" + "="*80)
    print("Saving Results")
    print("="*80)

    # Save metadata
    metadata_df = pd.DataFrame(proteins)
    metadata_file = output_dir / "protein_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    print(f"  Saved metadata: {metadata_file}")

    # Save ESM-C results
    esmc_emb_file = output_dir / "esmc_embeddings.npy"
    np.save(esmc_emb_file, esmc_embeddings)
    print(f"  Saved ESM-C embeddings: {esmc_emb_file}")

    esmc_dist_file = output_dir / "esmc_distance_matrix.npy"
    np.save(esmc_dist_file, esmc_distances)
    print(f"  Saved ESM-C distances: {esmc_dist_file}")

    # Save gLM2 results
    glm2_emb_file = output_dir / "glm2_protein_only_embeddings.npy"
    np.save(glm2_emb_file, glm2_embeddings)
    print(f"  Saved gLM2 embeddings: {glm2_emb_file}")

    glm2_dist_file = output_dir / "glm2_protein_only_distance_matrix.npy"
    np.save(glm2_dist_file, glm2_distances)
    print(f"  Saved gLM2 distances: {glm2_dist_file}")

    # Summary statistics
    print("\n" + "="*80)
    print("Summary Statistics")
    print("="*80)
    print(f"  Number of proteins: {len(proteins)}")
    print(f"\n  ESM-C:")
    print(f"    Embedding dimension: {esmc_embeddings.shape[1]}")
    print(f"    Mean pairwise distance: {esmc_distances[np.triu_indices_from(esmc_distances, k=1)].mean():.4f}")
    print(f"    Std pairwise distance: {esmc_distances[np.triu_indices_from(esmc_distances, k=1)].std():.4f}")

    print(f"\n  gLM2 (protein-only):")
    print(f"    Embedding dimension: {glm2_embeddings.shape[1]}")
    print(f"    Mean pairwise distance: {glm2_distances[np.triu_indices_from(glm2_distances, k=1)].mean():.4f}")
    print(f"    Std pairwise distance: {glm2_distances[np.triu_indices_from(glm2_distances, k=1)].std():.4f}")


def main():
    args = parse_args()

    print("="*80)
    print("Compute Protein-Only Embeddings: ESM-C vs gLM2")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Protein directory: {args.protein_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Device: {args.device}")

    # Check device
    if args.device == "cuda":
        if not torch.cuda.is_available():
            print("\n  WARNING: CUDA not available, falling back to CPU")
            args.device = "cpu"
        else:
            print(f"\n  GPU: {torch.cuda.get_device_name(0)}")
            print(f"  GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    # Load proteins
    protein_file = Path(args.protein_dir) / "proteins.faa"
    print(f"\nLoading proteins from {protein_file}...")
    proteins = load_proteins(protein_file)
    print(f"  Loaded {len(proteins)} proteins")
    print(f"  Avg length: {np.mean([p['length'] for p in proteins]):.1f} aa")

    # Compute ESM-C embeddings
    esmc_embeddings = compute_esmc_embeddings(
        proteins,
        batch_size=args.batch_size,
        device=args.device,
    )

    # Compute ESM-C distances
    esmc_distances = compute_distance_matrix(esmc_embeddings, metric="cosine")

    # Compute gLM2 embeddings (protein-only)
    glm2_embeddings = compute_glm2_embeddings(
        proteins,
        batch_size=args.batch_size,
        device=args.device,
    )

    # Compute gLM2 distances
    glm2_distances = compute_distance_matrix(glm2_embeddings, metric="cosine")

    # Save results
    save_results(
        proteins,
        esmc_embeddings,
        glm2_embeddings,
        esmc_distances,
        glm2_distances,
        Path(args.output_dir),
    )

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)
    print(f"\nResults saved to: {args.output_dir}")
    print("\nNext steps:")
    print("1. Compute gLM2 embeddings WITH genomic context")
    print("2. Compare distance matrices (protein-only vs with context)")
    print("3. Analyze correlation between ESM-C and gLM2 distances")


if __name__ == "__main__":
    main()
