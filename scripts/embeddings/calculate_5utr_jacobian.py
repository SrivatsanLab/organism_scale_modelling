#!/usr/bin/env python3
"""
Calculate categorical Jacobian for 5' UTR sequences using gLM2.

This script computes how mutations at each position in the 5' UTR affect
the embedding at every other position in the 5' UTR.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO


def load_sequence_data(protein_dir: Path, num_proteins: int = 5) -> List[Tuple[str, str]]:
    """
    Load 5' UTR sequences from FASTA file.

    Returns:
        List of tuples (seq_id, dna_5_sequence)
    """
    data = []

    # Load from flanking_5prime.fna file
    fasta_file = protein_dir / "flanking_5prime.fna"

    if not fasta_file.exists():
        return data

    records = list(SeqIO.parse(fasta_file, "fasta"))[:num_proteins]

    for record in records:
        seq_id = record.id
        dna_5_seq = str(record.seq)
        data.append((seq_id, dna_5_seq))

    return data


def calculate_glm2_5utr_jacobian(model, tokenizer, dna_sequence: str, device: str) -> np.ndarray:
    """
    Calculate categorical Jacobian for 5' UTR sequence.

    For each position j, mutate to each of 4 nucleotides and measure
    the change in embedding at each position i.

    Args:
        model: gLM2 model
        tokenizer: gLM2 tokenizer
        dna_sequence: 5' UTR DNA sequence
        device: Device to run on

    Returns:
        Jacobian matrix (L x L) where L is sequence length
    """
    L = len(dna_sequence)
    D = 1280  # gLM2 embedding dimension

    # Tokenize and get baseline embedding
    dna_input = f"<|>{dna_sequence}"
    baseline_tokens = tokenizer(dna_input, return_tensors="pt")
    # Remove token_type_ids if present (gLM2 doesn't use it)
    baseline_tokens = {k: v.to(device) for k, v in baseline_tokens.items() if k != 'token_type_ids'}

    with torch.no_grad():
        baseline_output = model(**baseline_tokens)

    # Convert BFloat16 to float32 before numpy
    baseline_emb = baseline_output.last_hidden_state[0].float().cpu().numpy()

    # Initialize Jacobian
    jacobian = np.zeros((L, L))

    nucleotides = 'ACGT'

    # For each position j (column - mutated position)
    for j in range(L):
        original_nt = dna_sequence[j]

        # Try each nucleotide mutation
        mutation_deltas = []

        for new_nt in nucleotides:
            if new_nt == original_nt:
                continue

            # Create mutated sequence
            mutated_seq = dna_sequence[:j] + new_nt + dna_sequence[j+1:]
            mutated_input = f"<|>{mutated_seq}"
            mutated_tokens = tokenizer(mutated_input, return_tensors="pt")
            # Remove token_type_ids if present
            mutated_tokens = {k: v.to(device) for k, v in mutated_tokens.items() if k != 'token_type_ids'}

            with torch.no_grad():
                mutated_output = model(**mutated_tokens)

            mutated_emb = mutated_output.last_hidden_state[0].float().cpu().numpy()

            # Calculate L2 norm of embedding change at each position
            delta = np.linalg.norm(mutated_emb - baseline_emb, axis=1)
            mutation_deltas.append(delta)

        # Average across all mutations at position j
        avg_delta = np.mean(mutation_deltas, axis=0)

        # Store in column j
        jacobian[:, j] = avg_delta

    return jacobian


def analyze_jacobian_statistics(jacobian: np.ndarray) -> dict:
    """Calculate summary statistics for a Jacobian matrix."""
    return {
        'mean': jacobian.mean(),
        'std': jacobian.std(),
        'max': jacobian.max(),
        'min': jacobian.min(),
        'median': np.median(jacobian),
    }


def plot_jacobian_heatmap(jacobian: np.ndarray,
                         output_path: Path,
                         title: str,
                         sequence: str):
    """
    Plot Jacobian as a heatmap with log10 scaling.

    Args:
        jacobian: The Jacobian matrix (L x L)
        output_path: Path to save the figure
        title: Title for the plot
        sequence: The DNA sequence (for length info)
    """
    L = len(sequence)

    # Apply log10 transform for visualization
    jacobian_log = np.log10(jacobian + 1e-10)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot heatmap
    sns.heatmap(jacobian_log, cmap='viridis', ax=ax,
                cbar_kws={'label': 'log₁₀(Coupling strength)'})

    ax.set_title(f'{title}\n5\' UTR Jacobian (gLM2)\n{L} nucleotides',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Position j (mutated)', fontsize=11)
    ax.set_ylabel('Position i (measured)', fontsize=11)

    # Add statistics text
    stats = analyze_jacobian_statistics(jacobian)
    stats_text = f"Mean: {stats['mean']:.4f}\nMax: {stats['max']:.4f}"
    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    # Also save the raw Jacobian
    npy_path = output_path.with_suffix('.npy')
    np.save(npy_path, jacobian)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate 5' UTR Jacobian using gLM2"
    )
    parser.add_argument(
        "--protein-dir",
        type=str,
        default="data/protein_samples",
        help="Directory containing FASTA files with sequences"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_5utr",
        help="Output directory"
    )
    parser.add_argument(
        "--num-proteins",
        type=int,
        default=5,
        help="Number of proteins to process"
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to use (cuda or cpu)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("5' UTR Categorical Jacobian Analysis")
    print("gLM2 Model")
    print("=" * 80)
    print(f"Protein directory: {args.protein_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of proteins: {args.num_proteins}")
    print(f"Device: {args.device}")
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load sequence data
    protein_dir = Path(args.protein_dir)
    print(f"Loading sequences from {protein_dir}...")
    data = load_sequence_data(protein_dir, args.num_proteins)
    print(f"  Loaded {len(data)} sequences")
    print()

    # Load gLM2 model
    print("Loading gLM2 model...")
    from transformers import AutoModel, AutoTokenizer
    glm2_model = AutoModel.from_pretrained(
        "tattabio/gLM2_650M",
        torch_dtype=torch.bfloat16,
        trust_remote_code=True,
    )
    glm2_model = glm2_model.to(args.device)
    glm2_model.eval()

    glm2_tokenizer = AutoTokenizer.from_pretrained(
        "tattabio/gLM2_650M",
        trust_remote_code=True,
    )
    print(f"  gLM2 loaded on {args.device}")
    print()

    # Process each sequence
    for seq_id, dna_5_seq in data:
        print(f"Processing: {seq_id}")
        print(f"  5' UTR: {len(dna_5_seq)} nt")

        # Calculate 5' UTR Jacobian
        print("  Computing 5' UTR Jacobian...")
        jacobian_5utr = calculate_glm2_5utr_jacobian(
            glm2_model, glm2_tokenizer, dna_5_seq, args.device
        )

        stats = analyze_jacobian_statistics(jacobian_5utr)

        # Plot and save
        safe_id = seq_id.replace("|", "_")
        output_path = output_dir / f"{safe_id}_5utr_jacobian.png"

        plot_jacobian_heatmap(
            jacobian_5utr,
            output_path,
            seq_id,
            dna_5_seq
        )

        print(f"  Saved: {output_path}")
        print(f"    Mean coupling: {stats['mean']:.4f}")
        print(f"    Max coupling: {stats['max']:.4f}")
        print()

    print("=" * 80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
