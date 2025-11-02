#!/usr/bin/env python
"""
Calculate categorical Jacobian for protein sequences using ESM-C and gLM2.

The categorical Jacobian measures positional dependencies in the model's learned
representations. For proteins, this can reveal:
- Coevolutionary constraints
- Structural dependencies
- Functional residue coupling

For each position i and j:
J[i,j] = how much the embedding at position i changes when position j is mutated

Args:
    Computes Jacobian by measuring embedding changes from single-position mutations.
"""

import argparse
import numpy as np
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate categorical Jacobian for protein sequences"
    )
    parser.add_argument(
        "--protein-dir",
        type=str,
        default="data/protein_samples",
        help="Directory containing protein FASTA file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_analysis",
        help="Output directory for results",
    )
    parser.add_argument(
        "--num-proteins",
        type=int,
        default=5,
        help="Number of proteins to analyze (default: 5)",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to use (cuda/cpu)",
    )
    return parser.parse_args()


def calculate_esmc_jacobian(model, sequence: str, device: str) -> np.ndarray:
    """
    Calculate categorical Jacobian for ESM-C.

    Computes how embedding at each position changes when other positions are mutated.
    J[i,j] = ||embedding(mutated at j)[i] - embedding(original)[i]||

    Args:
        model: ESM-C model
        sequence: Protein sequence string
        device: Device for computation

    Returns:
        L x L Jacobian matrix where L is sequence length
    """
    from esm.sdk.api import ESMProtein, LogitsConfig

    L = len(sequence)
    jacobian = np.zeros((L, L))

    # Get baseline embeddings
    protein = ESMProtein(sequence=sequence)
    protein_tensor = model.encode(protein)
    baseline_output = model.logits(
        protein_tensor,
        LogitsConfig(sequence=True, return_embeddings=True)
    )
    baseline_emb = baseline_output.embeddings.detach().cpu().numpy()[0]  # (seq_len, D)

    # ESM-C may add special tokens, extract only protein positions
    # Typically: [BOS, ...protein..., EOS] or similar
    # Use the middle L positions
    emb_len = baseline_emb.shape[0]
    if emb_len > L:
        # Extract middle L positions (skip special tokens)
        offset = (emb_len - L) // 2
        baseline_emb = baseline_emb[offset:offset+L]
    elif emb_len < L:
        raise ValueError(f"Embedding length {emb_len} < sequence length {L}")

    # Amino acid vocabulary (20 standard amino acids)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # For each position j, try mutations
    for j in tqdm(range(L), desc="Computing Jacobian", leave=False):
        original_aa = sequence[j]

        mutation_deltas = []

        # Try each possible amino acid at position j
        for new_aa in amino_acids:
            if new_aa == original_aa:
                continue

            # Create mutated sequence
            mutated_seq = sequence[:j] + new_aa + sequence[j+1:]

            # Get embeddings for mutated sequence
            mut_protein = ESMProtein(sequence=mutated_seq)
            mut_tensor = model.encode(mut_protein)
            mut_output = model.logits(
                mut_tensor,
                LogitsConfig(sequence=True, return_embeddings=True)
            )
            mut_emb = mut_output.embeddings.detach().cpu().numpy()[0]  # (seq_len, D)

            # Extract protein positions (same as baseline)
            if mut_emb.shape[0] > L:
                offset = (mut_emb.shape[0] - L) // 2
                mut_emb = mut_emb[offset:offset+L]

            # Compute change in embedding at each position i
            delta = np.linalg.norm(mut_emb - baseline_emb, axis=1)  # (L,)
            mutation_deltas.append(delta)

        # Average across all mutations at position j
        if mutation_deltas:
            avg_delta = np.mean(mutation_deltas, axis=0)  # (L,)
            jacobian[:, j] = avg_delta

    return jacobian


def calculate_glm2_jacobian(model, tokenizer, sequence: str, device: str) -> np.ndarray:
    """
    Calculate categorical Jacobian for gLM2.

    Args:
        model: gLM2 model
        tokenizer: gLM2 tokenizer
        sequence: Protein sequence string
        device: Device for computation

    Returns:
        L x L Jacobian matrix where L is sequence length
    """
    L = len(sequence)
    jacobian = np.zeros((L, L))

    # Get baseline embeddings
    baseline_seq = f"<+>{sequence}"
    baseline_inputs = tokenizer(
        baseline_seq,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=1024,
    )
    # Remove token_type_ids if present
    if 'token_type_ids' in baseline_inputs:
        del baseline_inputs['token_type_ids']
    baseline_inputs = {k: v.to(device) for k, v in baseline_inputs.items()}

    with torch.no_grad():
        baseline_output = model(**baseline_inputs)
        baseline_emb = baseline_output.last_hidden_state[0].float().cpu().numpy()  # (seq_len, D)

    # Amino acid vocabulary
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # For each position j in the protein sequence (accounting for special tokens)
    # Note: gLM2 tokenizer may split the sequence differently, so we need to be careful
    # For simplicity, we'll assume 1:1 mapping after the special tokens
    token_offset = 1  # Offset for <+> token

    for j in tqdm(range(L), desc="Computing Jacobian", leave=False):
        original_aa = sequence[j]

        mutation_deltas = []

        # Try each possible amino acid at position j
        for new_aa in amino_acids:
            if new_aa == original_aa:
                continue

            # Create mutated sequence
            mutated_seq = f"<+>{sequence[:j]}{new_aa}{sequence[j+1:]}"

            # Get embeddings for mutated sequence
            mut_inputs = tokenizer(
                mutated_seq,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=1024,
            )
            if 'token_type_ids' in mut_inputs:
                del mut_inputs['token_type_ids']
            mut_inputs = {k: v.to(device) for k, v in mut_inputs.items()}

            with torch.no_grad():
                mut_output = model(**mut_inputs)
                mut_emb = mut_output.last_hidden_state[0].float().cpu().numpy()  # (seq_len, D)

            # Compute change in embedding at each position
            # Handle potential length mismatches
            min_len = min(baseline_emb.shape[0], mut_emb.shape[0])
            delta = np.linalg.norm(mut_emb[:min_len] - baseline_emb[:min_len], axis=1)

            # Pad if necessary
            if len(delta) < baseline_emb.shape[0]:
                delta = np.pad(delta, (0, baseline_emb.shape[0] - len(delta)))

            mutation_deltas.append(delta)

        # Average across all mutations at position j
        if mutation_deltas:
            avg_delta = np.mean(mutation_deltas, axis=0)
            # Map back to protein positions (skip special tokens)
            if len(avg_delta) > token_offset:
                jacobian[:, j] = avg_delta[token_offset:token_offset+L]

    return jacobian


def plot_jacobian_heatmap(jacobian: np.ndarray, output_path: Path, title: str, sequence: str):
    """Plot Jacobian matrix as heatmap with log10 scale."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Apply log10 transform, adding small epsilon to avoid log(0)
    jacobian_log = np.log10(jacobian + 1e-10)

    sns.heatmap(
        jacobian_log,
        cmap='viridis',
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'log10(Embedding Change Magnitude)'},
        ax=ax
    )

    ax.set_xlabel('Mutated Position j', fontsize=12)
    ax.set_ylabel('Affected Position i', fontsize=12)
    ax.set_title(f'{title}\nSequence length: {len(sequence)} aa', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def analyze_jacobian_statistics(jacobian: np.ndarray) -> dict:
    """Compute summary statistics for Jacobian matrix."""
    # Remove diagonal (self-dependencies)
    off_diagonal = jacobian.copy()
    np.fill_diagonal(off_diagonal, 0)

    stats = {
        'mean_coupling': np.mean(off_diagonal),
        'max_coupling': np.max(off_diagonal),
        'median_coupling': np.median(off_diagonal),
        'std_coupling': np.std(off_diagonal),
        'mean_diagonal': np.mean(np.diag(jacobian)),
        'frobenius_norm': np.linalg.norm(jacobian, 'fro'),
    }

    return stats


def main():
    args = parse_args()

    print("="*80)
    print("Categorical Jacobian Analysis for Protein Sequences")
    print("="*80)
    print(f"Protein directory: {args.protein_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of proteins: {args.num_proteins}")
    print(f"Device: {args.device}")
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load proteins
    protein_file = Path(args.protein_dir) / "proteins.faa"
    print(f"Loading proteins from {protein_file}...")
    proteins = []
    for i, record in enumerate(SeqIO.parse(protein_file, "fasta")):
        if i >= args.num_proteins:
            break
        proteins.append((record.id, str(record.seq)))
    print(f"  Loaded {len(proteins)} proteins")
    print()

    # Load ESM-C model
    print("Loading ESM-C model...")
    from esm.models.esmc import ESMC
    device_obj = torch.device(args.device)
    esmc_model = ESMC.from_pretrained("esmc_300m", device=device_obj)
    print(f"  ESM-C loaded on {args.device}")
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

    # Process each protein
    for seq_id, sequence in proteins:
        print(f"Processing: {seq_id} ({len(sequence)} aa)")

        # Calculate ESM-C Jacobian
        print("  Computing ESM-C Jacobian...")
        esmc_jacobian = calculate_esmc_jacobian(esmc_model, sequence, args.device)
        esmc_stats = analyze_jacobian_statistics(esmc_jacobian)

        # Plot ESM-C Jacobian
        safe_id = seq_id.replace("|", "_").replace("/", "_")
        plot_jacobian_heatmap(
            esmc_jacobian,
            output_dir / f"{safe_id}_esmc_jacobian.png",
            f"ESM-C Categorical Jacobian: {seq_id}",
            sequence
        )

        # Save ESM-C Jacobian
        np.save(output_dir / f"{safe_id}_esmc_jacobian.npy", esmc_jacobian)

        print(f"    Mean coupling: {esmc_stats['mean_coupling']:.4f}")
        print(f"    Max coupling: {esmc_stats['max_coupling']:.4f}")
        print()

        # Calculate gLM2 Jacobian
        print("  Computing gLM2 Jacobian...")
        glm2_jacobian = calculate_glm2_jacobian(glm2_model, glm2_tokenizer, sequence, args.device)
        glm2_stats = analyze_jacobian_statistics(glm2_jacobian)

        # Plot gLM2 Jacobian
        plot_jacobian_heatmap(
            glm2_jacobian,
            output_dir / f"{safe_id}_glm2_jacobian.png",
            f"gLM2 Categorical Jacobian: {seq_id}",
            sequence
        )

        # Save gLM2 Jacobian
        np.save(output_dir / f"{safe_id}_glm2_jacobian.npy", glm2_jacobian)

        print(f"    Mean coupling: {glm2_stats['mean_coupling']:.4f}")
        print(f"    Max coupling: {glm2_stats['max_coupling']:.4f}")
        print()

    print("="*80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
