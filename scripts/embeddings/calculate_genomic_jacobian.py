#!/usr/bin/env python
"""
Calculate categorical Jacobian for genomic sequences (DNA+protein) using gLM2.

For gLM2, the model can handle mixed DNA and protein sequences. This script computes
the Jacobian where:
- DNA positions are mutated to other DNA nucleotides (A, C, G, T)
- Protein positions are mutated to other amino acids (20 standard AAs)

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
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate categorical Jacobian for genomic sequences with gLM2"
    )
    parser.add_argument(
        "--protein-dir",
        type=str,
        default="data/protein_samples",
        help="Directory containing protein and genomic context files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_genomic",
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


def load_sequence_data(protein_dir: Path, num_proteins: int):
    """
    Load proteins and their corresponding genomic context sequences.

    Returns:
        List of tuples (seq_id, dna_sequence, protein_sequence, cds_start, cds_end)
    """
    # Load metadata to get CDS positions
    metadata = pd.read_csv(protein_dir / "protein_metadata.csv")
    metadata = metadata.head(num_proteins)

    # Load genomic context sequences (DNA)
    genomic_seqs = {}
    for record in SeqIO.parse(protein_dir / "genomic_context.fna", "fasta"):
        genomic_seqs[record.id] = str(record.seq)

    # Load protein sequences
    protein_seqs = {}
    for record in SeqIO.parse(protein_dir / "proteins.faa", "fasta"):
        protein_seqs[record.id] = str(record.seq)

    # Load CDS sequences to determine positions
    cds_seqs = {}
    for record in SeqIO.parse(protein_dir / "cds_sequences.fna", "fasta"):
        cds_seqs[record.id] = str(record.seq)

    # Combine data
    data = []
    for _, row in metadata.iterrows():
        seq_id = f"{row['genome_id']}|{row['gene_id']}"

        if seq_id in genomic_seqs and seq_id in protein_seqs and seq_id in cds_seqs:
            dna_seq = genomic_seqs[seq_id]
            protein_seq = protein_seqs[seq_id]
            cds_seq = cds_seqs[seq_id]

            # Find CDS start position in genomic context
            # The genomic context is: 5' flanking + CDS + 3' flanking
            cds_start = dna_seq.find(cds_seq)
            if cds_start == -1:
                # Try reverse complement for minus strand
                from Bio.Seq import Seq
                cds_seq_rc = str(Seq(cds_seq).reverse_complement())
                cds_start = dna_seq.find(cds_seq_rc)
                if cds_start != -1:
                    cds_seq = cds_seq_rc

            if cds_start != -1:
                cds_end = cds_start + len(cds_seq)
                data.append((seq_id, dna_seq, protein_seq, cds_seq, cds_start, cds_end))

    return data


def calculate_glm2_genomic_jacobian(
    model,
    tokenizer,
    dna_sequence: str,
    protein_sequence: str,
    cds_sequence: str,
    cds_start: int,
    cds_end: int,
    device: str
) -> np.ndarray:
    """
    Calculate categorical Jacobian for gLM2 with genomic context.

    Args:
        model: gLM2 model
        tokenizer: gLM2 tokenizer
        dna_sequence: Full genomic context (DNA)
        protein_sequence: Protein sequence
        cds_sequence: CDS sequence
        cds_start: Start position of CDS in genomic context
        cds_end: End position of CDS in genomic context
        device: Device for computation

    Returns:
        L x L Jacobian matrix where L is total sequence length (DNA + protein)
    """
    # For gLM2 with genomic context, we use format: <|>dna_5prime<+>protein<|>dna_3prime
    dna_5prime = dna_sequence[:cds_start]
    dna_3prime = dna_sequence[cds_end:]

    # Total sequence for gLM2
    full_seq = f"<|>{dna_5prime}<+>{protein_sequence}<|>{dna_3prime}"

    # Calculate total length (excluding special tokens for now)
    L_dna_5 = len(dna_5prime)
    L_protein = len(protein_sequence)
    L_dna_3 = len(dna_3prime)
    L_total = L_dna_5 + L_protein + L_dna_3

    jacobian = np.zeros((L_total, L_total))

    # Get baseline embeddings
    baseline_inputs = tokenizer(
        full_seq,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=2048,
    )
    if 'token_type_ids' in baseline_inputs:
        del baseline_inputs['token_type_ids']
    baseline_inputs = {k: v.to(device) for k, v in baseline_inputs.items()}

    with torch.no_grad():
        baseline_output = model(**baseline_inputs)
        baseline_emb = baseline_output.last_hidden_state[0].float().cpu().numpy()

    # Define mutation vocabularies
    dna_nucleotides = 'ACGT'
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # Iterate through each position in the full sequence
    for j in tqdm(range(L_total), desc="Computing Jacobian", leave=False):
        # Determine if this is DNA or protein position
        if j < L_dna_5:
            # 5' DNA flanking region
            pos_in_dna5 = j
            original_char = dna_5prime[pos_in_dna5]
            vocabulary = dna_nucleotides
            seq_type = "dna5"
        elif j < L_dna_5 + L_protein:
            # Protein coding region
            pos_in_protein = j - L_dna_5
            original_char = protein_sequence[pos_in_protein]
            vocabulary = amino_acids
            seq_type = "protein"
        else:
            # 3' DNA flanking region
            pos_in_dna3 = j - L_dna_5 - L_protein
            original_char = dna_3prime[pos_in_dna3]
            vocabulary = dna_nucleotides
            seq_type = "dna3"

        mutation_deltas = []

        # Try each possible character from the appropriate vocabulary
        for new_char in vocabulary:
            if new_char == original_char:
                continue

            # Create mutated sequence
            if seq_type == "dna5":
                mut_dna5 = dna_5prime[:pos_in_dna5] + new_char + dna_5prime[pos_in_dna5+1:]
                mutated_seq = f"<|>{mut_dna5}<+>{protein_sequence}<|>{dna_3prime}"
            elif seq_type == "protein":
                mut_protein = protein_sequence[:pos_in_protein] + new_char + protein_sequence[pos_in_protein+1:]
                mutated_seq = f"<|>{dna_5prime}<+>{mut_protein}<|>{dna_3prime}"
            else:  # dna3
                mut_dna3 = dna_3prime[:pos_in_dna3] + new_char + dna_3prime[pos_in_dna3+1:]
                mutated_seq = f"<|>{dna_5prime}<+>{protein_sequence}<|>{mut_dna3}"

            # Get embeddings for mutated sequence
            mut_inputs = tokenizer(
                mutated_seq,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=2048,
            )
            if 'token_type_ids' in mut_inputs:
                del mut_inputs['token_type_ids']
            mut_inputs = {k: v.to(device) for k, v in mut_inputs.items()}

            with torch.no_grad():
                mut_output = model(**mut_inputs)
                mut_emb = mut_output.last_hidden_state[0].float().cpu().numpy()

            # Compute change in embedding at each position
            min_len = min(baseline_emb.shape[0], mut_emb.shape[0])
            delta = np.linalg.norm(mut_emb[:min_len] - baseline_emb[:min_len], axis=1)

            # Pad if necessary
            if len(delta) < baseline_emb.shape[0]:
                delta = np.pad(delta, (0, baseline_emb.shape[0] - len(delta)))

            mutation_deltas.append(delta)

        # Average across all mutations at position j
        if mutation_deltas:
            avg_delta = np.mean(mutation_deltas, axis=0)

            # Map back to sequence positions (accounting for special tokens)
            # The tokenization is: <|> dna_5prime <+> protein <|> dna_3prime
            # Special tokens: <|>, <+>, <|> = 3 tokens
            # We need to extract the middle L_total positions

            # Calculate token offset - gLM2 uses special tokens
            # Typically: <|> (1 token) + sequence + <+> (1 token) + sequence + <|> (1 token)
            # For safety, extract positions based on expected length
            expected_len = baseline_emb.shape[0]
            if expected_len >= L_total:
                # Extract sequence positions (skip first special token)
                start_offset = 1  # After <|>
                jacobian[:, j] = avg_delta[start_offset:start_offset+L_total]

    return jacobian


def plot_jacobian_heatmap(
    jacobian: np.ndarray,
    output_path: Path,
    title: str,
    L_dna_5: int,
    L_protein: int,
    L_dna_3: int
):
    """Plot Jacobian matrix as heatmap with regions marked and log10 scale."""
    fig, ax = plt.subplots(figsize=(14, 12))

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

    # Add lines to demarcate regions
    L_total = L_dna_5 + L_protein + L_dna_3

    # Vertical lines
    ax.axvline(x=L_dna_5, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax.axvline(x=L_dna_5+L_protein, color='red', linestyle='--', linewidth=2, alpha=0.7)

    # Horizontal lines
    ax.axhline(y=L_dna_5, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax.axhline(y=L_dna_5+L_protein, color='red', linestyle='--', linewidth=2, alpha=0.7)

    # Add region labels
    ax.text(L_dna_5/2, -20, "5' DNA", ha='center', fontsize=10, fontweight='bold')
    ax.text(L_dna_5 + L_protein/2, -20, "Protein", ha='center', fontsize=10, fontweight='bold')
    ax.text(L_dna_5 + L_protein + L_dna_3/2, -20, "3' DNA", ha='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('Mutated Position j', fontsize=12)
    ax.set_ylabel('Affected Position i', fontsize=12)
    ax.set_title(
        f'{title}\\n5\' DNA: {L_dna_5} nt | Protein: {L_protein} aa | 3\' DNA: {L_dna_3} nt',
        fontsize=14,
        fontweight='bold'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def analyze_jacobian_statistics(
    jacobian: np.ndarray,
    L_dna_5: int,
    L_protein: int,
    L_dna_3: int
) -> dict:
    """Compute summary statistics for Jacobian matrix by region."""
    # Define regions
    dna5_end = L_dna_5
    protein_end = L_dna_5 + L_protein

    # Extract sub-matrices
    dna5_dna5 = jacobian[:dna5_end, :dna5_end]
    dna5_protein = jacobian[:dna5_end, dna5_end:protein_end]
    dna5_dna3 = jacobian[:dna5_end, protein_end:]

    protein_dna5 = jacobian[dna5_end:protein_end, :dna5_end]
    protein_protein = jacobian[dna5_end:protein_end, dna5_end:protein_end]
    protein_dna3 = jacobian[dna5_end:protein_end, protein_end:]

    dna3_dna5 = jacobian[protein_end:, :dna5_end]
    dna3_protein = jacobian[protein_end:, dna5_end:protein_end]
    dna3_dna3 = jacobian[protein_end:, protein_end:]

    stats = {
        # Within-region coupling (off-diagonal)
        'dna5_dna5_mean': np.mean(dna5_dna5 - np.diag(np.diag(dna5_dna5))) if dna5_dna5.size > 0 else 0,
        'protein_protein_mean': np.mean(protein_protein - np.diag(np.diag(protein_protein))),
        'dna3_dna3_mean': np.mean(dna3_dna3 - np.diag(np.diag(dna3_dna3))) if dna3_dna3.size > 0 else 0,

        # Cross-region coupling
        'dna5_protein_mean': np.mean(dna5_protein) if dna5_protein.size > 0 else 0,
        'protein_dna5_mean': np.mean(protein_dna5) if protein_dna5.size > 0 else 0,
        'protein_dna3_mean': np.mean(protein_dna3) if protein_dna3.size > 0 else 0,
        'dna3_protein_mean': np.mean(dna3_protein) if dna3_protein.size > 0 else 0,

        # Overall statistics
        'mean_coupling': np.mean(jacobian),
        'max_coupling': np.max(jacobian),
        'frobenius_norm': np.linalg.norm(jacobian, 'fro'),
    }

    return stats


def main():
    args = parse_args()

    print("="*80)
    print("Categorical Jacobian Analysis for Genomic Sequences (DNA+Protein)")
    print("gLM2 Model")
    print("="*80)
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
    for seq_id, dna_seq, protein_seq, cds_seq, cds_start, cds_end in data:
        L_dna_5 = cds_start
        L_protein = len(protein_seq)
        L_dna_3 = len(dna_seq) - cds_end
        L_total = L_dna_5 + L_protein + L_dna_3

        print(f"Processing: {seq_id}")
        print(f"  5' DNA: {L_dna_5} nt | Protein: {L_protein} aa | 3' DNA: {L_dna_3} nt | Total: {L_total}")

        # Calculate gLM2 Jacobian
        print("  Computing gLM2 genomic Jacobian...")
        jacobian = calculate_glm2_genomic_jacobian(
            glm2_model, glm2_tokenizer,
            dna_seq, protein_seq, cds_seq,
            cds_start, cds_end,
            args.device
        )

        # Analyze statistics
        stats = analyze_jacobian_statistics(jacobian, L_dna_5, L_protein, L_dna_3)

        # Plot Jacobian
        safe_id = seq_id.replace("|", "_").replace("/", "_")
        plot_jacobian_heatmap(
            jacobian,
            output_dir / f"{safe_id}_glm2_genomic_jacobian.png",
            f"gLM2 Genomic Jacobian: {seq_id}",
            L_dna_5, L_protein, L_dna_3
        )

        # Save Jacobian
        np.save(output_dir / f"{safe_id}_glm2_genomic_jacobian.npy", jacobian)

        print(f"    Overall mean coupling: {stats['mean_coupling']:.4f}")
        print(f"    Overall max coupling: {stats['max_coupling']:.4f}")
        print(f"    Protein-protein coupling: {stats['protein_protein_mean']:.4f}")
        print(f"    DNA5-protein coupling: {stats['dna5_protein_mean']:.4f}")
        print(f"    Protein-DNA3 coupling: {stats['protein_dna3_mean']:.4f}")
        print()

    print("="*80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
