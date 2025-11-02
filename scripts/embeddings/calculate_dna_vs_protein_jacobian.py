#!/usr/bin/env python
"""
Compare categorical Jacobians for DNA-encoded vs protein-encoded sequences using gLM2.

For the same protein, compute Jacobians for:
1. CDS (DNA) representation - mutate codons to synonymous/non-synonymous variants
2. Protein representation - mutate amino acids to other amino acids

This reveals how gLM2's learned representations differ between DNA and protein encodings.

Args:
    Computes Jacobians by measuring embedding changes from single-position mutations.
"""

import argparse
import numpy as np
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare DNA vs Protein Jacobians using gLM2"
    )
    parser.add_argument(
        "--protein-dir",
        type=str,
        default="data/protein_samples",
        help="Directory containing protein and CDS files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_dna_vs_protein",
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
    Load proteins and their corresponding CDS sequences.

    Returns:
        List of tuples (seq_id, protein_sequence, cds_sequence)
    """
    # Load protein sequences
    protein_seqs = {}
    for record in SeqIO.parse(protein_dir / "proteins.faa", "fasta"):
        protein_seqs[record.id] = str(record.seq)

    # Load CDS sequences
    cds_seqs = {}
    for record in SeqIO.parse(protein_dir / "cds_sequences.fna", "fasta"):
        cds_seqs[record.id] = str(record.seq)

    # Combine data
    data = []
    for seq_id in list(protein_seqs.keys())[:num_proteins]:
        if seq_id in cds_seqs:
            data.append((seq_id, protein_seqs[seq_id], cds_seqs[seq_id]))

    return data


def calculate_glm2_protein_jacobian(
    model,
    tokenizer,
    protein_sequence: str,
    device: str
) -> np.ndarray:
    """
    Calculate categorical Jacobian for protein sequence.

    Args:
        model: gLM2 model
        tokenizer: gLM2 tokenizer
        protein_sequence: Protein sequence string
        device: Device for computation

    Returns:
        L x L Jacobian matrix where L is protein length
    """
    L = len(protein_sequence)
    jacobian = np.zeros((L, L))

    # Get baseline embeddings
    baseline_seq = f"<+>{protein_sequence}"
    baseline_inputs = tokenizer(
        baseline_seq,
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

    # Amino acid vocabulary
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # For each position j in the protein sequence
    token_offset = 1  # Offset for <+> token

    for j in tqdm(range(L), desc="Computing Protein Jacobian", leave=False):
        original_aa = protein_sequence[j]

        mutation_deltas = []

        # Try each possible amino acid at position j
        for new_aa in amino_acids:
            if new_aa == original_aa:
                continue

            # Create mutated sequence
            mutated_seq = f"<+>{protein_sequence[:j]}{new_aa}{protein_sequence[j+1:]}"

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
            # Map back to protein positions (skip special tokens)
            if len(avg_delta) > token_offset:
                jacobian[:, j] = avg_delta[token_offset:token_offset+L]

    return jacobian


def calculate_glm2_dna_jacobian(
    model,
    tokenizer,
    cds_sequence: str,
    protein_sequence: str,
    device: str
) -> np.ndarray:
    """
    Calculate categorical Jacobian for CDS (DNA) sequence, translated to amino acid level.

    For each amino acid position, we:
    1. Mutate the corresponding codon to all possible amino acid changes
    2. Average the embedding changes across all synonymous codons for each AA change
    3. This makes it directly comparable to the protein Jacobian

    Args:
        model: gLM2 model
        tokenizer: gLM2 tokenizer
        cds_sequence: CDS nucleotide sequence string
        protein_sequence: Protein sequence (for translation verification)
        device: Device for computation

    Returns:
        L_protein x L_protein Jacobian matrix at amino acid level
    """
    L_protein = len(protein_sequence)
    jacobian = np.zeros((L_protein, L_protein))

    # Get baseline embeddings for DNA
    baseline_seq = f"<|>{cds_sequence}"
    baseline_inputs = tokenizer(
        baseline_seq,
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

    # Extract nucleotide-level embeddings (skip special token)
    token_offset = 1
    L_nt = len(cds_sequence)
    if baseline_emb.shape[0] > token_offset + L_nt:
        baseline_nt_emb = baseline_emb[token_offset:token_offset+L_nt]
    else:
        baseline_nt_emb = baseline_emb[token_offset:]

    # Amino acid vocabulary
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    # For each amino acid position
    for aa_idx in tqdm(range(L_protein), desc="Computing DNA Jacobian", leave=False):
        # Get the corresponding codon
        codon_start = aa_idx * 3
        original_codon = cds_sequence[codon_start:codon_start+3]
        original_aa = str(Seq(original_codon).translate())

        # Verify translation matches
        if original_aa != protein_sequence[aa_idx]:
            print(f"  Warning: Translation mismatch at position {aa_idx}: {original_aa} != {protein_sequence[aa_idx]}")

        # For each target amino acid
        aa_mutation_deltas = {}

        # Generate all possible single-nucleotide mutations of this codon
        for pos_in_codon in range(3):
            nt_pos = codon_start + pos_in_codon
            original_nt = cds_sequence[nt_pos]

            for new_nt in 'ACGT':
                if new_nt == original_nt:
                    continue

                # Create mutated codon
                mutated_codon = list(original_codon)
                mutated_codon[pos_in_codon] = new_nt
                mutated_codon = ''.join(mutated_codon)

                # Translate to amino acid
                try:
                    mutated_aa = str(Seq(mutated_codon).translate())
                    if mutated_aa == '*':  # Stop codon
                        continue
                except:
                    continue

                # Create full mutated sequence
                mutated_seq = f"<|>{cds_sequence[:nt_pos]}{new_nt}{cds_sequence[nt_pos+1:]}"

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

                # Extract nucleotide-level embeddings
                if mut_emb.shape[0] > token_offset + L_nt:
                    mut_nt_emb = mut_emb[token_offset:token_offset+L_nt]
                else:
                    mut_nt_emb = mut_emb[token_offset:]

                # Compute change at nucleotide level
                min_len = min(baseline_nt_emb.shape[0], mut_nt_emb.shape[0])
                nt_delta = np.linalg.norm(mut_nt_emb[:min_len] - baseline_nt_emb[:min_len], axis=1)

                # Pad if necessary
                if len(nt_delta) < L_nt:
                    nt_delta = np.pad(nt_delta, (0, L_nt - len(nt_delta)))

                # Group by resulting amino acid
                if mutated_aa not in aa_mutation_deltas:
                    aa_mutation_deltas[mutated_aa] = []
                aa_mutation_deltas[mutated_aa].append(nt_delta)

        # Average across all mutations for each target amino acid
        # Then average across all target amino acids
        all_aa_deltas = []
        for target_aa in amino_acids:
            if target_aa == original_aa:
                continue
            if target_aa in aa_mutation_deltas:
                # Average all synonymous codon changes that result in this AA
                avg_for_aa = np.mean(aa_mutation_deltas[target_aa], axis=0)
                all_aa_deltas.append(avg_for_aa)

        # Average across all amino acid targets
        if all_aa_deltas:
            avg_delta = np.mean(all_aa_deltas, axis=0)

            # Convert nucleotide-level deltas to amino acid-level (average every 3 nt)
            aa_deltas = np.array([
                np.mean(avg_delta[i*3:min((i+1)*3, len(avg_delta))])
                for i in range(L_protein)
            ])

            jacobian[:, aa_idx] = aa_deltas

    return jacobian


def plot_comparison_heatmap(
    protein_jacobian: np.ndarray,
    dna_jacobian: np.ndarray,
    output_path: Path,
    title: str,
    protein_sequence: str
):
    """Plot side-by-side comparison of DNA and protein Jacobians with log10 scale."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 9))

    # Apply log10 transform
    protein_log = np.log10(protein_jacobian + 1e-10)
    dna_log = np.log10(dna_jacobian + 1e-10)

    # Find common color scale for both
    vmin = min(protein_log.min(), dna_log.min())
    vmax = max(protein_log.max(), dna_log.max())

    # Plot protein Jacobian
    sns.heatmap(
        protein_log,
        cmap='viridis',
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'log10(Embedding Change)'},
        ax=axes[0],
        vmin=vmin,
        vmax=vmax
    )
    axes[0].set_xlabel('Mutated Position j (AA)', fontsize=12)
    axes[0].set_ylabel('Affected Position i (AA)', fontsize=12)
    axes[0].set_title(f'Protein Encoding\n({len(protein_sequence)} AA)', fontsize=14, fontweight='bold')

    # Plot DNA Jacobian
    sns.heatmap(
        dna_log,
        cmap='viridis',
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'log10(Embedding Change)'},
        ax=axes[1],
        vmin=vmin,
        vmax=vmax
    )
    axes[1].set_xlabel('Mutated Codon j', fontsize=12)
    axes[1].set_ylabel('Affected Codon i', fontsize=12)
    axes[1].set_title(f'DNA Encoding (CDS)\n({len(protein_sequence)} codons)', fontsize=14, fontweight='bold')

    plt.suptitle(f'{title}\ngLM2: DNA vs Protein Jacobian Comparison',
                 fontsize=16, fontweight='bold', y=0.98)

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
    print("DNA vs Protein Jacobian Comparison")
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
    for seq_id, protein_seq, cds_seq in data:
        print(f"Processing: {seq_id} ({len(protein_seq)} aa, {len(cds_seq)} nt)")

        # Calculate Protein Jacobian
        print("  Computing Protein Jacobian...")
        protein_jacobian = calculate_glm2_protein_jacobian(
            glm2_model, glm2_tokenizer, protein_seq, args.device
        )
        protein_stats = analyze_jacobian_statistics(protein_jacobian)

        # Calculate DNA Jacobian (translated to amino acid level)
        print("  Computing DNA (CDS) Jacobian (translating to AA level)...")
        dna_jacobian = calculate_glm2_dna_jacobian(
            glm2_model, glm2_tokenizer, cds_seq, protein_seq, args.device
        )
        dna_stats = analyze_jacobian_statistics(dna_jacobian)

        # Plot comparison
        safe_id = seq_id.replace("|", "_").replace("/", "_")
        plot_comparison_heatmap(
            protein_jacobian,
            dna_jacobian,
            output_dir / f"{safe_id}_dna_vs_protein_comparison.png",
            seq_id,
            protein_seq
        )

        # Save Jacobians
        np.save(output_dir / f"{safe_id}_protein_jacobian.npy", protein_jacobian)
        np.save(output_dir / f"{safe_id}_dna_jacobian.npy", dna_jacobian)

        print(f"    Protein encoding - Mean: {protein_stats['mean_coupling']:.4f}, Max: {protein_stats['max_coupling']:.4f}")
        print(f"    DNA encoding     - Mean: {dna_stats['mean_coupling']:.4f}, Max: {dna_stats['max_coupling']:.4f}")
        print()

    print("="*80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
