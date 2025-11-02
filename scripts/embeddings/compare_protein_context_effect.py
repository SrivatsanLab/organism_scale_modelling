#!/usr/bin/env python3
"""
Compare protein Jacobians with and without genomic context.

This script extracts the protein-protein sub-matrix from the genomic Jacobian
and compares it side-by-side with the protein-only Jacobian to see if flanking
DNA context affects internal protein coupling patterns.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple


def load_jacobian(filepath: Path) -> np.ndarray:
    """Load a Jacobian matrix from a numpy file."""
    return np.load(filepath)


def extract_protein_submatrix(genomic_jacobian: np.ndarray,
                              L_dna_5: int,
                              L_protein: int) -> np.ndarray:
    """
    Extract the protein-protein sub-matrix from a genomic Jacobian.

    Args:
        genomic_jacobian: Full genomic Jacobian (DNA5 + Protein + DNA3)
        L_dna_5: Length of 5' DNA region
        L_protein: Length of protein region

    Returns:
        Protein-protein sub-matrix (L_protein x L_protein)
    """
    start = L_dna_5
    end = L_dna_5 + L_protein
    return genomic_jacobian[start:end, start:end]


def plot_comparison(protein_only: np.ndarray,
                   protein_from_genomic: np.ndarray,
                   seq_id: str,
                   output_path: Path):
    """
    Plot side-by-side comparison of protein Jacobians with and without context.

    Args:
        protein_only: Jacobian from protein-only encoding
        protein_from_genomic: Protein sub-matrix from genomic encoding
        seq_id: Sequence identifier
        output_path: Path to save the figure
    """
    L = protein_only.shape[0]

    # Apply log10 transform
    protein_only_log = np.log10(protein_only + 1e-10)
    protein_genomic_log = np.log10(protein_from_genomic + 1e-10)

    # Calculate difference
    difference = protein_from_genomic - protein_only
    difference_log = protein_genomic_log - protein_only_log

    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Determine shared color scale for log-transformed heatmaps
    vmin = min(protein_only_log.min(), protein_genomic_log.min())
    vmax = max(protein_only_log.max(), protein_genomic_log.max())

    # Plot 1: Protein-only Jacobian
    ax1 = axes[0, 0]
    sns.heatmap(protein_only_log, cmap='viridis', ax=ax1,
                vmin=vmin, vmax=vmax, cbar_kws={'label': 'log₁₀(Coupling)'})
    ax1.set_title('Protein-Only Encoding\n(gLM2)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Position j (mutated)', fontsize=11)
    ax1.set_ylabel('Position i (measured)', fontsize=11)

    # Add statistics
    mean_val = protein_only.mean()
    max_val = protein_only.max()
    ax1.text(0.02, 0.98, f'Mean: {mean_val:.2f}\nMax: {max_val:.2f}',
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             fontsize=10)

    # Plot 2: Protein from genomic context
    ax2 = axes[0, 1]
    sns.heatmap(protein_genomic_log, cmap='viridis', ax=ax2,
                vmin=vmin, vmax=vmax, cbar_kws={'label': 'log₁₀(Coupling)'})
    ax2.set_title('Protein with Genomic Context\n(gLM2: 5\'UTR-Protein-3\'UTR)',
                  fontsize=14, fontweight='bold')
    ax2.set_xlabel('Position j (mutated)', fontsize=11)
    ax2.set_ylabel('Position i (measured)', fontsize=11)

    # Add statistics
    mean_val = protein_from_genomic.mean()
    max_val = protein_from_genomic.max()
    ax2.text(0.02, 0.98, f'Mean: {mean_val:.2f}\nMax: {max_val:.2f}',
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             fontsize=10)

    # Plot 3: Absolute difference
    ax3 = axes[1, 0]
    max_abs = max(abs(difference.min()), abs(difference.max()))
    sns.heatmap(difference, cmap='RdBu_r', ax=ax3, center=0,
                vmin=-max_abs, vmax=max_abs,
                cbar_kws={'label': 'Difference (Genomic - Protein)'})
    ax3.set_title('Absolute Difference\n(Genomic Context - Protein Only)',
                  fontsize=14, fontweight='bold')
    ax3.set_xlabel('Position j (mutated)', fontsize=11)
    ax3.set_ylabel('Position i (measured)', fontsize=11)

    # Add statistics
    mean_diff = difference.mean()
    ax3.text(0.02, 0.98, f'Mean Δ: {mean_diff:.2f}',
             transform=ax3.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             fontsize=10)

    # Plot 4: Log-ratio
    ax4 = axes[1, 1]
    max_abs_log = max(abs(difference_log.min()), abs(difference_log.max()))
    sns.heatmap(difference_log, cmap='RdBu_r', ax=ax4, center=0,
                vmin=-max_abs_log, vmax=max_abs_log,
                cbar_kws={'label': 'log₁₀(Genomic / Protein)'})
    ax4.set_title('Log-Ratio\nlog₁₀(Genomic / Protein Only)',
                  fontsize=14, fontweight='bold')
    ax4.set_xlabel('Position j (mutated)', fontsize=11)
    ax4.set_ylabel('Position i (measured)', fontsize=11)

    # Add statistics
    mean_log_ratio = difference_log.mean()
    ax4.text(0.02, 0.98, f'Mean log-ratio: {mean_log_ratio:.3f}',
             transform=ax4.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             fontsize=10)

    # Main title
    fig.suptitle(f'Protein Co-evolution: Context Effect Comparison\n{seq_id}',
                 fontsize=16, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")
    print(f"    Protein-only     - Mean: {protein_only.mean():.4f}, Max: {protein_only.max():.4f}")
    print(f"    With context     - Mean: {protein_from_genomic.mean():.4f}, Max: {protein_from_genomic.max():.4f}")
    print(f"    Absolute diff    - Mean: {difference.mean():.4f}")
    print(f"    Relative change  - Mean: {(difference.mean() / protein_only.mean() * 100):.2f}%")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Compare protein Jacobians with and without genomic context"
    )
    parser.add_argument(
        "--protein-jacobian-dir",
        type=str,
        default="results/jacobian_analysis",
        help="Directory containing protein-only Jacobians"
    )
    parser.add_argument(
        "--genomic-jacobian-dir",
        type=str,
        default="results/jacobian_genomic",
        help="Directory containing genomic Jacobians"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_context_comparison",
        help="Output directory for comparison plots"
    )
    parser.add_argument(
        "--dna-flank-length",
        type=int,
        default=500,
        help="Length of DNA flanking regions in genomic analysis"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("Protein Jacobian Context Effect Analysis")
    print("Comparing protein co-evolution with and without genomic context")
    print("=" * 80)
    print(f"Protein-only directory: {args.protein_jacobian_dir}")
    print(f"Genomic directory: {args.genomic_jacobian_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"DNA flank length: {args.dna_flank_length} nt")
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    protein_dir = Path(args.protein_jacobian_dir)
    genomic_dir = Path(args.genomic_jacobian_dir)

    # Find all gLM2 protein Jacobian files
    protein_files = sorted(protein_dir.glob("*_glm2_jacobian.npy"))

    if not protein_files:
        print(f"ERROR: No protein Jacobian files found in {protein_dir}")
        return

    print(f"Found {len(protein_files)} protein sequences to compare")
    print()

    # Process each protein
    for protein_file in protein_files:
        # Extract sequence ID from filename
        # Format: GCF_024496245.1_1_2691_glm2_jacobian.npy
        seq_id = protein_file.stem.replace("_glm2_jacobian", "")

        print(f"Processing: {seq_id}")

        # Load protein-only Jacobian
        protein_jacobian = load_jacobian(protein_file)
        L_protein = protein_jacobian.shape[0]
        print(f"  Protein length: {L_protein} aa")

        # Find corresponding genomic Jacobian
        genomic_file = genomic_dir / f"{seq_id}_glm2_genomic_jacobian.npy"

        if not genomic_file.exists():
            print(f"  WARNING: Genomic Jacobian not found: {genomic_file}")
            print()
            continue

        # Load genomic Jacobian
        genomic_jacobian = load_jacobian(genomic_file)
        print(f"  Genomic total length: {genomic_jacobian.shape[0]} (DNA5 + Protein + DNA3)")

        # Extract protein sub-matrix from genomic Jacobian
        protein_from_genomic = extract_protein_submatrix(
            genomic_jacobian, args.dna_flank_length, L_protein
        )

        # Verify dimensions match
        if protein_jacobian.shape != protein_from_genomic.shape:
            print(f"  ERROR: Shape mismatch!")
            print(f"    Protein-only: {protein_jacobian.shape}")
            print(f"    From genomic: {protein_from_genomic.shape}")
            print()
            continue

        # Create comparison plot
        output_path = output_dir / f"{seq_id}_context_comparison.png"
        plot_comparison(protein_jacobian, protein_from_genomic, seq_id, output_path)

    print("=" * 80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
