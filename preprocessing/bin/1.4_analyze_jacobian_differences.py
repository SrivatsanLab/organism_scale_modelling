#!/usr/bin/env python
"""
Analyze the differences between DNA-encoded and protein-encoded Jacobians.

After computing Jacobians for both DNA (CDS) and protein representations,
this script analyzes:
1. The difference matrix (DNA - Protein)
2. Positions where DNA encoding shows stronger coupling
3. Positions where protein encoding shows stronger coupling
4. Correlation between the two Jacobians

Args:
    Loads pre-computed Jacobian matrices and performs differential analysis.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze differences between DNA and Protein Jacobians"
    )
    parser.add_argument(
        "--jacobian-dir",
        type=str,
        default="results/jacobian_dna_vs_protein",
        help="Directory containing Jacobian .npy files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/jacobian_differences",
        help="Output directory for analysis results",
    )
    return parser.parse_args()


def plot_difference_heatmap(
    difference: np.ndarray,
    output_path: Path,
    title: str,
    seq_id: str
):
    """Plot the difference between DNA and protein Jacobians."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Use diverging colormap centered at 0
    # Positive values = DNA encoding shows stronger coupling
    # Negative values = Protein encoding shows stronger coupling
    max_abs = np.abs(difference).max()

    sns.heatmap(
        difference,
        cmap='RdBu_r',
        center=0,
        vmin=-max_abs,
        vmax=max_abs,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'Coupling Difference\n(DNA - Protein)'},
        ax=ax
    )

    ax.set_xlabel('Position j', fontsize=12)
    ax.set_ylabel('Position i', fontsize=12)
    ax.set_title(
        f'{title}\nDifference: DNA - Protein Jacobian\n' +
        f'Red = DNA stronger, Blue = Protein stronger',
        fontsize=14,
        fontweight='bold'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def plot_log_difference_heatmap(
    dna_jacobian: np.ndarray,
    protein_jacobian: np.ndarray,
    output_path: Path,
    title: str,
    seq_id: str
):
    """Plot the log-ratio of DNA/Protein Jacobians."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Compute log ratio: log10(DNA / Protein)
    # Add small epsilon to avoid division by zero
    log_ratio = np.log10((dna_jacobian + 1e-10) / (protein_jacobian + 1e-10))

    # Clip extreme values for visualization
    vmax = np.percentile(np.abs(log_ratio), 95)

    sns.heatmap(
        log_ratio,
        cmap='RdBu_r',
        center=0,
        vmin=-vmax,
        vmax=vmax,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'log10(DNA / Protein)'},
        ax=ax
    )

    ax.set_xlabel('Position j', fontsize=12)
    ax.set_ylabel('Position i', fontsize=12)
    ax.set_title(
        f'{title}\nLog-ratio: DNA / Protein Jacobian\n' +
        f'Red = DNA stronger, Blue = Protein stronger',
        fontsize=14,
        fontweight='bold'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def plot_correlation_scatter(
    dna_jacobian: np.ndarray,
    protein_jacobian: np.ndarray,
    output_path: Path,
    title: str
):
    """Plot scatter plot of DNA vs Protein Jacobian values."""
    fig, ax = plt.subplots(figsize=(10, 10))

    # Flatten matrices and remove diagonal
    mask = ~np.eye(dna_jacobian.shape[0], dtype=bool)
    dna_flat = dna_jacobian[mask]
    protein_flat = protein_jacobian[mask]

    # Compute correlation
    pearson_r, pearson_p = stats.pearsonr(dna_flat, protein_flat)
    spearman_r, spearman_p = stats.spearmanr(dna_flat, protein_flat)

    # Hexbin plot for large datasets
    hexbin = ax.hexbin(
        protein_flat,
        dna_flat,
        gridsize=50,
        cmap='viridis',
        mincnt=1,
        bins='log'
    )

    # Add diagonal line
    max_val = max(dna_flat.max(), protein_flat.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, alpha=0.5, label='y=x')

    ax.set_xlabel('Protein Encoding Coupling', fontsize=12)
    ax.set_ylabel('DNA Encoding Coupling', fontsize=12)
    ax.set_title(
        f'{title}\nPearson r={pearson_r:.3f} (p={pearson_p:.2e})\n' +
        f'Spearman ρ={spearman_r:.3f} (p={spearman_p:.2e})',
        fontsize=14,
        fontweight='bold'
    )

    plt.colorbar(hexbin, ax=ax, label='log10(Count)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")

    return pearson_r, spearman_r


def analyze_top_differences(
    dna_jacobian: np.ndarray,
    protein_jacobian: np.ndarray,
    difference: np.ndarray,
    seq_id: str,
    top_n: int = 20
):
    """Analyze top positions with largest differences."""
    # Flatten and get positions
    L = dna_jacobian.shape[0]

    # Remove diagonal
    mask = ~np.eye(L, dtype=bool)
    diff_flat = difference[mask]

    # Get indices of top differences (both positive and negative)
    top_dna_indices = np.argsort(diff_flat)[-top_n:][::-1]  # Largest positive (DNA stronger)
    top_protein_indices = np.argsort(diff_flat)[:top_n]     # Largest negative (Protein stronger)

    # Convert flat indices back to 2D
    positions = np.where(mask)

    results = {
        'dna_stronger': [],
        'protein_stronger': []
    }

    # DNA stronger positions
    for idx in top_dna_indices:
        i, j = positions[0][idx], positions[1][idx]
        results['dna_stronger'].append({
            'position_i': int(i),
            'position_j': int(j),
            'dna_coupling': float(dna_jacobian[i, j]),
            'protein_coupling': float(protein_jacobian[i, j]),
            'difference': float(difference[i, j]),
            'fold_change': float(dna_jacobian[i, j] / (protein_jacobian[i, j] + 1e-10))
        })

    # Protein stronger positions
    for idx in top_protein_indices:
        i, j = positions[0][idx], positions[1][idx]
        results['protein_stronger'].append({
            'position_i': int(i),
            'position_j': int(j),
            'dna_coupling': float(dna_jacobian[i, j]),
            'protein_coupling': float(protein_jacobian[i, j]),
            'difference': float(difference[i, j]),
            'fold_change': float(protein_jacobian[i, j] / (dna_jacobian[i, j] + 1e-10))
        })

    return results


def compute_summary_statistics(
    dna_jacobian: np.ndarray,
    protein_jacobian: np.ndarray,
    difference: np.ndarray
):
    """Compute summary statistics for the Jacobians."""
    # Remove diagonal
    mask = ~np.eye(dna_jacobian.shape[0], dtype=bool)

    dna_off_diag = dna_jacobian[mask]
    protein_off_diag = protein_jacobian[mask]
    diff_off_diag = difference[mask]

    stats_dict = {
        'dna_mean': float(np.mean(dna_off_diag)),
        'dna_std': float(np.std(dna_off_diag)),
        'dna_median': float(np.median(dna_off_diag)),
        'protein_mean': float(np.mean(protein_off_diag)),
        'protein_std': float(np.std(protein_off_diag)),
        'protein_median': float(np.median(protein_off_diag)),
        'diff_mean': float(np.mean(diff_off_diag)),
        'diff_std': float(np.std(diff_off_diag)),
        'diff_median': float(np.median(diff_off_diag)),
        'fraction_dna_stronger': float(np.sum(diff_off_diag > 0) / len(diff_off_diag)),
        'fraction_protein_stronger': float(np.sum(diff_off_diag < 0) / len(diff_off_diag)),
    }

    return stats_dict


def main():
    args = parse_args()

    print("="*80)
    print("DNA vs Protein Jacobian Difference Analysis")
    print("="*80)
    print(f"Input directory: {args.jacobian_dir}")
    print(f"Output directory: {args.output_dir}")
    print()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Find all protein/DNA Jacobian pairs
    jacobian_dir = Path(args.jacobian_dir)
    protein_files = sorted(jacobian_dir.glob("*_protein_jacobian.npy"))

    all_stats = []

    for protein_file in protein_files:
        # Get corresponding DNA file
        base_name = protein_file.stem.replace("_protein_jacobian", "")
        dna_file = jacobian_dir / f"{base_name}_dna_jacobian.npy"

        if not dna_file.exists():
            print(f"Warning: No DNA Jacobian found for {base_name}")
            continue

        seq_id = base_name.replace("_", "|", 1)  # Convert back to original ID format
        print(f"Analyzing: {seq_id}")

        # Load Jacobians
        protein_jacobian = np.load(protein_file)
        dna_jacobian = np.load(dna_file)

        print(f"  Protein Jacobian shape: {protein_jacobian.shape}")
        print(f"  DNA Jacobian shape: {dna_jacobian.shape}")

        # Compute difference
        difference = dna_jacobian - protein_jacobian

        # Plot difference heatmap
        plot_difference_heatmap(
            difference,
            output_dir / f"{base_name}_difference.png",
            seq_id,
            base_name
        )

        # Plot log-ratio heatmap
        plot_log_difference_heatmap(
            dna_jacobian,
            protein_jacobian,
            output_dir / f"{base_name}_log_ratio.png",
            seq_id,
            base_name
        )

        # Plot correlation
        pearson_r, spearman_r = plot_correlation_scatter(
            dna_jacobian,
            protein_jacobian,
            output_dir / f"{base_name}_correlation.png",
            seq_id
        )

        # Compute summary statistics
        stats = compute_summary_statistics(dna_jacobian, protein_jacobian, difference)
        stats['seq_id'] = seq_id
        stats['pearson_r'] = pearson_r
        stats['spearman_r'] = spearman_r
        all_stats.append(stats)

        print(f"  DNA mean coupling: {stats['dna_mean']:.4f}")
        print(f"  Protein mean coupling: {stats['protein_mean']:.4f}")
        print(f"  Difference mean: {stats['diff_mean']:.4f}")
        print(f"  Fraction DNA stronger: {stats['fraction_dna_stronger']:.2%}")
        print(f"  Correlation (Pearson): {pearson_r:.3f}")
        print(f"  Correlation (Spearman): {spearman_r:.3f}")

        # Analyze top differences
        top_diffs = analyze_top_differences(
            dna_jacobian,
            protein_jacobian,
            difference,
            seq_id,
            top_n=20
        )

        # Save top differences
        import json
        with open(output_dir / f"{base_name}_top_differences.json", 'w') as f:
            json.dump(top_diffs, f, indent=2)

        print(f"  Saved top differences to JSON")
        print()

    # Save summary statistics
    df_stats = pd.DataFrame(all_stats)
    df_stats.to_csv(output_dir / "summary_statistics.csv", index=False)
    print(f"Saved summary statistics to: {output_dir / 'summary_statistics.csv'}")

    # Create summary plot
    if len(all_stats) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # Plot 1: Mean coupling comparison
        ax = axes[0, 0]
        x = np.arange(len(all_stats))
        width = 0.35
        ax.bar(x - width/2, [s['dna_mean'] for s in all_stats], width, label='DNA', alpha=0.8)
        ax.bar(x + width/2, [s['protein_mean'] for s in all_stats], width, label='Protein', alpha=0.8)
        ax.set_xlabel('Protein')
        ax.set_ylabel('Mean Coupling')
        ax.set_title('Mean Coupling: DNA vs Protein')
        ax.set_xticks(x)
        ax.set_xticklabels([s['seq_id'].split('|')[1] for s in all_stats], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 2: Correlation comparison
        ax = axes[0, 1]
        ax.scatter([s['pearson_r'] for s in all_stats],
                  [s['spearman_r'] for s in all_stats], s=100, alpha=0.6)
        for i, s in enumerate(all_stats):
            ax.annotate(s['seq_id'].split('|')[1],
                       ([s['pearson_r'] for s in all_stats][i],
                        [s['spearman_r'] for s in all_stats][i]),
                       fontsize=8)
        ax.set_xlabel('Pearson r')
        ax.set_ylabel('Spearman ρ')
        ax.set_title('Correlation between DNA and Protein Jacobians')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)

        # Plot 3: Fraction DNA vs Protein stronger
        ax = axes[1, 0]
        dna_fractions = [s['fraction_dna_stronger'] for s in all_stats]
        protein_fractions = [s['fraction_protein_stronger'] for s in all_stats]
        ax.bar(x - width/2, dna_fractions, width, label='DNA stronger', alpha=0.8)
        ax.bar(x + width/2, protein_fractions, width, label='Protein stronger', alpha=0.8)
        ax.set_xlabel('Protein')
        ax.set_ylabel('Fraction of position pairs')
        ax.set_title('Fraction of positions with stronger coupling')
        ax.set_xticks(x)
        ax.set_xticklabels([s['seq_id'].split('|')[1] for s in all_stats], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 4: Difference distribution
        ax = axes[1, 1]
        for i, s in enumerate(all_stats):
            jacobian_file = jacobian_dir / f"{s['seq_id'].replace('|', '_', 1)}_protein_jacobian.npy"
            dna_file = jacobian_dir / f"{s['seq_id'].replace('|', '_', 1)}_dna_jacobian.npy"
            protein_jac = np.load(jacobian_file)
            dna_jac = np.load(dna_file)
            diff = dna_jac - protein_jac
            mask = ~np.eye(diff.shape[0], dtype=bool)
            ax.hist(diff[mask], bins=50, alpha=0.5, label=s['seq_id'].split('|')[1])
        ax.set_xlabel('Difference (DNA - Protein)')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of coupling differences')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', linewidth=2)

        plt.tight_layout()
        plt.savefig(output_dir / "summary_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved summary plot to: {output_dir / 'summary_analysis.png'}")

    print()
    print("="*80)
    print("✓ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
