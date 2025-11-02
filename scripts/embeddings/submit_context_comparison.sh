#!/bin/bash
#SBATCH --job-name=context_comparison
#SBATCH --output=logs/context_comparison_%j.out
#SBATCH --error=logs/context_comparison_%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Print job information
echo "========================================"
echo "Protein Context Effect Comparison"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Activate micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate glm2_env

echo "========================================"
echo "Comparing Protein Jacobians"
echo "========================================"
echo ""

# Run the analysis
python scripts/embeddings/compare_protein_context_effect.py \
    --protein-jacobian-dir results/jacobian_analysis \
    --genomic-jacobian-dir results/jacobian_genomic \
    --output-dir results/jacobian_context_comparison \
    --dna-flank-length 500

echo ""
echo "========================================"
echo "Job Complete!"
echo "========================================"
