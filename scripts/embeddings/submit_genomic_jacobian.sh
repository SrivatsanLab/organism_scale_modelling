#!/bin/bash
#SBATCH --job-name=genomic_jac
#SBATCH --output=logs/genomic_jacobian_%j.out
#SBATCH --error=logs/genomic_jacobian_%j.err
#SBATCH --partition=chorus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00

# Print job info
echo "========================================"
echo "Genomic Jacobian Analysis"
echo "gLM2 with DNA+Protein Context"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Check GPU
echo "GPU Info:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo ""

# Load environment
source ~/.bashrc
micromamba activate glm2_env

# Run analysis
echo "========================================"
echo "Computing Jacobians"
echo "========================================"
echo ""

python scripts/embeddings/calculate_genomic_jacobian.py \
    --protein-dir data/protein_samples \
    --output-dir results/jacobian_genomic \
    --num-proteins 5 \
    --device cuda

echo ""
echo "========================================"
echo "Job Complete!"
echo "========================================"
