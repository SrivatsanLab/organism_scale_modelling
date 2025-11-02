#!/bin/bash
#SBATCH --job-name=jacobian_analysis
#SBATCH --partition=chorus
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/jacobian_analysis_%j.out
#SBATCH --error=logs/jacobian_analysis_%j.err

echo "========================================"
echo "Categorical Jacobian Analysis"
echo "ESM-C vs gLM2"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load CUDA module
module load CUDA/11.7.0

# Load environment
source ~/.bashrc
micromamba activate glm2_env

# Check GPU
echo "GPU Info:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
echo ""

# Run Jacobian analysis
echo "========================================"
echo "Computing Jacobians"
echo "========================================"
echo ""

python scripts/embeddings/calculate_protein_jacobian.py \
    --protein-dir data/protein_samples \
    --output-dir results/jacobian_analysis \
    --num-proteins 5 \
    --device cuda

echo ""
echo "========================================"
echo "Job Complete!"
echo "========================================"
