#!/bin/bash
#SBATCH --job-name=5utr_jacobian
#SBATCH --output=logs/5utr_jacobian_%j.out
#SBATCH --error=logs/5utr_jacobian_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=chorus
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Print job information
echo "========================================"
echo "5' UTR Jacobian Analysis"
echo "gLM2 Model"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Check GPU
echo "GPU Info:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo ""

# Activate micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate glm2_env

echo "========================================"
echo "Computing 5' UTR Jacobians"
echo "========================================"
echo ""

# Run the analysis
python scripts/embeddings/calculate_5utr_jacobian.py \
    --protein-dir data/protein_samples \
    --output-dir results/jacobian_5utr \
    --num-proteins 5 \
    --device cuda

echo ""
echo "========================================"
echo "Job Complete!"
echo "========================================"
