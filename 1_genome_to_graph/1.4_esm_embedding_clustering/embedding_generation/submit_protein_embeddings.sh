#!/bin/bash
#SBATCH --job-name=protein_embeddings
#SBATCH --partition=chorus
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/protein_embeddings_%j.out
#SBATCH --error=logs/protein_embeddings_%j.err

echo "========================================"
echo "Compute Protein-Only Embeddings"
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

# Check PyTorch CUDA
echo "PyTorch CUDA Status:"
python -c "
import torch
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA version: {torch.version.cuda}')
    print(f'GPU: {torch.cuda.get_device_name(0)}')
"
echo ""

# Run embedding computation
echo "========================================"
echo "Computing Embeddings"
echo "========================================"
echo ""

python scripts/embeddings/compute_protein_embeddings.py \
    --protein-dir data/protein_samples \
    --output-dir data/embeddings/protein_only \
    --batch-size 32 \
    --device cuda

echo ""
echo "========================================"
echo "Job Complete!"
echo "========================================"
