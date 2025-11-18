#!/bin/bash
#SBATCH --job-name=gpu_test
#SBATCH --output=logs/gpu_test_%j.out
#SBATCH --error=logs/gpu_test_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=30:00
#SBATCH --partition=campus-new
#SBATCH --gres=gpu:1

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
TEST_SCRIPT="${PROJECT_DIR}/scripts/embeddings/test_gpu.py"

mkdir -p "${PROJECT_DIR}/logs"

echo "=========================================="
echo "GPU Environment Test"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "GPU: ${SLURM_JOB_GPUS}"
echo "Node: ${SLURMD_NODENAME}"
echo "=========================================="
echo ""

# Check available modules
echo "Checking for CUDA modules..."
module avail 2>&1 | grep -i "CUDA/11" | head -5
echo ""

# Try loading CUDA module
echo "Loading CUDA module..."
module load CUDA/11.7.0 2>/dev/null || module load CUDA/11.8.0 2>/dev/null || echo "No CUDA module loaded"
echo ""

# Check nvidia-smi
echo "Running nvidia-smi..."
nvidia-smi
echo ""

# Test with existing Python environment
echo "Testing with esm3_env..."
/home/dmullane/micromamba/envs/esm3_env/bin/python "${TEST_SCRIPT}"
echo ""

# Try installing cuml via micromamba
echo "=========================================="
echo "Attempting to install cuML via micromamba"
echo "=========================================="

/home/dmullane/.local/bin/micromamba install -n esm3_env -c rapidsai -c conda-forge -c nvidia \
    cuml=23.10 python=3.10 cudatoolkit=11.8 -y

echo ""
echo "Testing again after cuML installation..."
/home/dmullane/micromamba/envs/esm3_env/bin/python "${TEST_SCRIPT}"

echo ""
echo "=========================================="
echo "GPU Test Complete"
echo "=========================================="
