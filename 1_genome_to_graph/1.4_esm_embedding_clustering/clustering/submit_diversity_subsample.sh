#!/bin/bash
#SBATCH --job-name=div_sample
#SBATCH --output=logs/diversity_subsample_%A.out
#SBATCH --error=logs/diversity_subsample_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --partition=largenode

# Diversity-based subsampling from full 29M gene dataset

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/diversity_subsample.py"
INPUT_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache_full.npz"  # 29M genes with PCA
OUTPUT_DIR="${PROJECT_DIR}/results/subsampling"

cd "${PROJECT_DIR}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Parameters
N_SAMPLES=1000000  # Target 1M diverse samples
N_CLUSTERS=10000   # Use 10k clusters for diversity
METHOD="hybrid"    # Cluster + farthest-point (best balance)

OUTPUT_FILE="${OUTPUT_DIR}/diverse_subsample_${N_SAMPLES}_${METHOD}.npz"

echo "=========================================="
echo "Diversity-Based Subsampling"
echo "=========================================="
echo "Input: Full 29M gene dataset"
echo "Method: ${METHOD}"
echo "Clusters: ${N_CLUSTERS}"
echo "Target samples: ${N_SAMPLES}"
echo "Output: $(basename ${OUTPUT_FILE})"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Run subsampling
echo "Starting diversity subsampling..."
echo "Start time: $(date)"
echo ""

/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --input "${INPUT_FILE}" \
    --output "${OUTPUT_FILE}" \
    --n-samples "${N_SAMPLES}" \
    --method "${METHOD}" \
    --n-clusters "${N_CLUSTERS}" \
    --use-pca

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: Diversity subsampling complete"
    echo "Output: ${OUTPUT_FILE}"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
