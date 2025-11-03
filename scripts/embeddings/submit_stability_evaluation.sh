#!/bin/bash
#SBATCH --job-name=stability_eval
#SBATCH --output=logs/stability_%A_%a.out
#SBATCH --error=logs/stability_%A_%a.err
#SBATCH --array=0-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=campus-new

# Stability evaluation for key clustering configurations
# Tests if gene pairs co-cluster consistently across multiple subsamples

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/evaluate_clustering_stability.py"
PCA_FILE="${PROJECT_DIR}/results/umap/pca_cache.npz"

cd "${PROJECT_DIR}"

mkdir -p logs
mkdir -p results/clustering/stability

# Define configurations to test
# Format: resolution,n_neighbors,cog_only
CONFIGS=(
    "1500,15,true"    # Best performer - COG-only
    "1500,15,false"   # Best performer - all genes
    "1000,15,true"    # Second best - COG-only
    "750,15,true"     # Third best - COG-only
)

CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"
IFS=',' read -r RES NN COG <<< "$CONFIG"

# Set output file
if [ "$COG" == "true" ]; then
    OUTPUT_FILE="${PROJECT_DIR}/results/clustering/stability/stability_res${RES}_nn${NN}_cogonly.npz"
    COG_FLAG="--cog-only"
else
    OUTPUT_FILE="${PROJECT_DIR}/results/clustering/stability/stability_res${RES}_nn${NN}_all.npz"
    COG_FLAG=""
fi

echo "=========================================="
echo "Stability Evaluation"
echo "Task: ${SLURM_ARRAY_TASK_ID}/4"
echo "=========================================="
echo "Configuration:"
echo "  Resolution: ${RES}"
echo "  N neighbors: ${NN}"
echo "  COG-only: ${COG}"
echo ""
echo "Output: ${OUTPUT_FILE}"
echo "Start time: $(date)"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Run stability evaluation
/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --pca "${PCA_FILE}" \
    --n-subsamples 10 \
    --subsample-size 100000 \
    --resolution ${RES} \
    --n-neighbors ${NN} \
    ${COG_FLAG} \
    --min-cooccur 3 \
    --output "${OUTPUT_FILE}"

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: res=${RES}, nn=${NN}, cog=${COG}"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
