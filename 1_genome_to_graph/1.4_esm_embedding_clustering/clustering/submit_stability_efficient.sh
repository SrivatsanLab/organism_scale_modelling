#!/bin/bash
#SBATCH --job-name=stability_eff
#SBATCH --output=logs/stability_eff_%A_%a.out
#SBATCH --error=logs/stability_eff_%A_%a.err
#SBATCH --array=0-49
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --partition=campus-new

# Comprehensive stability evaluation across parameter space
# Tests: 5 resolutions × 5 n_neighbors × 2 datasets = 50 configs

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/evaluate_clustering_stability_efficient.py"
PCA_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"

cd "${PROJECT_DIR}"

mkdir -p logs
mkdir -p results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability

# Define parameter grid
RESOLUTIONS=(300 500 750 1000 1500)
N_NEIGHBORS=(15 25 50 100 200)
DATASETS=("true" "false")  # COG-only, all genes

# Calculate indices for this task
TASK_ID=${SLURM_ARRAY_TASK_ID}

# Total: 5 res × 5 nn × 2 datasets = 50
# Layout: [res0_nn0_cog, res0_nn0_all, res0_nn1_cog, res0_nn1_all, ...]

# Unpack task ID
DATASET_IDX=$((TASK_ID % 2))
NN_IDX=$(((TASK_ID / 2) % 5))
RES_IDX=$(((TASK_ID / 10)))

RES=${RESOLUTIONS[$RES_IDX]}
NN=${N_NEIGHBORS[$NN_IDX]}
COG=${DATASETS[$DATASET_IDX]}

# Set output file
if [ "$COG" == "true" ]; then
    OUTPUT_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability/stability_eff_res${RES}_nn${NN}_cogonly.npz"
    COG_FLAG="--cog-only"
else
    OUTPUT_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/stability/stability_eff_res${RES}_nn${NN}_all.npz"
    COG_FLAG=""
fi

echo "=========================================="
echo "Efficient Stability Evaluation"
echo "Task: ${SLURM_ARRAY_TASK_ID}/50"
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
echo "Starting Python script..."

/home/dmullane/micromamba/envs/esm3_env/bin/python -u "${SCRIPT}" \
    --pca "${PCA_FILE}" \
    --n-subsamples 10 \
    --subsample-size 100000 \
    --resolution ${RES} \
    --n-neighbors ${NN} \
    ${COG_FLAG} \
    --output "${OUTPUT_FILE}"

EXIT_CODE=$?
echo "Python script completed with exit code: ${EXIT_CODE}"

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
