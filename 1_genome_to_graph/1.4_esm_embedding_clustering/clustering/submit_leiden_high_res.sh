#!/bin/bash
#SBATCH --job-name=leiden_highres
#SBATCH --output=logs/leiden_highres_%A_%a.out
#SBATCH --error=logs/leiden_highres_%A_%a.err
#SBATCH --array=0-19
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --partition=campus-new

# Higher resolution Leiden clustering
# - Resolutions: 750, 1000
# - N neighbors: 15, 25, 50, 100, 200
# - Versions: all genes, COG-annotated only
# Total: 2 × 5 × 2 = 20 jobs

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/cluster_leiden_comprehensive.py"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering"
PCA_CACHE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"

# Create output directory
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Parameter grid
RESOLUTIONS=(750 1000)
N_NEIGHBORS=(15 25 50 100 200)
VERSIONS=("all" "cog_only")

# Total configs: 2 × 5 × 2 = 20
N_RES=${#RESOLUTIONS[@]}
N_NEIGHBORS_VALS=${#N_NEIGHBORS[@]}
N_VERSIONS=${#VERSIONS[@]}

# Calculate which config this task is
RES_IDX=$((SLURM_ARRAY_TASK_ID / (N_NEIGHBORS_VALS * N_VERSIONS)))
NEIGHBORS_IDX=$(((SLURM_ARRAY_TASK_ID / N_VERSIONS) % N_NEIGHBORS_VALS))
VERSION_IDX=$((SLURM_ARRAY_TASK_ID % N_VERSIONS))

RESOLUTION=${RESOLUTIONS[$RES_IDX]}
N_NEIGHBOR=${N_NEIGHBORS[$NEIGHBORS_IDX]}
VERSION=${VERSIONS[$VERSION_IDX]}

# Build output filename
if [ "${VERSION}" = "cog_only" ]; then
    SUFFIX="leiden_res${RESOLUTION}_nn${N_NEIGHBOR}_cogonly"
    COG_FLAG="--cog-only"
else
    SUFFIX="leiden_res${RESOLUTION}_nn${N_NEIGHBOR}_all"
    COG_FLAG=""
fi

OUTPUT_FILE="${RESULTS_DIR}/clusters_${SUFFIX}.npz"

echo "=========================================="
echo "High Resolution Leiden Clustering"
echo "Task: ${SLURM_ARRAY_TASK_ID}/20"
echo "=========================================="
echo "Resolution: ${RESOLUTION}"
echo "N neighbors: ${N_NEIGHBOR}"
echo "Version: ${VERSION}"
echo "Output: ${SUFFIX}"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Run clustering with evaluation
echo "Starting clustering..."
echo "Start time: $(date)"
echo ""

/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --pca-cache "${PCA_CACHE}" \
    --output "${OUTPUT_FILE}" \
    --resolution "${RESOLUTION}" \
    --n-neighbors "${N_NEIGHBOR}" \
    ${COG_FLAG} \
    --evaluate

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: ${SUFFIX}"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
