#!/bin/bash
#SBATCH --job-name=umap_rapids
#SBATCH --output=logs/umap_%A_%a.out
#SBATCH --error=logs/umap_%A_%a.err
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=campus-new
#SBATCH --gres=gpu:1

# UMAP with RAPIDS Container
# Uses Singularity container with cuML pre-installed

# Configuration
N_NEIGHBORS_VALUES=(15 25 50 100 200)
N_PCS=50
MIN_DIST=0.1
SUBSAMPLE=100000

# Paths
PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
EMBEDDINGS_DIR="${PROJECT_DIR}/data/esm_embeddings"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/compute_umap_array.py"

# RAPIDS container (will be pulled on first run)
CONTAINER_DIR="${PROJECT_DIR}/containers"
RAPIDS_CONTAINER="${CONTAINER_DIR}/rapids_24.10.sif"

# Create directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${CONTAINER_DIR}"

# Pull RAPIDS container if it doesn't exist
if [ ! -f "${RAPIDS_CONTAINER}" ]; then
    echo "Pulling RAPIDS container (one-time setup, ~5GB download)..."
    singularity pull "${RAPIDS_CONTAINER}" docker://rapidsai/rapidsai:24.10-cuda12.0-runtime-ubuntu22.04-py3.10
fi

# Get parameters
N_NEIGHBORS=${N_NEIGHBORS_VALUES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "RAPIDS Container UMAP"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "n_neighbors: ${N_NEIGHBORS}"
echo "=========================================="

# Output file
if [ -n "$SUBSAMPLE" ]; then
    OUTPUT_FILE="${RESULTS_DIR}/umap_n${N_NEIGHBORS}_subsample${SUBSAMPLE}.npz"
else
    OUTPUT_FILE="${RESULTS_DIR}/umap_n${N_NEIGHBORS}_full.npz"
fi

# Build command
CMD="python ${SCRIPT}"
CMD="${CMD} --embeddings-dir ${EMBEDDINGS_DIR}"
CMD="${CMD} --output ${OUTPUT_FILE}"
CMD="${CMD} --n-neighbors ${N_NEIGHBORS}"
CMD="${CMD} --n-pcs ${N_PCS}"
CMD="${CMD} --min-dist ${MIN_DIST}"
CMD="${CMD} --use-gpu"  # Container has cuML

if [ -n "$SUBSAMPLE" ]; then
    CMD="${CMD} --subsample ${SUBSAMPLE}"
fi

# Run in Singularity container with GPU
echo "Running: ${CMD}"
singularity exec --nv "${RAPIDS_CONTAINER}" ${CMD}

EXIT_CODE=$?

echo "=========================================="
echo "Exit code: ${EXIT_CODE}"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Successfully computed UMAP"
else
    echo "✗ Error computing UMAP"
fi
echo "=========================================="

exit ${EXIT_CODE}
