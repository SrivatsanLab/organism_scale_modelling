#!/bin/bash
#SBATCH --job-name=plot_leiden
#SBATCH --output=logs/plot_leiden_%A_%a.out
#SBATCH --error=logs/plot_leiden_%A_%a.err
#SBATCH --array=0-644
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new

# Plot all Leiden clusterings × UMAPs
# 129 clusterings × 5 UMAPs = 645 jobs

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/plot_leiden_clustering.py"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering"
PLOTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/plots"

cd "${PROJECT_DIR}"

# Create output directories
mkdir -p "${PLOTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Get list of all Leiden clustering files (sorted)
readarray -t CLUSTERING_FILES < <(ls -1 "${RESULTS_DIR}"/clusters_leiden_*.npz | sort)
N_CLUSTERINGS=${#CLUSTERING_FILES[@]}

# UMAP n_neighbors values
N_NEIGHBORS=(15 25 50 100 200)
N_UMAPS=${#N_NEIGHBORS[@]}

# Calculate which clustering and UMAP this task is
CLUSTERING_IDX=$((SLURM_ARRAY_TASK_ID / N_UMAPS))
UMAP_IDX=$((SLURM_ARRAY_TASK_ID % N_UMAPS))

# Get the actual files/values
CLUSTERING_FILE="${CLUSTERING_FILES[$CLUSTERING_IDX]}"
UMAP_N="${N_NEIGHBORS[$UMAP_IDX]}"

# Extract clustering name from file
CLUSTERING_NAME=$(basename "${CLUSTERING_FILE}" .npz | sed 's/clusters_//')

# Output file path
OUTPUT_DIR="${PLOTS_DIR}/umap_n${UMAP_N}"
OUTPUT_FILE="${OUTPUT_DIR}/${CLUSTERING_NAME}.png"

echo "=========================================="
echo "Plot Leiden Clustering Array"
echo "Task: ${SLURM_ARRAY_TASK_ID}/645"
echo "=========================================="
echo "Clustering: ${CLUSTERING_NAME}"
echo "UMAP n_neighbors: ${UMAP_N}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Run plotting
echo "Starting plot generation..."
echo "Start time: $(date)"
echo ""

/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --clustering "${CLUSTERING_FILE}" \
    --umap-n "${UMAP_N}" \
    --output "${OUTPUT_FILE}"

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: ${CLUSTERING_NAME} × UMAP n=${UMAP_N}"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
