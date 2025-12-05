#!/bin/bash
#SBATCH --job-name=umap_cogonly
#SBATCH --output=logs/umap_cogonly_%A_%a.out
#SBATCH --error=logs/umap_cogonly_%A_%a.err
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=campus-new

# Compute UMAP for COG-annotated genes only
# For n_neighbors: 15, 25, 50, 100, 200
# Total: 5 jobs

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/compute_umap_cogonly.py"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap"
PCA_CACHE="${RESULTS_DIR}/pca_cache.npz"

cd "${PROJECT_DIR}"

# Create output directory
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# N neighbors values
N_NEIGHBORS=(15 25 50 100 200)
N_NEIGHBOR=${N_NEIGHBORS[$SLURM_ARRAY_TASK_ID]}

# Output file
OUTPUT_FILE="${RESULTS_DIR}/umap_n${N_NEIGHBOR}_subsample1000000_cogonly.npz"

echo "=========================================="
echo "UMAP for COG-Annotated Genes"
echo "Task: ${SLURM_ARRAY_TASK_ID}/5"
echo "=========================================="
echo "N neighbors: ${N_NEIGHBOR}"
echo "Output: $(basename ${OUTPUT_FILE})"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Run UMAP
echo "Starting UMAP computation..."
echo "Start time: $(date)"
echo ""

/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --pca-cache "${PCA_CACHE}" \
    --output "${OUTPUT_FILE}" \
    --n-neighbors "${N_NEIGHBOR}" \
    --min-dist 0.1

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: UMAP n=${N_NEIGHBOR} (COG-only)"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
