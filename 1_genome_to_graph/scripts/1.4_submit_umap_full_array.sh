#!/bin/bash
#SBATCH --job-name=umap_full
#SBATCH --output=logs/umap_full_%A_%a.out
#SBATCH --error=logs/umap_full_%A_%a.err
#SBATCH --array=0-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --partition=campus-new

# UMAP Array Job for Full Dataset (~29M genes)
# Each job computes UMAP with a different n_neighbors on ALL genes

# Configuration
N_NEIGHBORS_VALUES=(15 25 50 100 200)  # Different n_neighbors to test
N_PCS=50                                 # Number of PCA components
MIN_DIST=0.1                            # UMAP min_dist parameter
SUBSAMPLE=""                            # Empty = use all genes (~29M)
USE_GPU=""                              # CPU mode (GPU likely won't have enough memory)

# Paths
PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
EMBEDDINGS_DIR="${PROJECT_DIR}/data/esm_embeddings"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap"
PCA_CACHE="${RESULTS_DIR}/pca_cache_full.npz"  # Different cache for full dataset
SCRIPT="${PROJECT_DIR}/scripts/embeddings/compute_umap_array.py"

# Create output directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Get n_neighbors value for this array task
N_NEIGHBORS=${N_NEIGHBORS_VALUES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "UMAP FULL DATASET (~29M genes)"
echo "=========================================="
echo "SLURM Job Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "n_neighbors: ${N_NEIGHBORS}"
echo "n_pcs: ${N_PCS}"
echo "Subsample: ALL GENES (no subsampling)"
echo "Memory: 256G"
echo "CPUs: 16"
echo "Time limit: 12 hours"
echo "=========================================="

# Output file for this n_neighbors value
OUTPUT_FILE="${RESULTS_DIR}/umap_n${N_NEIGHBORS}_full.npz"

# Build command
CMD="/home/dmullane/micromamba/envs/esm3_env/bin/python ${SCRIPT}"
CMD="${CMD} --embeddings-dir ${EMBEDDINGS_DIR}"
CMD="${CMD} --output ${OUTPUT_FILE}"
CMD="${CMD} --n-neighbors ${N_NEIGHBORS}"
CMD="${CMD} --n-pcs ${N_PCS}"
CMD="${CMD} --min-dist ${MIN_DIST}"

# Don't add subsample flag - process all genes

# Use cached PCA if it exists
# First job (task 0) creates the cache, others use it
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    # First job: compute and save PCA
    CMD="${CMD} --save-pca ${PCA_CACHE}"
    echo "First job - will compute and save PCA to cache"
    echo "This will take significant time for ~29M genes"
else
    # Other jobs: wait for PCA cache to be created, then use it
    echo "Waiting for PCA cache to be created by task 0..."

    # Wait up to 4 hours for PCA cache (PCA on 29M genes will take a while)
    for i in {1..1440}; do
        if [ -f "${PCA_CACHE}" ]; then
            echo "PCA cache found! Using cached PCA."
            CMD="${CMD} --load-cached-pca ${PCA_CACHE}"
            break
        fi

        if [ $i -eq 1440 ]; then
            echo "ERROR: PCA cache not found after 4 hours"
            exit 1
        fi

        sleep 10
    done
fi

echo ""
echo "Command:"
echo "${CMD}"
echo ""
echo "Starting at $(date)"
echo ""

# Run the job
${CMD}

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Completed at $(date)"
echo "Job completed with exit code: ${EXIT_CODE}"
echo "Output: ${OUTPUT_FILE}"
echo "=========================================="

exit ${EXIT_CODE}
