#!/bin/bash
#SBATCH --job-name=cluster_hires
#SBATCH --output=logs/cluster_hires_%A_%a.out
#SBATCH --error=logs/cluster_hires_%A_%a.err
#SBATCH --array=0-15
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --partition=campus-new

# High-Resolution Clustering Parameter Sweep
# Focus: Finding THOUSANDS of clusters using high-resolution methods
# Methods: Leiden (high resolution), K-means (high k), MiniBatch K-means

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/cluster_umap_array.py"
UMAP_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample1000000.npz"
PCA_CACHE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Define all clustering configurations
# Format: METHOD|PARAMS_JSON|USE_PCA|OUTPUT_SUFFIX
# FOCUS: High resolution clustering on 50D PCA space

CONFIGS=(
    # Leiden on PCA - VERY HIGH resolutions for thousands of clusters
    'leiden|{"resolution":5.0,"n_neighbors":30}|true|leiden_pca_res5.0'
    'leiden|{"resolution":10.0,"n_neighbors":30}|true|leiden_pca_res10.0'
    'leiden|{"resolution":20.0,"n_neighbors":30}|true|leiden_pca_res20.0'
    'leiden|{"resolution":50.0,"n_neighbors":30}|true|leiden_pca_res50.0'
    'leiden|{"resolution":100.0,"n_neighbors":30}|true|leiden_pca_res100.0'
    'leiden|{"resolution":200.0,"n_neighbors":30}|true|leiden_pca_res200.0'

    # K-means on PCA - high k values
    'kmeans|{"n_clusters":500}|true|kmeans_pca_k500'
    'kmeans|{"n_clusters":1000}|true|kmeans_pca_k1000'
    'kmeans|{"n_clusters":2000}|true|kmeans_pca_k2000'
    'kmeans|{"n_clusters":5000}|true|kmeans_pca_k5000'

    # MiniBatch K-means on PCA - very high k (faster than regular k-means)
    'minibatch_kmeans|{"n_clusters":1000}|true|minibatch_pca_k1000'
    'minibatch_kmeans|{"n_clusters":2000}|true|minibatch_pca_k2000'
    'minibatch_kmeans|{"n_clusters":5000}|true|minibatch_pca_k5000'
    'minibatch_kmeans|{"n_clusters":10000}|true|minibatch_pca_k10000'

    # HDBSCAN on PCA - very small min_cluster_size for more clusters
    'hdbscan|{"min_cluster_size":10}|true|hdbscan_pca_minclust10'
    'hdbscan|{"min_cluster_size":25}|true|hdbscan_pca_minclust25'
)

# Get configuration for this task
CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

# Parse configuration
IFS='|' read -r METHOD PARAMS USE_PCA OUTPUT_SUFFIX <<< "$CONFIG"

echo "=========================================="
echo "High-Resolution Clustering Sweep"
echo "Task: ${SLURM_ARRAY_TASK_ID} / $((${#CONFIGS[@]} - 1))"
echo "=========================================="
echo "Method: ${METHOD}"
echo "Parameters: ${PARAMS}"
echo "Use PCA: ${USE_PCA}"
echo "Output: ${RESULTS_DIR}/clusters_${OUTPUT_SUFFIX}.npz"
echo "=========================================="
echo ""

# Build command args
OUTPUT_PATH="${RESULTS_DIR}/clusters_${OUTPUT_SUFFIX}.npz"

echo "Running:"
echo "/home/dmullane/micromamba/envs/esm3_env/bin/python ${SCRIPT}"
echo "  --umap ${UMAP_FILE}"
echo "  --output ${OUTPUT_PATH}"
echo "  --method ${METHOD}"
echo "  --params ${PARAMS}"

if [ "${USE_PCA}" = "true" ]; then
    echo "  --use-pca"
    echo "  --pca-cache ${PCA_CACHE}"
fi
echo ""

# Run clustering
if [ "${USE_PCA}" = "true" ]; then
    /home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
        --umap "${UMAP_FILE}" \
        --output "${OUTPUT_PATH}" \
        --method "${METHOD}" \
        --params "${PARAMS}" \
        --use-pca \
        --pca-cache "${PCA_CACHE}"
else
    /home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
        --umap "${UMAP_FILE}" \
        --output "${OUTPUT_PATH}" \
        --method "${METHOD}" \
        --params "${PARAMS}"
fi

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Success"
else
    echo "✗ Error"
fi
echo "Exit code: ${EXIT_CODE}"
echo "=========================================="

exit ${EXIT_CODE}
