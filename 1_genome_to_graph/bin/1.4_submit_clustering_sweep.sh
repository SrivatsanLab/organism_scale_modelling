#!/bin/bash
#SBATCH --job-name=cluster_sweep
#SBATCH --output=logs/cluster_%A_%a.out
#SBATCH --error=logs/cluster_%A_%a.err
#SBATCH --array=0-25
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new

# Clustering parameter sweep
# Tests multiple methods and parameters

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/cluster_umap_array.py"
UMAP_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample1000000.npz"
PCA_CACHE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Activate environment
source /home/dmullane/.local/share/mamba/etc/profile.d/micromamba.sh
micromamba activate esm3_env

# Define all clustering configurations
# Format: METHOD|PARAMS_JSON|USE_PCA|OUTPUT_SUFFIX
# FOCUS: Clustering on 50D PCA space (91.8% variance from 1152D ESM embeddings)

CONFIGS=(
    # HDBSCAN on PCA - varying min_cluster_size (most important parameter)
    'hdbscan|{"min_cluster_size":50}|true|hdbscan_pca_minclust50'
    'hdbscan|{"min_cluster_size":100}|true|hdbscan_pca_minclust100'
    'hdbscan|{"min_cluster_size":200}|true|hdbscan_pca_minclust200'
    'hdbscan|{"min_cluster_size":500}|true|hdbscan_pca_minclust500'
    'hdbscan|{"min_cluster_size":1000}|true|hdbscan_pca_minclust1000'
    'hdbscan|{"min_cluster_size":2000}|true|hdbscan_pca_minclust2000'

    # HDBSCAN on PCA - varying cluster_selection_epsilon
    'hdbscan|{"min_cluster_size":500,"cluster_selection_epsilon":0.1}|true|hdbscan_pca_eps0.1'
    'hdbscan|{"min_cluster_size":500,"cluster_selection_epsilon":0.5}|true|hdbscan_pca_eps0.5'

    # K-means on PCA - varying k
    'kmeans|{"n_clusters":25}|true|kmeans_pca_k25'
    'kmeans|{"n_clusters":50}|true|kmeans_pca_k50'
    'kmeans|{"n_clusters":100}|true|kmeans_pca_k100'
    'kmeans|{"n_clusters":200}|true|kmeans_pca_k200'
    'kmeans|{"n_clusters":500}|true|kmeans_pca_k500'

    # Leiden on PCA - varying resolution
    'leiden|{"resolution":0.5,"n_neighbors":30}|true|leiden_pca_res0.5'
    'leiden|{"resolution":1.0,"n_neighbors":30}|true|leiden_pca_res1.0'
    'leiden|{"resolution":2.0,"n_neighbors":30}|true|leiden_pca_res2.0'
    'leiden|{"resolution":5.0,"n_neighbors":30}|true|leiden_pca_res5.0'

    # DBSCAN on PCA (with appropriate eps for PCA space)
    'dbscan|{"eps":3.0,"min_samples":50}|true|dbscan_pca_eps3.0'
    'dbscan|{"eps":5.0,"min_samples":50}|true|dbscan_pca_eps5.0'
    'dbscan|{"eps":10.0,"min_samples":50}|true|dbscan_pca_eps10.0'

    # Gaussian Mixture on PCA
    'gaussian_mixture|{"n_components":50}|true|gmm_pca_k50'
    'gaussian_mixture|{"n_components":100}|true|gmm_pca_k100'
    'gaussian_mixture|{"n_components":200}|true|gmm_pca_k200'

    # MiniBatch K-means on PCA (faster alternative for quick testing)
    'minibatch_kmeans|{"n_clusters":50}|true|minibatch_pca_k50'
    'minibatch_kmeans|{"n_clusters":100}|true|minibatch_pca_k100'
    'minibatch_kmeans|{"n_clusters":200}|true|minibatch_pca_k200'
)

# Get configuration for this task
CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

IFS='|' read -r METHOD PARAMS USE_PCA SUFFIX <<< "$CONFIG"

OUTPUT_FILE="${RESULTS_DIR}/clusters_${SUFFIX}.npz"

echo "=========================================="
echo "Clustering Sweep"
echo "Task: ${SLURM_ARRAY_TASK_ID} / ${#CONFIGS[@]}"
echo "Method: ${METHOD}"
echo "Parameters: ${PARAMS}"
echo "Use PCA: ${USE_PCA}"
echo "Output: ${OUTPUT_FILE}"
echo "=========================================="

# Build command
CMD="/home/dmullane/micromamba/envs/esm3_env/bin/python ${SCRIPT}"
CMD="${CMD} --umap ${UMAP_FILE}"
CMD="${CMD} --output ${OUTPUT_FILE}"
CMD="${CMD} --method ${METHOD}"
CMD="${CMD} --params '${PARAMS}'"

if [ "${USE_PCA}" = "true" ]; then
    CMD="${CMD} --use-pca --pca-cache ${PCA_CACHE}"
fi

echo ""
echo "Command:"
echo "${CMD}"
echo ""

# Run clustering
eval ${CMD}

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Exit code: ${EXIT_CODE}"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Clustering complete"
    echo "Output: ${OUTPUT_FILE}"
else
    echo "✗ Error"
fi
echo "=========================================="

exit ${EXIT_CODE}
