#!/bin/bash
#SBATCH --job-name=gpu_embed
#SBATCH --output=logs/gpu_embed_%A_%a.out
#SBATCH --error=logs/gpu_embed_%A_%a.err
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new
#SBATCH --gres=gpu:1

# GPU-Accelerated Embeddings using RAPIDS cuML
# Runs UMAP, t-SNE, and clustering on GPU

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/compute_gpu_embeddings.py"
PCA_CACHE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"  # 1M genes, 50D PCA
OUTPUT_DIR="${PROJECT_DIR}/results/gpu_embeddings"
CONTAINER_DIR="${PROJECT_DIR}/containers"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${CONTAINER_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# RAPIDS container
RAPIDS_CONTAINER="docker://rapidsai/rapidsai:23.10-cuda11.8-runtime-ubuntu22.04-py3.10"
CONTAINER_IMAGE="${CONTAINER_DIR}/rapids_23.10.sif"

# Pull container if not exists
if [ ! -f "${CONTAINER_IMAGE}" ]; then
    echo "Pulling RAPIDS container..."
    singularity pull "${CONTAINER_IMAGE}" "${RAPIDS_CONTAINER}"
fi

# Configurations: METHOD|PARAM_NAME|PARAM_VALUE|CLUSTER|CLUSTER_PARAM|SUFFIX
CONFIGS=(
    # GPU UMAP with different n_neighbors
    "umap|n-neighbors|15|none|0|gpu_umap_n15"
    "umap|n-neighbors|50|none|0|gpu_umap_n50"
    "umap|n-neighbors|100|none|0|gpu_umap_n100"
    
    # GPU t-SNE
    "tsne|perplexity|30|none|0|gpu_tsne_p30"
    
    # GPU UMAP + HDBSCAN clustering
    "umap|n-neighbors|50|hdbscan|100|gpu_umap_n50_hdbscan100"
    "umap|n-neighbors|50|hdbscan|500|gpu_umap_n50_hdbscan500"
    
    # GPU UMAP + K-means clustering
    "umap|n-neighbors|50|kmeans|1000|gpu_umap_n50_kmeans1000"
    "umap|n-neighbors|50|kmeans|5000|gpu_umap_n50_kmeans5000"
)

CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"
IFS='|' read -r METHOD PARAM_NAME PARAM_VALUE CLUSTER CLUSTER_PARAM SUFFIX <<< "$CONFIG"

OUTPUT_FILE="${OUTPUT_DIR}/${SUFFIX}.npz"

echo "=========================================="
echo "GPU Embedding + Clustering"
echo "=========================================="
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Method: ${METHOD}"
echo "Clustering: ${CLUSTER}"
echo "Output: ${SUFFIX}"
echo "GPU: ${SLURM_JOB_GPUS}"
echo "=========================================="
echo ""

# Build command
CMD="python ${SCRIPT}"
CMD="${CMD} --input ${PCA_CACHE}"
CMD="${CMD} --output ${OUTPUT_FILE}"
CMD="${CMD} --method ${METHOD}"

if [ "${METHOD}" = "umap" ]; then
    CMD="${CMD} --n-neighbors ${PARAM_VALUE}"
elif [ "${METHOD}" = "tsne" ]; then
    CMD="${CMD} --perplexity ${PARAM_VALUE}"
fi

if [ "${CLUSTER}" = "hdbscan" ]; then
    CMD="${CMD} --cluster hdbscan --min-cluster-size ${CLUSTER_PARAM}"
elif [ "${CLUSTER}" = "kmeans" ]; then
    CMD="${CMD} --cluster kmeans --n-clusters ${CLUSTER_PARAM}"
fi

echo "Command:"
echo "${CMD}"
echo ""

# Run in RAPIDS container with GPU
singularity exec --nv \
    --bind ${PROJECT_DIR}:${PROJECT_DIR} \
    "${CONTAINER_IMAGE}" \
    ${CMD}

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Success"
else
    echo "✗ Failed"
fi
echo "=========================================="

exit ${EXIT_CODE}
