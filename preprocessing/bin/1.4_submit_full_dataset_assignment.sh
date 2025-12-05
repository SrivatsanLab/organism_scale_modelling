#!/bin/bash
#SBATCH --job-name=full_assign
#SBATCH --output=logs/full_assign_%j.out
#SBATCH --error=logs/full_assign_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=campus-new

# Assign all 29M genes to clusters based on best subsample clustering
# This is a critical step for getting final cluster assignments

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/assign_full_dataset_to_clusters.py"
CLUSTERING_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/clusters_leiden_res1500_nn15_cogonly.npz"
PCA_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"
EMBEDDING_DIR="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/data/refseq_esm_embeddings"
OUTPUT_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/full_dataset_assignments_res1500_nn15_cogonly.npz"

cd "${PROJECT_DIR}"

mkdir -p logs
mkdir -p results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering

echo "=========================================="
echo "Full Dataset Cluster Assignment (29M genes)"
echo "=========================================="
echo "Clustering: ${CLUSTERING_FILE}"
echo "PCA cache: ${PCA_FILE}"
echo "Embeddings: ${EMBEDDING_DIR}"
echo "Output: ${OUTPUT_FILE}"
echo "Start time: $(date)"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "WARNING: Output file already exists: ${OUTPUT_FILE}"
    echo "Remove it first if you want to regenerate."
    # exit 0
fi

# Run assignment
echo "Starting Python script..."

/home/dmullane/micromamba/envs/esm3_env/bin/python -u "${SCRIPT}" \
    --clustering "${CLUSTERING_FILE}" \
    --pca "${PCA_FILE}" \
    --embedding-dir "${EMBEDDING_DIR}" \
    --output "${OUTPUT_FILE}" \
    --compute-statistics \
    --batch-size 100000

EXIT_CODE=$?
echo "Python script completed with exit code: ${EXIT_CODE}"

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success"
    echo "=========================================="
    echo ""
    echo "Output file: ${OUTPUT_FILE}"
    echo ""
    echo "To load results in Python:"
    echo "  import numpy as np"
    echo "  data = np.load('${OUTPUT_FILE}')"
    echo "  gene_ids = data['gene_ids']"
    echo "  cluster_assignments = data['cluster_assignments']"
    echo "  cluster_means = data['cluster_means']"
    echo "  cluster_stds = data['cluster_stds']"
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
