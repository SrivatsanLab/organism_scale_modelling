#!/bin/bash
#SBATCH --job-name=quality_eval
#SBATCH --output=logs/quality_eval_%j.out
#SBATCH --error=logs/quality_eval_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=campus-new

# Comprehensive quality evaluation for all Leiden clusterings
# Computes ARI, AMI, silhouette, Davies-Bouldin, size distribution

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/evaluate_clustering_quality.py"
PCA_FILE="${PROJECT_DIR}/results/umap/pca_cache.npz"
OUTPUT_FILE="${PROJECT_DIR}/results/clustering/quality_metrics_comprehensive.csv"

cd "${PROJECT_DIR}"

mkdir -p logs
mkdir -p results/clustering

echo "=========================================="
echo "Comprehensive Quality Evaluation"
echo "=========================================="
echo "PCA file: ${PCA_FILE}"
echo "Output: ${OUTPUT_FILE}"
echo "Start time: $(date)"
echo ""

/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}" \
    --pca "${PCA_FILE}" \
    --clustering results/clustering \
    --output "${OUTPUT_FILE}" \
    --sample-size 50000 \
    --batch

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
