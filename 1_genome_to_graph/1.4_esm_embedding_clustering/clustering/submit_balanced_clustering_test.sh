#!/bin/bash
#SBATCH --job-name=balanced_cluster
#SBATCH --output=logs/balanced_clustering_%j.out
#SBATCH --error=logs/balanced_clustering_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=campus-new

# Test Leiden clustering stability on balanced sample

echo "========================================="
echo "Balanced Sample Clustering Test"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "========================================="
echo

# Activate conda environment
source /home/dmullane/.bashrc
micromamba activate esm3_env

# Set parameters
GENE_IDS="data/balanced_sample_gene_ids.txt"
PCA_CACHE="results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"
OUTPUT="results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/balanced_stability_test.npz"

# Clustering parameters
RESOLUTION=750
N_NEIGHBORS=15
N_SUBSAMPLES=10
SUBSAMPLE_SIZE=4000

echo "Parameters:"
echo "  Balanced sample: $GENE_IDS"
echo "  PCA cache: $PCA_CACHE"
echo "  Resolution: $RESOLUTION"
echo "  n_neighbors: $N_NEIGHBORS"
echo "  Subsamples: $N_SUBSAMPLES"
echo "  Subsample size: $SUBSAMPLE_SIZE"
echo "  Output: $OUTPUT"
echo

# Check input files
if [ ! -f "$GENE_IDS" ]; then
    echo "ERROR: Gene IDs file not found: $GENE_IDS"
    exit 1
fi

if [ ! -f "$PCA_CACHE" ]; then
    echo "ERROR: PCA cache not found: $PCA_CACHE"
    exit 1
fi

echo "Input validation:"
echo "  Gene IDs: $(wc -l < $GENE_IDS) genes"
echo

# Create output directory
mkdir -p results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering

# Run stability evaluation
echo "========================================="
echo "Running stability evaluation..."
echo "========================================="
echo

python scripts/analysis/evaluate_balanced_clustering.py \
    --gene-ids $GENE_IDS \
    --pca-cache $PCA_CACHE \
    --resolution $RESOLUTION \
    --n-neighbors $N_NEIGHBORS \
    --n-subsamples $N_SUBSAMPLES \
    --subsample-size $SUBSAMPLE_SIZE \
    --output $OUTPUT \
    --seed 42

EXIT_CODE=$?

echo
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Evaluation complete!"
    echo "Results saved to: $OUTPUT"
else
    echo "❌ Evaluation failed with exit code: $EXIT_CODE"
fi
echo "========================================="
echo

exit $EXIT_CODE
