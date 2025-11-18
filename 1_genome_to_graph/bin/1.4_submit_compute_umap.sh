#!/bin/bash
#SBATCH --job-name=compute_umap
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/compute_umap_%j.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/compute_umap_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --partition=campus-new

# Compute UMAP embeddings from PCA cache for all 11.8M proteins

set -e

echo "==========================================="
echo "Compute UMAP Embeddings"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Input: PCA cache with 11.8M proteins
PCA_CACHE="results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz"

# Output: UMAP embeddings
OUTPUT="results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_full_n15.npz"

echo "Input PCA cache: $PCA_CACHE"
echo "Output UMAP file: $OUTPUT"
echo ""

# Run UMAP computation
$PYTHON_ENV 1_genome_to_graph/1.4_esm_embedding_clustering/clustering/compute_umap_array.py \
    --load-cached-pca "$PCA_CACHE" \
    --output "$OUTPUT" \
    --n-neighbors 15 \
    --min-dist 0.1 \
    --n-pcs 50 \
    --seed 42

echo ""
echo "==========================================="
echo "UMAP Computation Complete"
echo "End time: $(date)"
echo "==========================================="
