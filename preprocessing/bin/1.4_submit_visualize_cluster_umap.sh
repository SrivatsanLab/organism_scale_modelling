#!/bin/bash
#SBATCH --job-name=cluster_umap
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/cluster_umap_%j.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/cluster_umap_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

# Create cluster-level UMAP visualizations

set -e

echo "==========================================="
echo "Cluster-Level UMAP Visualization"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

$PYTHON_ENV -u 1_genome_to_graph/1.4_esm_embedding_clustering/clustering/visualize_cluster_umap.py \
    --umap-file results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_full_n15.npz \
    --cluster-stats results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/mmseqs_cluster_statistics.csv \
    --pca-cache results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
    --output-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/figures

echo ""
echo "==========================================="
echo "Visualization Complete"
echo "End time: $(date)"
echo "==========================================="
