#!/bin/bash
#SBATCH --job-name=cluster_tightness
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/cluster_tightness_%j.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/cluster_tightness_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new

# Analyze MMseqs2 cluster tightness in ESM embedding space

set -e

echo "==========================================="
echo "MMseqs2 Cluster Tightness Analysis"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Step 1: Analyze cluster tightness
echo "Step 1: Computing cluster statistics..."
$PYTHON_ENV 1_genome_to_graph/1.4_esm_embedding_clustering/clustering/analyze_mmseqs_cluster_tightness.py \
    --pca-cache results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
    --mmseqs-clusters results/1_genome_to_graph/1.3_msa/mmseqs_seqid_0p7/clusters.tsv \
    --output-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis \
    --min-cluster-size 2

echo ""
echo "Step 2: Creating visualizations..."
$PYTHON_ENV 1_genome_to_graph/1.4_esm_embedding_clustering/clustering/visualize_cluster_tightness.py \
    --stats-file results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/mmseqs_cluster_statistics.csv \
    --per-dim-file results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/mmseqs_cluster_per_dimension_stats.csv \
    --output-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/cluster_analysis/figures

echo ""
echo "==========================================="
echo "Analysis Complete"
echo "End time: $(date)"
echo "==========================================="
