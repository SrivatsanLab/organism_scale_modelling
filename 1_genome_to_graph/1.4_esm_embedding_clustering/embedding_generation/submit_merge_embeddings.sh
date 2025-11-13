#!/bin/bash
#SBATCH --job-name=merge_embeddings
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/merge_embeddings_%j.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/merge_embeddings_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new

# Merge embedding batches and create PCA cache

set -e

echo "==========================================="
echo "Merge Embedding Batches"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Run merge script
$PYTHON_ENV 1_genome_to_graph/1.4_esm_embedding_clustering/embedding_generation/merge_embedding_batches.py

echo ""
echo "==========================================="
echo "Merge Complete"
echo "End time: $(date)"
echo "==========================================="
