#!/bin/bash
#SBATCH --job-name=cluster_emb
#SBATCH --output=../logs/1.3_cluster_embeddings_%j.out
#SBATCH --error=../logs/1.3_cluster_embeddings_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new

#
# submit_cluster_embeddings.sh
# Create mean ESM embeddings for each MMseqs2 protein cluster
#
# This job:
# - Loads ~30M protein cluster assignments
# - Loads ESM embeddings from 7,664 genome files
# - Computes mean 1152D embedding for each of ~6M clusters
# - Saves cluster × embedding matrix
#
# Usage: sbatch scripts/1.3_submit_cluster_embeddings.sh
#

set -e

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.3_create_cluster_embeddings.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Creating cluster-level ESM embeddings"
echo "=========================================="
echo ""

# Run the script with the ESM environment Python
$PYTHON_ENV "$SCRIPT" \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --embeddings data/esm_embeddings \
    --output-dir 1_genome_to_graph/graph_outputs/protein_graph

echo ""
echo "Job completed successfully!"
