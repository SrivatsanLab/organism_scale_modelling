#!/bin/bash
#SBATCH --job-name=toy_cluster_emb
#SBATCH --output=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.3_toy_cluster_embeddings_%j.out
#SBATCH --error=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.3_toy_cluster_embeddings_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

#
# submit_toy_cluster_embeddings.sh
# Create toy protein graph directly from embeddings (GitHub-friendly version)
#
# This job creates a compressed version with:
# - High-confidence clusters only (size >= 10)
# - PCA-reduced embeddings (1152D → 50D)
# - Total size < 100 MB
#
# Usage: sbatch scripts/1.3_submit_toy_cluster_embeddings.sh
#

set -e

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.3_create_toy_cluster_embeddings_direct.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Creating toy protein graph"
echo "=========================================="
echo ""

# Run the script with the ESM environment Python
$PYTHON_ENV "$SCRIPT" \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --embeddings data/esm_embeddings \
    --min-size 10 \
    --n-components 50 \
    --output-dir 1_genome_to_graph/graph_outputs/protein_graph_toy

echo ""
echo "Job completed successfully!"
