#!/bin/bash
#
# Run script to create toy protein graph (GitHub-friendly version)
#
# This is a lightweight script that can run on the login node or via SLURM.
# It filters to high-confidence clusters and applies PCA for compression.
#
# Usage:
#   bash scripts/1.3_run_toy_cluster_embeddings.sh
#   OR
#   sbatch scripts/1.3_run_toy_cluster_embeddings.sh  (if SLURM needed)
#

set -e

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.3_create_toy_cluster_embeddings.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Creating Toy Protein Graph"
echo "=========================================="
echo ""

# This script requires the full cluster embeddings to exist first
# Check if they exist
if [ ! -f "1_genome_to_graph/graph_outputs/protein_graph/cluster_embeddings.npz" ]; then
    echo "ERROR: Full cluster embeddings not found!"
    echo "Please run 1.3_submit_cluster_embeddings.sh first"
    exit 1
fi

# Run the script
$PYTHON_ENV "$SCRIPT" \
    --embeddings 1_genome_to_graph/graph_outputs/protein_graph/cluster_embeddings.npz \
    --index 1_genome_to_graph/graph_outputs/protein_graph/cluster_index.csv \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --min-seq-id 0.9 \
    --min-size 10 \
    --n-components 50 \
    --output-dir 1_genome_to_graph/graph_outputs/protein_graph_toy

echo ""
echo "Toy graph created successfully!"
echo "Output: 1_genome_to_graph/graph_outputs/protein_graph_toy/"
