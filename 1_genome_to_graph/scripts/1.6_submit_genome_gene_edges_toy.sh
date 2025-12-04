#!/bin/bash
#SBATCH --job-name=genome_gene_toy
#SBATCH --output=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.6_genome_gene_edges_toy_%j.out
#SBATCH --error=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.6_genome_gene_edges_toy_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

#
# submit_genome_gene_edges.sh
# Create genome-gene presence/absence matrix for heterogeneous graph
#
# This job:
# - Reads MMseqs2 cluster assignments
# - Creates binary matrix of which clusters (genes) are in which genomes
# - Saves as sparse matrix and edge list
#
# Usage: sbatch scripts/1.6_submit_genome_gene_edges.sh
#

set -euo pipefail

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.6_create_genome_gene_edges.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Creating genome-gene edge matrix"
echo "=========================================="
echo ""

# Verify input files exist
echo "Checking input files..."
if [ ! -f "1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv" ]; then
    echo "ERROR: Cluster assignments file not found"
    exit 1
fi

if [ ! -f "1_genome_to_graph/graph_outputs/protein_graph_toy/cluster_index.csv" ]; then
    echo "ERROR: Toy cluster index file not found"
    exit 1
fi

if [ ! -f "data/genome_metadata.csv" ]; then
    echo "ERROR: Genome metadata file not found"
    exit 1
fi

echo "All input files found"
echo ""

# Run the script
$PYTHON_ENV "$SCRIPT" \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --cluster-stats 1_genome_to_graph/graph_outputs/protein_graph_toy/cluster_index.csv \
    --genome-list data/genome_metadata.csv \
    --output 1_genome_to_graph/graph_outputs/protein_graph_toy

echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
