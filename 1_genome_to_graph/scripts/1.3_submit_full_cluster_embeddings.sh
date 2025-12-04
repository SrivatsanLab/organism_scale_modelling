#!/bin/bash
#SBATCH --job-name=full_cluster_emb
#SBATCH --output=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.3_full_cluster_embeddings_%j.out
#SBATCH --error=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.3_full_cluster_embeddings_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new

# Full cluster embeddings creation (1152D for 388k clusters)
# Uses direct approach from embedding files (like toy graph)

set -euo pipefail

SCRIPT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/1_genome_to_graph/scripts/1.3_create_full_cluster_embeddings.py"

# Python environment with sklearn
export PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Change to the project root directory
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling

echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="
echo ""

# Verify files exist
echo "Checking input files..."
if [ ! -f "1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv" ]; then
    echo "ERROR: Cluster statistics file not found"
    exit 1
fi

if [ ! -f "1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv" ]; then
    echo "ERROR: Clusters TSV file not found"
    exit 1
fi

if [ ! -d "data/esm_embeddings" ]; then
    echo "ERROR: ESM embeddings directory not found"
    exit 1
fi

echo "All input files found"
echo ""

echo "=========================================="
echo "Creating full cluster embeddings (1152D)"
echo "=========================================="
echo ""

# Run the script
$PYTHON_ENV "$SCRIPT" \
    --cluster-stats 1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --embeddings data/esm_embeddings \
    --output 1_genome_to_graph/graph_outputs/protein_graph

echo ""
echo "=========================================="
echo "Job finished at: $(date)"
echo "=========================================="
