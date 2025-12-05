#!/bin/bash
#SBATCH --job-name=gene_gene_edges
#SBATCH --output=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.5_gene_gene_edges_%j.out
#SBATCH --error=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.5_gene_gene_edges_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

#
# submit_gene_gene_edges.sh
# Create gene-gene (cluster-cluster) adjacency matrix from operon predictions
#
# This job:
# - Runs Prodigal to predict genes
# - Predicts operons based on gene proximity and strand
# - Creates cluster-cluster edges when genes co-occur in operons
#
# Usage: sbatch scripts/1.5_submit_gene_gene_edges.sh
#

set -euo pipefail

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.5_create_gene_gene_edges.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Load Prodigal module
module load prodigal/2.6.3-GCCcore-11.2.0

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Creating gene-gene edges from operons"
echo "=========================================="
echo ""

# Verify input files exist
echo "Checking input files..."
if [ ! -f "1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv" ]; then
    echo "ERROR: Cluster assignments file not found"
    exit 1
fi

if [ ! -f "1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv" ]; then
    echo "ERROR: Cluster statistics file not found"
    exit 1
fi

if [ ! -f "data/genome_metadata_with_files.csv" ]; then
    echo "ERROR: Genome metadata file not found"
    exit 1
fi

echo "All input files found"
echo ""

# Run the script
$PYTHON_ENV "$SCRIPT" \
    --genomes-dir data/refseq_genomes \
    --genome-list data/genome_metadata_with_files.csv \
    --clusters 1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
    --cluster-stats 1_genome_to_graph/intermediate/protein/esm_embeddings/cluster_analysis/mmseqs_cluster_statistics.csv \
    --max-distance 150 \
    --min-cooccurrence 1 \
    --output-dir 1_genome_to_graph/graph_outputs/protein_graph \
    --intermediate-dir 1_genome_to_graph/intermediate/protein/operons \
    --prodigal-dir 1_genome_to_graph/intermediate/protein/gene_predictions

echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
