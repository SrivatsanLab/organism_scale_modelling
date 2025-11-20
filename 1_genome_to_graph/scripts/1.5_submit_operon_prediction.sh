#!/bin/bash
#SBATCH --job-name=operon_pred
#SBATCH --output=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.5_operon_prediction_%j.out
#SBATCH --error=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/logs/1.5_operon_prediction_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

#
# submit_operon_prediction.sh
# Predict operons from bacterial genomes to create protein-protein adjacency edges
#
# This job:
# - Runs Prodigal to predict genes and get coordinates (if not already done)
# - Predicts operons based on gene proximity and strand
# - Creates protein-protein edges for graph construction
#
# Usage: sbatch scripts/1.5_submit_operon_prediction.sh
#

set -e

PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_ROOT}/1_genome_to_graph/scripts/1.5_predict_operons.py"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Load Prodigal module if needed
module load prodigal/2.6.3-GCCcore-11.2.0

cd "$PROJECT_ROOT"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Predicting operons"
echo "=========================================="
echo ""

# Run the script
$PYTHON_ENV "$SCRIPT" \
    --genomes-dir data/refseq_genomes \
    --genome-list data/genome_metadata_with_files.csv \
    --max-distance 150 \
    --output-dir 1_genome_to_graph/graph_outputs/protein_graph \
    --intermediate-dir 1_genome_to_graph/intermediate/protein/operons \
    --prodigal-dir 1_genome_to_graph/intermediate/protein/gene_predictions

echo ""
echo "Job completed successfully!"
