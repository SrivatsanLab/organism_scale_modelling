#!/bin/bash
#SBATCH --job-name=batch_embeddings
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/batch_embeddings_%A_%a.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/batch_embeddings_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new
#SBATCH --gres=gpu:1
#SBATCH --array=0-1183%50

# Array job to generate ESM embeddings in batches
# With batch size of 10,000 proteins, 11,837,414 proteins requires 1,184 batches
# Limit to 50 concurrent jobs to avoid overwhelming the cluster

set -e

echo "==========================================="
echo "Batch Embedding Generation"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Configuration
FASTA_FILE="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/data/proteins_for_embedding.faa"
BATCH_SIZE=10000
OUTPUT_DIR="/fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches"
PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

# Create output directory
mkdir -p $OUTPUT_DIR

# Run batch embedding generation
echo "Processing batch $SLURM_ARRAY_TASK_ID..."
echo "  Batch size: $BATCH_SIZE"
echo "  Output: $OUTPUT_DIR"
echo ""

$PYTHON_ENV 1_genome_to_graph/1.4_esm_embedding_clustering/embedding_generation/batch_generate_embeddings.py \
    --fasta $FASTA_FILE \
    --batch-idx $SLURM_ARRAY_TASK_ID \
    --batch-size $BATCH_SIZE \
    --output-dir $OUTPUT_DIR \
    --device cuda

echo ""
echo "==========================================="
echo "Batch $SLURM_ARRAY_TASK_ID Complete"
echo "End time: $(date)"
echo "==========================================="
