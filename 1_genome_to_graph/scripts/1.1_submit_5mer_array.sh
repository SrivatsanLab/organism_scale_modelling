#!/bin/bash
#SBATCH --job-name=5mer_profile
#SBATCH --output=../../logs/1.1_kmer_profile_%A_%a.out
#SBATCH --error=../../logs/1.1_kmer_profile_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new
#SBATCH --array=1-7664

#
# submit_5mer_array.sh
# SLURM array job to generate 5-mer profiles for all 7,664 bacterial genomes
#
# Usage: sbatch submit_5mer_array.sh
#

set -e

# Directories
PROJECT_ROOT="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling"
GENOME_DIR="$PROJECT_ROOT/data/refseq_genomes"
OUTPUT_DIR="$PROJECT_ROOT/results/1_genome_to_graph/1.1_kmer_profiling"
SCRIPT_DIR="$PROJECT_ROOT/1_genome_to_graph/1.1_kmer_profiling"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get genome file for this array task
GENOME_FILE=$(ls "$GENOME_DIR"/*.fasta | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$GENOME_FILE" ]; then
    echo "Error: Could not find genome file for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

GENOME_ID=$(basename "$GENOME_FILE" | sed 's/\.fasta$//')

echo "=========================================="
echo "SLURM Array Job: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Genome: $GENOME_ID"
echo "File: $GENOME_FILE"
echo "=========================================="
echo ""

# Check if output already exists (for resuming interrupted jobs)
DUMP_FILE="${OUTPUT_DIR}/${GENOME_ID}_5mer.txt"
if [ -f "$DUMP_FILE" ]; then
    echo "Output already exists, skipping: $DUMP_FILE"
    exit 0
fi

# Run 5-mer profiling
bash "$SCRIPT_DIR/generate_5mer_profiles.sh" "$GENOME_FILE" "$OUTPUT_DIR"

echo ""
echo "=========================================="
echo "Task completed successfully!"
echo "=========================================="
