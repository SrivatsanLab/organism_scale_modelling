#!/bin/bash
#SBATCH --job-name=mmseqs_test
#SBATCH --output=logs/mmseqs_test_%j.out
#SBATCH --error=logs/mmseqs_test_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=campus-new

# MMseqs2 clustering test on downloaded protein sequences

echo "========================================="
echo "MMseqs2 Protein Clustering Test"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "========================================="
echo

# Load modules
module purge
module load MMseqs2/13-45111-gompi-2021b

# Check MMseqs2
echo "MMseqs2 version:"
mmseqs version
echo

# Activate conda environment
source /home/dmullane/.bashrc
micromamba activate esm3_env

# Set parameters
INPUT_FASTA="data/test_proteins_50genomes.faa"
OUTPUT_DIR="data/mmseqs_test"
MIN_SEQ_ID=0.5
COVERAGE=0.8
THREADS=8

echo "Parameters:"
echo "  Input: $INPUT_FASTA"
echo "  Output directory: $OUTPUT_DIR"
echo "  Min sequence identity: $MIN_SEQ_ID"
echo "  Min coverage: $COVERAGE"
echo "  Threads: $THREADS"
echo

# Check input file
if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: Input file not found: $INPUT_FASTA"
    exit 1
fi

echo "Input file size: $(wc -l < $INPUT_FASTA) lines"
echo "Number of sequences: $(grep -c '^>' $INPUT_FASTA)"
echo

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/tmp

# Run MMseqs2 clustering
echo "========================================="
echo "Running MMseqs2 clustering"
echo "========================================="
echo

# Create sequence database
echo "[1/4] Creating sequence database..."
mmseqs createdb $INPUT_FASTA $OUTPUT_DIR/seq_db

echo "[2/4] Running clustering..."
mmseqs cluster $OUTPUT_DIR/seq_db $OUTPUT_DIR/cluster_db $OUTPUT_DIR/tmp \
    --min-seq-id $MIN_SEQ_ID \
    -c $COVERAGE \
    --cluster-mode 0 \
    --threads $THREADS \
    -v 3

echo "[3/4] Converting to TSV..."
mmseqs createtsv $OUTPUT_DIR/seq_db $OUTPUT_DIR/seq_db \
    $OUTPUT_DIR/cluster_db $OUTPUT_DIR/clusters.tsv

echo "[4/4] Analyzing results..."
python scripts/analysis/analyze_mmseqs_clusters.py \
    --clusters $OUTPUT_DIR/clusters.tsv \
    --output $OUTPUT_DIR/cluster_summary.csv

echo
echo "========================================="
echo "Clustering complete!"
echo "========================================="
echo "Results saved to $OUTPUT_DIR/"
echo

# Clean up
echo "Cleaning up temporary files..."
rm -rf $OUTPUT_DIR/tmp
rm -f $OUTPUT_DIR/seq_db*
rm -f $OUTPUT_DIR/cluster_db*

echo "Done!"
