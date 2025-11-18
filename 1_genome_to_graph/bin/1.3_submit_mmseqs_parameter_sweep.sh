#!/bin/bash
#SBATCH --job-name=mmseqs_sweep
#SBATCH --output=logs/mmseqs_sweep_%A_%a.out
#SBATCH --error=logs/mmseqs_sweep_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --partition=campus-new

# Parameter sweep for MMseqs2 clustering
# Testing different sequence identity thresholds to find optimal cluster sizes

echo "========================================="
echo "MMseqs2 Parameter Sweep"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "========================================="
echo

# Activate conda environment
source /home/dmullane/.bashrc
micromamba activate esm3_env

# Load MMseqs2 module
module load MMseqs2/13-45111-gompi-2021b

# Parameter arrays
SEQ_IDS=(0.3 0.4 0.5 0.6 0.7)
COVERAGE=0.8
CLUSTER_MODE=0
THREADS=$SLURM_CPUS_PER_TASK

# Get parameters for this task
SEQ_ID=${SEQ_IDS[$SLURM_ARRAY_TASK_ID]}

# Paths
PROTEIN_FASTA="data/all_proteins.faa"
OUTPUT_DIR="data/mmseqs_seqid_${SEQ_ID/./p}"  # e.g., mmseqs_seqid_0p3

echo "Parameters for this task:"
echo "  Sequence identity: $SEQ_ID"
echo "  Coverage: $COVERAGE"
echo "  Cluster mode: $CLUSTER_MODE"
echo "  Threads: $THREADS"
echo "  Output directory: $OUTPUT_DIR"
echo

# Check input file
if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "ERROR: Protein FASTA not found: $PROTEIN_FASTA"
    exit 1
fi

echo "Input validation:"
echo "  Protein FASTA: $PROTEIN_FASTA"
echo "  File size: $(du -h $PROTEIN_FASTA | cut -f1)"
echo "  Protein count: $(grep -c '^>' $PROTEIN_FASTA)"
echo

# Run MMseqs2 clustering
echo "========================================="
echo "Running MMseqs2 clustering..."
echo "========================================="
echo

python scripts/analysis/cluster_proteins_mmseqs.py \
    --input "$PROTEIN_FASTA" \
    --output-dir "$OUTPUT_DIR" \
    --min-seq-id $SEQ_ID \
    --coverage $COVERAGE \
    --cluster-mode $CLUSTER_MODE \
    --threads $THREADS

EXIT_CODE=$?

echo
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ MMseqs2 clustering complete!"
    echo "Results in: $OUTPUT_DIR"
    echo
    echo "Quick summary:"
    echo "  Sequence identity: $SEQ_ID"
    tail -5 $OUTPUT_DIR/cluster_summary.csv | head -3
else
    echo "❌ MMseqs2 clustering failed with exit code: $EXIT_CODE"
fi
echo "========================================="
echo

exit $EXIT_CODE
