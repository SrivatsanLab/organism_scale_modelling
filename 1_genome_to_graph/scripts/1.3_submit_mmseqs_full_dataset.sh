#!/bin/bash
#SBATCH --job-name=mmseqs_full
#SBATCH --output=logs/mmseqs_full_%j.out
#SBATCH --error=logs/mmseqs_full_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --partition=campus-new

# Run MMseqs2 clustering on all ~29M proteins from 7,664 genomes

echo "========================================="
echo "MMseqs2 Full Dataset Clustering"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "========================================="
echo

# Activate conda environment
source /home/dmullane/.bashrc
micromamba activate esm3_env

# Load MMseqs2 module
module load MMseqs2/13-45111-gompi-2021b

# Parameters
INPUT_DIR="/fh/working/srivatsan_s/ShotgunDomestication/data/refseq_genomes/gene_annotations"
OUTPUT_DIR="data/mmseqs_full_dataset"
PROTEIN_FASTA="data/all_proteins.faa"

# MMseqs2 parameters
MIN_SEQ_ID=0.5
COVERAGE=0.8
CLUSTER_MODE=0  # Greedy clustering
THREADS=$SLURM_CPUS_PER_TASK

echo "Parameters:"
echo "  Input directory: $INPUT_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Protein FASTA: $PROTEIN_FASTA"
echo "  Min sequence identity: $MIN_SEQ_ID"
echo "  Coverage: $COVERAGE"
echo "  Cluster mode: $CLUSTER_MODE"
echo "  Threads: $THREADS"
echo

# Step 1: Concatenate all protein sequences
echo "========================================="
echo "Step 1: Concatenating protein sequences"
echo "========================================="
echo

if [ -f "$PROTEIN_FASTA" ]; then
    echo "Protein FASTA already exists: $PROTEIN_FASTA"
    echo "Skipping concatenation step."
    echo "File size: $(du -h $PROTEIN_FASTA | cut -f1)"
    echo "Protein count: $(grep -c '^>' $PROTEIN_FASTA)"
else
    python scripts/analysis/concatenate_all_proteins.py \
        --input-dir "$INPUT_DIR" \
        --output "$PROTEIN_FASTA" \
        --pattern "*_prodigal_proteins.faa"

    if [ $? -ne 0 ]; then
        echo "ERROR: Protein concatenation failed!"
        exit 1
    fi
fi

echo

# Step 2: Run MMseqs2 clustering
echo "========================================="
echo "Step 2: Running MMseqs2 clustering"
echo "========================================="
echo

python scripts/analysis/cluster_proteins_mmseqs.py \
    --input "$PROTEIN_FASTA" \
    --output-dir "$OUTPUT_DIR" \
    --min-seq-id $MIN_SEQ_ID \
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
    echo "Cluster summary:"
    wc -l $OUTPUT_DIR/clusters.tsv
    head -5 $OUTPUT_DIR/cluster_summary.csv
else
    echo "❌ MMseqs2 clustering failed with exit code: $EXIT_CODE"
fi
echo "========================================="
echo

exit $EXIT_CODE
