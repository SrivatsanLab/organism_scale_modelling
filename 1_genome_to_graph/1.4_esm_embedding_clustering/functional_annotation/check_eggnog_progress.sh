#!/bin/bash
# Quick script to check eggNOG-mapper annotation progress

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"

echo "========================================"
echo "eggNOG-mapper Progress Check"
echo "========================================"
echo ""

# Count completed genomes
COMPLETED=$(ls ${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation/*/GC*.emapper.annotations 2>/dev/null | wc -l)
TOTAL=7664

echo "Completed: ${COMPLETED} / ${TOTAL} genomes"
echo "Progress: $(echo "scale=1; 100*${COMPLETED}/${TOTAL}" | bc)%"
echo ""

# Check running jobs
echo "Running jobs:"
squeue -u $USER -n eggnog -o "%.10i %.8T %.10M" | head -10
RUNNING=$(squeue -u $USER -n eggnog -h | wc -l)
echo "Total running/pending: ${RUNNING}"
echo ""

# Show recent completions
echo "Most recent completions:"
ls -lt ${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation/*/GC*.emapper.annotations 2>/dev/null | head -5 | awk '{print $9}'
echo ""

# Estimate time remaining
if [ ${COMPLETED} -gt 0 ]; then
    FIRST_LOG=$(ls -t ${PROJECT_DIR}/logs/eggnog_*.out | tail -1)
    if [ -f "${FIRST_LOG}" ]; then
        # Get start time from log
        echo "Estimated time remaining will be calculated after more genomes complete..."
    fi
fi

echo "========================================"
