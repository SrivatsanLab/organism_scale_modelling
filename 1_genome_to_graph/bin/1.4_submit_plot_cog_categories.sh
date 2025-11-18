#!/bin/bash
#SBATCH --job-name=plot_cogs
#SBATCH --output=logs/plot_cogs_%j.out
#SBATCH --error=logs/plot_cogs_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new

# Plot individual COG categories on UMAP embeddings
# Creates one plot per COG category per n_neighbors value

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/plot_individual_cogs.py"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cog_categories"

echo "=========================================="
echo "Plotting Individual COG Categories"
echo "=========================================="
echo ""
echo "Script: ${SCRIPT}"
echo "Output: ${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cog_categories/"
echo ""
echo "Starting at $(date)"
echo "=========================================="
echo ""

cd "${PROJECT_DIR}"

# Run plotting script
/home/dmullane/micromamba/envs/esm3_env/bin/python "${SCRIPT}"

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Completed at $(date)"
echo "Exit code: ${EXIT_CODE}"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ All COG category plots generated successfully"
    echo "Output directory: ${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cog_categories/"
else
    echo "✗ Error generating plots"
fi

echo "=========================================="

exit ${EXIT_CODE}
