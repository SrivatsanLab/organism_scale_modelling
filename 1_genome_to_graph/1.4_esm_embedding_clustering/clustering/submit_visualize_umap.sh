#!/bin/bash
#SBATCH --job-name=umap_viz
#SBATCH --output=/fh/working/srivatsan_s/dmullane_organism_scale/logs/umap_viz_%j.out
#SBATCH --error=/fh/working/srivatsan_s/dmullane_organism_scale/logs/umap_viz_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=campus-new

# Create UMAP visualizations

set -e

echo "==========================================="
echo "UMAP Visualization"
echo "==========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

PYTHON_ENV="/home/dmullane/micromamba/envs/esm3_env/bin/python"

$PYTHON_ENV -u 1_genome_to_graph/1.4_esm_embedding_clustering/clustering/visualize_umap.py \
    --umap-file results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_full_n15.npz \
    --mmseqs-clusters results/1_genome_to_graph/1.3_msa/mmseqs_seqid_0p7/clusters.tsv \
    --output-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/figures \
    --sample-size 100000

echo ""
echo "==========================================="
echo "Visualization Complete"
echo "End time: $(date)"
echo "==========================================="
