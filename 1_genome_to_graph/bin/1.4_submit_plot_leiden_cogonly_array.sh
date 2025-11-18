#!/bin/bash
#SBATCH --job-name=plot_leiden_cog
#SBATCH --output=logs/plot_leiden_cog_%A_%a.out
#SBATCH --error=logs/plot_leiden_cog_%A_%a.err
#SBATCH --array=0-274
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new

# Plot COG-only Leiden clusterings using COG-only UMAPs
# ~55 COG-only clusterings × 5 UMAPs = 275 jobs

PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
SCRIPT="${PROJECT_DIR}/scripts/embeddings/plot_leiden_clustering.py"
RESULTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering"
PLOTS_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/plots/cogonly"

cd "${PROJECT_DIR}"

# Create output directories
mkdir -p "${PLOTS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Get list of COG-only Leiden clustering files (sorted)
readarray -t CLUSTERING_FILES < <(ls -1 "${RESULTS_DIR}"/clusters_leiden_*_cogonly.npz | sort)
N_CLUSTERINGS=${#CLUSTERING_FILES[@]}

# UMAP n_neighbors values
N_NEIGHBORS=(15 25 50 100 200)
N_UMAPS=${#N_NEIGHBORS[@]}

# Calculate which clustering and UMAP this task is
CLUSTERING_IDX=$((SLURM_ARRAY_TASK_ID / N_UMAPS))
UMAP_IDX=$((SLURM_ARRAY_TASK_ID % N_UMAPS))

# Check if indices are valid
if [ ${CLUSTERING_IDX} -ge ${N_CLUSTERINGS} ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID}: Clustering index ${CLUSTERING_IDX} >= ${N_CLUSTERINGS}, exiting"
    exit 0
fi

# Get the actual files/values
CLUSTERING_FILE="${CLUSTERING_FILES[$CLUSTERING_IDX]}"
UMAP_N="${N_NEIGHBORS[$UMAP_IDX]}"

# Extract clustering name from file
CLUSTERING_NAME=$(basename "${CLUSTERING_FILE}" .npz | sed 's/clusters_//')

# Output file path (in separate cogonly subdirectory)
OUTPUT_DIR="${PLOTS_DIR}/umap_n${UMAP_N}"
OUTPUT_FILE="${OUTPUT_DIR}/${CLUSTERING_NAME}.png"

echo "=========================================="
echo "Plot COG-Only Leiden Clustering"
echo "Task: ${SLURM_ARRAY_TASK_ID}/${N_CLUSTERINGS}×${N_UMAPS}"
echo "=========================================="
echo "Clustering: ${CLUSTERING_NAME}"
echo "UMAP n_neighbors: ${UMAP_N} (COG-only)"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Check if already done
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Already exists: ${OUTPUT_FILE}"
    echo "Skipping..."
    exit 0
fi

# Check if COG-only UMAP exists
UMAP_FILE="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n${UMAP_N}_subsample1000000_cogonly.npz"
if [ ! -f "${UMAP_FILE}" ]; then
    echo "ERROR: COG-only UMAP not found: ${UMAP_FILE}"
    echo "Please run compute_umap_cogonly.py first"
    exit 1
fi

# Run plotting with COG-only UMAP
echo "Starting plot generation..."
echo "Start time: $(date)"
echo ""

# Temporarily modify the script to use COG-only UMAP
# Use Python to plot with custom UMAP path
/home/dmullane/micromamba/envs/esm3_env/bin/python << EOF
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '${PROJECT_DIR}')
from scripts.embeddings.plot_all_clusterings import COG_COLORS

# Load clustering
clustering_data = np.load('${CLUSTERING_FILE}', allow_pickle=True)
cluster_labels = clustering_data['labels']
cluster_gene_ids = clustering_data['gene_ids']

# Load COG-only UMAP
umap_data = np.load('${UMAP_FILE}', allow_pickle=True)
umap_coords = umap_data['umap_embedding']
umap_gene_ids = umap_data['gene_ids']

# Since both are COG-only, gene IDs should match exactly
assert len(cluster_gene_ids) == len(umap_gene_ids), "Gene count mismatch!"

# Verify gene IDs match (they should be identical)
cluster_set = set(cluster_gene_ids)
umap_set = set(umap_gene_ids)
overlap = cluster_set & umap_set
print(f"Gene overlap: {len(overlap):,} / {len(cluster_gene_ids):,} ({100*len(overlap)/len(cluster_gene_ids):.1f}%)")

# Map cluster labels to UMAP order
gene_to_cluster = dict(zip(cluster_gene_ids, cluster_labels))
aligned_labels = np.array([gene_to_cluster.get(gid, -1) for gid in umap_gene_ids])

# For COG-only, all genes are annotated, so create dummy COG categories
# (We'll just show that they're all annotated)
cog_categories = np.array(['COG' for _ in umap_gene_ids])

# Plot
fig, axes = plt.subplots(1, 2, figsize=(24, 10))

# LEFT: Clusters
ax = axes[0]
mask_noise = aligned_labels == -1
mask_clustered = ~mask_noise

if np.any(mask_noise):
    ax.scatter(umap_coords[mask_noise, 0], umap_coords[mask_noise, 1],
               c='lightgray', s=0.1, alpha=0.3, rasterized=True)
if np.any(mask_clustered):
    ax.scatter(umap_coords[mask_clustered, 0], umap_coords[mask_clustered, 1],
               c=aligned_labels[mask_clustered], s=0.5, alpha=0.5,
               cmap='tab20', rasterized=True)

n_clusters = len(set(aligned_labels)) - (1 if -1 in aligned_labels else 0)
n_noise = np.sum(mask_noise)
ax.set_title(f'${CLUSTERING_NAME}\\\n{n_clusters:,} clusters (COG-annotated genes only)',
             fontsize=14, fontweight='bold')
ax.set_xlabel(f'UMAP 1 (n_neighbors=${UMAP_N}, COG-only)', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.grid(True, alpha=0.2)

# RIGHT: Show it's all annotated
ax = axes[1]
ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
           c='green', s=0.5, alpha=0.3, rasterized=True, label='COG annotated')
ax.set_title(f'COG-Annotated Genes Only\\\n100% annotated ({len(umap_gene_ids):,} genes)',
             fontsize=14, fontweight='bold')
ax.set_xlabel(f'UMAP 1 (n_neighbors=${UMAP_N}, COG-only)', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.grid(True, alpha=0.2)
ax.legend(fontsize=12)

plt.tight_layout()

output_path = Path('${OUTPUT_FILE}')
output_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"Saved: {output_path}")
EOF

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "=========================================="
    echo "✓ Success: ${CLUSTERING_NAME} × UMAP n=${UMAP_N} (COG-only)"
    echo "=========================================="
else
    echo "=========================================="
    echo "✗ Failed with exit code ${EXIT_CODE}"
    echo "=========================================="
fi

exit ${EXIT_CODE}
