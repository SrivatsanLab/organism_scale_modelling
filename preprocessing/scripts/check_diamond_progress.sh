#!/bin/bash
# Check progress of diamond COG annotation job

OUTPUT_DIR="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/1_genome_to_graph/intermediate/protein/cog_annotations"
TOTAL_GENOMES=7664

echo "======================================"
echo "Diamond COG Annotation Progress"
echo "======================================"
echo ""

# Count completed genomes
COMPLETED=$(find "$OUTPUT_DIR" -name "*_cog_hits.tsv" 2>/dev/null | wc -l)
PERCENT=$(awk "BEGIN {printf \"%.1f\", $COMPLETED/$TOTAL_GENOMES*100}")

echo "Completed: $COMPLETED / $TOTAL_GENOMES genomes ($PERCENT%)"
echo ""

# Check running jobs
echo "SLURM job status:"
squeue -u dmullane -n diamond_cog | head -10

echo ""
echo "Recent completions (last 5):"
find "$OUTPUT_DIR" -name "*_cog_hits.tsv" -printf '%T+ %p\n' 2>/dev/null | sort -r | head -5 | while read -r line; do
    timestamp=$(echo "$line" | cut -d' ' -f1)
    file=$(basename "$(echo "$line" | cut -d' ' -f2)" _cog_hits.tsv)
    echo "  $timestamp - $file"
done

echo ""
echo "Recent errors (last 5):"
tail -5 logs/diamond_cog_43078588_*.err 2>/dev/null | head -20

echo ""
echo "To regenerate metadata after completion, run:"
echo "python scripts/create_cluster_metadata.py \\"
echo "  --cluster-index graph_outputs/protein_graph/cluster_index.csv \\"
echo "  --clusters intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \\"
echo "  --cog-dir intermediate/protein/cog_annotations \\"
echo "  --output graph_outputs/protein_graph/metadata/cluster_metadata.csv"
