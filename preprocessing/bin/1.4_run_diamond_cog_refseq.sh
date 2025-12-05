#!/bin/bash
#SBATCH --job-name=diamond_cog
#SBATCH --output=logs/diamond_cog_%A_%a.out
#SBATCH --error=logs/diamond_cog_%A_%a.err
#SBATCH --array=0-7663%50
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=campus-new

# DIAMOND COG Annotation for RefSeq genomes
# Faster alternative to eggNOG-mapper

# Configuration
PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
GENE_DIR="${PROJECT_DIR}/data/refseq_gene_annotations"
OUTPUT_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation"
COG_DB="/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/cog"

# Set temp directory to /tmp
export TMPDIR="/tmp/${USER}/diamond_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Load DIAMOND module
module purge
module load DIAMOND/2.0.13-GCC-11.2.0

# Get list of protein files
PROTEIN_FILES=(${GENE_DIR}/*_prodigal_proteins.faa)
TOTAL_FILES=${#PROTEIN_FILES[@]}

# Check if array task ID is valid
if [ ${SLURM_ARRAY_TASK_ID} -ge ${TOTAL_FILES} ]; then
    echo "Array task ID ${SLURM_ARRAY_TASK_ID} >= total files ${TOTAL_FILES}, exiting"
    exit 0
fi

# Get protein file for this task
PROTEIN_FILE="${PROTEIN_FILES[$SLURM_ARRAY_TASK_ID]}"
GENOME_ID=$(basename "${PROTEIN_FILE}" _prodigal_proteins.faa)

echo "=========================================="
echo "DIAMOND COG Annotation"
echo "Task: ${SLURM_ARRAY_TASK_ID}/${TOTAL_FILES}"
echo "Genome: ${GENOME_ID}"
echo "=========================================="

# Output directory for this genome
SAMPLE_OUTPUT="${OUTPUT_DIR}/${GENOME_ID}_cog_diamond"
mkdir -p "${SAMPLE_OUTPUT}"

# Check if already processed
if [ -f "${SAMPLE_OUTPUT}/${GENOME_ID}_cog_hits.tsv" ]; then
    echo "Already processed: ${GENOME_ID}"
    exit 0
fi

# Count proteins
N_PROTEINS=$(grep -c "^>" "${PROTEIN_FILE}" || echo 0)
echo "Number of proteins: ${N_PROTEINS}"

if [ "${N_PROTEINS}" -eq 0 ]; then
    echo "Error: No proteins found in input file"
    exit 1
fi

echo ""
echo "Running DIAMOND search against COG database..."
echo "Start time: $(date)"

START_TIME=$(date +%s)

diamond blastp \
  --db "${COG_DB}/cog-20" \
  --query "${PROTEIN_FILE}" \
  --out "${SAMPLE_OUTPUT}/${GENOME_ID}_cog_hits.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads ${SLURM_CPUS_PER_TASK} \
  --sensitive \
  --tmpdir "${TMPDIR}"

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "End time: $(date)"
echo "Elapsed time: ${ELAPSED} seconds"

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "✗ Error: DIAMOND failed with exit code ${EXIT_CODE}"
    rm -rf "${TMPDIR}"
    exit ${EXIT_CODE}
fi

# Count hits
HITS_FILE="${SAMPLE_OUTPUT}/${GENOME_ID}_cog_hits.tsv"
N_HITS=$(wc -l < "${HITS_FILE}" || echo 0)
echo "DIAMOND hits: ${N_HITS}"

# Parse COG categories
if [ ${N_HITS} -gt 0 ]; then
    echo ""
    echo "Parsing COG categories..."

    /home/dmullane/micromamba/envs/esm3_env/bin/python3 << 'PYEOF'
import sys
from collections import Counter
from pathlib import Path

# Get variables from bash
n_proteins = ${N_PROTEINS}
sample_name = "${GENOME_ID}"
sample_output = "${SAMPLE_OUTPUT}"
hits_file = "${HITS_FILE}"
cog_db_dir = "${COG_DB}"

# Load protein-to-COG mapping
print("Loading COG mapping...")
prot_to_cog = {}
with open(f"{cog_db_dir}/cog-20.cog.csv", 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        fields = line.strip().split(',')
        if len(fields) >= 7:
            prot_id = fields[2]
            cog_id = fields[6]
            if prot_id not in prot_to_cog:
                prot_to_cog[prot_id] = cog_id

# Load COG definitions
print("Loading COG definitions...")
cog_definitions = {}
with open(f"{cog_db_dir}/cog-20.def.tab", 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) >= 3:
            cog_id = fields[0]
            cog_category = fields[1]
            cog_definitions[cog_id] = cog_category

# Map DIAMOND hits to COG IDs
protein_cogs = {}
with open(hits_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        query_id = fields[0]
        subject_id = fields[1]

        # Convert underscore to dot for matching
        if '_' in subject_id:
            parts = subject_id.rsplit('_', 1)
            subject_id_dot = f"{parts[0]}.{parts[1]}"
        else:
            subject_id_dot = subject_id

        if subject_id_dot in prot_to_cog:
            cog_id = prot_to_cog[subject_id_dot]
            if query_id not in protein_cogs:
                protein_cogs[query_id] = cog_id

# Count COG categories
category_counts = Counter()
proteins_with_cog = 0

for cog_id in protein_cogs.values():
    proteins_with_cog += 1
    if cog_id in cog_definitions:
        categories = cog_definitions[cog_id]
        for cat in categories:
            if cat != '-':
                category_counts[cat] += 1

# Write summary
with open(f"{sample_output}/{sample_name}_cog_summary.txt", 'w') as out:
    out.write(f"Sample: {sample_name}\n")
    out.write(f"Total proteins: {n_proteins}\n")
    out.write(f"Proteins with COG hits: {proteins_with_cog}\n")
    out.write(f"Annotation rate: {proteins_with_cog/n_proteins*100:.2f}%\n")
    out.write(f"\nCOG category counts:\n")
    for cat in sorted(category_counts.keys()):
        out.write(f"  {cat}: {category_counts[cat]}\n")

print(f"\nProteins with COG: {proteins_with_cog} ({proteins_with_cog/n_proteins*100:.1f}%)")
print(f"COG categories: {len(category_counts)}")
PYEOF

fi

echo ""
echo "=========================================="
if [ ${N_HITS} -gt 0 ]; then
    echo "✓ Successfully annotated ${GENOME_ID}"
else
    echo "⚠ No COG hits found for ${GENOME_ID}"
fi
echo "=========================================="

# Clean up temp directory
if [ -d "${TMPDIR}" ]; then
    rm -rf "${TMPDIR}"
fi

exit 0
