#!/bin/bash
#SBATCH --job-name=eggnog
#SBATCH --output=logs/eggnog_%A_%a.out
#SBATCH --error=logs/eggnog_%A_%a.err
#SBATCH --array=0-7663%50
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=campus-new

# Run eggNOG-mapper on RefSeq gene annotations
# This annotates proteins with COG functional categories

# Configuration
PROJECT_DIR="/home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling"
GENE_DIR="${PROJECT_DIR}/data/refseq_gene_annotations"
OUTPUT_DIR="${PROJECT_DIR}/results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation"
EGGNOG_DB="/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/eggnog"

# Use existing eggNOG database
export EGGNOG_DATA_DIR="${EGGNOG_DB}"

# Set temp directory to /tmp instead of project directory
# This prevents eggNOG from creating emappertmp_* directories in the project
export TMPDIR="/tmp/${USER}/eggnog_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# Load eggNOG-mapper module
module purge
module load eggnog-mapper/2.1.7-foss-2021b

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
echo "eggNOG-mapper Annotation"
echo "Task: ${SLURM_ARRAY_TASK_ID}/${TOTAL_FILES}"
echo "Genome: ${GENOME_ID}"
echo "=========================================="

# Output directory for this genome
SAMPLE_OUTPUT="${OUTPUT_DIR}/${GENOME_ID}_eggnog"
mkdir -p "${SAMPLE_OUTPUT}"

# Check if already processed
if [ -f "${SAMPLE_OUTPUT}/${GENOME_ID}.emapper.annotations" ]; then
    echo "Already processed: ${GENOME_ID}"
    echo "Delete ${SAMPLE_OUTPUT} to reprocess"
    exit 0
fi

# Run eggNOG-mapper
echo ""
echo "Running eggNOG-mapper..."
echo "Input: ${PROTEIN_FILE}"
echo "Output: ${SAMPLE_OUTPUT}"
echo ""

emapper.py \
    -i "${PROTEIN_FILE}" \
    --output "${GENOME_ID}" \
    --output_dir "${SAMPLE_OUTPUT}" \
    --data_dir "${EGGNOG_DB}" \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -m diamond \
    --dmnd_db "${EGGNOG_DB}/eggnog_proteins.dmnd" \
    --override

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Exit code: ${EXIT_CODE}"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Successfully annotated ${GENOME_ID}"

    # Print summary
    if [ -f "${SAMPLE_OUTPUT}/${GENOME_ID}.emapper.annotations" ]; then
        TOTAL=$(grep -v '^#' "${SAMPLE_OUTPUT}/${GENOME_ID}.emapper.annotations" | wc -l)
        ANNOTATED=$(grep -v '^#' "${SAMPLE_OUTPUT}/${GENOME_ID}.emapper.annotations" | awk -F'\t' '$7 != "-"' | wc -l)
        echo "  Total proteins: ${TOTAL}"
        echo "  Annotated: ${ANNOTATED}"
        if [ ${TOTAL} -gt 0 ]; then
            PCT=$((100 * ANNOTATED / TOTAL))
            echo "  Annotation rate: ${PCT}%"
        fi
    fi
else
    echo "✗ Error annotating ${GENOME_ID}"
fi

echo "=========================================="

# Clean up temp directory
if [ -d "${TMPDIR}" ]; then
    rm -rf "${TMPDIR}"
    echo "Cleaned up temp directory: ${TMPDIR}"
fi

exit ${EXIT_CODE}
