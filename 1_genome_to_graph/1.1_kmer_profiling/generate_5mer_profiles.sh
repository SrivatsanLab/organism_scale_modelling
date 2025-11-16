#!/bin/bash
#
# generate_5mer_profiles.sh
# Generate 5-mer frequency profiles for bacterial genomes using Jellyfish
#
# Usage: ./generate_5mer_profiles.sh <genome.fasta> <output_dir>
#

set -e  # Exit on error

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome.fasta> <output_dir>"
    echo "Example: $0 genome.fasta results/kmer_profiles/"
    exit 1
fi

GENOME_FILE="$1"
OUTPUT_DIR="$2"

# Extract genome ID from filename
GENOME_ID=$(basename "$GENOME_FILE" | sed 's/\.fasta$//' | sed 's/\.fa$//' | sed 's/\.fna$//')

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Output files
JF_FILE="${OUTPUT_DIR}/${GENOME_ID}_5mer.jf"
DUMP_FILE="${OUTPUT_DIR}/${GENOME_ID}_5mer.txt"
STATS_FILE="${OUTPUT_DIR}/${GENOME_ID}_5mer_stats.txt"

echo "Processing: $GENOME_ID"
echo "Input: $GENOME_FILE"
echo "Output directory: $OUTPUT_DIR"

# Load Jellyfish module
module load Jellyfish/2.3.0-GCC-11.2.0

# Step 1: Count 5-mers
# - Use canonical k-mers (count both strands together)
# - Initial hash size: 100M (sufficient for bacterial genomes)
# - For 5-mers, there are 4^5 = 1024 possible k-mers
echo "Counting 5-mers..."
jellyfish count \
    -m 5 \
    -s 100M \
    -t 4 \
    -C \
    -o "$JF_FILE" \
    "$GENOME_FILE"

# Step 2: Export counts to text format
echo "Exporting k-mer counts..."
jellyfish dump \
    -c \
    -t \
    -o "$DUMP_FILE" \
    "$JF_FILE"

# Step 3: Generate statistics
echo "Generating statistics..."
{
    echo "Genome: $GENOME_ID"
    echo "Date: $(date)"
    echo "---"
    echo "Total unique 5-mers:"
    wc -l < "$DUMP_FILE"
    echo ""
    echo "Total 5-mer observations:"
    awk '{sum+=$2} END {print sum}' "$DUMP_FILE"
    echo ""
    echo "Most frequent 5-mers (top 10):"
    sort -k2,2nr "$DUMP_FILE" | head -10
    echo ""
    echo "Least frequent 5-mers (bottom 10):"
    sort -k2,2n "$DUMP_FILE" | head -10
} > "$STATS_FILE"

# Clean up binary file (keep text dump)
rm "$JF_FILE"

echo "Complete! Results in:"
echo "  Counts: $DUMP_FILE"
echo "  Stats:  $STATS_FILE"
