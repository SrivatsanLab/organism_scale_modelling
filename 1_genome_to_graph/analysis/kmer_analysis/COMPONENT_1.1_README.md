# 5-mer Profiling Pipeline

Generate 5-mer frequency profiles for all 7,664 bacterial genomes using Jellyfish.

## Overview

This pipeline:
1. Counts 5-mers in each genome using Jellyfish (fast, memory-efficient)
2. Generates individual k-mer profiles (text files)
3. Combines profiles into a single matrix for downstream analysis

## Quick Start

### 1. Test on Single Genome

```bash
# Test the script on one genome
bash scripts/kmer_profiling/generate_5mer_profiles.sh \
    data/refseq_genomes/GCF_000006985.1_Chlorobaculum_tepidum_TLS.fasta \
    results/kmer_profiles/test/

# Check output
ls results/kmer_profiles/test/
cat results/kmer_profiles/test/*_stats.txt
```

### 2. Run on All Genomes (SLURM Array)

```bash
# Create logs directory
mkdir -p logs

# Submit array job for all 7,664 genomes
sbatch scripts/kmer_profiling/submit_5mer_array.sh

# Monitor progress
squeue -u $USER
watch "ls results/kmer_profiles/5mer/*.txt | wc -l"

# Check for errors
grep -i error logs/5mer_profile_*.err
```

**Runtime**: ~30 min total (parallel processing)
- Each genome: ~2-3 minutes
- 7,664 genomes run in parallel (array job)
- Actual walltime: depends on cluster queue

### 3. Create Combined Matrix

After all individual profiles complete:

```bash
# Create sparse matrix (recommended)
python scripts/kmer_profiling/create_kmer_matrix.py \
    --input-dir results/kmer_profiles/5mer \
    --output results/kmer_profiles/5mer_matrix.npz \
    --normalize frequency

# Or create CSV (only for small subsets!)
python scripts/kmer_profiling/create_kmer_matrix.py \
    --input-dir results/kmer_profiles/test \
    --output results/kmer_profiles/test_matrix.csv \
    --format csv \
    --normalize frequency
```

## Pipeline Scripts

### `generate_5mer_profiles.sh`

Counts 5-mers for a single genome.

**Key parameters**:
- `-m 5`: 5-mer length
- `-C`: Canonical k-mers (count both strands)
- `-s 100M`: Hash size (sufficient for bacteria)
- `-t 4`: 4 CPU threads

**Output per genome**:
- `{genome_id}_5mer.txt`: K-mer counts (tab-separated)
- `{genome_id}_5mer_stats.txt`: Statistics summary

### `submit_5mer_array.sh`

SLURM array job wrapper.

**SLURM parameters**:
- `--array=1-7664`: One task per genome
- `--cpus-per-task=4`: Parallelization per genome
- `--mem=8G`: Sufficient for bacterial genomes
- `--time=00:30:00`: Conservative time limit

**Features**:
- Skip completed genomes (resume-friendly)
- Individual logs per genome
- Error checking

### `create_kmer_matrix.py`

Combines individual profiles into a matrix.

**Options**:
- `--normalize frequency`: K-mer frequencies (default)
- `--normalize presence`: Binary (present/absent)
- `--normalize none`: Raw counts
- `--format npz`: Sparse matrix (efficient)
- `--format csv`: Dense matrix (for small datasets)

## Output Files

### Individual Profiles

```
results/kmer_profiles/5mer/
├── GCF_000006985.1_Chlorobaculum_tepidum_TLS_5mer.txt
├── GCF_000006985.1_Chlorobaculum_tepidum_TLS_5mer_stats.txt
├── GCF_000007085.1_..._5mer.txt
├── GCF_000007085.1_..._5mer_stats.txt
...
```

**5mer.txt format**:
```
AAAAA   1234
AAAAC   567
AAAAG   890
...
```

### Combined Matrix

```
results/kmer_profiles/
├── 5mer_matrix.npz         # Sparse matrix (7664 × 1024)
└── 5mer_matrix.meta.npz    # Metadata (genome IDs, k-mer list)
```

## Using the K-mer Matrix

### Load in Python

```python
import numpy as np
from scipy.sparse import load_npz

# Load matrix and metadata
matrix = load_npz('results/kmer_profiles/5mer_matrix.npz')
meta = np.load('results/kmer_profiles/5mer_matrix.meta.npz', allow_pickle=True)

genome_ids = meta['genome_ids']
kmer_list = meta['kmer_list']

print(f"Matrix shape: {matrix.shape}")
print(f"Genomes: {len(genome_ids)}")
print(f"K-mers: {len(kmer_list)}")

# Convert to dense array (if needed)
matrix_dense = matrix.toarray()

# Get k-mer profile for first genome
genome_0_profile = matrix[0, :].toarray().flatten()
```

### Downstream Analysis Ideas

1. **Phylogenetic clustering**: Cluster genomes by k-mer similarity
2. **Dimensionality reduction**: PCA/UMAP on k-mer profiles
3. **Genome comparison**: Cosine similarity between k-mer vectors
4. **GC content**: Compute from k-mer frequencies
5. **Sequence motifs**: Identify enriched k-mers

## Expected Results

### K-mer Statistics

For 5-mers:
- **Possible k-mers**: 4^5 = 1,024 (canonical: 512)
- **Observed per genome**: ~200-512 (depends on genome size/diversity)
- **Most frequent**: AT-rich or GC-rich motifs
- **Rare k-mers**: May indicate horizontal gene transfer or unique regions

### Matrix Size

- **Dimensions**: 7,664 genomes × 1,024 k-mers
- **Sparse format**: ~10-50 MB (most k-mers absent in most genomes)
- **Dense format**: ~32 MB (7,664 × 1,024 × 4 bytes)

## Troubleshooting

### Job Failures

```bash
# Check which genomes failed
for i in {1..7664}; do
    genome=$(ls data/refseq_genomes/*.fasta | sed -n "${i}p")
    genome_id=$(basename "$genome" .fasta)
    if [ ! -f "results/kmer_profiles/5mer/${genome_id}_5mer.txt" ]; then
        echo "Missing: $genome_id (task $i)"
    fi
done

# Resubmit specific array tasks
sbatch --array=42,100,523 scripts/kmer_profiling/submit_5mer_array.sh
```

### Memory Issues

If jobs run out of memory (unlikely for bacteria):
```bash
# Increase memory allocation
#SBATCH --mem=16G
```

### Jellyfish Not Found

```bash
# Load module manually
module load Jellyfish/2.3.0-GCC-11.2.0
module list
```

## Performance

### Individual Genome

- **Small genome** (2 Mbp): ~30 seconds
- **Medium genome** (4 Mbp): ~1-2 minutes
- **Large genome** (8 Mbp): ~3-5 minutes

### Full Dataset

- **Serial processing**: ~50-150 hours
- **Parallel (array)**: ~30 min (walltime)
- **Matrix creation**: ~5-10 minutes

## Advanced Usage

### Custom K-mer Length

To use different k-mer lengths (e.g., 4-mers or 6-mers):

```bash
# Edit generate_5mer_profiles.sh
# Change: -m 5  to  -m 6

# Then run normally
sbatch scripts/kmer_profiling/submit_5mer_array.sh
```

### Non-canonical K-mers

To count forward and reverse separately (remove `-C` flag):

```bash
# In generate_5mer_profiles.sh, remove the -C flag
jellyfish count -m 5 -s 100M -t 4 -o "$JF_FILE" "$GENOME_FILE"
```

This doubles the k-mer space: 1,024 → 2,048 for 5-mers.

## Integration with Pipeline

The k-mer matrix can complement protein embeddings:

```bash
# 1. Generate k-mer profiles (genome-level)
sbatch scripts/kmer_profiling/submit_5mer_array.sh

# 2. Generate protein embeddings (protein-level)
sbatch scripts/embeddings/generate_esm_embeddings.sh

# 3. Combine both for multi-scale analysis
python scripts/analysis/combine_genome_protein_features.py \
    --kmer-matrix results/kmer_profiles/5mer_matrix.npz \
    --protein-embeddings data/refseq_esm_embeddings/ \
    --output results/combined_features/
```

## References

- **Jellyfish**: Marçais & Kingsford (2011) "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers"
- **K-mer analysis**: Compeau et al. (2011) "How to apply de Bruijn graphs to genome assembly"

## Contact

For issues or questions, check SLURM logs in `logs/` directory.
