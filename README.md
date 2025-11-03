# Organism Scale Modelling

Pipeline for extracting genes from bacterial genomes and generating protein embeddings using ESM-C (Evolutionary Scale Modeling for proteins).

## Project Structure

```
.
├── data/
│   ├── refseq_genomes/          # Full bacterial genome sequences from refseq
│   └── refseq_gene_annotations/ # Annotations from refseq
├── scripts/
│   ├── gene_prediction/         # Tools for extracting genes from genomes
│   └── embeddings/              # Tools for generating ESM-C embeddings
├── docs/                        # Documentation
├── notebooks/                   # Jupyter notebooks
└── environment.yml             # Conda/micromamba environment
```

## Quick Start

### 1. Setup Environment

```bash
# Install environment
micromamba env create -f environment.yml
micromamba activate esm3_env

# Or update existing environment
micromamba env update -f environment.yml -n esm3_env
```

### 2. Test Installation

```bash
# Test gene prediction tools
python scripts/gene_prediction/test_gene_prediction.py

# Test ESM-C embedding generation
python scripts/embeddings/test_embeddings.py
```

### 3. Run Pipeline

**Option A: Use existing RefSeq annotations**

```bash
# Generate embeddings from pre-computed gene predictions
python scripts/embeddings/get_esm_embeddings.py \
    --gene-dir data/refseq_gene_annotations/ \
    --output-dir embeddings/
```

**Option B: Start from full genomes**

```bash
# Step 1: Predict genes from genomes
python scripts/gene_prediction/predict_genes.py \
    --genome-dir data/refseq_genomes/ \
    --output-dir predicted_genes/

# Step 2: Generate embeddings from predicted proteins
python scripts/embeddings/get_esm_embeddings.py \
    --gene-dir predicted_genes/ \
    --output-dir embeddings/
```

## Documentation

📚 **[Complete Documentation Index](docs/index.md)**

### Core Pipeline
- **[Pipeline Overview](docs/README.md)** - Complete RefSeq clustering pipeline (29M proteins)
- **[Gene Prediction Guide](docs/README_gene_prediction.md)** - Extract protein sequences from genomes
- **[ESM-C Embeddings Guide](docs/README_embeddings.md)** - Generate protein embeddings

### Advanced Analysis
- **[Evaluation Framework](docs/evaluation.md)** - Multi-metric clustering evaluation (ARI, AMI, stability)
- **[UMAP Array System](docs/umap_array.md)** - Parallelized dimensionality reduction
- **[Clustering Guide](docs/clustering.md)** - Usage examples and visualization

### Alternative Models
- **[gLM2 Integration](docs/README_glm2.md)** - Genomic language model comparison

## Scripts

### Gene Prediction (`scripts/gene_prediction/`)

| Script | Description |
|--------|-------------|
| `predict_genes.py` | Professional gene prediction using Prodigal |
| `extract_simple_orfs.py` | Lightweight ORF extraction (no dependencies) |
| `test_gene_prediction.py` | Test suite for gene prediction tools |

### ESM-C Embeddings (`scripts/embeddings/`)

| Script | Description |
|--------|-------------|
| `get_esm_embeddings.py` | Generate ESM-C embeddings for proteins |
| `test_embeddings.py` | Test ESM-C model and pipeline |
| `inspect_embeddings.py` | Analyze and summarize embedding results |

## Workflow Overview

```
┌─────────────────────┐
│ Bacterial Genome    │
│ (FASTA)             │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Gene Prediction     │
│ (Prodigal)          │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Protein Sequences   │
│ (.faa files)        │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ ESM-C Embedding     │
│ Generation          │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Gene Embeddings     │
│ (HDF5/NPZ)          │
└─────────────────────┘
```

## Data

- **RefSeq Genomes**: 7,664 bacterial genome assemblies
- **Gene Annotations**: Pre-computed Prodigal predictions for all RefSeq genomes
- Symlinks in `data/` point to shared storage

## Environment

Key packages:
- **Python 3.10**
- **ESM** - Protein language models
- **Prodigal/Pyrodigal** - Gene prediction
- **Biopython** - Sequence manipulation
- **PyTorch** - Deep learning (with CUDA support)
- **HDF5/NumPy** - Data storage

See `environment.yml` for full details.

## Hardware Requirements

- **CPU**: Any modern CPU (gene prediction)
- **GPU**: Highly recommended for embedding generation (CUDA-compatible)
  - ESM-C 600M model benefits significantly from GPU acceleration
- **Memory**: 16+ GB RAM recommended
- **Storage**:
  - Genomes: ~20 GB
  - Embeddings: ~1-5 GB (depends on number of genes)

## Output Formats

### Gene Prediction Output

For each genome:
- `{genome_id}_proteins.faa` - Amino acid sequences
- `{genome_id}_genes.fna` - Nucleotide sequences
- `{genome_id}_genes.gff` - Gene coordinates
- `{genome_id}_stats.txt` - Statistics

### Embedding Output

**HDF5 format** (default, single file):
```python
import h5py
with h5py.File("embeddings/esmc_embeddings.h5", "r") as f:
    embeddings = f[genome_id]["embeddings"][:]  # (n_genes, 960)
    gene_ids = f[genome_id]["gene_ids"][:]
    lengths = f[genome_id]["lengths"][:]
```

**NPZ format** (one file per genome):
```python
import numpy as np
data = np.load(f"embeddings/{genome_id}_embeddings.npz")
embeddings = data["embeddings"]  # (n_genes, 960)
```

## Tips

1. **Always test first**: Use `--num-genomes 10` to test on a small subset
2. **Use existing annotations**: RefSeq annotations are already high-quality
3. **GPU acceleration**: Run embeddings on a GPU node for 10-100x speedup
4. **Batch processing**: Adjust `--batch-size` based on available GPU memory
5. **Monitor progress**: All scripts use tqdm for progress bars

## Troubleshooting

**Import errors**: Make sure environment is activated
```bash
micromamba activate esm3_env
```

**Prodigal not found**: Install dependencies
```bash
micromamba env update -f environment.yml -n esm3_env
```

**CUDA out of memory**: Reduce batch size
```bash
python scripts/embeddings/get_esm_embeddings.py --batch-size 16
```

**Symlinks broken**: Check data directory paths
```bash
ls -la data/
```

## Citation

If you use this pipeline, please cite:

- **ESM**: [Evolutionary-scale prediction of atomic-level protein structure with a language model](https://doi.org/10.1126/science.ade2574)
- **Prodigal**: [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://doi.org/10.1186/1471-2105-11-119)

## Contact

See lab documentation for support and questions.
