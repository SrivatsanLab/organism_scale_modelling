# Component 1: Genome to Graph Pipeline

This directory contains the complete pipeline for converting bacterial genomes into graph representations suitable for graph neural network training.

## Directory Structure

```
1_genome_to_graph/
├── scripts/              # Essential pipeline scripts (production)
├── bin/                  # Exploratory and analysis scripts
├── graph_outputs/        # Final graph data (git-tracked)
│   ├── genome_graph/    # Genome-level graph (k-mer based)
│   └── protein_graph/   # Protein-level graph (ESM embedding based)
├── analysis/            # Analysis results and documentation
│   ├── kmer_analysis/   # K-mer profiling analysis and plots
│   └── protein_clustering/  # Protein clustering analysis
└── intermediate/        # Intermediate data (NOT git-tracked)
    ├── kmer/           # K-mer profiles and matrices
    └── protein/        # MSA clusters and ESM embeddings
```

## Pipeline Components

### Component 1.1: K-mer Profiling
**Purpose**: Generate k-mer frequency profiles from bacterial genomes and construct genome similarity graphs

**Input**: Bacterial genome FASTA files (from `data/genomes/`)

**Output**:
- Genome × k-mer frequency matrix (compressed)
- k-NN genome similarity graph (adjacency matrix)
- Pairwise cosine distance matrix
- Genome and k-mer indices with taxonomy

**Scripts** (in `scripts/`):
- `1.1_generate_5mer_profiles.sh` - Generate 5-mer profiles using Jellyfish
- `1.1_generate_6mer_profiles.sh` - Generate 6-mer profiles using Jellyfish
- `1.1_submit_5mer_array.sh` - SLURM array job for 5-mer profiling
- `1.1_submit_6mer_array.sh` - SLURM array job for 6-mer profiling
- `1.1_create_kmer_matrix.py` - Consolidate profiles into genome × k-mer matrix
  - Input: Individual k-mer profile files from `intermediate/kmer/{5,6}mer/`
  - Output: Compressed k-mer matrix to `intermediate/kmer/`
- `1.1_create_final_outputs.py` - Generate final graph outputs
  - Input: K-mer matrix from previous step
  - Output: Adjacency matrix, indices, metadata to `graph_outputs/genome_graph/`
- `1.1_create_distance_matrix.py` - Compute pairwise genome distances
  - Input: K-mer matrix
  - Output: Full distance matrix to `intermediate/kmer/`
- `1.1_compress_distance_matrix.py` - Compress distance matrix for storage
  - Input: Full distance matrix
  - Output: Compressed distance matrix to `graph_outputs/genome_graph/`

**Key Parameters**:
- K-mer size: 6 (selected based on phylogenetic signal)
- k-NN graph: k=10 neighbors (optimized for genus-level purity: 38.4%)
- Distance metric: Cosine distance

---

### Component 1.2: Genome Parser
**Purpose**: Predict protein-coding genes from bacterial genome sequences

**Input**: Bacterial genome FASTA files (from `data/genomes/`)

**Output**:
- Predicted protein sequences (amino acid FASTA)
- Gene annotations (GFF format)
- Gene nucleotide sequences

**Scripts** (in `scripts/`):
- `1.2_predict_genes.py` - Main gene prediction script
  - Tool: Prodigal (bacterial gene caller)
  - Input: Genome FASTA file(s)
  - Output: Protein sequences to `intermediate/protein/predicted_genes/`
  - Usage: `python 1.2_predict_genes.py --genome-dir data/genomes/ --output-dir intermediate/protein/predicted_genes/`

---

### Component 1.3: Multiple Sequence Alignment (MSA)
**Purpose**: Cluster homologous proteins using sequence similarity to reduce redundancy

**Input**: Protein sequences from Component 1.2

**Output**:
- Protein cluster assignments (TSV)
- Cluster statistics and analysis
- Representative sequences per cluster

**Scripts** (in `scripts/`):
- `1.3_concatenate_all_proteins.py` - Combine all predicted proteins into single FASTA
  - Input: Individual protein FASTA files from `intermediate/protein/predicted_genes/`
  - Output: Concatenated FASTA to `intermediate/protein/all_proteins.faa`
- `1.3_cluster_proteins_mmseqs.py` - Run MMseqs2 clustering
  - Tool: MMseqs2
  - Input: Concatenated protein FASTA
  - Output: Cluster TSV to `intermediate/protein/msa_clusters/`
  - Parameters: Sequence identity threshold (default: 0.5), coverage (default: 0.8)
- `1.3_submit_mmseqs_full_dataset.sh` - SLURM submission for full dataset clustering
  - Launches MMseqs2 clustering job on HPC

**Key Parameters**:
- Minimum sequence identity: 50% (determined by parameter sweep)
- Minimum coverage: 80%
- Clustering mode: Greedy set cover

---

### Component 1.4: ESM Embedding + Clustering
**Purpose**: Generate ESM protein embeddings and cluster proteins by functional/structural similarity

**Input**: Protein sequences (filtered through MSA clusters)

**Output**:
- ESM-2 embeddings (1152-dimensional vectors)
- UMAP-reduced embeddings (2D or 50D)
- Protein cluster assignments
- Functional annotations (COG categories)

**Scripts** (in `scripts/`):
- `1.4_compute_gpu_embeddings.py` - Generate ESM embeddings using GPU
  - Model: ESM-2 (650M parameters)
  - Input: Protein sequences from MSA clusters
  - Output: Embedding NPZ files to `intermediate/protein/esm_embeddings/`
  - Usage: Run via SLURM GPU nodes
- `1.4_submit_gpu_array.sh` - SLURM GPU array job for embedding generation
  - Parallelizes embedding computation across multiple GPUs
- `1.4_compute_umap_array.py` - Compute UMAP dimensionality reduction
  - Input: ESM embeddings from previous step
  - Output: UMAP coordinates to `intermediate/protein/esm_embeddings/umap/`
  - Dimensions: 2D for visualization, 50D for clustering
- `1.4_submit_umap_full_array.sh` - SLURM submission for UMAP computation
  - Uses RAPIDS cuML for GPU-accelerated UMAP

**Subdirectories in `intermediate/protein/esm_embeddings/`**:
- `embedding_generation/` - Raw ESM-2 embeddings
- `clustering/` - Leiden clustering results
- `umap/` - UMAP-reduced embeddings
- `functional_annotation/` - COG/eggNOG annotations

**Key Parameters**:
- ESM model: ESM-2 (esm2_t33_650M_UR50D)
- Embedding dimension: 1152
- UMAP neighbors: 15
- UMAP min_dist: 0.1

---

### Component 1.5: Graph Assembly
**Status**: Pending implementation

**Purpose**: Construct final graph structures
- Genome-to-genome graph (using k-mer similarity)
- Protein-to-protein graph (using ESM embedding similarity)
- Combined multi-level graph

---

## Running the Pipeline

### Prerequisites
- Prodigal (gene prediction)
- Jellyfish (k-mer counting)
- MMseqs2 (sequence clustering)
- PyTorch + ESM (protein embeddings)
- SLURM cluster access (for large-scale processing)

### Execution Order

1. **K-mer Profiling** (Component 1.1)
   ```bash
   # Generate k-mer profiles
   sbatch scripts/1.1_submit_6mer_array.sh

   # Create k-mer matrix
   python scripts/1.1_create_kmer_matrix.py --kmer-size 6 --input-dir intermediate/kmer/6mer/

   # Generate final outputs
   python scripts/1.1_create_final_outputs.py --input intermediate/kmer/6mer_matrix.npz

   # Create distance matrices
   python scripts/1.1_create_distance_matrix.py
   python scripts/1.1_compress_distance_matrix.py
   ```

2. **Gene Prediction** (Component 1.2)
   ```bash
   python scripts/1.2_predict_genes.py --genome-dir data/genomes/ --output-dir intermediate/protein/predicted_genes/
   ```

3. **Protein Clustering** (Component 1.3)
   ```bash
   # Concatenate proteins
   python scripts/1.3_concatenate_all_proteins.py

   # Run MMseqs2 clustering
   sbatch scripts/1.3_submit_mmseqs_full_dataset.sh
   ```

4. **ESM Embeddings** (Component 1.4)
   ```bash
   # Generate embeddings
   sbatch scripts/1.4_submit_gpu_array.sh

   # Compute UMAP
   sbatch scripts/1.4_submit_umap_full_array.sh
   ```

### Output Files

Final outputs are in `graph_outputs/`:

**Genome Graph** (`graph_outputs/genome_graph/`):
- `6mer_matrix_compressed.npz` - Genome × k-mer frequency matrix (7,664 × 2,080)
- `adjacency_matrix.npz` - k-NN genome similarity graph (sparse, k=10)
- `distance_matrix_compressed.npz` - Pairwise genome distances (compressed)
- `genome_index.csv` - Genome IDs with taxonomy (genus, species)
- `kmer_index.csv` - K-mer sequences (canonical form)
- `metadata.txt` - Pipeline parameters and metrics
- `README.md` - Detailed output documentation

**Protein Graph** (`graph_outputs/protein_graph/`):
- Status: Pending (Component 1.5)

---

## Exploratory Analysis

The `bin/` directory contains 83+ analysis and visualization scripts:
- K-mer phylogenetic analysis
- MMseqs parameter testing
- ESM embedding visualization
- Clustering quality evaluation
- COG functional annotation
- UMAP visualizations

All scripts are prefixed by component (1.1_, 1.2_, 1.3_, 1.4_) for easy identification.

---

## Documentation

- Component 1.1: `analysis/kmer_analysis/COMPONENT_1.1_README.md`
- Component 1.4: `analysis/protein_clustering/docs/`
  - `CLUSTERING_README.md`
  - `EMBEDDING_README.md`
  - `ANNOTATION_README.md`

---

## Current Status (November 2024)

- ✅ Component 1.1: Complete (7,664 bacterial genomes processed)
- ✅ Component 1.2: Complete (gene prediction functional)
- ✅ Component 1.3: Complete (MMseqs2 clustering optimized)
- 🔄 Component 1.4: In progress (embeddings generated, clustering ongoing)
- ⏳ Component 1.5: Pending
