# Masked Graph Learning for Bacterial Genomes

Project: Masked graph learning for bacterial genomics
Goals:
1. Define gene nodes (mostly done)
2. Represent genomes as graphs
3. Implement [GraphMAE](https://arxiv.org/abs/2205.10803)

![Proposed Graph Structure](graph_concept.png)
 
The data:
(https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/)
 
For defining gene nodes:
[ESM Atlas](https://esmatlas.com/)
[ESM Cambrian](https://www.evolutionaryscale.ai/blog/esm-cambrian)

## Current Status (2025-11-12)

### Component 1: Genome to Graph - ✅ Major Progress

**Completed:**
- ✅ Generated ESM-C embeddings for **11.8M proteins** (1,184 GPU batches)
- ✅ Created PCA cache: 50D embeddings with 89.2% variance explained (4.3GB)
- ✅ MMseqs2 clustering: 30M proteins clustered at 70% sequence identity
- ✅ Optimized storage: Moved 17GB to cheaper `/working` storage

**In Progress:**
- 🔄 UMAP computation for full 11.8M protein dataset
- 🔄 Analyzing cluster tightness in embedding space to determine if MMseqs2 clusters are sufficient

**Next:**
- Decide on clustering strategy (MMseqs2 vs Leiden) based on embedding space analysis
- Generate comprehensive visualizations
- Create gene node definitions for graph construction

See [EMBEDDING_PIPELINE.md](docs/EMBEDDING_PIPELINE.md) for detailed progress.

## Project Structure

```
.
├── 1_genome_to_graph/
│   ├── 1.3_msa/                          # MMseqs2 clustering results
│   └── 1.4_esm_embedding_clustering/     # ESM embedding generation & analysis
├── data/
│   ├── all_proteins.faa                  # 30M proteins from 7,664 genomes
│   ├── refseq_genomes/                   # Genome sequences (symlink)
│   └── refseq_gene_annotations/          # Gene annotations (symlink)
├── results/
│   └── 1_genome_to_graph/
│       ├── 1.3_msa/mmseqs_seqid_0p7/    # 70% ID clusters (30M assignments)
│       └── 1.4_esm_embedding_clustering/ # PCA cache, UMAP, analysis
├── docs/                                 # Documentation
└── environment.yml                       # Conda/micromamba environment
```
