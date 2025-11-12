# Component 1: Genome to Graph Pipeline

This directory contains the complete pipeline for converting bacterial genomes into graph representations suitable for graph neural network training.

## Subcomponents

### 1.1 K-mer Profiling (Not yet implemented)
- Purpose: Generate k-mer profiles (5 or 6-mer counts) as genome node attributes
- Status: Pending implementation

### 1.2 Genome Parser
Location: `1.2_genome_parser/`
- Purpose: Parse protein-coding genes from bacterial genomes
- Tool: Prodigal
- Input: Raw genome sequences (FASTA)
- Output: Protein sequences with unique IDs

### 1.3 Multiple Sequence Alignment (MSA)
Location: `1.3_msa/`
- Purpose: Group similar sequences using sequence alignment to reduce redundancy
- Tool: MMseqs2
- Input: Protein sequences from 1.2
- Output: Clusters of homologous proteins
- Key scripts:
  - `cluster_proteins_mmseqs.py` - Run MMseqs2 clustering
  - `submit_mmseqs_parameter_sweep.sh` - Test multiple identity thresholds
  - `concatenate_all_proteins.py` - Prepare input for MMseqs2
  - `analyze_mmseqs_clusters.py` - Analyze cluster distributions

### 1.4 ESM Embedding + Clustering
Location: `1.4_esm_embedding_clustering/`
- Purpose: Generate ESM embeddings and cluster proteins by functional similarity
- Input: Protein sequences (filtered through MSA)
- Output: ESM embedding vectors and protein clusters
- Subdirectories:
  - `embedding_generation/` - ESM-C embedding pipeline (currently running)
  - `clustering/` - Leiden clustering and evaluation
  - `visualization/` - UMAP and cluster visualization
  - `functional_annotation/` - COG functional annotations

### 1.5 Graph Assembly (Not yet implemented)
- Purpose: Construct genomic graphs (gene-to-gene or genome-to-genome)
- Status: Pending implementation
