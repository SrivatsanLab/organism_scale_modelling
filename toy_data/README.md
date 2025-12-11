# Overview:
This directory contains the files necessary to assemble the full heterogeneous graph represention of **NCBI Refseq Prokaryotic**. Raw `.fasta` files for each of the 7,664 high quality prokaryotic genome assemblies were downloaded from refseq, and processed on the Fred Hutch HPC cluster with the scripts in `preprocessing/` (scripts built for any SLURM system). These data are the most compressed versions, with some information loss (e.g. 1152 dimensional ESM embeddings for gene node attributes are PCA reduced to 50 dims), intended to fit on `github`. 

### `genome_node_index.csv`:
Genome index mapping matrix rows to genome IDs
- Columns:
  - `genome_index`: Row index in the matrix (0-7663)
  - `genome_id`: RefSeq genome identifier
  - `genus`: Taxonomic genus
  - `species`: Species name
  - `full_name`: Genus_species combined
- 7,664 genomes across 2,140 unique genera

### `genome_node_attributes.npz`:
Compressed genome × k-mer matrix
- Format: Scipy CSR sparse matrix (float32)
- Shape: 7,664 genomes × 2,080 canonical 6-mers
- Load with: `scipy.sparse.load_npz()`

### `genome_node_attribute_index.csv`:
K-mer index mapping matrix columns to canonical 6-mer sequences
- Columns:
  - `kmer_index`: Column index in the matrix (0-2079)
  - `kmer`: 6-mer sequence (canonical representation)
- 2,080 canonical 6-mers (out of 4^6 = 4,096 possible, considering reverse complements)

### `genome_genome_edges.npz`:
Genome-genome adjacency matrix
- Format: Scipy sparse matrix (int64, undirected)
- Shape: 7,664 × 7,664
- Edges: 50,303 (undirected)
- Average degree: 13.1 neighbors per genome
- Construction method:
  - Distance metric: **Cosine similarity** in full 6-mer space (2080D)
  - k-neighbors: **k=10** (optimized for genus-level structure)
  - Graph is symmetrized (undirected) with self-loops removed
- Load with: `scipy.sparse.load_npz('genome_genome_edges.npz')`

### `genome_distance_matrix.npz`: 
Pairwise cosine distance matrix (genome × genome)
- Format: Custom compressed format (upper triangle only, float16)
- Full matrix: 7,664 × 7,664 (symmetric with zeros on diagonal)
- Compression: 4.13x smaller than full matrix
- Precision: Max error 0.000244, mean error 0.000060
- Distance range: [0.0, 0.992], Mean: 0.349, Median: 0.270

### `gene_gene_edges.npz`:
Gene-gene (cluster-cluster) edge matrix from operon predictions
- Format: Numpy arrays (sparse matrix components)
- Shape: (365,211, 365,211)
- Arrays: `data` (edge weights, float32), `indices` (column indices), `indptr` (row pointers), `shape`, `cluster_ids`
- Edges: Gene clusters predicted to be in operons (adjacent on genome, same strand)
- Load with: `np.load('gene_gene_edges.npz')`

### `gene_node_attributes.npz`:
Gene Node attributes
- Format: Numpy arrays
- Shape: (365,211 clusters × 50 PCs)
- Arrays: `embeddings` (float32), `cluster_ids`
- PCA reduced ESM-2 600M protein language model embeddings (top 50 pcs, mean-pooled per cluster)
- Load with: `data = np.load('gene_node_attributes.npz'); embeddings = data['embeddings']`

### `gene_node_metadata.csv`:
Gene cluster functional annotations from COG-2020 database
- Columns:
  - `cluster_id`: Gene cluster identifier
  - `cluster_size`: Number of proteins in cluster
  - `top_cog_category`: Most common COG category (single letter code, e.g., 'J', 'E')
  - `top_cog_name`: Full name of top COG category (e.g., 'Translation, ribosomal structure')
  - `cog_counts`: JSON dict of all COG categories in cluster with counts
  - `annotation_rate`: Fraction of cluster members with COG annotations (0.0-1.0)
  - `representative_genome`, `representative_contig`, `representative_gene`: Cluster representative
- 365,211 clusters with functional annotations
- Annotation method: Diamond BLASTp against COG-2020 database
- Top COG categories:
  - E (Amino acid transport/metabolism): 22,156 clusters (6.1%)
  - K (Transcription): 21,606 clusters (5.9%)
  - S (Function unknown): 19,868 clusters (5.4%)
- Load with: `pd.read_csv()`

### `gene_node_index.csv`:
Gene Node Index 
- Load with: `pd.read_csv()`

### `genome_gene_edges.npz`:
Genome-gene bipartite graph edge matrix
- Format: Numpy arrays (sparse matrix components)
- Shape: (7,664, 365,211)
- Arrays: `data` (presence indicator, int8), `indices` (column indices), `indptr` (row pointers), `shape`, `genome_ids`, `cluster$
- File size: 35 MB
- Edges: Which gene clusters are present in which genomes
- Load with: `np.load('genome_gene_edges.npz')`