# Cluster Metadata

This directory contains metadata files for gene clusters in the protein graph.

## Files

### `cluster_metadata_simple.csv`
Basic cluster metadata with parsed cluster ID components.

**Columns:**
- `cluster_id`: Cluster representative ID (e.g., `GCF_000006985.1_NC_002932.3_1196`)
- `cluster_size`: Number of genes in cluster
- `representative_genome`: Genome accession (e.g., `GCF_000006985.1`)
- `representative_contig`: Contig/chromosome ID (e.g., `NC_002932.3`)
- `representative_gene`: Gene number from Prodigal (e.g., `1196`)

**Size:** ~2 MB
**Rows:** 27,654 clusters
**Genomes:** 5,815 unique genomes represented

### `cluster_metadata.csv`
Full cluster metadata including COG functional annotations (when available).

**Additional columns beyond simple metadata:**
- `top_cog_category`: Most common COG functional category in cluster (single letter)
- `top_cog_name`: Full name of top COG category
- `cog_counts`: JSON dict of all COG categories and their counts within cluster
- `annotation_rate`: Fraction of cluster members with COG annotation (0.0-1.0)

**Note:** COG columns will be empty unless diamond/eggNOG annotations have been run.

## Cluster Naming Convention

Cluster IDs follow the format:
```
{genome_accession}_{contig_id}_{gene_number}
```

**Examples:**
- `GCF_000006985.1_NC_002932.3_1196`
  - Genome: GCF_000006985.1 (Chlorobaculum tepidum TLS)
  - Contig: NC_002932.3
  - Gene: 1196 (gene prediction #1196 from Prodigal)

- `GCF_000008505.1_NC_005957.1_2572`
  - Genome: GCF_000008505.1
  - Contig: NC_005957.1
  - Gene: 2572

The cluster ID represents the **cluster representative** (seed gene chosen by MMseqs2 clustering). The cluster contains this gene plus all similar genes from across the entire dataset.

## Usage Examples

### Load metadata
```python
import pandas as pd

# Load simple metadata
metadata = pd.read_csv('graph_outputs/protein_graph/metadata/cluster_metadata_simple.csv')

# Load full metadata with COG
metadata_full = pd.read_csv('graph_outputs/protein_graph/metadata/cluster_metadata.csv')
```

### Merge with cluster embeddings
```python
# Merge with any cluster-level data
cluster_embeddings = pd.read_csv('cluster_embeddings.csv')
annotated = cluster_embeddings.merge(metadata, on='cluster_id', how='left')
```

### Filter by genome
```python
# Get all clusters from a specific genome
genome_clusters = metadata[metadata['representative_genome'] == 'GCF_000006985.1']
```

### Filter by COG category (if available)
```python
# Get translation/ribosomal clusters (COG category J)
translation_clusters = metadata_full[metadata_full['top_cog_category'] == 'J']

# Get energy production clusters (COG category C)
energy_clusters = metadata_full[metadata_full['top_cog_category'] == 'C']
```

## COG Functional Categories

When COG annotations are available, the `top_cog_category` uses single letter codes:

| Code | Category |
|------|----------|
| J | Translation, ribosomal structure and biogenesis |
| K | Transcription |
| L | Replication, recombination and repair |
| D | Cell cycle control, cell division |
| V | Defense mechanisms |
| T | Signal transduction mechanisms |
| M | Cell wall/membrane/envelope biogenesis |
| N | Cell motility |
| U | Intracellular trafficking, secretion |
| O | Posttranslational modification, protein turnover |
| C | Energy production and conversion |
| G | Carbohydrate transport and metabolism |
| E | Amino acid transport and metabolism |
| F | Nucleotide transport and metabolism |
| H | Coenzyme transport and metabolism |
| I | Lipid transport and metabolism |
| P | Inorganic ion transport and metabolism |
| Q | Secondary metabolites biosynthesis |
| R | General function prediction only |
| S | Function unknown |
| X | Mobilome (prophages, transposons) |

## Adding COG Annotations

To add COG functional annotations:

1. Run diamond against COG database:
```bash
sbatch 1_genome_to_graph/bin/1.4_run_diamond_cog_refseq.sh
```

2. Regenerate metadata with COG info:
```bash
python scripts/create_cluster_metadata.py \
  --cluster-index graph_outputs/protein_graph/cluster_index.csv \
  --clusters intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv \
  --cog-dir results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation \
  --output graph_outputs/protein_graph/metadata/cluster_metadata.csv
```

## Data Source

- **Cluster assignments:** MMseqs2 clustering of ESM-2 protein embeddings
- **Cluster sizes:** Number of genes per cluster (min=10 for protein_graph)
- **Representative selection:** Chosen by MMseqs2 based on sequence similarity
- **COG annotations:** Diamond BLASTp against COG-2020 database (when available)
- **Gene predictions:** Prodigal v2.6.3 gene calling

## Files by Graph Type

### protein_graph (stringent clustering)
- 27,654 clusters
- Minimum cluster size: 10 genes
- 5,815 unique genomes represented

### protein_graph_toy (permissive clustering)
- 365,211 clusters
- Used for full genome-gene graph visualization
- 7,659 unique genomes represented
