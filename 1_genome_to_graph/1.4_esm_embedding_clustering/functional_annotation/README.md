# COG Functional Annotation Workflow

Annotate gene embeddings with COG (Clusters of Orthologous Genes) functional categories using eggNOG-mapper, then visualize functional organization in UMAP space.

## Overview

This workflow assigns functional categories to genes based on homology to known protein families, then overlays these annotations onto UMAP embeddings to reveal functional organization across genomes.

**COG Categories** are 26 functional groups (J, K, L, D, V, T, M, N, etc.) representing biological processes like:
- **J**: Translation, ribosomal structure
- **K**: Transcription
- **C**: Energy production
- **E**: Amino acid metabolism
- **S**: Function unknown
- ... and 21 more

## Workflow Steps

```
RefSeq Proteins → eggNOG-mapper → Parse to Metadata → Merge with Any Analysis → Visualize
    (.faa)           (DIAMOND)      (per-genome TSV)    (UMAP, clusters, etc.)   (.png)
```

## Key Design: Reusable Metadata Files

Unlike traditional workflows that merge annotations with a specific analysis, this workflow creates **per-genome COG metadata files** that can be reused across any analysis:

- One TSV file per genome (e.g., `GCF_000005845_cog.tsv`)
- Indexed by gene_id for easy merging
- Works with any UMAP, clustering, differential expression, etc.
- No need to re-annotate for different analyses

## Files

- **`run_eggnog_refseq.sh`** - SLURM array job for eggNOG-mapper annotation
- **`parse_eggnog_to_metadata.py`** - Parse eggNOG results to per-genome metadata files
- **`load_cog_metadata.py`** - Utility functions to load and merge metadata
- **`plot_cog_umap.py`** - Visualize UMAP colored by COG categories
- ~~`annotate_umap_with_cogs.py`~~ - Deprecated (use metadata approach instead)

## Quick Start

### 1. Run eggNOG-mapper Annotation

Annotate all RefSeq proteins with COG categories:

```bash
cd /home/dmullane/SrivatsanLab/Dustin/organism_scale_modelling

# Submit SLURM array job (processes 100 genomes, 10 at a time)
sbatch scripts/functional_annotation/run_eggnog_refseq.sh
```

This will:
- Process each genome's protein sequences (`*_prodigal_proteins.faa`)
- Run DIAMOND search against eggNOG database
- Output COG annotations to `results/functional_annotation/{genome}_eggnog/`

**Expected time**: ~1-2 hours per genome with 8 CPUs

### 2. Monitor Progress

```bash
# Check job status
squeue -u $USER

# View logs
tail -f logs/eggnog_*.out

# Check how many genomes completed
ls results/functional_annotation/*/GC*.emapper.annotations | wc -l
```

### 3. Parse eggNOG Results to Metadata

After eggNOG-mapper completes, parse results into per-genome metadata files:

```bash
python scripts/functional_annotation/parse_eggnog_to_metadata.py \
    --annotations results/functional_annotation \
    --output results/metadata/cog_annotations
```

**Output**: Per-genome TSV files with columns:
- `gene_id` - Gene identifier
- `cog_category` - Single letter COG category (J, K, L, ...)
- `cog_categories_all` - All COG categories if multiple (e.g., "KL")
- `cog_description` - Gene functional description
- `category_name` - Full category description

**Files created**:
```
results/metadata/cog_annotations/
├── GCF_000005845_cog.tsv
├── GCF_000006945_cog.tsv
├── ... (one per genome)
└── annotation_summary.tsv  # Overall stats
```

### 4. Visualize UMAP with COGs

Create publication-quality plots (automatically loads and merges metadata):

```bash
python scripts/functional_annotation/plot_cog_umap.py \
    --umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample100000.npz \
    --cog-metadata results/metadata/cog_annotations \
    --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_cog_visualization.png
```

This creates:
- **Main plot**: UMAP colored by COG categories
- **Bar chart**: Category distribution
- **Summary stats**: Annotation rates
- **Individual plots**: Each major COG category highlighted separately

### 5. Use Metadata in Your Own Analyses

```python
# In a notebook or script
from load_cog_metadata import load_all_cog_metadata, merge_with_cog

# Load all COG metadata
cog_df = load_all_cog_metadata('results/metadata/cog_annotations')

# Merge with any analysis
my_analysis_df = pd.DataFrame({
    'gene_id': [...],
    'cluster': [...],
    'expression': [...]
})

annotated_df = merge_with_cog(my_analysis_df, cog_df)
# Now you have COG annotations for your analysis!
```

## Detailed Usage

### eggNOG-mapper Options

Edit `run_eggnog_refseq.sh` to customize:

```bash
# Number of genomes to process
#SBATCH --array=0-99%10    # Process 100 genomes, 10 concurrent

# Resources per job
#SBATCH --cpus-per-task=8  # CPU cores
#SBATCH --mem=32G          # Memory
#SBATCH --time=2:00:00     # Time limit

# Database location (using Sanjay's existing database)
EGGNOG_DB="/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/eggnog"
```

**Database size**: ~50-100GB (already downloaded, no need to redownload!)

### Parse to Metadata Options

```bash
python scripts/functional_annotation/parse_eggnog_to_metadata.py \
    --annotations <annotation_dir>    # Directory with eggNOG results
    --output <output_dir>             # Directory for metadata TSV files
```

This creates one TSV file per genome with COG annotations indexed by gene_id.

### Load Metadata (Python API)

```python
from load_cog_metadata import (
    load_all_cog_metadata,      # Load all genomes
    load_cog_for_genomes,       # Load specific genomes
    merge_with_cog,             # Merge with any DataFrame
    filter_by_cog_category,     # Filter to specific categories
    get_cog_summary             # Get category distribution
)

# Load all metadata
cog_df = load_all_cog_metadata('results/metadata/cog_annotations')

# Load specific genomes only
cog_df = load_cog_for_genomes(
    'results/metadata/cog_annotations',
    ['GCF_000005845', 'GCF_000006945']
)

# Merge with your data
annotated_df = merge_with_cog(my_df, cog_df, on='gene_id', how='left')

# Filter to specific categories
translation_genes = filter_by_cog_category(annotated_df, 'J')
metabolism_genes = filter_by_cog_category(annotated_df, ['C', 'E', 'G'])

# Get category distribution
summary = get_cog_summary(annotated_df)
```

### Visualization Options

```bash
# Option 1: Load UMAP and merge with metadata (recommended)
python scripts/functional_annotation/plot_cog_umap.py \
    --umap <umap_file.npz>          # UMAP from compute_umap_array.py
    --cog-metadata <metadata_dir>   # COG metadata directory
    --output <plot.png>             # Output plot file
    --max-points 100000             # Subsample for plotting (optional)
    --dpi 300                       # Plot resolution (default: 300)
    --point-size 0.5                # Point size (default: 0.5)
    --alpha 0.3                     # Transparency (default: 0.3)

# Option 2: Use pre-merged CSV (if you already merged)
python scripts/functional_annotation/plot_cog_umap.py \
    --annotated-csv <merged.csv>    # Pre-merged CSV file
    --output <plot.png>
```

## Expected Results

### Annotation Rates

Typical annotation rates for bacterial genomes:
- **Annotated**: 70-85% of genes
- **Unannotated**: 15-30% (novel or poorly characterized genes)

### COG Category Distribution

Most abundant categories (typical):
1. **R**: General function prediction (~15%)
2. **S**: Function unknown (~10%)
3. **E**: Amino acid metabolism (~8%)
4. **K**: Transcription (~7%)
5. **J**: Translation (~6%)

## Output Files

```
results/functional_annotation/          # eggNOG-mapper raw output
├── GCF_000005845_eggnog/
│   ├── GCF_000005845.emapper.annotations    # Main annotation file
│   ├── GCF_000005845.emapper.seed_orthologs
│   └── GCF_000005845.emapper.hits
├── GCF_000006945_eggnog/
│   └── ...
└── ... (one directory per genome)

results/metadata/cog_annotations/       # Parsed metadata (reusable!)
├── GCF_000005845_cog.tsv               # Per-genome COG metadata
├── GCF_000006945_cog.tsv
├── ... (one TSV per genome)
└── annotation_summary.tsv              # Overall statistics

results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/                            # Visualizations
├── umap_cog_visualization.png          # Main visualization
└── cog_individual/
    ├── umap_cog_J.png                  # Translation genes
    ├── umap_cog_K.png                  # Transcription genes
    ├── umap_cog_C.png                  # Energy metabolism genes
    └── ... (one plot per major category)
```

**Key Benefit**: The `results/metadata/cog_annotations/` directory contains reusable metadata files that can be merged with any future analysis (different UMAP, clustering, differential expression, etc.) without re-running eggNOG-mapper!

## COG Category Reference

### Information Storage and Processing
- **J**: Translation, ribosomal structure and biogenesis
- **A**: RNA processing and modification
- **K**: Transcription
- **L**: Replication, recombination and repair
- **B**: Chromatin structure and dynamics

### Cellular Processes and Signaling
- **D**: Cell cycle control, cell division
- **Y**: Nuclear structure
- **V**: Defense mechanisms
- **T**: Signal transduction mechanisms
- **M**: Cell wall/membrane/envelope biogenesis
- **N**: Cell motility
- **Z**: Cytoskeleton
- **W**: Extracellular structures
- **U**: Intracellular trafficking, secretion
- **O**: Posttranslational modification, protein turnover

### Metabolism
- **C**: Energy production and conversion
- **G**: Carbohydrate transport and metabolism
- **E**: Amino acid transport and metabolism
- **F**: Nucleotide transport and metabolism
- **H**: Coenzyme transport and metabolism
- **I**: Lipid transport and metabolism
- **P**: Inorganic ion transport and metabolism
- **Q**: Secondary metabolites biosynthesis

### Poorly Characterized
- **R**: General function prediction only
- **S**: Function unknown

### Mobile Elements
- **X**: Mobilome (prophages, transposons)

## Integration with UMAP Array Workflow

This COG annotation workflow is designed to work with UMAP results from the job array system:

```bash
# 1. Compute UMAP embeddings (from scripts/embeddings/)
sbatch scripts/embeddings/submit_umap_array.sh

# 2. Run eggNOG-mapper annotation (parallel with UMAP)
sbatch scripts/functional_annotation/run_eggnog_refseq.sh

# 3. Parse eggNOG results to reusable metadata files
python scripts/functional_annotation/parse_eggnog_to_metadata.py \
    --annotations results/functional_annotation \
    --output results/metadata/cog_annotations

# 4. Visualize any UMAP with COGs (no pre-merging needed!)
python scripts/functional_annotation/plot_cog_umap.py \
    --umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_subsample100000.npz \
    --cog-metadata results/metadata/cog_annotations \
    --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_cog_visualization.png

# 5. Try different UMAP parameters - just reuse the same metadata!
python scripts/functional_annotation/plot_cog_umap.py \
    --umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n50_subsample100000.npz \
    --cog-metadata results/metadata/cog_annotations \
    --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n50_cog_visualization.png
```

**Flexibility**: Because metadata is separate from UMAP, you can:
- Visualize any UMAP parameter sweep
- Combine with clustering results
- Annotate subsets of genes
- Use in differential expression analysis
- All without re-running eggNOG-mapper!

## Troubleshooting

### eggNOG-mapper fails with "Database not found"

```bash
# Check database path
ls /fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/eggnog/

# Should contain:
# - eggnog.db
# - eggnog_proteins.dmnd
# - og2level.txt.gz
# - etc.
```

If missing, ask Sanjay or download from http://eggnog5.embl.de/download/

### Annotation rate is very low (<50%)

Possible causes:
- Poor quality genome assembly
- Non-bacterial sequences
- Truncated protein sequences

Check input protein files:
```bash
head data/refseq_gene_annotations/GCF_*.faa
```

### Visualization runs out of memory

Subsample to fewer points:
```bash
python scripts/functional_annotation/plot_cog_umap.py \
    --annotated-umap results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15_with_cogs.csv \
    --output results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_cog_visualization.png \
    --max-points 50000  # Reduce from default
```

### UMAP and annotations don't match

Check that gene IDs match between UMAP and eggNOG files:
```bash
# From UMAP
python -c "import numpy as np; d=np.load('results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/umap_n15.npz', allow_pickle=True); print(d['gene_ids'][0])"

# From eggNOG
head -20 results/functional_annotation/*/GC*.emapper.annotations | grep -v "^#"
```

Gene IDs should be identical (e.g., `GCF_000005845_00001`, `GCF_000005845_00002`, etc.)

## Performance

### eggNOG-mapper
- **Speed**: ~1-2 hours per genome (8 CPUs)
- **Memory**: ~8-16GB per job
- **Total time**: ~10-20 hours for 100 genomes (parallel)

### Annotation merging
- **Speed**: ~1-2 minutes for 100k genes
- **Memory**: ~2-4GB

### Visualization
- **Speed**: ~30 seconds for 100k points
- **Memory**: ~4-8GB
- **Tip**: Use `--max-points` to subsample for faster plotting

## References

- **eggNOG 5.0**: http://eggnog5.embl.de/
- **eggNOG-mapper**: https://github.com/eggnogdb/eggnog-mapper
- **COG categories**: https://www.ncbi.nlm.nih.gov/COG/

## Citation

If you use this workflow, please cite:

- eggNOG-mapper: Huerta-Cepas et al. (2019) Mol Biol Evol. 36(10):2284-2289
- eggNOG 5.0: Huerta-Cepas et al. (2019) Nucleic Acids Res. 47(D1):D309-D314
