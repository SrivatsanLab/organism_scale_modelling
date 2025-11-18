# ESM Embedding Generation

This directory contains scripts for generating ESM-C embeddings for all proteins in the dataset.

## Current Pipeline (Active)

**Job ID**: 41725169 (1,184 batch jobs)
**Status**: Running
**Purpose**: Generate embeddings for 11.8M proteins in MMseqs2-filtered clusters

### Key Scripts

- `batch_generate_embeddings.py` - Generate embeddings for a batch of proteins
- `submit_batch_embeddings.sh` - SLURM array job for batch processing
- `identify_proteins_needing_embeddings.py` - Find proteins without embeddings
- `extract_sequences_for_embedding.py` - Extract FASTA sequences for target proteins
- `merge_embedding_batches.py` - Merge batch results into single cache

### Legacy Scripts (Original ESM Pipeline)

- `compute_gpu_embeddings.py` - Original GPU embedding computation
- `get_esm_embeddings.py` - Original ESM embedding script
- `compute_protein_embeddings.py` - Alternative embedding computation

## Data Storage

- Working directory: `/fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches/`
- Batch outputs: `embeddings_batch_XXXXX.npz` (1,184 files)
- Final cache: `results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache_full.npz` (after merging)

## Progress

Check current status:
```bash
# Count completed batches
ls /fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches/*.npz | wc -l

# Check job queue
squeue -j 41725169 | head -20

# Check a sample batch output
ls -lh /fh/working/srivatsan_s/dmullane_organism_scale/embeddings/batches/ | head -10

# View job progress
sacct -j 41725169 --format=JobID,State,Elapsed,MaxRSS | head -30
```

## Technical Notes

### BFloat16 Compatibility Fix

The ESM-C 600M model uses BFloat16 precision by default, which is not supported by older GPU architectures (Pascal/Turing: GTX 1080 Ti, RTX 2080 Ti). The script converts the model to Float32 for compatibility:

```python
client = ESMC.from_pretrained("esmc_600m", device=device)
client = client.float()  # Convert BFloat16 -> Float32
```

This conversion allows the embeddings to run on all available GPUs in the campus-new partition.

### Performance

- Processing speed: ~10-15 proteins/second per GPU
- Time per batch (10,000 proteins): ~15-20 minutes
- Parallel jobs: 50 concurrent GPUs
- Total completion time: ~6-8 hours

### Troubleshooting

**Error: "Got unsupported ScalarType BFloat16"**
- Cause: Model using BFloat16 on incompatible GPU
- Solution: Already fixed with `.float()` conversion in batch_generate_embeddings.py

**Empty output files**
- Check logs: `/fh/working/srivatsan_s/dmullane_organism_scale/logs/batch_embeddings_*`
- Verify model conversion is working
- Check GPU availability: `sinfo -o "%P %G %N" | grep gpu`
