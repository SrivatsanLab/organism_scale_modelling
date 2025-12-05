# 5-mer Profiling Quick Start

## Summary

Complete 5-mer profiling pipeline for 7,664 bacterial genomes using Jellyfish.

**Status**: ✅ Tested and ready to run
**Location**: `1_genome_to_graph/1.1_kmer_profiling/`

## One-Command Setup

```bash
# From project root directory
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling

# 1. Test on single genome (verify installation)
bash 1_genome_to_graph/1.1_kmer_profiling/generate_5mer_profiles.sh \
    data/refseq_genomes/GCF_000006985.1_Chlorobaculum_tepidum_TLS.fasta \
    results/1_genome_to_graph/1.1_kmer_profiling/test/

# 2. Submit full job (all 7,664 genomes)
cd 1_genome_to_graph/1.1_kmer_profiling
sbatch submit_5mer_array.sh

# 3. Monitor progress
watch "ls ../../results/1_genome_to_graph/1.1_kmer_profiling/*_5mer.txt | wc -l"

# 4. Create combined matrix (after all genomes complete)
cd ../..
/home/dmullane/micromamba/envs/esm3_env/bin/python \
    1_genome_to_graph/1.1_kmer_profiling/create_kmer_matrix.py \
    --input-dir results/1_genome_to_graph/1.1_kmer_profiling \
    --output results/1_genome_to_graph/1.1_kmer_profiling/5mer_matrix.npz \
    --normalize frequency
```

## Expected Output

### Individual Profiles (per genome)
Location: `results/1_genome_to_graph/1.1_kmer_profiling/`
- `{genome_id}_5mer.txt`: K-mer counts (512 canonical 5-mers)
- `{genome_id}_5mer_stats.txt`: Statistics summary

### Combined Matrix
Location: `results/1_genome_to_graph/1.1_kmer_profiling/`
- `5mer_matrix.npz`: Sparse matrix (7664 × 512)
- `5mer_matrix.meta.npz`: Metadata (genome IDs, k-mer list)

## Performance

- **Per genome**: ~1-3 minutes
- **Total (parallel)**: ~30 minutes walltime
- **Matrix creation**: ~5-10 minutes

## Test Results

✅ Successfully tested on 3 genomes:
- GCF_000006985.1_Chlorobaculum_tepidum_TLS
- GCF_000007085.1_Caldanaerobacter_subterraneus_subsp__tengcongensis
- GCF_000007485.1_Tropheryma_whipplei_str__Twist

Matrix verified: 3 × 512 (100% density expected for 5-mers)

## Next Steps

After matrix creation, use for:
- Genome-genome distance calculations (component 1.6)
- Genome node attributes in heterogeneous graph (component 1.7)
- Phylogenetic clustering
- Integration with protein embeddings

## Full Documentation

See `1_genome_to_graph/1.1_kmer_profiling/README.md` for complete details.
