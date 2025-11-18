# Component 1.1: K-mer Profiling - Status Report

## Overview
Generate 5-mer frequency profiles for all 7,664 bacterial genomes to serve as genome node attributes in the heterogeneous graph.

## Implementation Status: ✅ COMPLETE

### What's Implemented

1. **Core Profiling Script** (`generate_5mer_profiles.sh`)
   - Uses Jellyfish 2.3.0 for fast k-mer counting
   - Canonical k-mers (both strands counted together)
   - Generates both counts and statistics per genome
   - Runtime: ~1-3 min per genome

2. **SLURM Array Job** (`submit_5mer_array.sh`)
   - Processes all 7,664 genomes in parallel
   - Resume capability (skips completed genomes)
   - Proper logging to project logs directory
   - Memory: 8GB, 4 CPUs per genome

3. **Matrix Creation** (`create_kmer_matrix.py`)
   - Combines individual profiles into sparse matrix
   - Supports multiple normalization methods
   - Saves metadata (genome IDs, k-mer list)
   - Optimized for large datasets

## Testing Status: ✅ VERIFIED

Successfully tested on 3 genomes:
- GCF_000006985.1_Chlorobaculum_tepidum_TLS
- GCF_000007085.1_Caldanaerobacter_subterraneus_subsp__tengcongensis
- GCF_000007485.1_Tropheryma_whipplei_str__Twist

Matrix creation verified: 3 × 512 (100% density)

## Files & Structure

```
1_genome_to_graph/1.1_kmer_profiling/
├── generate_5mer_profiles.sh    # Core processing script
├── submit_5mer_array.sh          # SLURM submission
├── create_kmer_matrix.py         # Matrix assembly
├── README.md                     # Full documentation
├── QUICKSTART.md                 # Quick reference
└── STATUS.md                     # This file

results/1_genome_to_graph/1.1_kmer_profiling/
├── {genome_id}_5mer.txt          # Individual k-mer counts
├── {genome_id}_5mer_stats.txt    # Statistics per genome
├── 5mer_matrix.npz               # Combined sparse matrix
└── 5mer_matrix.meta.npz          # Metadata
```

## Ready to Run

**Command:**
```bash
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/organism_scale_modelling/1_genome_to_graph/1.1_kmer_profiling
sbatch submit_5mer_array.sh
```

**Expected runtime:** ~30 minutes (walltime, parallel)

## Output Specifications

### K-mer Profiles
- **Format**: Tab-separated (k-mer, count)
- **Size**: 512 canonical 5-mers per genome
- **Total observations**: ~2M k-mers per genome (avg)

### Combined Matrix
- **Dimensions**: 7,664 genomes × 512 k-mers
- **Format**: Sparse CSR matrix (.npz)
- **Density**: 100% (all k-mers expected)
- **Size**: ~10-50 MB

## Integration Points

This component provides input for:

1. **Component 1.6** (Genome-genome graph assembly)
   - K-mer distance calculations
   - Phylogenetic similarity metrics

2. **Component 1.7** (Full graph assembly)
   - Genome node attributes
   - 512-dimensional feature vectors per genome

## Dependencies

- **Software**: Jellyfish/2.3.0-GCC-11.2.0 (module)
- **Python**: /home/dmullane/micromamba/envs/esm3_env/bin/python
- **Packages**: numpy, scipy, pandas, tqdm

## Notes

- Symlink to genome data verified and working
- All paths follow component directory structure
- Scripts tested and ready for production run
- Logs will be saved to: `logs/1.1_kmer_profile_*.{out,err}`

## Next Steps

1. Submit full array job (7,664 genomes)
2. Monitor completion
3. Create combined matrix
4. Use matrix for genome-genome distance calculations (1.6)
