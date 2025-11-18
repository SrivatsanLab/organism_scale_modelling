# Component 1.1 Final Outputs

This directory contains the final compressed outputs from the k-mer profiling analysis of 7,664 bacterial genomes.

## Files

### 1. **6mer_matrix_compressed.npz** (47.9 MB)
Compressed genome × k-mer matrix
- Format: Scipy CSR sparse matrix (float32)
- Shape: 7,664 genomes × 2,080 canonical 6-mers
- Sparsity: 0.9999 (highly sparse)
- Load with: `scipy.sparse.load_npz()`

### 2. **genome_index.csv** (727 KB)
Genome index mapping matrix rows to genome IDs
- Columns:
  - `genome_index`: Row index in the matrix (0-7663)
  - `genome_id`: RefSeq genome identifier
  - `genus`: Taxonomic genus
  - `species`: Species name
  - `full_name`: Genus_species combined
- 7,664 genomes across 2,140 unique genera

### 3. **kmer_index.csv** (24 KB)
K-mer index mapping matrix columns to canonical 6-mer sequences
- Columns:
  - `kmer_index`: Column index in the matrix (0-2079)
  - `kmer`: 6-mer sequence (canonical representation)
- 2,080 canonical 6-mers (out of 4^6 = 4,096 possible, considering reverse complements)

### 4. **adjacency_matrix.npz** (198 KB)
K-nearest neighbor graph adjacency matrix
- Format: Scipy CSR sparse matrix (binary, undirected)
- Shape: 7,664 × 7,664
- Edges: 50,303 (undirected)
- Average degree: 13.1 neighbors per genome
- Construction method:
  - Distance metric: **Cosine similarity** in full 6-mer space (2080D)
  - k-neighbors: **k=10** (optimized for genus-level structure)
  - Graph is symmetrized (undirected) with self-loops removed
- **Genus purity: 0.384** (38.4% of neighbors share the same genus)
- Load with: `scipy.sparse.load_npz()`

### 5. **distance_matrix_compressed.npz** (48.2 MB) ⭐ Recommended
Pairwise cosine distance matrix (genome × genome) - **Git-friendly compressed format**
- Format: Upper triangle only, float16 precision
- Full matrix: 7,664 × 7,664 (symmetric with zeros on diagonal)
- Compression: 4.13x smaller than full matrix
- Precision: Max error 0.000244, mean error 0.000060
- Distance range: [0.0, 0.992], Mean: 0.349, Median: 0.270
- Load with helper function (see usage examples below)

### 5b. **distance_matrix.npz** (198.8 MB) - Full version (not in git)
Full pairwise distance matrix (uncompressed)
- Format: Full symmetric matrix (float32)
- Available locally but not tracked in git due to size (>100 MB GitHub limit)
- Identical data to compressed version, just uncompressed

### 6. **k_optimization.csv** (220 B)
Results from k-parameter optimization
- Shows genus purity for different k values (10, 15, 20, ..., 50)
- k=10 achieves highest genus purity (0.384)

### 7. **k_optimization.png**
Visualization of k-optimization showing genus purity vs k

### 8. **metadata.txt**
Summary metadata for all outputs

## Design Choices

### Why cosine distance in full 6-mer space?
- **Full 6-mer space vs PCA**: The full 2080D 6-mer space was chosen over 50D PCA space because it achieved slightly better genus purity (0.283 vs 0.279 at k=20)
- **Cosine similarity**: Most appropriate for frequency-based vectors, as it measures angular distance independent of magnitude
- This captures compositional similarity rather than absolute count differences

### Why k=10?
The k-nearest neighbors parameter was optimized to maximize genus-level taxonomic purity:
- k=10: 38.4% genus purity ← **OPTIMAL**
- k=15: 32.4% genus purity
- k=20: 28.3% genus purity
- k=30: 23.0% genus purity
- k=50: 17.3% genus purity

Smaller k values preserve local genus-level structure better. k=10 represents a good balance between:
1. Capturing meaningful genus-level relationships (38.4% of neighbors are same genus)
2. Maintaining sufficient connectivity for graph-based analyses
3. Average degree of ~13 neighbors per genome

### Matrix compression
- Scipy sparse matrices only support float32 as minimum float type (float16 not supported)
- CSR format chosen for efficient row slicing (genome-wise access)
- Already near-optimal compression at 47.9 MB for 7,664 × 2,080 matrix

## Usage Examples

### Load k-mer matrix
```python
from scipy.sparse import load_npz
import pandas as pd

# Load matrix and indices
matrix = load_npz('6mer_matrix_compressed.npz')
genomes = pd.read_csv('genome_index.csv')
kmers = pd.read_csv('kmer_index.csv')

# Get 6-mer profile for first genome
genome_idx = 0
profile = matrix[genome_idx].toarray().flatten()
print(f"Genome: {genomes.iloc[genome_idx]['genome_id']}")
print(f"Profile shape: {profile.shape}")
```

### Load adjacency graph
```python
from scipy.sparse import load_npz
import pandas as pd

# Load adjacency matrix and genome index
adj = load_npz('adjacency_matrix.npz')
genomes = pd.read_csv('genome_index.csv')

# Get neighbors for a genome
genome_idx = 0
neighbors_idx = adj[genome_idx].nonzero()[1]  # Column indices
neighbor_genomes = genomes.iloc[neighbors_idx]

print(f"Genome {genomes.iloc[genome_idx]['genome_id']} has {len(neighbors_idx)} neighbors")
print(neighbor_genomes[['genome_id', 'genus']])
```

### Compute genus purity for a genome
```python
genome_idx = 0
my_genus = genomes.iloc[genome_idx]['genus']
neighbors_idx = adj[genome_idx].nonzero()[1]
neighbor_genera = genomes.iloc[neighbors_idx]['genus']
purity = (neighbor_genera == my_genus).mean()
print(f"Genus purity: {purity:.3f}")
```

### Load and use distance matrix (compressed version)
```python
import numpy as np
import pandas as pd

# Method 1: Use helper function (recommended)
import sys
sys.path.append('../../1_genome_to_graph/1.1_kmer_profiling')
from compress_distance_matrix import decompress_distance_matrix

dist_matrix = decompress_distance_matrix('distance_matrix_compressed.npz')
genomes = pd.read_csv('genome_index.csv')

# Method 2: Manual decompression
# data = np.load('distance_matrix_compressed.npz')
# upper = data['upper_triangle'].astype(np.float32)
# n = int(data['n'])
# dist_matrix = np.zeros((n, n), dtype=np.float32)
# triu_idx = np.triu_indices(n, k=1)
# dist_matrix[triu_idx] = upper
# dist_matrix = dist_matrix + dist_matrix.T  # Symmetrize

# Get distances from one genome to all others
genome_idx = 0
distances = dist_matrix[genome_idx]

# Find k nearest neighbors
k = 10
nearest_indices = np.argsort(distances)[1:k+1]  # Skip self (distance=0)
nearest_distances = distances[nearest_indices]

print(f"Nearest neighbors for {genomes.iloc[genome_idx]['genome_id']}:")
for i, (idx, dist) in enumerate(zip(nearest_indices, nearest_distances)):
    print(f"{i+1}. {genomes.iloc[idx]['genome_id']} (distance={dist:.4f})")
```

### Get pairwise distance between two genomes
```python
genome_idx1 = 0
genome_idx2 = 100
distance = dist_matrix[genome_idx1, genome_idx2]
print(f"Distance: {distance:.4f}")
```

### Load full uncompressed version (if available locally)
```python
# Only if you have the full version locally (not in git)
dist_matrix = np.load('distance_matrix.npz')['distance_matrix']
```

## Notes

- The adjacency matrix is **undirected** (symmetric) and has **no self-loops**
- Both matrices use 0-based indexing aligned with the index files
- All distances/similarities computed in normalized 6-mer frequency space
- Genus purity of 0.384 indicates moderate phylogenetic signal in 6-mer profiles
- For graph algorithms, consider the adjacency matrix as an unweighted graph
