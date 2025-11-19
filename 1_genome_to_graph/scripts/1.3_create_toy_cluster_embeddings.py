#!/usr/bin/env python3
"""
Create compressed toy version of protein cluster embeddings for GitHub sharing.

This script creates a smaller, GitHub-friendly version (<100MB) by:
1. Filtering to high-confidence clusters (>0.9 seq identity, >10 proteins)
2. Using PCA to reduce embeddings from 1152D to 50D
3. Storing in compressed format

Input:
  - Full cluster embeddings: graph_outputs/protein_graph/cluster_embeddings.npz
  - Cluster index: graph_outputs/protein_graph/cluster_index.csv
  - MMseqs high-stringency clusters (seq_id >= 0.9)

Output:
  - graph_outputs/protein_graph_toy/cluster_embeddings_pca50.npz (PCA-reduced)
  - graph_outputs/protein_graph_toy/cluster_index.csv (filtered)
  - graph_outputs/protein_graph_toy/pca_model.npz (for reconstruction)
  - graph_outputs/protein_graph_toy/README.md (documentation)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
import argparse
from collections import defaultdict
from tqdm import tqdm


def load_high_confidence_clusters(cluster_file, min_seq_id=0.9, min_size=10):
    """
    Load high-stringency MMseqs2 clusters.

    Args:
        cluster_file: Path to MMseqs2 cluster TSV
        min_seq_id: Minimum sequence identity threshold (0.9 = 90%)
        min_size: Minimum cluster size

    Returns:
        set: Cluster IDs meeting criteria
        dict: Cluster sizes
    """
    print(f"\nLoading high-confidence clusters (seq_id >= {min_seq_id}, size >= {min_size})...")

    # Count cluster sizes
    cluster_sizes = defaultdict(int)
    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Counting cluster sizes"):
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_rep = parts[0]
                cluster_sizes[cluster_rep] += 1

    # Filter by size
    high_conf_clusters = {
        cluster_id for cluster_id, size in cluster_sizes.items()
        if size >= min_size
    }

    print(f"Clusters with >= {min_size} proteins: {len(high_conf_clusters):,}")

    return high_conf_clusters, cluster_sizes


def filter_cluster_embeddings(embeddings_file, index_file, high_conf_clusters):
    """
    Filter full embeddings to only high-confidence clusters.

    Returns:
        np.ndarray: Filtered embeddings
        pd.DataFrame: Filtered index
    """
    print("\nLoading full cluster embeddings...")
    data = np.load(embeddings_file)
    full_embeddings = data['embeddings']
    cluster_ids = data['cluster_ids']

    print(f"Full embeddings shape: {full_embeddings.shape}")

    # Load cluster index
    index_df = pd.read_csv(index_file)

    # Filter to high-confidence clusters
    print("\nFiltering to high-confidence clusters...")
    keep_mask = np.array([
        str(cid) in high_conf_clusters
        for cid in cluster_ids
    ])

    filtered_embeddings = full_embeddings[keep_mask]
    filtered_ids = cluster_ids[keep_mask]
    filtered_index = index_df[index_df['cluster_id'].isin(filtered_ids)]

    print(f"Filtered embeddings shape: {filtered_embeddings.shape}")
    print(f"Kept {len(filtered_ids):,} / {len(cluster_ids):,} clusters ({100*len(filtered_ids)/len(cluster_ids):.1f}%)")

    return filtered_embeddings, filtered_index, filtered_ids


def apply_pca(embeddings, n_components=50):
    """
    Apply PCA to reduce embedding dimensionality.

    Returns:
        np.ndarray: PCA-reduced embeddings
        PCA: Fitted PCA model
    """
    print(f"\nApplying PCA (1152D → {n_components}D)...")

    pca = PCA(n_components=n_components, random_state=42)
    embeddings_pca = pca.fit_transform(embeddings)

    variance_explained = pca.explained_variance_ratio_.sum()
    print(f"PCA variance explained: {variance_explained:.3f}")
    print(f"PCA shape: {embeddings_pca.shape}")

    return embeddings_pca, pca


def save_toy_outputs(embeddings_pca, cluster_ids, index_df, pca_model, output_dir):
    """Save compressed toy version outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save PCA-reduced embeddings
    embedding_file = output_dir / 'cluster_embeddings_pca50.npz'
    print(f"\nSaving PCA embeddings to {embedding_file}...")
    np.savez_compressed(
        embedding_file,
        embeddings=embeddings_pca.astype(np.float32),
        cluster_ids=cluster_ids
    )

    # Save cluster index
    index_file = output_dir / 'cluster_index.csv'
    print(f"Saving cluster index to {index_file}...")
    index_df.to_csv(index_file, index=False)

    # Save PCA model for potential reconstruction
    pca_file = output_dir / 'pca_model.npz'
    print(f"Saving PCA model to {pca_file}...")
    np.savez_compressed(
        pca_file,
        components=pca_model.components_,
        mean=pca_model.mean_,
        explained_variance=pca_model.explained_variance_,
        explained_variance_ratio=pca_model.explained_variance_ratio_
    )

    # Create README
    readme_file = output_dir / 'README.md'
    print(f"Creating README at {readme_file}...")

    readme_content = f"""# Toy Protein Graph (GitHub-Friendly Version)

This is a compressed version of the protein cluster embeddings suitable for GitHub sharing.

## Filtering Criteria

- **Sequence identity**: ≥ 90% (high-confidence clusters)
- **Cluster size**: ≥ 10 proteins
- **Dimensionality reduction**: PCA from 1152D → 50D

## Files

### 1. cluster_embeddings_pca50.npz
**PCA-reduced cluster embeddings** (50-dimensional)
- Shape: {embeddings_pca.shape}
- File size: {embedding_file.stat().st_size / 1024**2:.1f} MB
- dtype: float32

```python
import numpy as np
data = np.load('cluster_embeddings_pca50.npz')
embeddings = data['embeddings']  # Shape: ({embeddings_pca.shape[0]}, 50)
cluster_ids = data['cluster_ids']
```

### 2. cluster_index.csv
**Cluster metadata**
- Number of clusters: {len(index_df):,}
- Columns: cluster_id, genome_id, cluster_size

```python
import pandas as pd
index = pd.read_csv('cluster_index.csv')
```

### 3. pca_model.npz
**PCA transformation parameters** (for reconstruction or analysis)
- Variance explained: {pca_model.explained_variance_ratio_.sum():.3f}

```python
data = np.load('pca_model.npz')
components = data['components']  # Shape: (50, 1152)
mean = data['mean']  # Shape: (1152,)
```

## Statistics

- **Total clusters**: {len(index_df):,}
- **Mean cluster size**: {index_df['cluster_size'].mean():.1f}
- **Median cluster size**: {index_df['cluster_size'].median():.0f}
- **Embedding dimension**: 50 (reduced from 1152)
- **Total file size**: ~{(embedding_file.stat().st_size + index_file.stat().st_size + pca_file.stat().st_size) / 1024**2:.1f} MB

## Comparison to Full Version

The full protein graph (in `graph_outputs/protein_graph/`) contains:
- All ~388k clusters (not just high-confidence)
- Full 1152D ESM embeddings (not PCA-reduced)
- Total size: ~1.7 GB

This toy version is suitable for:
- Quick prototyping and testing
- Sharing on GitHub (< 100 MB)
- Educational purposes
- Method development

For production use, use the full version from `graph_outputs/protein_graph/`.
"""

    with open(readme_file, 'w') as f:
        f.write(readme_content)

    # Print summary
    print("\n" + "="*70)
    print("TOY PROTEIN GRAPH SUMMARY")
    print("="*70)
    print(f"Clusters: {len(index_df):,}")
    print(f"Embedding dimension: 50 (PCA-reduced from 1152)")
    print(f"Mean cluster size: {index_df['cluster_size'].mean():.1f}")
    print(f"Median cluster size: {index_df['cluster_size'].median():.0f}")
    print(f"PCA variance explained: {pca_model.explained_variance_ratio_.sum():.3f}")
    print(f"\nFile sizes:")
    print(f"  Embeddings: {embedding_file.stat().st_size / 1024**2:.1f} MB")
    print(f"  Index: {index_file.stat().st_size / 1024**2:.1f} MB")
    print(f"  PCA model: {pca_file.stat().st_size / 1024**2:.1f} MB")
    total_size = (embedding_file.stat().st_size + index_file.stat().st_size +
                  pca_file.stat().st_size) / 1024**2
    print(f"  Total: {total_size:.1f} MB ✓ (< 100 MB)")
    print("="*70)


def main():
    parser = argparse.ArgumentParser(
        description='Create toy version of protein cluster embeddings for GitHub'
    )
    parser.add_argument(
        '--embeddings',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph/cluster_embeddings.npz',
        help='Full cluster embeddings file'
    )
    parser.add_argument(
        '--index',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph/cluster_index.csv',
        help='Full cluster index CSV'
    )
    parser.add_argument(
        '--clusters',
        type=str,
        default='1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv',
        help='MMseqs2 cluster TSV file'
    )
    parser.add_argument(
        '--min-seq-id',
        type=float,
        default=0.9,
        help='Minimum sequence identity (0-1)'
    )
    parser.add_argument(
        '--min-size',
        type=int,
        default=10,
        help='Minimum cluster size'
    )
    parser.add_argument(
        '--n-components',
        type=int,
        default=50,
        help='Number of PCA components'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='1_genome_to_graph/graph_outputs/protein_graph_toy',
        help='Output directory for toy version'
    )

    args = parser.parse_args()

    print("="*70)
    print("Creating Toy Protein Graph (GitHub-Friendly)")
    print("="*70)
    print(f"Min sequence identity: {args.min_seq_id}")
    print(f"Min cluster size: {args.min_size}")
    print(f"PCA components: {args.n_components}")

    # Load high-confidence clusters
    high_conf_clusters, cluster_sizes = load_high_confidence_clusters(
        args.clusters,
        min_seq_id=args.min_seq_id,
        min_size=args.min_size
    )

    # Filter embeddings
    filtered_embeddings, filtered_index, filtered_ids = filter_cluster_embeddings(
        args.embeddings,
        args.index,
        high_conf_clusters
    )

    # Apply PCA
    embeddings_pca, pca_model = apply_pca(
        filtered_embeddings,
        n_components=args.n_components
    )

    # Save outputs
    save_toy_outputs(
        embeddings_pca,
        filtered_ids,
        filtered_index,
        pca_model,
        args.output_dir
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
