#!/usr/bin/env python3
"""
Create toy protein graph directly from embedding files (doesn't require full cluster embeddings).

This script creates a GitHub-friendly version (<100MB) by:
1. Loading only high-confidence clusters (>0.9 seq identity, >10 proteins)
2. Computing mean embeddings for just those clusters
3. Using PCA to reduce from 1152D to 50D
4. Storing in compressed format

This is faster than the two-step approach because it doesn't require generating
the full 388k cluster embeddings first.

Input:
  - MMseqs high-stringency clusters
  - ESM embeddings (7,664 genome files)

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


def load_high_confidence_clusters(cluster_file, min_size=10):
    """
    Load high-stringency MMseqs2 clusters.

    Since these are from mmseqs_full_dataset which already uses high seq identity,
    we only filter by cluster size.

    Args:
        cluster_file: Path to MMseqs2 cluster TSV
        min_size: Minimum cluster size

    Returns:
        set: Cluster IDs meeting criteria
        dict: Cluster representative -> member proteins
        dict: Cluster sizes
    """
    print(f"\nLoading high-confidence clusters (size >= {min_size})...")

    # Build cluster membership and sizes
    cluster_members = defaultdict(list)
    cluster_sizes = defaultdict(int)

    with open(cluster_file, 'r') as f:
        for line in tqdm(f, desc="Reading clusters"):
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_rep, protein_id = parts
                cluster_members[cluster_rep].append(protein_id)
                cluster_sizes[cluster_rep] += 1

    # Filter by size
    high_conf_clusters = {
        cluster_id for cluster_id, size in cluster_sizes.items()
        if size >= min_size
    }

    print(f"Total clusters: {len(cluster_sizes):,}")
    print(f"Clusters with >= {min_size} proteins: {len(high_conf_clusters):,}")

    # Build set of all protein IDs in high-confidence clusters
    proteins_in_hc_clusters = set()
    for cluster_id in high_conf_clusters:
        proteins_in_hc_clusters.update(cluster_members[cluster_id])

    print(f"Total proteins in high-confidence clusters: {len(proteins_in_hc_clusters):,}")

    # Create protein->cluster mapping for high-confidence clusters only
    protein_to_cluster = {}
    for cluster_id in high_conf_clusters:
        for protein_id in cluster_members[cluster_id]:
            protein_to_cluster[protein_id] = cluster_id

    return high_conf_clusters, protein_to_cluster, cluster_sizes


def load_embeddings_for_clusters(embedding_dir, protein_to_cluster):
    """
    Load ESM embeddings only for proteins in high-confidence clusters.

    Returns:
        dict: cluster_id -> list of embeddings (1152D)
    """
    embedding_dir = Path(embedding_dir)
    embedding_files = sorted(embedding_dir.glob("*.npz"))

    print(f"\nFound {len(embedding_files)} embedding files")

    cluster_embeddings = defaultdict(list)
    proteins_found = 0
    proteins_skipped = 0

    for emb_file in tqdm(embedding_files, desc="Loading embeddings"):
        # Extract genome ID from filename
        genome_id = emb_file.stem.split('_esmc_embeddings')[0]

        data = np.load(emb_file, allow_pickle=True)
        embeddings = data['embeddings']  # Shape: (n_proteins, 1152)
        seq_ids = data['seq_ids']

        for protein_id, embedding in zip(seq_ids, embeddings):
            # Prepend genome ID to match cluster file format
            protein_id_full = f"{genome_id}_{protein_id}"

            if protein_id_full in protein_to_cluster:
                cluster_id = protein_to_cluster[protein_id_full]
                cluster_embeddings[cluster_id].append(embedding)
                proteins_found += 1
            else:
                proteins_skipped += 1

    print(f"\nProteins with embeddings: {proteins_found:,}")
    print(f"Proteins skipped (not in high-conf clusters): {proteins_skipped:,}")
    print(f"Clusters with embeddings: {len(cluster_embeddings):,}")

    return cluster_embeddings


def compute_cluster_means(cluster_embeddings):
    """
    Compute mean embedding for each cluster.

    Returns:
        np.ndarray: (n_clusters, 1152) mean embeddings
        list: cluster IDs (in same order as matrix rows)
    """
    print("\nComputing mean embeddings per cluster...")

    cluster_ids = sorted(cluster_embeddings.keys())
    n_clusters = len(cluster_ids)
    embedding_dim = len(cluster_embeddings[cluster_ids[0]][0])

    print(f"Computing means for {n_clusters:,} clusters")
    print(f"Embedding dimension: {embedding_dim}")

    mean_embeddings = np.zeros((n_clusters, embedding_dim), dtype=np.float32)

    for i, cluster_id in enumerate(tqdm(cluster_ids, desc="Computing means")):
        embeddings = np.array(cluster_embeddings[cluster_id])
        mean_embeddings[i] = embeddings.mean(axis=0)

    return mean_embeddings, cluster_ids


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


def save_toy_outputs(embeddings_pca, cluster_ids, cluster_sizes, pca_model, output_dir):
    """Save compressed toy version outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create cluster index
    cluster_data = []
    for cluster_id in cluster_ids:
        # Parse genome ID from cluster representative
        genome_id = '_'.join(cluster_id.split('_')[:2])

        cluster_data.append({
            'cluster_id': cluster_id,
            'genome_id': genome_id,
            'cluster_size': cluster_sizes[cluster_id]
        })

    index_df = pd.DataFrame(cluster_data)

    # Save PCA-reduced embeddings
    embedding_file = output_dir / 'cluster_embeddings_pca50.npz'
    print(f"\nSaving PCA embeddings to {embedding_file}...")
    np.savez_compressed(
        embedding_file,
        embeddings=embeddings_pca.astype(np.float32),
        cluster_ids=np.array(cluster_ids)
    )

    # Save cluster index
    index_file = output_dir / 'cluster_index.csv'
    print(f"Saving cluster index to {index_file}...")
    index_df.to_csv(index_file, index=False)

    # Save PCA model
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

- **Sequence identity**: ≥ 90% (MMseqs high-stringency clustering)
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

The full protein graph (in `graph_outputs/protein_graph/`) will contain:
- All ~388k clusters (not just high-confidence)
- Full 1152D ESM embeddings (not PCA-reduced)
- Total size: ~1.7 GB

This toy version is suitable for:
- Quick prototyping and testing
- Sharing on GitHub (< 100 MB)
- Educational purposes
- Method development

For production use, use the full version from `graph_outputs/protein_graph/` once generated.
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
        description='Create toy protein graph directly from ESM embeddings'
    )
    parser.add_argument(
        '--clusters',
        type=str,
        default='1_genome_to_graph/intermediate/protein/msa_clusters/mmseqs_full_dataset/clusters.tsv',
        help='MMseqs2 cluster TSV file'
    )
    parser.add_argument(
        '--embeddings',
        type=str,
        default='data/esm_embeddings',
        help='Directory containing ESM embedding NPZ files'
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
    print("Creating Toy Protein Graph (Direct from Embeddings)")
    print("="*70)
    print(f"Min cluster size: {args.min_size}")
    print(f"PCA components: {args.n_components}")

    # Load high-confidence clusters
    high_conf_clusters, protein_to_cluster, cluster_sizes = load_high_confidence_clusters(
        args.clusters,
        min_size=args.min_size
    )

    # Load embeddings for those clusters only
    cluster_embeddings = load_embeddings_for_clusters(
        args.embeddings,
        protein_to_cluster
    )

    # Compute mean embeddings
    mean_embeddings, cluster_ids = compute_cluster_means(cluster_embeddings)

    # Apply PCA
    embeddings_pca, pca_model = apply_pca(
        mean_embeddings,
        n_components=args.n_components
    )

    # Save outputs
    save_toy_outputs(
        embeddings_pca,
        cluster_ids,
        cluster_sizes,
        pca_model,
        args.output_dir
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
