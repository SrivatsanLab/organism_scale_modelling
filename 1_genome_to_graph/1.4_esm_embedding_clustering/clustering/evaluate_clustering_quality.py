#!/usr/bin/env python
"""
Comprehensive clustering quality evaluation.

Implements multiple metrics beyond simple homogeneity:
1. Adjusted Rand Index (ARI) - similarity to COG annotations
2. Adjusted Mutual Information (AMI) - information-theoretic similarity
3. Silhouette Score - cluster separation quality
4. Davies-Bouldin Index - cluster compactness
5. Cluster size distribution analysis

Usage:
    python evaluate_clustering_quality.py \
        --clustering results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/clusters_leiden_res1500_nn15_cogonly.npz \
        --pca results/1_genome_to_graph/1.4_esm_embedding_clustering/umap/pca_cache.npz \
        --output results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/quality_metrics.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from sklearn.metrics import (
    adjusted_rand_score,
    adjusted_mutual_info_score,
    silhouette_score,
    davies_bouldin_score,
    silhouette_samples
)
from collections import Counter
import sys

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def load_cog_annotations(gene_ids, genome_ids):
    """Load COG category annotations for genes."""
    diamond_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation')
    cog_db_dir = Path('/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/cog')
    cog_csv_file = cog_db_dir / 'cog-20.cog.csv'
    cog_def_file = cog_db_dir / 'cog-20.def.tab'

    # Load protein → COG mapping
    prot_to_cog = {}
    with open(cog_csv_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            fields = line.strip().split(',')
            if len(fields) >= 7:
                prot_id = fields[2]
                cog_id = fields[6]
                if prot_id not in prot_to_cog:
                    prot_to_cog[prot_id] = cog_id

    # Load COG → category mapping
    cog_def = pd.read_csv(cog_def_file, sep='\t', header=None, encoding='latin-1',
                         names=['cog_id', 'category', 'description', 'gene', 'pathway', 'extra1', 'extra2'])
    cog_to_category = dict(zip(cog_def['cog_id'], cog_def['category']))

    # Process DIAMOND hits to get gene → COG category mapping
    cog_lookup = {}
    unique_genomes = set(genome_ids)
    gene_id_set = set(gene_ids)

    for genome_id in unique_genomes:
        hit_file = diamond_dir / f'{genome_id}_cog_diamond' / f'{genome_id}_cog_hits.tsv'
        if not hit_file.exists():
            continue

        try:
            hits_df = pd.read_csv(
                hit_file, sep='\t',
                names=['gene_id', 'protein_id', 'pident', 'length', 'mismatch',
                       'gapopen', 'qstart', 'qend', 'sstart', 'send',
                       'evalue', 'bitscore', 'qcovhsp', 'scovhsp']
            )

            hits_df = hits_df[hits_df['gene_id'].isin(gene_id_set)]

            for gene_id, group in hits_df.groupby('gene_id'):
                best_hit = group.iloc[0]
                protein_id = best_hit['protein_id']

                if pd.notna(protein_id):
                    # Convert underscore to dot
                    if '_' in protein_id:
                        parts = protein_id.rsplit('_', 1)
                        protein_id_dot = f"{parts[0]}.{parts[1]}"
                    else:
                        protein_id_dot = protein_id

                    # Map protein → COG → category
                    if protein_id_dot in prot_to_cog:
                        cog_id = prot_to_cog[protein_id_dot]
                        if cog_id in cog_to_category:
                            categories = cog_to_category[cog_id]
                            if len(categories) > 0 and categories[0] != '-':
                                cog_lookup[gene_id] = categories[0]
        except Exception as e:
            continue

    return cog_lookup


def evaluate_against_cog(cluster_labels, gene_ids, genome_ids):
    """
    Evaluate clustering quality against COG annotations.

    Returns:
        dict with ARI, AMI scores
    """
    print("Loading COG annotations...")
    cog_lookup = load_cog_annotations(gene_ids, genome_ids)

    # Create COG label array aligned with clustering
    cog_labels = []
    valid_indices = []

    for i, gene_id in enumerate(gene_ids):
        if gene_id in cog_lookup:
            cog_labels.append(cog_lookup[gene_id])
            valid_indices.append(i)

    if len(valid_indices) == 0:
        print("Warning: No COG annotations found!")
        return {
            'ari_vs_cog': None,
            'ami_vs_cog': None,
            'n_annotated': 0,
            'annotation_rate': 0.0
        }

    # Filter to annotated genes only
    cluster_labels_filtered = cluster_labels[valid_indices]

    # Compute metrics
    ari = adjusted_rand_score(cog_labels, cluster_labels_filtered)
    ami = adjusted_mutual_info_score(cog_labels, cluster_labels_filtered)

    print(f"  Annotated genes: {len(valid_indices):,} / {len(gene_ids):,} "
          f"({100*len(valid_indices)/len(gene_ids):.1f}%)")
    print(f"  ARI vs COG: {ari:.4f}")
    print(f"  AMI vs COG: {ami:.4f}")

    return {
        'ari_vs_cog': ari,
        'ami_vs_cog': ami,
        'n_annotated': len(valid_indices),
        'annotation_rate': len(valid_indices) / len(gene_ids)
    }


def evaluate_embedding_metrics(embeddings, cluster_labels, sample_size=50000):
    """
    Evaluate clustering quality using embedding-based metrics.

    Uses sampling to make silhouette score tractable for large datasets.

    Returns:
        dict with silhouette_score, davies_bouldin_score
    """
    print("Computing embedding-based metrics...")

    # Filter out noise points (label -1)
    valid_mask = cluster_labels >= 0
    embeddings_valid = embeddings[valid_mask]
    labels_valid = cluster_labels[valid_mask]

    n_valid = len(labels_valid)
    n_clusters = len(np.unique(labels_valid))

    print(f"  Valid points: {n_valid:,}")
    print(f"  Clusters: {n_clusters:,}")

    if n_clusters < 2:
        print("  Warning: Need at least 2 clusters for metrics")
        return {
            'silhouette_score': None,
            'davies_bouldin_score': None,
            'silhouette_score_sampled': None
        }

    # Davies-Bouldin (relatively fast even for large data)
    print("  Computing Davies-Bouldin score...")
    db_score = davies_bouldin_score(embeddings_valid, labels_valid)
    print(f"    Davies-Bouldin: {db_score:.4f}")

    # Silhouette on sample (too slow for full dataset)
    silhouette_sampled = None
    if n_valid > sample_size:
        print(f"  Computing silhouette score on sample of {sample_size:,}...")
        sample_indices = np.random.choice(n_valid, sample_size, replace=False)
        embeddings_sample = embeddings_valid[sample_indices]
        labels_sample = labels_valid[sample_indices]
        silhouette_sampled = silhouette_score(embeddings_sample, labels_sample)
        print(f"    Silhouette (sampled): {silhouette_sampled:.4f}")
    else:
        print(f"  Computing silhouette score on all {n_valid:,} points...")
        silhouette_sampled = silhouette_score(embeddings_valid, labels_valid)
        print(f"    Silhouette: {silhouette_sampled:.4f}")

    return {
        'silhouette_score_sampled': silhouette_sampled,
        'davies_bouldin_score': db_score,
        'n_valid_points': n_valid,
        'n_clusters_nonzero': n_clusters
    }


def analyze_cluster_size_distribution(cluster_labels):
    """
    Analyze the distribution of cluster sizes.

    Returns:
        dict with size distribution statistics
    """
    print("Analyzing cluster size distribution...")

    # Count cluster sizes
    cluster_counts = Counter(cluster_labels[cluster_labels >= 0])
    sizes = np.array(list(cluster_counts.values()))

    if len(sizes) == 0:
        return {
            'n_clusters': 0,
            'mean_size': 0,
            'median_size': 0,
            'min_size': 0,
            'max_size': 0,
            'size_std': 0,
            'size_gini': 0
        }

    # Compute statistics
    n_clusters = len(sizes)
    mean_size = np.mean(sizes)
    median_size = np.median(sizes)
    min_size = np.min(sizes)
    max_size = np.max(sizes)
    size_std = np.std(sizes)

    # Gini coefficient (inequality measure)
    # Perfect equality = 0, perfect inequality = 1
    sorted_sizes = np.sort(sizes)
    n = len(sorted_sizes)
    index = np.arange(1, n + 1)
    gini = (2 * np.sum(index * sorted_sizes)) / (n * np.sum(sorted_sizes)) - (n + 1) / n

    print(f"  Clusters: {n_clusters:,}")
    print(f"  Size - Mean: {mean_size:.1f}, Median: {median_size:.0f}")
    print(f"  Size - Min: {min_size}, Max: {max_size:,}")
    print(f"  Size - Std: {size_std:.1f}")
    print(f"  Size - Gini: {gini:.4f}")

    # Percentiles
    percentiles = {}
    for p in [10, 25, 50, 75, 90, 95, 99]:
        percentiles[f'size_p{p}'] = np.percentile(sizes, p)

    return {
        'n_clusters': n_clusters,
        'mean_size': float(mean_size),
        'median_size': float(median_size),
        'min_size': int(min_size),
        'max_size': int(max_size),
        'size_std': float(size_std),
        'size_gini': float(gini),
        **{k: float(v) for k, v in percentiles.items()}
    }


def evaluate_clustering(clustering_file, pca_file, sample_size=50000):
    """
    Comprehensive evaluation of a single clustering.

    Args:
        clustering_file: Path to clustering .npz file
        pca_file: Path to PCA embeddings .npz file
        sample_size: Number of points to sample for silhouette score

    Returns:
        dict with all metrics
    """
    print(f"\n{'='*80}")
    print(f"Evaluating: {clustering_file.name}")
    print(f"{'='*80}\n")

    # Load clustering
    print("Loading clustering...")
    clustering_data = np.load(clustering_file, allow_pickle=True)
    cluster_labels = clustering_data['labels']
    gene_ids = clustering_data['gene_ids']
    genome_ids = clustering_data['genome_ids']

    # Load PCA embeddings
    print("Loading PCA embeddings...")
    pca_data = np.load(pca_file, allow_pickle=True)
    pca_embeddings = pca_data['embeddings_pca']
    pca_gene_ids = pca_data['gene_ids']

    # Align embeddings to clustering genes
    print("Aligning embeddings to clustering genes...")
    gene_to_idx = {gid: i for i, gid in enumerate(pca_gene_ids)}
    embedding_indices = [gene_to_idx[gid] for gid in gene_ids if gid in gene_to_idx]

    if len(embedding_indices) != len(gene_ids):
        print(f"  Warning: Only {len(embedding_indices)}/{len(gene_ids)} genes found in PCA")

    embeddings_aligned = pca_embeddings[embedding_indices]

    # Extract metadata
    metadata = {
        'file': clustering_file.name,
        'resolution': float(clustering_data.get('resolution', 0)),
        'n_neighbors': int(clustering_data.get('n_neighbors', 0)),
        'cog_only': bool(clustering_data.get('cog_only', False)),
        'n_genes': len(gene_ids)
    }

    # Run evaluations
    results = {**metadata}

    # 1. COG comparison metrics
    cog_metrics = evaluate_against_cog(cluster_labels, gene_ids, genome_ids)
    results.update(cog_metrics)

    # 2. Embedding-based metrics
    embedding_metrics = evaluate_embedding_metrics(
        embeddings_aligned, cluster_labels, sample_size
    )
    results.update(embedding_metrics)

    # 3. Size distribution
    size_metrics = analyze_cluster_size_distribution(cluster_labels)
    results.update(size_metrics)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive clustering quality evaluation'
    )
    parser.add_argument('--clustering', type=str, required=True,
                        help='Path to clustering .npz file')
    parser.add_argument('--pca', type=str, required=True,
                        help='Path to PCA embeddings .npz file')
    parser.add_argument('--output', type=str,
                        default='results/1_genome_to_graph/1.4_esm_embedding_clustering/clustering/quality_metrics.csv',
                        help='Output CSV file')
    parser.add_argument('--sample-size', type=int, default=50000,
                        help='Sample size for silhouette score')
    parser.add_argument('--batch', action='store_true',
                        help='Evaluate all clusterings in directory')

    args = parser.parse_args()

    pca_file = Path(args.pca)
    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if args.batch:
        # Evaluate all clusterings
        clustering_dir = Path(args.clustering)
        clustering_files = sorted(clustering_dir.glob('clusters_leiden_*.npz'))

        print(f"Found {len(clustering_files)} clustering files")

        results = []
        for clustering_file in clustering_files:
            try:
                result = evaluate_clustering(clustering_file, pca_file, args.sample_size)
                results.append(result)
            except Exception as e:
                print(f"Error evaluating {clustering_file.name}: {e}")
                continue

        # Save results
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        print(f"\n{'='*80}")
        print(f"Saved results to {output_file}")
        print(f"  Total evaluations: {len(df)}")
        print(f"{'='*80}")

    else:
        # Evaluate single clustering
        clustering_file = Path(args.clustering)
        result = evaluate_clustering(clustering_file, pca_file, args.sample_size)

        # Save as single-row CSV
        df = pd.DataFrame([result])
        df.to_csv(output_file, index=False)
        print(f"\nSaved results to {output_file}")


if __name__ == '__main__':
    main()
