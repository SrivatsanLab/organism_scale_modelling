#!/usr/bin/env python
"""
Comprehensive Leiden clustering with COG filtering and evaluation.

Features:
- Multiple resolution parameters
- Multiple n_neighbors for graph construction
- Optional filtering to COG-annotated genes only
- Cluster quality evaluation (COG homogeneity)
"""

import argparse
import numpy as np
from pathlib import Path
import time
import sys

def load_pca(pca_file):
    """Load PCA embeddings."""
    print(f"Loading PCA from {pca_file}...")
    data = np.load(pca_file, allow_pickle=True)

    pca_embedding = data['embeddings_pca']
    gene_ids = data['gene_ids']
    genome_ids = data['genome_ids']

    print(f"  PCA shape: {pca_embedding.shape}")
    print(f"  Genes: {len(gene_ids):,}")

    return pca_embedding, gene_ids, genome_ids


def load_cog_annotations(gene_ids, genome_ids):
    """Load COG annotations for genes."""
    import pandas as pd

    print("Loading COG annotations...")

    cog_db_dir = Path('/fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/databases/cog')
    cog_csv_file = cog_db_dir / 'cog-20.cog.csv'
    cog_def_file = cog_db_dir / 'cog-20.def.tab'

    # Load protein → COG mapping
    print(f"  Loading protein→COG mapping...")
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
    print(f"  Loading COG→category mapping...")
    cog_def = pd.read_csv(cog_def_file, sep='\t', header=None, encoding='latin-1',
                         names=['cog_id', 'category', 'description', 'gene', 'pathway', 'extra1', 'extra2'])
    cog_to_category = dict(zip(cog_def['cog_id'], cog_def['category']))

    # Load DIAMOND hits
    unique_genomes = set(genome_ids)
    gene_id_set = set(gene_ids)
    diamond_dir = Path('results/1_genome_to_graph/1.4_esm_embedding_clustering/functional_annotation')
    cog_lookup = {}

    print(f"  Loading annotations from {len(unique_genomes)} genomes...")
    for i, genome_id in enumerate(unique_genomes):
        if (i + 1) % 100 == 0:
            print(f"    Progress: {i+1}/{len(unique_genomes)} genomes")

        hit_file = diamond_dir / f'{genome_id}_cog_diamond' / f'{genome_id}_cog_hits.tsv'
        if not hit_file.exists():
            continue

        try:
            hits_df = pd.read_csv(hit_file, sep='\t',
                                 names=['gene_id', 'protein_id', 'pident', 'length', 'mismatch',
                                       'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                       'evalue', 'bitscore', 'qcovhsp', 'scovhsp'])

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

    n_annotated = len(cog_lookup)
    pct_annotated = 100 * n_annotated / len(gene_ids)
    print(f"  Annotated: {n_annotated:,} / {len(gene_ids):,} ({pct_annotated:.1f}%)")

    return cog_lookup


def filter_to_annotated(embeddings, gene_ids, genome_ids, cog_lookup):
    """Filter to only COG-annotated genes."""
    print("\nFiltering to COG-annotated genes...")

    # Create mask for annotated genes
    mask = np.array([gid in cog_lookup for gid in gene_ids])

    filtered_embeddings = embeddings[mask]
    filtered_gene_ids = gene_ids[mask]
    filtered_genome_ids = genome_ids[mask]

    print(f"  Kept {len(filtered_gene_ids):,} / {len(gene_ids):,} genes ({100*len(filtered_gene_ids)/len(gene_ids):.1f}%)")

    return filtered_embeddings, filtered_gene_ids, filtered_genome_ids


def cluster_leiden(embeddings, resolution=1.0, n_neighbors=30, random_state=42):
    """
    Cluster using Leiden algorithm.

    Args:
        embeddings: Input embeddings (N x D)
        resolution: Resolution parameter (higher = more clusters)
        n_neighbors: Number of neighbors for graph construction
        random_state: Random seed

    Returns:
        labels: Cluster labels
    """
    import igraph as ig
    import leidenalg
    from sklearn.neighbors import kneighbors_graph

    print(f"\nClustering with Leiden (resolution={resolution}, n_neighbors={n_neighbors})...")

    # Build k-NN graph
    print("  Building k-NN graph...")
    start = time.time()
    knn_graph = kneighbors_graph(embeddings, n_neighbors=n_neighbors, mode='connectivity', include_self=False)
    print(f"    Time: {time.time() - start:.1f}s")

    # Convert to igraph
    print("  Converting to igraph...")
    start = time.time()
    sources, targets = knn_graph.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))

    g = ig.Graph(n=embeddings.shape[0], edges=edgelist, directed=False)
    g.simplify(multiple=True, loops=True)
    print(f"    Time: {time.time() - start:.1f}s")
    print(f"    Vertices: {g.vcount():,}, Edges: {g.ecount():,}")

    # Run Leiden
    print(f"  Running Leiden algorithm...")
    start = time.time()
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        seed=random_state,
        n_iterations=-1
    )
    print(f"    Time: {time.time() - start:.1f}s")

    labels = np.array(partition.membership)
    n_clusters = len(set(labels))
    print(f"  Found {n_clusters:,} clusters")

    return labels


def evaluate_clusters(labels, gene_ids, cog_lookup):
    """
    Evaluate cluster quality based on COG homogeneity.

    Metrics:
    - Per-cluster annotation rate
    - Per-cluster COG homogeneity (most common COG / total)
    - Overall weighted homogeneity
    """
    from collections import Counter

    print("\nEvaluating cluster quality...")

    # Get COG labels for all genes
    gene_cogs = np.array([cog_lookup.get(gid, None) for gid in gene_ids])

    cluster_metrics = []

    for cluster_id in sorted(set(labels)):
        mask = labels == cluster_id
        cluster_genes = gene_ids[mask]
        cluster_cogs = gene_cogs[mask]

        # Annotation rate
        annotated_mask = cluster_cogs != None
        n_annotated = annotated_mask.sum()
        annotation_rate = n_annotated / len(cluster_genes)

        # COG homogeneity (only among annotated genes)
        if n_annotated > 0:
            annotated_cogs = cluster_cogs[annotated_mask]
            cog_counts = Counter(annotated_cogs)
            most_common_cog, most_common_count = cog_counts.most_common(1)[0]
            homogeneity = most_common_count / n_annotated
        else:
            most_common_cog = None
            homogeneity = 0.0

        cluster_metrics.append({
            'cluster_id': cluster_id,
            'size': len(cluster_genes),
            'n_annotated': n_annotated,
            'annotation_rate': annotation_rate,
            'most_common_cog': most_common_cog,
            'homogeneity': homogeneity
        })

    # Overall metrics (weighted by cluster size)
    total_genes = len(gene_ids)
    total_annotated = (gene_cogs != None).sum()

    weighted_homogeneity = sum(m['homogeneity'] * m['size'] for m in cluster_metrics) / total_genes
    mean_homogeneity = np.mean([m['homogeneity'] for m in cluster_metrics])

    print(f"  Total clusters: {len(cluster_metrics):,}")
    print(f"  Overall annotation rate: {100*total_annotated/total_genes:.1f}%")
    print(f"  Mean COG homogeneity: {mean_homogeneity:.3f}")
    print(f"  Weighted COG homogeneity: {weighted_homogeneity:.3f}")

    return {
        'cluster_metrics': cluster_metrics,
        'mean_homogeneity': mean_homogeneity,
        'weighted_homogeneity': weighted_homogeneity,
        'annotation_rate': total_annotated / total_genes
    }


def main():
    parser = argparse.ArgumentParser(description='Comprehensive Leiden clustering')
    parser.add_argument('--pca-cache', required=True, help='PCA cache file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--resolution', type=float, default=1.0, help='Leiden resolution')
    parser.add_argument('--n-neighbors', type=int, default=30, help='Number of neighbors for graph')
    parser.add_argument('--cog-only', action='store_true', help='Filter to COG-annotated genes only')
    parser.add_argument('--evaluate', action='store_true', help='Evaluate cluster quality')

    args = parser.parse_args()

    print("=" * 80)
    print("Leiden Clustering with COG Evaluation")
    print("=" * 80)
    print(f"Resolution: {args.resolution}")
    print(f"N neighbors: {args.n_neighbors}")
    print(f"COG filter: {args.cog_only}")
    print(f"Evaluate: {args.evaluate}")
    print()

    # Load data
    embeddings, gene_ids, genome_ids = load_pca(args.pca_cache)

    # Load COG annotations if needed
    if args.cog_only or args.evaluate:
        cog_lookup = load_cog_annotations(gene_ids, genome_ids)

    # Filter to annotated genes if requested
    if args.cog_only:
        embeddings, gene_ids, genome_ids = filter_to_annotated(embeddings, gene_ids, genome_ids, cog_lookup)

    # Cluster
    labels = cluster_leiden(embeddings, resolution=args.resolution, n_neighbors=args.n_neighbors)

    # Evaluate if requested
    if args.evaluate:
        evaluation = evaluate_clusters(labels, gene_ids, cog_lookup)
    else:
        evaluation = None

    # Save results
    print(f"\nSaving to {args.output}...")
    output_data = {
        'labels': labels,
        'gene_ids': gene_ids,
        'genome_ids': genome_ids,
        'resolution': args.resolution,
        'n_neighbors': args.n_neighbors,
        'cog_only': args.cog_only,
        'n_clusters': len(set(labels))
    }

    if evaluation:
        output_data['evaluation'] = evaluation

    np.savez_compressed(args.output, **output_data)

    print("\n" + "=" * 80)
    print("✓ Clustering complete!")
    print(f"Clusters: {len(set(labels)):,}")
    if evaluation:
        print(f"Mean COG homogeneity: {evaluation['mean_homogeneity']:.3f}")
        print(f"Weighted COG homogeneity: {evaluation['weighted_homogeneity']:.3f}")
    print("=" * 80)


if __name__ == '__main__':
    main()
