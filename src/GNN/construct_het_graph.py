import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import torch
import torch_geometric as tg

def _unpack_upper_strict_rowmajor(flat):
    m = len(flat)
    # Solve m = n(n-1)/2 → n = (1 + sqrt(1+8m))/2
    n = int((np.sqrt(1 + 8*m) + 1) / 2)
    if n * (n - 1) // 2 != m:
        raise ValueError("length is not n(n-1)/2 for any integer n")

    M = np.zeros((n, n), dtype=np.array(flat).dtype)
    r, c = np.triu_indices(n, k=1)     # k=1 → strictly upper triangular
    M[r, c] = flat
    return M

def construct_het_graph(
    genome_features: str,
    gene_features: str,
    genome_genome_edges: str,
    genome_genome_edge_features: str,
    genome_gene_edges: str,
    gene_gene_edges: str,
    save_path: str | None = None,
    return_graph: bool = True,
):
    data = tg.data.HeteroData()

    # Load node features
    data['genome'].x = torch.tensor(load_npz(genome_features).toarray(), dtype=torch.float32)
    data['gene'].x = torch.tensor(np.load(gene_features)['embeddings'], dtype=torch.float32)

    # Load edge indices
    genome_genome_adj = load_npz(genome_genome_edges)
    data['genome', 'similarity', 'genome'].edge_index = tg.utils.from_scipy_sparse_matrix(
        genome_genome_adj
    )[0]

    genome_gene_edge_data = np.load(genome_gene_edges)
    n_genome_gene_edges = len(genome_gene_edge_data['indptr']) - 1
    n_genes_per_genome = genome_gene_edge_data['indptr'][1:] - genome_gene_edge_data['indptr'][:-1]
    genome_gene_src = torch.repeat_interleave(torch.arange(n_genome_gene_edges), torch.tensor(n_genes_per_genome))
    genome_gene_dest = torch.tensor(genome_gene_edge_data['indices'], dtype=torch.long)
    data['genome', 'present', 'gene'].edge_index = torch.stack([genome_gene_src, genome_gene_dest], dim=0)

    gene_gene_edge_data = np.load(gene_gene_edges)
    n_gene_gene_edges = len(gene_gene_edge_data['indptr']) - 1
    n_genes_per_gene = gene_gene_edge_data['indptr'][1:] - gene_gene_edge_data['indptr'][:-1]
    gene_gene_src = torch.repeat_interleave(torch.arange(n_gene_gene_edges), torch.tensor(n_genes_per_gene))
    gene_gene_dest = torch.tensor(gene_gene_edge_data['indices'], dtype=torch.long)
    data['gene', 'interacts', 'gene'].edge_index = torch.stack([gene_gene_src, gene_gene_dest], dim=0)

    # Load edge features, only for genome-genome edges
    genome_genome_edge_feats = np.load(genome_genome_edge_features)['upper_triangle'].astype(np.float32)
    genome_genome_edge_feats = _unpack_upper_strict_rowmajor(genome_genome_edge_feats)
    genome_genome_edge_feats = genome_genome_edge_feats + genome_genome_edge_feats.T  # Symmetrize
    row, col = genome_genome_adj.nonzero()
    genome_genome_edge_feats = genome_genome_edge_feats[row, col]
    data['genome', 'similarity', 'genome'].edge_attr = torch.tensor(
        genome_genome_edge_feats, dtype=torch.float32
    )

    # Set num nodes
    data['genome'].num_nodes = data['genome'].x.shape[0]
    data['gene'].num_nodes = data['gene'].x.shape[0]
    
    data.validate(raise_on_error=True)
    if save_path is not None:
        save_het_graph(data, save_path)

    if return_graph:
        return data

def save_het_graph(
    data: tg.data.HeteroData,
    save_path: str,
):
    """Save a `torch_geometric.data.HeteroData` graph to disk.

    Parameters:
        data: The heterogeneous graph to save.
        save_path: Path to save the constructed graph (e.g. `graph.pt` or directory).
    """
    torch.save(data.to_dict(), save_path)

def load_het_graph(
    load_path: str,
) -> tg.data.HeteroData:
    """Load a `torch_geometric.data.HeteroData` graph from disk.

    Parameters:
        load_path: Path to load the constructed graph (e.g. `graph.pt` or directory).

    Returns:
        torch_geometric.data.HeteroData: The loaded heterogeneous graph.
    """
    return tg.data.HeteroData.from_dict(torch.load(load_path, weights_only=False))
