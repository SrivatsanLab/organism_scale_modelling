
import random
from typing import List, Optional


def sample_mask_indices(
    num_edges: int,
    mask_ratio: float,
    seed: Optional[int] = None,
) -> List[int]:
    """
    Sample indices of edges to mask, given a total number of edges and a mask ratio.

    Parameters
    ----------
    num_edges:
        Total number of edges E in the graph. If E <= 0, returns an empty list.
    mask_ratio:
        Fraction of edges to mask (0.0–1.0). Values outside this range are
        clipped into [0.0, 1.0].
    seed:
        Optional random seed for reproducibility.

    Returns
    -------
    mask_indices:
        A sorted list of unique edge indices in [0, E-1] that have been selected
        to be masked.
    """
    if num_edges <= 0:
        return []

    # Clip mask_ratio to [0, 1]
    ratio = max(0.0, min(1.0, float(mask_ratio)))

    if ratio == 0.0:
        return []

    num_to_mask = int(round(ratio * num_edges))
    if num_to_mask <= 0:
        return []

    rng = random.Random(seed)
    indices = list(range(num_edges))
    rng.shuffle(indices)

    chosen = indices[:num_to_mask]
    chosen.sort()
    return chosen

def sample_mask_node_indices(
    num_nodes: int,
    mask_ratio: float,
    seed: Optional[int] = None,
):
    if num_nodes <= 0:
        return []
    r = max(0.0, min(1.0, float(mask_ratio)))
    if r == 0:
        return []
    n = int(round(num_nodes * r))
    if n <= 0:
        return []
    rng = random.Random(seed)
    idx = list(range(num_nodes))
    rng.shuffle(idx)
    selected = idx[:n]
    selected.sort()
    return selected


import torch
from torch_geometric.data import HeteroData


def _mask_one_type(
    data: HeteroData,
    node_type: str,
    mask_ratio: float,
    seed: Optional[int],
):
    """
    Internal helper for node masking:
    returns (kept_mask, kept_idx, new_index, masked_idx).
    """
    n = data[node_type].num_nodes

    # select nodes to mask
    masked_idx = torch.tensor(
        sample_mask_node_indices(n, mask_ratio, seed),
        dtype=torch.long,
    )

    kept_mask = torch.ones(n, dtype=torch.bool)
    kept_mask[masked_idx] = False
    kept_idx = kept_mask.nonzero(as_tuple=False).view(-1)

    new_index = -torch.ones(n, dtype=torch.long)
    new_index[kept_idx] = torch.arange(kept_idx.size(0))

    return kept_mask, kept_idx, new_index, masked_idx

def mask_nodes(
    data: HeteroData,
    mask_gene_ratio: float,
    mask_genome_ratio: float = 0.0,
    seed: Optional[int] = None,
):
    """
    Randomly mask (delete) gene nodes (and optionally genome nodes) from a
    heterogeneous graph. All edges touching masked nodes are removed.

    Returns a new HeteroData and metadata describing which nodes were masked.
    """
    masked = HeteroData()

    # --- Mask gene nodes ---
    k_gene, idx_gene, map_gene, m_gene = _mask_one_type(
        data, 'gene', mask_gene_ratio, seed
    )
    masked['gene'].x = data['gene'].x[idx_gene]
    masked['gene'].num_nodes = idx_gene.size(0)

    # --- Mask genome nodes ---
    if 'genome' in data.node_types:
        k_genome, idx_genome, map_genome, m_genome = _mask_one_type(
            data, 'genome', mask_genome_ratio, seed
        )
        masked['genome'].x = data['genome'].x[idx_genome]
        masked['genome'].num_nodes = idx_genome.size(0)
    else:
        k_genome = torch.tensor([], dtype=torch.bool)
        m_genome = torch.tensor([], dtype=torch.long)
        map_genome = None

    # --- Rebuild edges ---
    for edge_type in data.edge_types:
        s_type, _, d_type = edge_type
        src = data[edge_type].edge_index[0]
        dst = data[edge_type].edge_index[1]

        # source
        if s_type == 'gene':
            s_kept = k_gene[src]
            s_new = map_gene[src[s_kept]]
        else:
            s_kept = k_genome[src]
            s_new = map_genome[src[s_kept]]

        # destination
        if d_type == 'gene':
            d_kept = k_gene[dst]
            d_new = map_gene[dst[d_kept]]
        else:
            d_kept = k_genome[dst]
            d_new = map_genome[dst[d_kept]]

        # both endpoints must be kept
        keep = s_kept & d_kept
        new_src = s_new[keep[s_kept]]
        new_dst = d_new[keep[d_kept]]

        masked[edge_type].edge_index = torch.stack([new_src, new_dst])

    return masked, dict(
        masked_gene_idx=m_gene,
        masked_genome_idx=m_genome,
        kept_gene_mask=k_gene,
        kept_genome_mask=k_genome,
    )


def mask_features(
    data: HeteroData,
    mask_gene_ratio: float,
    mask_genome_ratio: float = 0.0,
    mask_token_gene: Optional[torch.Tensor] = None,
    mask_token_genome: Optional[torch.Tensor] = None,
    seed: Optional[int] = None,
):
    """
    Mask node features (GraphMAE-style) without removing nodes from the graph.
    Graph structure remains intact - only features are masked with mask tokens.

    Parameters
    ----------
    data:
        HeteroData graph
    mask_gene_ratio:
        Fraction of gene nodes to mask (0.0-1.0)
    mask_genome_ratio:
        Fraction of genome nodes to mask (0.0-1.0)
    mask_token_gene:
        Learnable mask token for gene features (shape: [1, gene_feat_dim])
        If None, uses zeros
    mask_token_genome:
        Learnable mask token for genome features (shape: [1, genome_feat_dim])
        If None, uses zeros
    seed:
        Optional random seed for reproducibility

    Returns
    -------
    masked_data:
        HeteroData with same structure but masked features
    mask_info:
        Dictionary with masked indices and masks for each node type
    """
    # Clone the data to avoid modifying original
    masked_data = data.clone()

    # Sample nodes to mask
    num_gene = data['gene'].num_nodes
    num_genome = data['genome'].num_nodes if 'genome' in data.node_types else 0

    masked_gene_idx = torch.tensor(
        sample_mask_node_indices(num_gene, mask_gene_ratio, seed),
        dtype=torch.long
    )

    if num_genome > 0:
        masked_genome_idx = torch.tensor(
            sample_mask_node_indices(num_genome, mask_genome_ratio, seed),
            dtype=torch.long
        )
    else:
        masked_genome_idx = torch.tensor([], dtype=torch.long)

    # Create masks
    gene_mask = torch.zeros(num_gene, dtype=torch.bool)
    if len(masked_gene_idx) > 0:
        gene_mask[masked_gene_idx] = True

    genome_mask = torch.zeros(num_genome, dtype=torch.bool) if num_genome > 0 else torch.tensor([], dtype=torch.bool)
    if len(masked_genome_idx) > 0:
        genome_mask[masked_genome_idx] = True

    # Apply mask tokens
    if len(masked_gene_idx) > 0:
        if mask_token_gene is not None:
            masked_data['gene'].x[masked_gene_idx] = mask_token_gene
        else:
            masked_data['gene'].x[masked_gene_idx] = 0.0

    if len(masked_genome_idx) > 0:
        if mask_token_genome is not None:
            masked_data['genome'].x[masked_genome_idx] = mask_token_genome
        else:
            masked_data['genome'].x[masked_genome_idx] = 0.0

    return masked_data, dict(
        masked_gene_idx=masked_gene_idx,
        masked_genome_idx=masked_genome_idx,
        gene_mask=gene_mask,
        genome_mask=genome_mask,
    )


def mask_edges(
    data: HeteroData,
    mask_ratio_dict: dict,
    seed: Optional[int] = None,
):
    """
    Mask edges by removing a fraction of edges from the graph.
    Nodes and their features remain intact.

    Parameters
    ----------
    data:
        HeteroData graph
    mask_ratio_dict:
        Dictionary mapping edge_type to mask ratio
        Example: {('gene', 'interacts', 'gene'): 0.2, ...}
    seed:
        Optional random seed for reproducibility

    Returns
    -------
    masked_data:
        HeteroData with masked edges removed
    mask_info:
        Dictionary with masked edge indices for each edge type
    """
    # Clone the data
    masked_data = data.clone()
    mask_info = {}

    for edge_type, mask_ratio in mask_ratio_dict.items():
        if edge_type not in data.edge_types:
            continue

        edge_index = data[edge_type].edge_index
        num_edges = edge_index.shape[1]

        # Sample edges to mask
        masked_edge_idx = torch.tensor(
            sample_mask_indices(num_edges, mask_ratio, seed),
            dtype=torch.long
        )

        # Create mask for keeping edges
        keep_mask = torch.ones(num_edges, dtype=torch.bool)
        if len(masked_edge_idx) > 0:
            keep_mask[masked_edge_idx] = False

        # Keep only non-masked edges
        masked_data[edge_type].edge_index = edge_index[:, keep_mask]

        # Store masked edge indices and original edges
        mask_info[edge_type] = {
            'masked_edge_idx': masked_edge_idx,
            'masked_edges': edge_index[:, masked_edge_idx] if len(masked_edge_idx) > 0 else torch.tensor([[], []], dtype=torch.long),
            'num_masked': len(masked_edge_idx),
            'num_total': num_edges,
        }

    return masked_data, mask_info

