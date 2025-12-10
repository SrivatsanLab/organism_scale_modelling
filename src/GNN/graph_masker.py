
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

