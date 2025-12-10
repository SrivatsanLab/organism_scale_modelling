import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import torch
import torch_geometric as tg

def _decompress_distance_matrix(matrix_fname: str) -> np.ndarray:
    """Decompress a distance matrix stored in a compressed format.

    Args:
        matrix_fname (str): Path to the compressed distance matrix file.
    Returns:
        np.ndarray: The decompressed distance matrix.

    The compressed format stores only the upper triangular part of the distance matrix
    to save space. This function reconstructs the full symmetric distance matrix.
    """
    data = np.load(matrix_fname)
    upper = data['upper_triangle'].astype(np.float32)
    n = int(data['n'])
    dist_matrix = np.zeros((n, n), dtype=np.float32)
    triu_idx = np.triu_indices(n, k=1)
    dist_matrix[triu_idx] = upper
    dist_matrix = dist_matrix + dist_matrix.T  # Symmetrize
    return dist_matrix

def construct_hom_graph(
    node_features_fname: str,
    adjacency_matrix_fname: str,
    distance_matrix_fname: str,
    save_path: str | None = None,
    return_graph: bool = True,
):
    """Construct a `torch_geometric.data.Data` graph for a single genome.

    Parameters:
        node_features_fname: Path to compressed sparse `.npz` file containing node feature matrix.
        adjacency_matrix_fname: Path to compressed sparse `.npz` file containing adjacency matrix.
        distance_matrix_fname: Path to compressed distance matrix file (upper triangle + n).
        save_path: Optional path to save the constructed graph (e.g. `graph.pt` or directory).
        return_graph: If True, return the in-memory graph object.

    Returns:
        torch_geometric.data.Data | None: The constructed graph if `return_graph` is True.
    """
    node_features = load_npz(node_features_fname)
    adjacency_matrix = load_npz(adjacency_matrix_fname)
    dist_matrix = _decompress_distance_matrix(distance_matrix_fname)

    edge_index = tg.utils.from_scipy_sparse_matrix(adjacency_matrix)[0]
    row, col = adjacency_matrix.nonzero()
    edge_features = dist_matrix[row, col].reshape(-1, 1)

    single_graph = tg.data.Data(
        x=torch.tensor(node_features.toarray(), dtype=torch.float32),
        edge_index=edge_index,
        edge_attr=torch.tensor(edge_features, dtype=torch.float32)
    )
    single_graph.num_nodes = node_features.shape[0]
    single_graph.validate(raise_on_error=True)

    if save_path is not None:
        tg.io.write_torch_geometric_data(single_graph, save_path)
    if return_graph:
        return single_graph