import torch
import torch.nn as nn
from torch_geometric.nn import Linear


class HeteroEdgeDecoder(nn.Module):
    """Edge decoder for predicting graph structure in heterogeneous graphs.

    Predicts which edges should exist between node pairs based on node embeddings.
    Supports multiple edge types with separate scoring functions.
    """

    def __init__(
        self,
        hidden_dim: int,
        edge_types: list,
        decoder_type: str = "dot",  # "dot", "bilinear", or "mlp"
        dropout: float = 0.1,
    ):
        """Initialize the heterogeneous edge decoder.

        Parameters:
            hidden_dim: Dimension of node embeddings
            edge_types: List of edge types (src_type, relation, dst_type)
            decoder_type: Type of edge scoring ("dot", "bilinear", or "mlp")
            dropout: Dropout probability
        """
        super().__init__()

        self.hidden_dim = hidden_dim
        self.edge_types = edge_types
        self.decoder_type = decoder_type

        # Create scoring functions for each edge type
        self.scorers = nn.ModuleDict()

        for edge_type in edge_types:
            src_type, rel, dst_type = edge_type
            edge_key = f"{src_type}__{rel}__{dst_type}"

            if decoder_type == "dot":
                # Simple dot product (no parameters)
                self.scorers[edge_key] = nn.Identity()

            elif decoder_type == "bilinear":
                # Bilinear scoring: src^T W dst
                self.scorers[edge_key] = nn.Bilinear(hidden_dim, hidden_dim, 1)

            elif decoder_type == "mlp":
                # MLP scoring: MLP([src; dst])
                self.scorers[edge_key] = nn.Sequential(
                    Linear(hidden_dim * 2, hidden_dim),
                    nn.PReLU(),
                    nn.Dropout(dropout),
                    Linear(hidden_dim, 1),
                )
            else:
                raise ValueError(f"Unknown decoder type: {decoder_type}")

    def forward(self, x_dict, edge_index_dict, edge_types_to_decode=None):
        """Predict edge scores for given node embeddings.

        Parameters:
            x_dict: Dictionary of node embeddings {node_type: embeddings}
            edge_index_dict: Dictionary of edge indices to score {edge_type: edge_index}
            edge_types_to_decode: Optional list of edge types to decode (defaults to all)

        Returns:
            Dictionary of edge scores {edge_type: scores} where scores are [num_edges]
        """
        if edge_types_to_decode is None:
            edge_types_to_decode = self.edge_types

        edge_scores = {}

        for edge_type in edge_types_to_decode:
            if edge_type not in edge_index_dict:
                continue

            src_type, rel, dst_type = edge_type
            edge_key = f"{src_type}__{rel}__{dst_type}"

            edge_index = edge_index_dict[edge_type]
            src_idx = edge_index[0]
            dst_idx = edge_index[1]

            src_emb = x_dict[src_type][src_idx]  # [num_edges, hidden_dim]
            dst_emb = x_dict[dst_type][dst_idx]  # [num_edges, hidden_dim]

            if self.decoder_type == "dot":
                # Dot product: sum(src * dst, dim=1)
                scores = (src_emb * dst_emb).sum(dim=1)  # [num_edges]

            elif self.decoder_type == "bilinear":
                # Bilinear: src^T W dst
                scores = self.scorers[edge_key](src_emb, dst_emb).squeeze(-1)  # [num_edges]

            elif self.decoder_type == "mlp":
                # MLP: MLP([src; dst])
                combined = torch.cat([src_emb, dst_emb], dim=1)  # [num_edges, 2*hidden_dim]
                scores = self.scorers[edge_key](combined).squeeze(-1)  # [num_edges]

            edge_scores[edge_type] = scores

        return edge_scores

    def decode_all_pairs(self, x_dict, edge_type, max_pairs=None):
        """Decode scores for all possible node pairs (expensive for large graphs).

        Parameters:
            x_dict: Dictionary of node embeddings
            edge_type: Edge type to decode
            max_pairs: Optional maximum number of pairs to sample (for efficiency)

        Returns:
            edge_index: [2, num_pairs] sampled node pairs
            scores: [num_pairs] edge scores
        """
        src_type, rel, dst_type = edge_type
        edge_key = f"{src_type}__{rel}__{dst_type}"

        num_src = x_dict[src_type].shape[0]
        num_dst = x_dict[dst_type].shape[0]
        total_pairs = num_src * num_dst

        # Sample pairs if too many
        if max_pairs is not None and total_pairs > max_pairs:
            # Random sampling
            sampled_indices = torch.randperm(total_pairs, device=x_dict[src_type].device)[:max_pairs]
            src_idx = sampled_indices // num_dst
            dst_idx = sampled_indices % num_dst
        else:
            # All pairs
            src_idx = torch.arange(num_src, device=x_dict[src_type].device).repeat_interleave(num_dst)
            dst_idx = torch.arange(num_dst, device=x_dict[dst_type].device).repeat(num_src)

        edge_index = torch.stack([src_idx, dst_idx])

        # Score the pairs
        src_emb = x_dict[src_type][src_idx]
        dst_emb = x_dict[dst_type][dst_idx]

        if self.decoder_type == "dot":
            scores = (src_emb * dst_emb).sum(dim=1)
        elif self.decoder_type == "bilinear":
            scores = self.scorers[edge_key](src_emb, dst_emb).squeeze(-1)
        elif self.decoder_type == "mlp":
            combined = torch.cat([src_emb, dst_emb], dim=1)
            scores = self.scorers[edge_key](combined).squeeze(-1)

        return edge_index, scores
