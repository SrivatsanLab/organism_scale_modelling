import torch
import torch.nn as nn
from torch_geometric.nn import HeteroConv, SAGEConv, GATConv, Linear


class HeteroGraphEncoder(nn.Module):
    """Heterogeneous graph encoder for genome-gene graphs (GraphMAE-style).

    Encodes a heterogeneous graph with:
    - Node types: 'genome', 'gene'
    - Edge types: ('genome', 'similarity', 'genome'),
                  ('genome', 'present', 'gene'),
                  ('gene', 'interacts', 'gene')

    Compatible with GraphMAE masking strategy.
    """

    def __init__(
        self,
        genome_in_dim: int,
        gene_in_dim: int,
        num_hidden: int = 256,
        num_layers: int = 2,
        dropout: float = 0.1,
        activation: str = "prelu",
        norm: str = "batchnorm",
        encoder_type: str = "sage",  # "sage" or "gat"
    ):
        """Initialize the heterogeneous graph encoder.

        Parameters:
            genome_in_dim: Input feature dimension for genome nodes
            gene_in_dim: Input feature dimension for gene nodes
            num_hidden: Hidden layer dimension
            num_layers: Number of message passing layers
            dropout: Dropout probability
            activation: Activation function ("prelu", "relu", "gelu")
            norm: Normalization type ("batchnorm", "layernorm", or None)
            encoder_type: Type of GNN layer ("sage" or "gat")
        """
        super().__init__()

        self.genome_in_dim = genome_in_dim
        self.gene_in_dim = gene_in_dim
        self.num_hidden = num_hidden
        self.num_layers = num_layers
        self.dropout = dropout
        self.encoder_type = encoder_type

        # Input projections for each node type
        self.genome_input_proj = Linear(genome_in_dim, num_hidden)
        self.gene_input_proj = Linear(gene_in_dim, num_hidden)

        # Activation
        if activation == "prelu":
            self.activation = nn.PReLU()
        elif activation == "relu":
            self.activation = nn.ReLU()
        elif activation == "gelu":
            self.activation = nn.GELU()
        else:
            raise ValueError(f"Unknown activation: {activation}")

        # Normalization
        self.norm_layers = nn.ModuleList()
        if norm == "batchnorm":
            for _ in range(num_layers):
                self.norm_layers.append(nn.ModuleDict({
                    'genome': nn.BatchNorm1d(num_hidden),
                    'gene': nn.BatchNorm1d(num_hidden),
                }))
        elif norm == "layernorm":
            for _ in range(num_layers):
                self.norm_layers.append(nn.ModuleDict({
                    'genome': nn.LayerNorm(num_hidden),
                    'gene': nn.LayerNorm(num_hidden),
                }))
        else:
            self.norm_layers = None

        # Heterogeneous convolution layers
        self.convs = nn.ModuleList()
        for i in range(num_layers):
            if encoder_type == "sage":
                conv = HeteroConv({
                    ('genome', 'similarity', 'genome'): SAGEConv(
                        num_hidden, num_hidden, normalize=True
                    ),
                    ('genome', 'present', 'gene'): SAGEConv(
                        (num_hidden, num_hidden), num_hidden, normalize=True
                    ),
                    ('gene', 'interacts', 'gene'): SAGEConv(
                        num_hidden, num_hidden, normalize=True
                    ),
                }, aggr='sum')
            elif encoder_type == "gat":
                conv = HeteroConv({
                    ('genome', 'similarity', 'genome'): GATConv(
                        num_hidden, num_hidden // 4, heads=4, concat=True, dropout=dropout
                    ),
                    ('genome', 'present', 'gene'): GATConv(
                        (num_hidden, num_hidden), num_hidden // 4, heads=4, concat=True, dropout=dropout
                    ),
                    ('gene', 'interacts', 'gene'): GATConv(
                        num_hidden, num_hidden // 4, heads=4, concat=True, dropout=dropout
                    ),
                }, aggr='sum')
            else:
                raise ValueError(f"Unknown encoder type: {encoder_type}")
            self.convs.append(conv)

        self.dropout_layer = nn.Dropout(dropout)

    def forward(self, x_dict, edge_index_dict, return_hidden=False):
        """Forward pass through the encoder.

        Parameters:
            x_dict: Dictionary of node features {node_type: features}
            edge_index_dict: Dictionary of edge indices {edge_type: edge_index}
            return_hidden: If True, return all hidden representations

        Returns:
            Final node embeddings dict, and optionally list of all hidden states
        """
        # Project inputs to hidden dimension
        x_dict = {
            'genome': self.genome_input_proj(x_dict['genome']),
            'gene': self.gene_input_proj(x_dict['gene']),
        }

        all_hidden = []

        # Apply message passing layers
        for i, conv in enumerate(self.convs):
            x_dict = conv(x_dict, edge_index_dict)

            # Apply normalization
            if self.norm_layers is not None:
                x_dict = {key: self.norm_layers[i][key](x) for key, x in x_dict.items()}

            # Apply activation and dropout
            x_dict = {key: self.activation(x) for key, x in x_dict.items()}
            x_dict = {key: self.dropout_layer(x) for key, x in x_dict.items()}

            if return_hidden:
                all_hidden.append(x_dict)

        if return_hidden:
            return x_dict, all_hidden
        return x_dict
