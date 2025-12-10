import torch
import torch.nn as nn
from torch_geometric.nn import HeteroConv, SAGEConv, GATConv, Linear


class HeteroGraphDecoder(nn.Module):
    """Heterogeneous graph decoder for masked node reconstruction (GraphMAE-style).

    Reconstructs original node features from encoded representations.
    Compatible with GraphMAE's re-masking strategy.
    """

    def __init__(
        self,
        genome_out_dim: int,
        gene_out_dim: int,
        num_hidden: int = 256,
        dropout: float = 0.1,
        activation: str = "prelu",
        decoder_type: str = "sage",  # "sage", "gat", "mlp", or "linear"
    ):
        """Initialize the heterogeneous graph decoder.

        Parameters:
            genome_out_dim: Output feature dimension for genome nodes (original dim)
            gene_out_dim: Output feature dimension for gene nodes (original dim)
            num_hidden: Hidden layer dimension (should match encoder)
            dropout: Dropout probability
            activation: Activation function ("prelu", "relu", "gelu")
            decoder_type: Type of decoder ("sage", "gat", "mlp", or "linear")
        """
        super().__init__()

        self.genome_out_dim = genome_out_dim
        self.gene_out_dim = gene_out_dim
        self.num_hidden = num_hidden
        self.decoder_type = decoder_type

        # Activation
        if activation == "prelu":
            self.activation = nn.PReLU()
        elif activation == "relu":
            self.activation = nn.ReLU()
        elif activation == "gelu":
            self.activation = nn.GELU()
        else:
            raise ValueError(f"Unknown activation: {activation}")

        if decoder_type == "sage":
            # Single-layer SAGE decoder
            self.conv = HeteroConv({
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

            # Output projections
            self.genome_output_proj = Linear(num_hidden, genome_out_dim)
            self.gene_output_proj = Linear(num_hidden, gene_out_dim)

        elif decoder_type == "gat":
            # Single-layer GAT decoder
            self.conv = HeteroConv({
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

            # Output projections
            self.genome_output_proj = Linear(num_hidden, genome_out_dim)
            self.gene_output_proj = Linear(num_hidden, gene_out_dim)

        elif decoder_type == "mlp":
            # MLP decoder (no message passing)
            self.genome_decoder = nn.Sequential(
                Linear(num_hidden, num_hidden),
                self.activation,
                nn.Dropout(dropout),
                Linear(num_hidden, genome_out_dim)
            )
            self.gene_decoder = nn.Sequential(
                Linear(num_hidden, num_hidden),
                self.activation,
                nn.Dropout(dropout),
                Linear(num_hidden, gene_out_dim)
            )

        elif decoder_type == "linear":
            # Linear decoder (simplest)
            self.genome_decoder = Linear(num_hidden, genome_out_dim)
            self.gene_decoder = Linear(num_hidden, gene_out_dim)

        else:
            raise ValueError(f"Unknown decoder type: {decoder_type}")

        self.dropout_layer = nn.Dropout(dropout)

    def forward(self, x_dict, edge_index_dict=None):
        """Forward pass through the decoder.

        Parameters:
            x_dict: Dictionary of node embeddings {node_type: embeddings}
            edge_index_dict: Dictionary of edge indices (only needed for GNN decoders)

        Returns:
            Dictionary of reconstructed node features {node_type: features}
        """
        if self.decoder_type in ["sage", "gat"]:
            # GNN-based decoder - requires graph structure
            if edge_index_dict is None:
                raise ValueError("GNN decoder requires edge_index_dict")

            # Apply convolution
            x_dict = self.conv(x_dict, edge_index_dict)
            x_dict = {key: self.activation(x) for key, x in x_dict.items()}
            x_dict = {key: self.dropout_layer(x) for key, x in x_dict.items()}

            # Project to output dimension
            recon_dict = {
                'genome': self.genome_output_proj(x_dict['genome']),
                'gene': self.gene_output_proj(x_dict['gene']),
            }

        elif self.decoder_type in ["mlp", "linear"]:
            # MLP/Linear decoder - no graph structure needed
            recon_dict = {
                'genome': self.genome_decoder(x_dict['genome']),
                'gene': self.gene_decoder(x_dict['gene']),
            }

        return recon_dict
