import torch
import torch.nn as nn
from typing import Optional

from .encoder import HeteroGraphEncoder
from .decoder import HeteroGraphDecoder
from .loss import hetero_sce_loss, hetero_mse_loss


class HeteroGraphMAE(nn.Module):
    """Heterogeneous Graph Masked Autoencoder (adapted from GraphMAE).

    Implements masked autoencoding for heterogeneous graphs with genome and gene nodes.
    Compatible with the mask_nodes() function from graph_masker.py.
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
        encoder_type: str = "sage",
        decoder_type: str = "sage",
        loss_fn: str = "sce",
        alpha_l: float = 3.0,
        concat_hidden: bool = False,
        remask_rate: float = 0.5,
    ):
        """Initialize HeteroGraphMAE.

        Parameters:
            genome_in_dim: Input feature dimension for genome nodes
            gene_in_dim: Input feature dimension for gene nodes
            num_hidden: Hidden layer dimension
            num_layers: Number of encoder layers
            dropout: Dropout probability
            activation: Activation function ("prelu", "relu", "gelu")
            norm: Normalization type ("batchnorm", "layernorm", or None)
            encoder_type: Type of encoder ("sage" or "gat")
            decoder_type: Type of decoder ("sage", "gat", "mlp", or "linear")
            loss_fn: Loss function ("sce" or "mse")
            alpha_l: Alpha parameter for SCE loss
            concat_hidden: Whether to concatenate all hidden layers
            remask_rate: Fraction of masked nodes to re-mask before decoder (0-1)
        """
        super().__init__()

        self.genome_in_dim = genome_in_dim
        self.gene_in_dim = gene_in_dim
        self.num_hidden = num_hidden
        self.num_layers = num_layers
        self.decoder_type = decoder_type
        self.loss_fn = loss_fn
        self.alpha_l = alpha_l
        self.concat_hidden = concat_hidden
        self.remask_rate = remask_rate

        # Build encoder
        self.encoder = HeteroGraphEncoder(
            genome_in_dim=genome_in_dim,
            gene_in_dim=gene_in_dim,
            num_hidden=num_hidden,
            num_layers=num_layers,
            dropout=dropout,
            activation=activation,
            norm=norm,
            encoder_type=encoder_type,
        )

        # Build decoder
        self.decoder = HeteroGraphDecoder(
            genome_out_dim=genome_in_dim,
            gene_out_dim=gene_in_dim,
            num_hidden=num_hidden,
            dropout=dropout,
            activation=activation,
            decoder_type=decoder_type,
        )

        # Encoder-to-decoder projection
        if concat_hidden:
            self.encoder_to_decoder = nn.Linear(num_hidden * num_layers, num_hidden, bias=False)
        else:
            self.encoder_to_decoder = nn.Linear(num_hidden, num_hidden, bias=False)

        # Learnable mask tokens for each node type
        self.genome_mask_token = nn.Parameter(torch.zeros(1, genome_in_dim))
        self.gene_mask_token = nn.Parameter(torch.zeros(1, gene_in_dim))

        # Setup loss function
        if loss_fn == "sce":
            self.criterion = lambda pred, target, mask: hetero_sce_loss(
                pred, target, mask, alpha=alpha_l
            )
        elif loss_fn == "mse":
            self.criterion = lambda pred, target, mask: hetero_mse_loss(pred, target, mask)
        else:
            raise ValueError(f"Unknown loss function: {loss_fn}")

    def mask_features(self, x_dict, masked_genome_idx, masked_gene_idx):
        """Apply mask tokens to masked nodes.

        Parameters:
            x_dict: Dictionary of node features {node_type: features}
            masked_genome_idx: Tensor of masked genome node indices
            masked_gene_idx: Tensor of masked gene node indices

        Returns:
            Dictionary of masked features
        """
        masked_x_dict = {
            'genome': x_dict['genome'].clone(),
            'gene': x_dict['gene'].clone(),
        }

        # Apply mask tokens
        if len(masked_genome_idx) > 0:
            masked_x_dict['genome'][masked_genome_idx] = self.genome_mask_token
        if len(masked_gene_idx) > 0:
            masked_x_dict['gene'][masked_gene_idx] = self.gene_mask_token

        return masked_x_dict

    def remask_encoded(self, enc_dict, masked_genome_idx, masked_gene_idx):
        """Re-mask encoded representations before decoder (GraphMAE strategy).

        Sets masked node embeddings to zero before passing to decoder.

        Parameters:
            enc_dict: Dictionary of encoded features
            masked_genome_idx: Tensor of masked genome node indices
            masked_gene_idx: Tensor of masked gene node indices

        Returns:
            Dictionary of re-masked encoded features
        """
        if self.remask_rate == 0.0:
            return enc_dict

        remasked_enc_dict = {
            'genome': enc_dict['genome'].clone(),
            'gene': enc_dict['gene'].clone(),
        }

        # Re-mask a subset of masked nodes
        if len(masked_genome_idx) > 0:
            num_remask = int(self.remask_rate * len(masked_genome_idx))
            if num_remask > 0:
                remask_genome_idx = masked_genome_idx[torch.randperm(len(masked_genome_idx))[:num_remask]]
                remasked_enc_dict['genome'][remask_genome_idx] = 0

        if len(masked_gene_idx) > 0:
            num_remask = int(self.remask_rate * len(masked_gene_idx))
            if num_remask > 0:
                remask_gene_idx = masked_gene_idx[torch.randperm(len(masked_gene_idx))[:num_remask]]
                remasked_enc_dict['gene'][remask_gene_idx] = 0

        return remasked_enc_dict

    def forward(self, data, x_dict, masked_genome_idx, masked_gene_idx):
        """Forward pass for training.

        Parameters:
            data: HeteroData graph (unmasked, for reconstruction targets)
            x_dict: Dictionary of node features (already masked by mask_nodes())
            masked_genome_idx: Tensor of masked genome node indices
            masked_gene_idx: Tensor of masked gene node indices

        Returns:
            loss: Scalar loss value
            loss_dict: Dictionary of loss components
        """
        # Store original features as targets
        target_x_dict = {
            'genome': data['genome'].x.clone(),
            'gene': data['gene'].x.clone(),
        }

        # Encode masked graph
        if self.concat_hidden:
            enc_rep, all_hidden = self.encoder(x_dict, data.edge_index_dict, return_hidden=True)
            # Concatenate all hidden layers
            enc_rep = {
                node_type: torch.cat([h[node_type] for h in all_hidden], dim=1)
                for node_type in enc_rep.keys()
            }
        else:
            enc_rep = self.encoder(x_dict, data.edge_index_dict, return_hidden=False)

        # Project encoder output to decoder input
        enc_rep = {
            node_type: self.encoder_to_decoder(feat)
            for node_type, feat in enc_rep.items()
        }

        # Re-mask encoded representations (GraphMAE strategy)
        if self.decoder_type not in ["mlp", "linear"]:
            enc_rep = self.remask_encoded(enc_rep, masked_genome_idx, masked_gene_idx)

        # Decode to reconstruct features
        if self.decoder_type in ["mlp", "linear"]:
            recon_x_dict = self.decoder(enc_rep)
        else:
            recon_x_dict = self.decoder(enc_rep, data.edge_index_dict)

        # Create mask dictionary for loss computation
        mask_dict = {
            'genome': torch.zeros(data['genome'].num_nodes, dtype=torch.bool, device=x_dict['genome'].device),
            'gene': torch.zeros(data['gene'].num_nodes, dtype=torch.bool, device=x_dict['gene'].device),
        }
        if len(masked_genome_idx) > 0:
            mask_dict['genome'][masked_genome_idx] = True
        if len(masked_gene_idx) > 0:
            mask_dict['gene'][masked_gene_idx] = True

        # Compute loss only on masked nodes
        loss = self.criterion(recon_x_dict, target_x_dict, mask_dict)

        loss_dict = {
            'loss': loss.item(),
            'num_masked_genome': len(masked_genome_idx),
            'num_masked_gene': len(masked_gene_idx),
        }

        return loss, loss_dict

    def embed(self, data, x_dict):
        """Generate embeddings for the graph (no masking).

        Parameters:
            data: HeteroData graph
            x_dict: Dictionary of node features

        Returns:
            Dictionary of node embeddings
        """
        self.eval()
        with torch.no_grad():
            enc_rep = self.encoder(x_dict, data.edge_index_dict, return_hidden=False)
        return enc_rep

    @property
    def enc_params(self):
        """Return encoder parameters."""
        return self.encoder.parameters()

    @property
    def dec_params(self):
        """Return decoder parameters."""
        from itertools import chain
        return chain(self.encoder_to_decoder.parameters(), self.decoder.parameters())
