import torch
import torch.nn as nn
from typing import Optional

from .encoder import HeteroGraphEncoder
from .decoder import HeteroGraphDecoder
from .edge_decoder import HeteroEdgeDecoder
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
        reconstruction_mode: str = "feature",  # "feature", "structure", or "joint"
        edge_decoder_type: str = "dot",  # "dot", "bilinear", or "mlp"
        edge_loss_weight: float = 1.0,  # Weight for edge reconstruction loss
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
            reconstruction_mode: Reconstruction task ("feature", "structure", or "joint")
            edge_decoder_type: Edge decoder type ("dot", "bilinear", or "mlp")
            edge_loss_weight: Weight for edge reconstruction loss in joint mode
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
        self.reconstruction_mode = reconstruction_mode
        self.edge_decoder_type = edge_decoder_type
        self.edge_loss_weight = edge_loss_weight

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

        # Build feature decoder
        self.decoder = HeteroGraphDecoder(
            genome_out_dim=genome_in_dim,
            gene_out_dim=gene_in_dim,
            num_hidden=num_hidden,
            dropout=dropout,
            activation=activation,
            decoder_type=decoder_type,
        )

        # Build edge decoder (for structure reconstruction)
        if reconstruction_mode in ["structure", "joint"]:
            # Define edge types for the heterogeneous graph
            edge_types = [
                ('genome', 'similarity', 'genome'),
                ('genome', 'present', 'gene'),
                ('gene', 'interacts', 'gene'),
            ]
            self.edge_decoder = HeteroEdgeDecoder(
                hidden_dim=num_hidden,
                edge_types=edge_types,
                decoder_type=edge_decoder_type,
                dropout=dropout,
            )
        else:
            self.edge_decoder = None

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
        # Handle both tensors and sparse matrices
        if hasattr(x_dict['genome'], 'clone'):
            genome_x = x_dict['genome'].clone()
        else:
            if hasattr(x_dict['genome'], 'toarray'):
                genome_x = torch.tensor(x_dict['genome'].toarray(), dtype=torch.float32)
            else:
                genome_x = torch.tensor(x_dict['genome'], dtype=torch.float32)

        if hasattr(x_dict['gene'], 'clone'):
            gene_x = x_dict['gene'].clone()
        else:
            if hasattr(x_dict['gene'], 'toarray'):
                gene_x = torch.tensor(x_dict['gene'].toarray(), dtype=torch.float32)
            else:
                gene_x = torch.tensor(x_dict['gene'], dtype=torch.float32)

        masked_x_dict = {
            'genome': genome_x,
            'gene': gene_x,
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

        # Handle both tensors and sparse matrices
        if hasattr(enc_dict['genome'], 'clone'):
            genome_enc = enc_dict['genome'].clone()
        else:
            if hasattr(enc_dict['genome'], 'toarray'):
                genome_enc = torch.tensor(enc_dict['genome'].toarray(), dtype=torch.float32)
            else:
                genome_enc = torch.tensor(enc_dict['genome'], dtype=torch.float32)

        if hasattr(enc_dict['gene'], 'clone'):
            gene_enc = enc_dict['gene'].clone()
        else:
            if hasattr(enc_dict['gene'], 'toarray'):
                gene_enc = torch.tensor(enc_dict['gene'].toarray(), dtype=torch.float32)
            else:
                gene_enc = torch.tensor(enc_dict['gene'], dtype=torch.float32)

        remasked_enc_dict = {
            'genome': genome_enc,
            'gene': gene_enc,
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

    def forward(self, data, x_dict, masked_genome_idx, masked_gene_idx, edge_mask_info=None):
        """Forward pass for training.

        Parameters:
            data: HeteroData graph (original unmasked, for reconstruction targets)
            x_dict: Dictionary of node features (already masked)
            masked_genome_idx: Tensor of masked genome node indices
            masked_gene_idx: Tensor of masked gene node indices
            edge_mask_info: Dictionary with masked edge information (for structure reconstruction)

        Returns:
            loss: Scalar loss value
            loss_dict: Dictionary of loss components
        """
        if self.reconstruction_mode == "feature":
            return self._forward_feature(data, x_dict, masked_genome_idx, masked_gene_idx)
        elif self.reconstruction_mode == "structure":
            return self._forward_structure(data, x_dict, edge_mask_info)
        elif self.reconstruction_mode == "joint":
            return self._forward_joint(data, x_dict, masked_genome_idx, masked_gene_idx, edge_mask_info)
        else:
            raise ValueError(f"Unknown reconstruction mode: {self.reconstruction_mode}")

    def _forward_feature(self, data, x_dict, masked_genome_idx, masked_gene_idx):
        """Forward pass for feature reconstruction only."""
        # Store original features as targets (handle both tensors and sparse matrices)
        if hasattr(data['genome'].x, 'clone'):
            genome_target = data['genome'].x.clone()
        else:
            if hasattr(data['genome'].x, 'toarray'):
                genome_target = torch.tensor(data['genome'].x.toarray(), dtype=torch.float32, device=x_dict['genome'].device)
            else:
                genome_target = torch.tensor(data['genome'].x, dtype=torch.float32, device=x_dict['genome'].device)

        if hasattr(data['gene'].x, 'clone'):
            gene_target = data['gene'].x.clone()
        else:
            if hasattr(data['gene'].x, 'toarray'):
                gene_target = torch.tensor(data['gene'].x.toarray(), dtype=torch.float32, device=x_dict['gene'].device)
            else:
                gene_target = torch.tensor(data['gene'].x, dtype=torch.float32, device=x_dict['gene'].device)

        target_x_dict = {
            'genome': genome_target,
            'gene': gene_target,
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

    def _forward_structure(self, data, x_dict, edge_mask_info):
        """Forward pass for structure reconstruction only."""
        if edge_mask_info is None:
            raise ValueError("edge_mask_info required for structure reconstruction")

        # Encode graph (using masked edges)
        if self.concat_hidden:
            enc_rep, all_hidden = self.encoder(x_dict, data.edge_index_dict, return_hidden=True)
            enc_rep = {
                node_type: torch.cat([h[node_type] for h in all_hidden], dim=1)
                for node_type in enc_rep.keys()
            }
        else:
            enc_rep = self.encoder(x_dict, data.edge_index_dict, return_hidden=False)

        # Project encoder output
        enc_rep = {
            node_type: self.encoder_to_decoder(feat)
            for node_type, feat in enc_rep.items()
        }

        # Extract masked edges from edge_mask_info for edge decoder
        masked_edges_dict = {
            edge_type: info['masked_edges']
            for edge_type, info in edge_mask_info.items()
            if info['num_masked'] > 0
        }

        # Predict edges using edge decoder
        edge_scores = self.edge_decoder(enc_rep, masked_edges_dict)

        # Compute edge reconstruction loss
        total_edge_loss = 0.0
        total_edges = 0
        edge_loss_components = {}

        for edge_type, info in edge_mask_info.items():
            masked_edges = info['masked_edges']  # [2, num_masked]
            num_masked = info['num_masked']

            if num_masked == 0:
                continue

            # Get scores for masked edges (positive samples)
            pos_scores = edge_scores[edge_type]

            # Create negative samples (randomly sampled non-edges)
            src_type, _, dst_type = edge_type
            num_neg = num_masked
            neg_src = torch.randint(0, data[src_type].num_nodes, (num_neg,), device=pos_scores.device)
            neg_dst = torch.randint(0, data[dst_type].num_nodes, (num_neg,), device=pos_scores.device)
            neg_edges = torch.stack([neg_src, neg_dst])

            # Score negative samples
            neg_scores = self.edge_decoder(enc_rep, {edge_type: neg_edges})[edge_type]

            # Binary cross-entropy loss
            pos_loss = -torch.log(torch.sigmoid(pos_scores) + 1e-10).mean()
            neg_loss = -torch.log(1 - torch.sigmoid(neg_scores) + 1e-10).mean()
            edge_loss = (pos_loss + neg_loss) / 2

            total_edge_loss += edge_loss * num_masked
            total_edges += num_masked
            edge_loss_components[str(edge_type)] = edge_loss.item()

        if total_edges > 0:
            total_edge_loss = total_edge_loss / total_edges
        else:
            total_edge_loss = torch.tensor(0.0, device=x_dict['genome'].device)

        loss_dict = {
            'loss': total_edge_loss.item(),
            'edge_loss': total_edge_loss.item(),
            **edge_loss_components,
        }

        return total_edge_loss, loss_dict

    def _forward_joint(self, data, x_dict, masked_genome_idx, masked_gene_idx, edge_mask_info):
        """Forward pass for joint feature + structure reconstruction."""
        # Feature reconstruction
        feat_loss, feat_loss_dict = self._forward_feature(data, x_dict, masked_genome_idx, masked_gene_idx)

        # Structure reconstruction
        if edge_mask_info is not None and self.edge_decoder is not None:
            edge_loss, edge_loss_dict = self._forward_structure(data, x_dict, edge_mask_info)

            # Combine losses
            total_loss = feat_loss + self.edge_loss_weight * edge_loss

            loss_dict = {
                'loss': total_loss.item(),
                'feature_loss': feat_loss.item(),
                'edge_loss': edge_loss.item(),
                'num_masked_genome': feat_loss_dict.get('num_masked_genome', 0),
                'num_masked_gene': feat_loss_dict.get('num_masked_gene', 0),
            }
        else:
            # Fall back to feature-only if no edge masking
            total_loss = feat_loss
            loss_dict = feat_loss_dict

        return total_loss, loss_dict

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
        params = [self.encoder_to_decoder.parameters(), self.decoder.parameters()]
        if self.edge_decoder is not None:
            params.append(self.edge_decoder.parameters())
        return chain(*params)
