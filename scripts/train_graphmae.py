#!/usr/bin/env python3
"""Training script for HeteroGraphMAE on genome-gene graphs."""

import sys
import os
import argparse
from pathlib import Path
import time

import torch
import torch.optim as optim
from tqdm import tqdm

# Add parent directory to path to import from src
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.GNN.construct_het_graph import construct_het_graph
from src.GNN.model import HeteroGraphMAE


def parse_args():
    parser = argparse.ArgumentParser(description='Train HeteroGraphMAE')

    # Data paths
    parser.add_argument('--data_dir', type=str, default='data',
                        help='Directory containing graph data files')

    # Model hyperparameters
    parser.add_argument('--num_hidden', type=int, default=256,
                        help='Hidden dimension')
    parser.add_argument('--num_layers', type=int, default=2,
                        help='Number of encoder layers')
    parser.add_argument('--dropout', type=float, default=0.1,
                        help='Dropout rate')
    parser.add_argument('--activation', type=str, default='prelu',
                        choices=['prelu', 'relu', 'gelu'],
                        help='Activation function')
    parser.add_argument('--norm', type=str, default='batchnorm',
                        choices=['batchnorm', 'layernorm', 'none'],
                        help='Normalization type')
    parser.add_argument('--encoder_type', type=str, default='sage',
                        choices=['sage', 'gat'],
                        help='Encoder architecture')
    parser.add_argument('--decoder_type', type=str, default='sage',
                        choices=['sage', 'gat', 'mlp', 'linear'],
                        help='Decoder architecture')

    # Loss function
    parser.add_argument('--loss_fn', type=str, default='sce',
                        choices=['sce', 'mse'],
                        help='Loss function')
    parser.add_argument('--alpha_l', type=float, default=3.0,
                        help='Alpha parameter for SCE loss')

    # Masking
    parser.add_argument('--mask_gene_ratio', type=float, default=0.3,
                        help='Ratio of gene nodes to mask')
    parser.add_argument('--mask_genome_ratio', type=float, default=0.0,
                        help='Ratio of genome nodes to mask')
    parser.add_argument('--remask_rate', type=float, default=0.5,
                        help='Re-masking rate before decoder')

    # Training
    parser.add_argument('--lr', type=float, default=0.001,
                        help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-5,
                        help='Weight decay')
    parser.add_argument('--max_epochs', type=int, default=100,
                        help='Maximum number of epochs')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--device', type=str, default='cpu',
                        help='Device (cpu or cuda:0)')

    # Checkpointing
    parser.add_argument('--save_every', type=int, default=10,
                        help='Save checkpoint every N epochs')
    parser.add_argument('--output_dir', type=str, default='checkpoints',
                        help='Directory to save checkpoints')

    return parser.parse_args()


def set_seed(seed):
    """Set random seeds for reproducibility."""
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def train_epoch(model, data, optimizer, args, device):
    """Train for one epoch."""
    model.train()

    # Move data to device
    data = data.to(device)

    # Sample random nodes to mask (GraphMAE-style: mask features, not nodes)
    num_genome = data['genome'].num_nodes
    num_gene = data['gene'].num_nodes

    # Sample masked node indices
    num_masked_genome = int(args.mask_genome_ratio * num_genome)
    num_masked_gene = int(args.mask_gene_ratio * num_gene)

    perm_genome = torch.randperm(num_genome, device=device)
    perm_gene = torch.randperm(num_gene, device=device)

    masked_genome_idx = perm_genome[:num_masked_genome]
    masked_gene_idx = perm_gene[:num_masked_gene]

    # Create masked features (mask features, keep graph structure)
    masked_x_dict = {
        'genome': data['genome'].x.clone(),
        'gene': data['gene'].x.clone(),
    }

    # Apply mask tokens using the model's mask_features method
    masked_x_dict = model.mask_features(masked_x_dict, masked_genome_idx, masked_gene_idx)

    # Forward pass (pass same graph, but with masked features)
    loss, loss_dict = model(data, masked_x_dict, masked_genome_idx, masked_gene_idx)

    # Backward pass
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    return loss.item(), loss_dict


def main():
    args = parse_args()

    # Set seed
    set_seed(args.seed)

    # Set device
    device = torch.device(args.device)
    print(f"Using device: {device}")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load graph
    print("Loading graph data...")
    data_dir = Path(args.data_dir)

    data = construct_het_graph(
        genome_features=str(data_dir / 'genome_node_attributes.npz'),
        gene_features=str(data_dir / 'gene_node_attributes.npz'),
        genome_genome_edges=str(data_dir / 'genome_genome_edges.npz'),
        genome_genome_edge_features=str(data_dir / 'genome_distance_matrix.npz'),
        genome_gene_edges=str(data_dir / 'genome_gene_edges.npz'),
        gene_gene_edges=str(data_dir / 'gene_gene_edges.npz'),
        save_path=None,
        return_graph=True,
    )

    print(f"Loaded graph:")
    print(f"  Genome nodes: {data['genome'].num_nodes}")
    print(f"  Gene nodes: {data['gene'].num_nodes}")
    print(f"  Genome features: {data['genome'].x.shape}")
    print(f"  Gene features: {data['gene'].x.shape}")

    # Get feature dimensions
    genome_in_dim = data['genome'].x.shape[1]
    gene_in_dim = data['gene'].x.shape[1]

    # Initialize model
    print("\nInitializing model...")
    model = HeteroGraphMAE(
        genome_in_dim=genome_in_dim,
        gene_in_dim=gene_in_dim,
        num_hidden=args.num_hidden,
        num_layers=args.num_layers,
        dropout=args.dropout,
        activation=args.activation,
        norm=args.norm if args.norm != 'none' else None,
        encoder_type=args.encoder_type,
        decoder_type=args.decoder_type,
        loss_fn=args.loss_fn,
        alpha_l=args.alpha_l,
        concat_hidden=False,
        remask_rate=args.remask_rate,
    )
    model = model.to(device)

    # Count parameters
    num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model has {num_params:,} trainable parameters")

    # Initialize optimizer
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    # Training loop
    print(f"\nStarting training for {args.max_epochs} epochs...")
    print(f"Masking: {args.mask_gene_ratio:.1%} gene nodes, {args.mask_genome_ratio:.1%} genome nodes")

    best_loss = float('inf')

    for epoch in range(1, args.max_epochs + 1):
        start_time = time.time()

        # Train
        loss, loss_dict = train_epoch(model, data, optimizer, args, device)

        epoch_time = time.time() - start_time

        # Log progress
        print(f"Epoch {epoch:3d} | Loss: {loss:.4f} | "
              f"Masked genes: {loss_dict['num_masked_gene']:5d} | "
              f"Masked genomes: {loss_dict['num_masked_genome']:3d} | "
              f"Time: {epoch_time:.2f}s")

        # Save checkpoint
        if epoch % args.save_every == 0:
            checkpoint_path = Path(args.output_dir) / f'checkpoint_epoch_{epoch}.pt'
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'loss': loss,
                'args': vars(args),
            }, checkpoint_path)
            print(f"  Saved checkpoint to {checkpoint_path}")

        # Save best model
        if loss < best_loss:
            best_loss = loss
            best_path = Path(args.output_dir) / 'best_model.pt'
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'loss': loss,
                'args': vars(args),
            }, best_path)

    print(f"\nTraining complete! Best loss: {best_loss:.4f}")
    print(f"Best model saved to {best_path}")


if __name__ == '__main__':
    main()
