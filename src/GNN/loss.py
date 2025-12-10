import torch
import torch.nn as nn
import torch.nn.functional as F


def sce_loss(x, y, alpha=3):
    """Scaled Cosine Error (SCE) loss from GraphMAE.

    Computes (1 - cosine_similarity(x, y))^alpha after L2 normalization.

    Parameters:
        x: Predicted features (N, D)
        y: Target features (N, D)
        alpha: Scaling exponent (default: 3)

    Returns:
        Scalar loss value
    """
    x = F.normalize(x, p=2, dim=-1)
    y = F.normalize(y, p=2, dim=-1)

    # Compute (1 - cosine_similarity)^alpha
    loss = (1 - (x * y).sum(dim=-1)).pow_(alpha)
    loss = loss.mean()

    return loss


def hetero_sce_loss(pred_dict, target_dict, mask_dict, alpha=3, node_type_weights=None):
    """Heterogeneous SCE loss for masked node reconstruction.

    Computes SCE loss only on masked nodes for each node type.

    Parameters:
        pred_dict: Dictionary of predicted features {node_type: features}
        target_dict: Dictionary of target features {node_type: features}
        mask_dict: Dictionary of boolean masks {node_type: mask} where True = masked
        alpha: Scaling exponent for SCE loss
        node_type_weights: Optional dict of weights for each node type

    Returns:
        Scalar loss value (weighted average across node types)
    """
    total_loss = 0.0
    total_count = 0

    for node_type in pred_dict.keys():
        if node_type not in target_dict or node_type not in mask_dict:
            continue

        # Get masked nodes for this type
        mask = mask_dict[node_type]
        if mask.sum() == 0:
            continue

        # Extract masked predictions and targets
        pred_masked = pred_dict[node_type][mask]
        target_masked = target_dict[node_type][mask]

        # Compute SCE loss for this node type
        loss = sce_loss(pred_masked, target_masked, alpha=alpha)

        # Apply node type weight if provided
        weight = 1.0
        if node_type_weights is not None and node_type in node_type_weights:
            weight = node_type_weights[node_type]

        total_loss += loss * weight * mask.sum().item()
        total_count += mask.sum().item()

    # Return weighted average
    if total_count > 0:
        return total_loss / total_count
    else:
        return torch.tensor(0.0, device=list(pred_dict.values())[0].device)


def mse_loss(x, y):
    """Mean squared error loss.

    Parameters:
        x: Predicted features (N, D)
        y: Target features (N, D)

    Returns:
        Scalar loss value
    """
    return F.mse_loss(x, y)


def hetero_mse_loss(pred_dict, target_dict, mask_dict, node_type_weights=None):
    """Heterogeneous MSE loss for masked node reconstruction.

    Computes MSE loss only on masked nodes for each node type.

    Parameters:
        pred_dict: Dictionary of predicted features {node_type: features}
        target_dict: Dictionary of target features {node_type: features}
        mask_dict: Dictionary of boolean masks {node_type: mask} where True = masked
        node_type_weights: Optional dict of weights for each node type

    Returns:
        Scalar loss value (weighted average across node types)
    """
    total_loss = 0.0
    total_count = 0

    for node_type in pred_dict.keys():
        if node_type not in target_dict or node_type not in mask_dict:
            continue

        # Get masked nodes for this type
        mask = mask_dict[node_type]
        if mask.sum() == 0:
            continue

        # Extract masked predictions and targets
        pred_masked = pred_dict[node_type][mask]
        target_masked = target_dict[node_type][mask]

        # Compute MSE loss for this node type
        loss = mse_loss(pred_masked, target_masked)

        # Apply node type weight if provided
        weight = 1.0
        if node_type_weights is not None and node_type in node_type_weights:
            weight = node_type_weights[node_type]

        total_loss += loss * weight * mask.sum().item()
        total_count += mask.sum().item()

    # Return weighted average
    if total_count > 0:
        return total_loss / total_count
    else:
        return torch.tensor(0.0, device=list(pred_dict.values())[0].device)
