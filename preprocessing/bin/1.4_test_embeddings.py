#!/usr/bin/env python
"""
Quick test script to verify ESM-C embedding generation works.

Tests on a small number of genomes before running the full pipeline.
"""

import sys
from pathlib import Path
import torch
from Bio import SeqIO
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein


def test_model_loading():
    """Test that the ESM-C model loads correctly."""
    print("Testing ESM-C model loading...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    try:
        model = ESMC.from_pretrained("esmc_600m").to(device)
        model.eval()
        print("✓ Model loaded successfully")
        return model, device
    except Exception as e:
        print(f"✗ Error loading model: {e}")
        sys.exit(1)


def test_single_sequence(model, device):
    """Test embedding generation on a single sequence."""
    print("\nTesting single sequence embedding...")

    # Test sequence
    test_seq = "MKTFIFLALLGAAVAFPVDDDDKIVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSG"

    try:
        protein = ESMProtein(sequence=test_seq)

        with torch.no_grad():
            embedding = model.encode(protein)

        # Mean pool to get sequence-level embedding
        pooled = embedding.mean(axis=0).cpu().numpy()

        print(f"✓ Generated embedding of shape: {pooled.shape}")
        print(f"  Embedding dimension: {pooled.shape[0]}")
        return True
    except Exception as e:
        print(f"✗ Error generating embedding: {e}")
        return False


def test_protein_file(model, device):
    """Test loading and processing a real protein file."""
    print("\nTesting real protein file processing...")

    # Find first protein file
    gene_dir = Path("data/refseq_gene_annotations")
    protein_files = sorted(gene_dir.glob("*_prodigal_proteins.faa"))

    if not protein_files:
        print("✗ No protein files found")
        return False

    test_file = protein_files[0]
    print(f"Testing with: {test_file.name}")

    try:
        # Load first 5 sequences
        sequences = []
        for i, record in enumerate(SeqIO.parse(test_file, "fasta")):
            if i >= 5:
                break
            seq = str(record.seq).replace("*", "")
            if len(seq) > 0 and len(seq) <= 1024:
                sequences.append({
                    "id": record.id,
                    "seq": seq,
                    "len": len(seq)
                })

        print(f"  Loaded {len(sequences)} sequences")

        # Process each sequence
        embeddings = []
        for seq_info in sequences:
            protein = ESMProtein(sequence=seq_info["seq"])
            with torch.no_grad():
                embedding = model.encode(protein)
            pooled = embedding.mean(axis=0).cpu().numpy()
            embeddings.append(pooled)
            print(f"    {seq_info['id']}: length={seq_info['len']}, "
                  f"embedding shape={pooled.shape}")

        print(f"✓ Successfully processed {len(embeddings)} sequences")
        return True

    except Exception as e:
        print(f"✗ Error processing file: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("=" * 60)
    print("ESM-C Embedding Generation Test Suite")
    print("=" * 60)

    # Test 1: Model loading
    model, device = test_model_loading()

    # Test 2: Single sequence
    if not test_single_sequence(model, device):
        sys.exit(1)

    # Test 3: Real protein file
    if not test_protein_file(model, device):
        sys.exit(1)

    print("\n" + "=" * 60)
    print("✓ All tests passed! Ready to run full embedding generation.")
    print("=" * 60)
    print("\nTo process a few genomes for testing:")
    print("  python get_esm_embeddings.py --num-genomes 10")
    print("\nTo process all genomes:")
    print("  python get_esm_embeddings.py")


if __name__ == "__main__":
    main()
