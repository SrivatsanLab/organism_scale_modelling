#!/usr/bin/env python
"""
Test gene prediction on a single genome to verify installation.

Usage:
    python test_gene_prediction.py
"""

import sys
from pathlib import Path
from Bio import SeqIO


def test_imports():
    """Test that required packages are available."""
    print("Testing imports...")

    packages = {
        "biopython": lambda: __import__("Bio"),
        "pandas": lambda: __import__("pandas"),
        "numpy": lambda: __import__("numpy"),
    }

    for name, import_func in packages.items():
        try:
            import_func()
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  ✗ {name} - NOT INSTALLED")
            return False

    return True


def test_pyrodigal():
    """Test if pyrodigal is available."""
    print("\nTesting pyrodigal...")
    try:
        import pyrodigal
        print(f"  ✓ pyrodigal {pyrodigal.__version__}")
        return True
    except ImportError:
        print("  ✗ pyrodigal not installed")
        print("    Install with: pip install pyrodigal")
        return False


def test_prodigal_cli():
    """Test if Prodigal command-line tool is available."""
    import subprocess

    print("\nTesting Prodigal CLI...")
    try:
        result = subprocess.run(
            ["prodigal", "-v"],
            capture_output=True,
            text=True
        )
        version = result.stderr.split()[1] if result.stderr else "unknown"
        print(f"  ✓ Prodigal {version}")
        return True
    except FileNotFoundError:
        print("  ✗ Prodigal not found in PATH")
        print("    Install with: conda install -c bioconda prodigal")
        return False


def test_gene_prediction():
    """Test gene prediction on a real genome."""
    print("\nTesting gene prediction on real data...")

    # Find a test genome
    genome_dir = Path("data/refseq_genomes")
    if not genome_dir.exists():
        print("  ✗ Genome directory not found: data/refseq_genomes")
        return False

    test_genomes = list(genome_dir.glob("*.fasta"))
    if not test_genomes:
        print("  ✗ No genome files found")
        return False

    test_genome = test_genomes[0]
    print(f"  Using test genome: {test_genome.name}")

    # Load genome
    records = list(SeqIO.parse(test_genome, "fasta"))
    print(f"  Genome has {len(records)} contig(s)")
    total_length = sum(len(r.seq) for r in records)
    print(f"  Total length: {total_length:,} bp")

    # Try gene prediction with pyrodigal
    try:
        import pyrodigal

        orf_finder = pyrodigal.GeneFinder()
        gene_count = 0

        for record in records:
            genes = orf_finder.find_genes(str(record.seq))
            gene_count += len(genes)

        print(f"  ✓ Predicted {gene_count} genes")
        print(f"  Average gene density: {gene_count / (total_length / 1000):.1f} genes/kb")

        # Show first few genes
        if gene_count > 0:
            first_contig_genes = orf_finder.find_genes(str(records[0].seq))
            print(f"\n  First 3 genes:")
            for i, gene in enumerate(list(first_contig_genes)[:3], start=1):
                protein = gene.translate()
                print(f"    {i}. Position {gene.begin}-{gene.end}, "
                      f"Strand {'+' if gene.strand == 1 else '-'}, "
                      f"{len(protein)} aa")

        return True

    except ImportError:
        print("  ⚠ Pyrodigal not available, skipping gene prediction test")
        return True
    except Exception as e:
        print(f"  ✗ Error during gene prediction: {e}")
        return False


def test_existing_annotations():
    """Test that existing gene annotations are accessible."""
    print("\nTesting existing gene annotations...")

    annot_dir = Path("data/refseq_gene_annotations")
    if not annot_dir.exists():
        print("  ⚠ Annotation directory not found (symlink may not be set up)")
        return True

    protein_files = list(annot_dir.glob("*_prodigal_proteins.faa"))
    print(f"  Found {len(protein_files)} protein annotation files")

    if protein_files:
        # Check first file
        test_file = protein_files[0]
        sequences = list(SeqIO.parse(test_file, "fasta"))
        print(f"  ✓ Example: {test_file.name} has {len(sequences)} proteins")
        if sequences:
            avg_len = sum(len(s.seq) for s in sequences) / len(sequences)
            print(f"    Average protein length: {avg_len:.1f} aa")
        return True
    else:
        print("  ⚠ No protein files found")
        return True


def main():
    """Run all tests."""
    print("=" * 70)
    print("Gene Prediction Test Suite")
    print("=" * 70)

    all_passed = True

    # Test 1: Imports
    if not test_imports():
        all_passed = False

    # Test 2: Pyrodigal
    has_pyrodigal = test_pyrodigal()

    # Test 3: Prodigal CLI
    has_prodigal = test_prodigal_cli()

    if not has_pyrodigal and not has_prodigal:
        print("\n" + "!" * 70)
        print("WARNING: Neither pyrodigal nor Prodigal CLI found!")
        print("Install at least one:")
        print("  pip install pyrodigal")
        print("  conda install -c bioconda prodigal")
        print("!" * 70)
        all_passed = False

    # Test 4: Gene prediction
    if has_pyrodigal or has_prodigal:
        if not test_gene_prediction():
            all_passed = False

    # Test 5: Existing annotations
    test_existing_annotations()

    # Summary
    print("\n" + "=" * 70)
    if all_passed:
        print("✓ All tests passed!")
        print("\nYou can now:")
        print("  1. Use existing annotations:")
        print("     python get_esm_embeddings.py --gene-dir data/refseq_gene_annotations/")
        print("\n  2. Or predict genes from new genomes:")
        print("     python predict_genes.py --genome genome.fasta --output-dir genes/")
    else:
        print("✗ Some tests failed. Check the output above.")
        sys.exit(1)
    print("=" * 70)


if __name__ == "__main__":
    main()
