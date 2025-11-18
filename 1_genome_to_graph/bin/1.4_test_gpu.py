#!/usr/bin/env python
"""Test GPU availability and try importing CUDA libraries."""

import sys

print("="*70)
print("GPU Test Script")
print("="*70)
print()

# Test 1: Check numpy/basic imports
print("Test 1: Basic imports")
try:
    import numpy as np
    print("  ✓ numpy:", np.__version__)
except Exception as e:
    print(f"  ✗ numpy failed: {e}")

# Test 2: Check if CUDA toolkit is accessible
print("\nTest 2: CUDA toolkit")
import subprocess
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    if result.returncode == 0:
        print("  ✓ nvidia-smi works")
        # Print GPU info
        lines = result.stdout.split('\n')
        for line in lines:
            if 'NVIDIA' in line or 'GPU' in line:
                print(f"    {line.strip()}")
    else:
        print(f"  ✗ nvidia-smi failed: {result.stderr}")
except Exception as e:
    print(f"  ✗ nvidia-smi not found: {e}")

# Test 3: Try PyTorch CUDA
print("\nTest 3: PyTorch CUDA")
try:
    import torch
    print(f"  ✓ PyTorch version: {torch.__version__}")
    print(f"  CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"  CUDA version: {torch.version.cuda}")
        print(f"  Device count: {torch.cuda.device_count()}")
        print(f"  Device name: {torch.cuda.get_device_name(0)}")
except ImportError:
    print("  ⚠ PyTorch not installed")
except Exception as e:
    print(f"  ✗ PyTorch error: {e}")

# Test 4: Try cuML (RAPIDS)
print("\nTest 4: RAPIDS cuML")
try:
    import cuml
    print(f"  ✓ cuML version: {cuml.__version__}")
    print("  ✓ RAPIDS cuML is available!")
except ImportError:
    print("  ⚠ cuML not installed (need to install RAPIDS)")
    print("  Installation: conda install -c rapidsai -c conda-forge cuml")
except Exception as e:
    print(f"  ✗ cuML error: {e}")

# Test 5: Try CuPy (alternative CUDA library)
print("\nTest 5: CuPy")
try:
    import cupy as cp
    print(f"  ✓ CuPy version: {cp.__version__}")
    print(f"  CUDA version: {cp.cuda.runtime.runtimeGetVersion()}")

    # Simple test
    x = cp.array([1, 2, 3])
    print(f"  ✓ Simple GPU array works: {x}")
except ImportError:
    print("  ⚠ CuPy not installed")
    print("  Installation: pip install cupy-cuda11x")
except Exception as e:
    print(f"  ✗ CuPy error: {e}")

# Test 6: Check available CUDA version
print("\nTest 6: CUDA version from nvcc")
try:
    result = subprocess.run(['nvcc', '--version'], capture_output=True, text=True)
    if result.returncode == 0:
        print("  ✓ nvcc found:")
        for line in result.stdout.split('\n'):
            if 'release' in line.lower():
                print(f"    {line.strip()}")
    else:
        print("  ⚠ nvcc not in PATH")
except Exception as e:
    print(f"  ⚠ nvcc not found: {e}")

print("\n" + "="*70)
print("GPU Test Complete")
print("="*70)
