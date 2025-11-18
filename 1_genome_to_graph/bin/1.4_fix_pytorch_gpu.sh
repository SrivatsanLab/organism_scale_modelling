#!/bin/bash
#SBATCH --job-name=fix_pytorch_gpu
#SBATCH --partition=chorus
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/fix_pytorch_gpu_%j.out
#SBATCH --error=logs/fix_pytorch_gpu_%j.err

echo "========================================"
echo "Fix PyTorch GPU Support in glm2_env"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load CUDA module
module load CUDA/11.7.0

# Load environment
source ~/.bashrc
micromamba activate glm2_env

# Check current PyTorch
echo "Current PyTorch installation:"
python -c "import torch; print(f'PyTorch version: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}')"
echo ""

# Check CUDA on system
echo "System CUDA version:"
nvidia-smi --query-gpu=driver_version --format=csv,noheader
nvcc --version 2>/dev/null || echo "nvcc not in PATH"
echo ""

# Uninstall CPU PyTorch
echo "Removing CPU-only PyTorch..."
pip uninstall -y torch torchvision torchaudio
micromamba remove -y pytorch pytorch-cuda --force
echo ""

# Install PyTorch with CUDA support using pip (more reliable than conda for GPU)
echo "Installing PyTorch with CUDA 11.8 support via pip..."
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
echo ""

# Verify GPU support
echo "========================================"
echo "Verification"
echo "========================================"
echo ""

python -c "
import torch
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')

if torch.cuda.is_available():
    print(f'CUDA version: {torch.version.cuda}')
    print(f'cuDNN version: {torch.backends.cudnn.version()}')
    print(f'GPU count: {torch.cuda.device_count()}')
    print(f'GPU 0: {torch.cuda.get_device_name(0)}')
    print(f'GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')

    # Test GPU computation
    print('')
    print('Testing GPU computation...')
    x = torch.randn(1000, 1000).cuda()
    y = torch.randn(1000, 1000).cuda()
    z = torch.matmul(x, y)
    print(f'✓ GPU computation successful!')
    print(f'  Result shape: {z.shape}')
    print(f'  Result device: {z.device}')
else:
    print('')
    print('ERROR: CUDA still not available!')
    print('Please check CUDA installation and drivers.')
"
echo ""

# Install missing package for gLM2
echo "Installing einops..."
pip install einops
echo ""

# Test gLM2 on GPU
echo "========================================"
echo "Testing gLM2 on GPU"
echo "========================================"
echo ""

python -c "
from transformers import AutoModel, AutoTokenizer
import torch

print('Loading gLM2_650M model...')
model = AutoModel.from_pretrained(
    'tattabio/gLM2_650M',
    torch_dtype=torch.bfloat16,
    trust_remote_code=True,
)

if torch.cuda.is_available():
    model = model.to('cuda')
    print('✓ Model loaded on GPU')
else:
    print('⚠ Model on CPU (GPU not available)')

tokenizer = AutoTokenizer.from_pretrained(
    'tattabio/gLM2_650M',
    trust_remote_code=True,
)

# Test inference
test_seq = '<+>MALTK<+>aattggcc<->MGKL'
print(f'Testing with: {test_seq}')

inputs = tokenizer(test_seq, return_tensors='pt')
if torch.cuda.is_available():
    inputs = {k: v.to('cuda') for k, v in inputs.items()}

with torch.no_grad():
    outputs = model(**inputs)
    embedding = outputs.last_hidden_state.mean(dim=1)

print(f'✓ Inference successful!')
print(f'  Embedding shape: {embedding.shape}')
print(f'  Embedding device: {embedding.device}')
print(f'  Expected dim: {model.config.hidden_size}')
"

echo ""
echo "========================================"
echo "Fix complete!"
echo "========================================"
echo ""
echo "To use the environment:"
echo "  micromamba activate glm2_env"
echo ""
echo "The environment now has:"
echo "  - PyTorch with CUDA 11.8 support"
echo "  - gLM2 model downloaded"
echo "  - All required dependencies"
