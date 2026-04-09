#!/bin/bash
# =============================================================================
# EquiBind Environment Setup for RunPod (A100 GPU)
# =============================================================================
# This script:
#   1. Installs Miniconda (if not present)
#   2. Creates the 'equibind' conda environment with all dependencies
#   3. Clones and sets up EquiBind
#   4. Registers the conda env as a Jupyter kernel
#
# Usage:
#   chmod +x setup_equibind_env.sh
#   ./setup_equibind_env.sh
#
# After running, restart the Jupyter kernel and select the 'equibind' kernel.
# =============================================================================

set -e  # Exit on any error

echo "================================================================"
echo " EquiBind Environment Setup"
echo "================================================================"

# ── 1. Ensure conda is available ─────────────────────────────────────────
# Use existing Anaconda/Miniconda installation
if command -v conda &>/dev/null; then
    CONDA_DIR="$(conda info --base)"
elif [ -d "$HOME/anaconda3" ]; then
    CONDA_DIR="$HOME/anaconda3"
elif [ -d "$HOME/miniconda3" ]; then
    CONDA_DIR="$HOME/miniconda3"
else
    echo "✗ No conda installation found. Please install Anaconda or Miniconda first."
    exit 1
fi

eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"
echo "Using existing conda at: $CONDA_DIR"
echo "Conda version: $(conda --version)"

# ── 2. Create equibind conda environment ──────────────────────────────────
ENV_NAME="equibind"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists"
else
    echo "Creating conda environment '$ENV_NAME' (Python 3.9)..."
    conda create -n "$ENV_NAME" --override-channels -c conda-forge python=3.9 -y -q
    echo "✓ Environment created"
fi

# Activate the environment
conda activate "$ENV_NAME"

echo "Python: $(python --version)"
echo "Installing dependencies..."

# ── 3. Install PyTorch with CUDA 12.4 support ────────────────────────────
echo "Installing PyTorch (CUDA 12.4)..."
pip install torch==2.4.1 --index-url https://download.pytorch.org/whl/cu124 -q

# ── 4. Install DGL with CUDA support ─────────────────────────────────────
echo "Installing DGL (CUDA 12.4)..."
pip install dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html -q

# ── 5. Install RDKit and scientific dependencies ─────────────────────────
echo "Installing RDKit and scientific packages..."
conda install -n "$ENV_NAME" --override-channels -c conda-forge rdkit numpy scipy pandas -y -q

# ── 6. Install remaining Python dependencies ─────────────────────────────
echo "Installing remaining pip dependencies..."
pip install -q \
    PyYAML \
    biopython \
    biopandas \
    POT \
    tqdm \
    e3nn \
    dgllife \
    joblib \
    tensorboard \
    jinja2

# ── 7. Clone and setup EquiBind ──────────────────────────────────────────
EQUIBIND_DIR="$HOME/tools/EquiBind"

if [ -d "$EQUIBIND_DIR" ]; then
    echo "✓ EquiBind already cloned at $EQUIBIND_DIR"
else
    echo "Cloning EquiBind..."
    git clone https://github.com/HannesStark/EquiBind.git "$EQUIBIND_DIR"
    echo "✓ EquiBind cloned"
fi

# Download pre-trained checkpoint if not present
CKPT_DIR="$EQUIBIND_DIR/runs/flexible_self_docking"
CKPT_FILE="$CKPT_DIR/best_checkpoint.pt"

if [ -f "$CKPT_FILE" ]; then
    echo "✓ EquiBind checkpoint found"
else
    echo "⚠  EquiBind checkpoint not found at $CKPT_FILE"
    echo "   You may need to download it manually or train the model."
    echo "   Check: https://github.com/HannesStark/EquiBind#pretrained-model"
    # Try the zenodo download if available
    mkdir -p "$CKPT_DIR"
    echo "   Attempting to download from GitHub releases..."
    # EquiBind ships the checkpoint in the repo itself under runs/
    cd "$EQUIBIND_DIR"
    git lfs pull 2>/dev/null || echo "   git-lfs not available, checkpoint may need manual download"
    cd -
fi

# ── 8. Install Jupyter kernel for the equibind env ───────────────────────
echo "Installing Jupyter kernel..."
pip install ipykernel -q
python -m ipykernel install --user --name equibind --display-name "Python (equibind)"
echo "✓ Jupyter kernel 'equibind' registered"

# ── 9. Verify installation ───────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Verification"
echo "================================================================"

python -c "
import torch
print(f'PyTorch:      {torch.__version__}')
print(f'CUDA avail:   {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU:          {torch.cuda.get_device_name(0)}')
    print(f'GPU memory:   {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')

import dgl
print(f'DGL:          {dgl.__version__}')

from rdkit import Chem
print(f'RDKit:        OK')

import numpy as np
print(f'NumPy:        {np.__version__}')

import pandas as pd
print(f'Pandas:       {pd.__version__}')

import yaml
print(f'PyYAML:       OK')

import scipy
print(f'SciPy:        {scipy.__version__}')
"

echo ""
echo "================================================================"
echo " Setup Complete"
echo "================================================================"
echo ""
echo "Next steps:"
echo "  1. Activate env:      conda activate equibind"
echo "  2. In Jupyter, select the 'Python (equibind)' kernel"
echo "  3. Run the notebook:  00_Equibind_Para_Runpod.ipynb"
echo ""
echo "Paths:"
echo "  EQUIBIND_DIR:     $HOME/tools/EquiBind"
echo "  EQUIBIND_DEVICE:  cuda"
echo ""
