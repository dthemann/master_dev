#!/bin/bash
# =============================================================================
# PoseBusters Environment Setup for RunPod
# =============================================================================
# This script:
#   1. Ensures conda is available
#   2. Creates the 'posebusters' conda environment with all dependencies
#   3. Installs Open Babel for PDBQT → SDF conversion
#   4. Registers the conda env as a Jupyter kernel
#
# Usage:
#   chmod +x setup_posebusters_env.sh
#   ./setup_posebusters_env.sh
#
# After running, restart the Jupyter kernel and select the 'posebusters' kernel.
# =============================================================================

set -e  # Exit on any error

echo "================================================================"
echo " PoseBusters Environment Setup"
echo "================================================================"

# ── 1. Ensure conda is available ─────────────────────────────────────────
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

# ── 2. Create posebusters conda environment ───────────────────────────────
ENV_NAME="posebusters"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists"
else
    echo "Creating conda environment '$ENV_NAME' (Python 3.10)..."
    conda create -n "$ENV_NAME" --override-channels -c conda-forge python=3.10 -y -q
    echo "✓ Environment created"
fi

# Activate the environment
conda activate "$ENV_NAME"

echo "Python: $(python --version)"
echo "Installing dependencies..."

# ── 3. Install RDKit via conda ────────────────────────────────────────────
echo "Installing RDKit..."
conda install -n "$ENV_NAME" --override-channels -c conda-forge rdkit -y -q

# ── 4. Install Open Babel (for PDBQT → SDF conversion) ───────────────────
echo "Installing Open Babel..."
conda install -n "$ENV_NAME" --override-channels -c conda-forge openbabel -y -q

# ── 5. Install Python dependencies via pip ────────────────────────────────
echo "Installing pip dependencies..."
pip install -q \
    numpy \
    pandas \
    matplotlib \
    posebusters

# ── 6. Install Jupyter kernel ────────────────────────────────────────────
echo "Installing Jupyter kernel..."
pip install ipykernel -q
python -m ipykernel install --user --name posebusters --display-name "Python (posebusters)"
echo "✓ Jupyter kernel 'posebusters' registered"

# ── 7. Verify installation ───────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Verification"
echo "================================================================"

python -c "
import numpy as np
print(f'NumPy:        {np.__version__}')

import pandas as pd
print(f'Pandas:       {pd.__version__}')

import matplotlib
print(f'Matplotlib:   {matplotlib.__version__}')

from rdkit import Chem
from rdkit.Chem import AllChem
print(f'RDKit:        OK')

from posebusters import PoseBusters
print(f'PoseBusters:  OK')
"

# Check obabel
if command -v obabel &>/dev/null; then
    echo "Open Babel:   $(obabel -V 2>&1 | head -1)"
else
    echo "⚠  obabel not found on PATH (PDBQT conversion will be unavailable)"
fi

echo ""
echo "================================================================"
echo " Setup Complete"
echo "================================================================"
echo ""
echo "Next steps:"
echo "  1. Activate env:      conda activate posebusters"
echo "  2. In Jupyter, select the 'Python (posebusters)' kernel"
echo "  3. Run the notebook:  pose_busters_para_refactored.ipynb"
echo ""
echo "Notes:"
echo "  - Open Babel (obabel) is only needed for AutoDock PDBQT conversion"
echo "  - For EquiBind/DiffDock SDF poses, only RDKit + PoseBusters are required"
echo ""
