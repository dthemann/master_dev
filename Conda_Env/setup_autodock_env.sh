#!/bin/bash
# =============================================================================
# AutoDock Environment Setup
# =============================================================================
# Creates a conda environment "autodock" with all dependencies needed to run
# the Autodock.ipynb notebook:
#   - RDKit  (ligand property computation)
#   - Open Babel  (PDBQT conversion via openbabel converter)
#   - Meeko  (PDBQT ligand/receptor preparation)
#   - numpy, pandas  (data handling)
#   - ipykernel  (Jupyter integration)
#
# AutoDock Vina itself is a standalone binary — the script checks for it but
# does not install it.  Set VINA_BIN in the notebook to point at your binary.
#
# Usage:
#   chmod +x setup_autodock_env.sh
#   ./setup_autodock_env.sh
#
# After running, select the "Python (autodock)" kernel in Jupyter / VS Code.
# =============================================================================

set -e  # Exit on any error

ENV_NAME="autodock"
PYTHON_VERSION="3.11"

echo "================================================================"
echo " AutoDock Environment Setup"
echo "================================================================"

# ── 1. Ensure conda is available ─────────────────────────────────────────
if command -v conda &>/dev/null; then
    CONDA_DIR="$(conda info --base)"
elif [ -d "$HOME/anaconda3" ]; then
    CONDA_DIR="$HOME/anaconda3"
elif [ -d "$HOME/miniconda3" ]; then
    CONDA_DIR="$HOME/miniconda3"
elif [ -d "/workspace/miniconda3" ]; then
    CONDA_DIR="/workspace/miniconda3"
else
    echo "✗ No conda installation found. Please install Anaconda or Miniconda first."
    exit 1
fi

eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"
echo "Using conda at: $CONDA_DIR  ($(conda --version))"

# ── 2. Create conda environment ──────────────────────────────────────────
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists — updating"
else
    echo "Creating conda environment '$ENV_NAME' (Python ${PYTHON_VERSION})..."
    conda create -n "$ENV_NAME" --override-channels -c conda-forge python="${PYTHON_VERSION}" -y -q
    echo "✓ Environment created"
fi

conda activate "$ENV_NAME"
echo "Python: $(python --version)"

# ── 3. Install conda packages ────────────────────────────────────────────
echo ""
echo "Installing conda packages (RDKit, Open Babel)..."
conda install -n "$ENV_NAME" --override-channels \
    -c conda-forge \
    rdkit \
    openbabel \
    numpy \
    -y -q
echo "✓ RDKit + Open Babel + NumPy installed"

# ── 4. Install pip packages ──────────────────────────────────────────────
echo ""
echo "Installing pip packages..."
pip install -q --upgrade \
    pandas \
    meeko \
    ipykernel

echo "✓ pandas, meeko, ipykernel installed"

# ── 5. Register Jupyter kernel ────────────────────────────────────────────
echo ""
echo "Registering Jupyter kernel..."
python -m ipykernel install --user --name "$ENV_NAME" --display-name "Python ($ENV_NAME)"
echo "✓ Jupyter kernel 'Python ($ENV_NAME)' registered"

# ── 6. Verify installation ───────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Verification"
echo "================================================================"

python -c "
import sys
print(f'Python:       {sys.version.split()[0]}')

import numpy as np
print(f'NumPy:        {np.__version__}')

import pandas as pd
print(f'Pandas:       {pd.__version__}')

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
print(f'RDKit:        OK')

import meeko
print(f'Meeko:        {meeko.__version__}')
"

# Check obabel
if command -v obabel &>/dev/null; then
    echo "Open Babel:   $(obabel -V 2>&1 | head -1)"
else
    echo "⚠  obabel not found on PATH after install"
fi

# Check mk_prepare_receptor.py (Meeko CLI)
if command -v mk_prepare_receptor.py &>/dev/null; then
    echo "Meeko CLI:    mk_prepare_receptor.py found"
else
    echo "⚠  mk_prepare_receptor.py not on PATH (may need: pip install meeko)"
fi

# Check for AutoDock Vina binary
VINA_PATHS=(
    "/home/manndo/AutoDock-Vina/build/linux/release/vina"
    "/workspace/miniconda3/envs/vina/bin/vina"
    "$(command -v vina 2>/dev/null || true)"
)
VINA_FOUND=""
for vp in "${VINA_PATHS[@]}"; do
    if [ -x "$vp" ] 2>/dev/null; then
        VINA_FOUND="$vp"
        break
    fi
done

if [ -n "$VINA_FOUND" ]; then
    echo "Vina binary:  $VINA_FOUND"
else
    echo "⚠  AutoDock Vina binary not found."
    echo "   Set VINA_BIN in the notebook to the correct path."
    echo "   Install from: https://github.com/ccsb-scripps/AutoDock-Vina"
fi

echo ""
echo "================================================================"
echo " Setup Complete"
echo "================================================================"
echo ""
echo "Next steps:"
echo "  1. conda activate $ENV_NAME"
echo "  2. In Jupyter / VS Code, select the 'Python ($ENV_NAME)' kernel"
echo "  3. Run:  Autodock.ipynb"
echo ""
