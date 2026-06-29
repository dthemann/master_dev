#!/bin/bash
# =============================================================================
# DiffDock (DiffDock-L) Environment Setup for RunPod (A100 / RTX 4090 GPU)
# =============================================================================
# This script:
#   1. Ensures conda is available
#   2. Creates the 'diffdock' conda environment from diffdock_runpod.yml
#      (Python 3.9 · torch 2.4.0+cu124 · PyG 2.6.1)
#   3. Clones DiffDock (gcorso/DiffDock — the DiffDock-L model)
#   4. Registers the conda env as a Jupyter kernel
#   5. Verifies the GPU + key imports
#
# Usage:
#   chmod +x setup_diffdock_env.sh
#   ./setup_diffdock_env.sh
#
# After running, restart the Jupyter kernel and select the 'diffdock' kernel.
# =============================================================================

set -e  # Exit on any error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/diffdock_runpod.yml"

echo "================================================================"
echo " DiffDock (DiffDock-L) Environment Setup"
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

# ── 2. Create the diffdock environment from the YAML ──────────────────────
ENV_NAME="diffdock"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists (skipping create)"
    echo "  To rebuild: conda env remove -n $ENV_NAME && re-run this script"
else
    if [ ! -f "$ENV_FILE" ]; then
        echo "✗ Env file not found: $ENV_FILE"
        exit 1
    fi
    echo "Creating conda environment '$ENV_NAME' from $(basename "$ENV_FILE")..."
    echo "  (downloads torch 2.4.0+cu124 + PyG wheels — a few minutes)"
    conda env create -f "$ENV_FILE"
    echo "✓ Environment created"
fi

# Activate the environment
conda activate "$ENV_NAME"
echo "Python: $(python --version)"

# ── 3. Clone DiffDock (DiffDock-L) ───────────────────────────────────────
# Mirrors the local convention: $HOME/docking_tools/DiffDock
DIFFDOCK_DIR="${DIFFDOCK_DIR:-$HOME/docking_tools/DiffDock}"

if [ -d "$DIFFDOCK_DIR/.git" ]; then
    echo "✓ DiffDock already cloned at $DIFFDOCK_DIR"
else
    echo "Cloning DiffDock into $DIFFDOCK_DIR ..."
    mkdir -p "$(dirname "$DIFFDOCK_DIR")"
    git clone https://github.com/gcorso/DiffDock.git "$DIFFDOCK_DIR"
    echo "✓ DiffDock cloned"
fi

# Model weights: DiffDock-L ships its config under workdir/ and downloads the
# .pt weights automatically on the first inference run. The SO(2)/SO(3) score
# lookup tables are also precomputed and cached on first run (a few minutes,
# one-time per pod). No manual checkpoint download is required.

# ── 4. Register Jupyter kernel ───────────────────────────────────────────
echo "Installing Jupyter kernel..."
pip install ipykernel -q
python -m ipykernel install --user --name diffdock --display-name "Python (diffdock)"
echo "✓ Jupyter kernel 'diffdock' registered"

# ── 5. Verify installation ───────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Verification"
echo "================================================================"

python -c "
import torch
print(f'PyTorch:        {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU:            {torch.cuda.get_device_name(0)}')
    print(f'Compute cap:    {torch.cuda.get_device_capability(0)}')
    print(f'GPU memory:     {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')

import torch_geometric, torch_scatter, torch_cluster
print(f'torch-geometric:{torch_geometric.__version__}')
print(f'torch-scatter:  {torch_scatter.__version__}')
print(f'torch-cluster:  {torch_cluster.__version__}')

import e3nn, esm, prody
from rdkit import Chem
import Bio
print(f'e3nn:           {e3nn.__version__}')
print('fair-esm:       OK')
print('rdkit / prody / biopython: OK')
"

echo ""
echo "================================================================"
echo " Setup Complete"
echo "================================================================"
echo ""
echo "Next steps:"
echo "  1. Activate env:   conda activate diffdock"
echo "  2. Smoke-test DiffDock on the bundled example:"
echo "       cd $DIFFDOCK_DIR"
echo "       python -m inference --config default_inference_args.yaml \\"
echo "         --protein_ligand_csv data/protein_ligand_example.csv \\"
echo "         --out_dir results/test"
echo "  3. For the batch pipeline, point run_diffdock.py at a pod config whose"
echo "     diffdock_dir / diffdock_python use the pod paths, e.g.:"
echo "       diffdock_dir:    $DIFFDOCK_DIR"
echo "       diffdock_python: $CONDA_DIR/envs/diffdock/bin/python"
echo "     then:  python run_diffdock.py -c <pod_config>.yaml --dry-run"
echo ""
echo "Notes:"
echo "  - First inference run precomputes SO(2)/SO(3) tables + downloads model"
echo "    weights (one-time, a few minutes)."
echo "  - GPU targets A100 SXM (sm_80) and RTX 4090 (sm_89) are both supported"
echo "    by torch 2.4.0+cu124."
echo ""
