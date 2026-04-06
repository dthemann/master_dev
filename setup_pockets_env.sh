#!/bin/bash
# =============================================================================
# Binding Site Detection Environment Setup for RunPod
# =============================================================================
# Installs:
#   1. Conda env 'pockets' with Python scientific stack
#   2. fpocket (pocket detection from protein geometry)
#   3. p2rank (machine-learning pocket prediction, requires Java)
#   4. Jupyter kernel for the 'pockets' env
#
# Usage:
#   chmod +x setup_pockets_env.sh
#   ./setup_pockets_env.sh
# =============================================================================

set -e

echo "================================================================"
echo " Binding Site Detection Environment Setup"
echo "================================================================"

# ── 1. Ensure conda is available ─────────────────────────────────────────
CONDA_DIR="$HOME/miniconda3"
CONDA_BIN="$CONDA_DIR/bin/conda"

if [ ! -f "$CONDA_BIN" ]; then
    echo "Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$CONDA_DIR"
    rm /tmp/miniconda.sh
fi

eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"
conda init bash 2>/dev/null || true
echo "Conda version: $(conda --version)"

# Accept TOS if needed (conda 26+)
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# ── 2. Create 'pockets' conda environment ────────────────────────────────
ENV_NAME="pockets"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists"
else
    echo "Creating conda environment '$ENV_NAME' (Python 3.11)..."
    conda create -n "$ENV_NAME" python=3.11 -y -q
    echo "✓ Environment created"
fi

conda activate "$ENV_NAME"
echo "Python: $(python --version)"

# ── 3. Install Python dependencies ───────────────────────────────────────
echo "Installing Python packages..."
conda install -n "$ENV_NAME" -c conda-forge \
    numpy scipy pandas matplotlib ipykernel -y -q

echo "✓ Python packages installed"

# ── 4. Install fpocket ───────────────────────────────────────────────────
FPOCKET_DIR="/workspace/fpocket"
FPOCKET_BIN="$FPOCKET_DIR/bin/fpocket"

if [ -f "$FPOCKET_BIN" ]; then
    echo "✓ fpocket already installed at $FPOCKET_BIN"
else
    echo "Installing fpocket..."

    # Install build dependencies
    apt-get update -qq && apt-get install -y -qq \
        build-essential git libnetcdf-dev > /dev/null 2>&1

    # Clone and build fpocket
    git clone https://github.com/Discngine/fpocket.git "$FPOCKET_DIR"
    cd "$FPOCKET_DIR"
    make -j$(nproc) 2>&1 | tail -5
    
    if [ -f "$FPOCKET_BIN" ]; then
        echo "✓ fpocket built successfully"
    else
        # Some fpocket builds put the binary in src/
        if [ -f "$FPOCKET_DIR/src/fpocket" ]; then
            mkdir -p "$FPOCKET_DIR/bin"
            cp "$FPOCKET_DIR/src/fpocket" "$FPOCKET_BIN"
            echo "✓ fpocket built (copied from src/)"
        else
            echo "✗ fpocket build failed, trying apt install..."
            apt-get install -y -qq fpocket 2>/dev/null || true
        fi
    fi
    cd /home/master_dev
fi

# ── 5. Install Java (required by p2rank) ─────────────────────────────────
if command -v java &>/dev/null; then
    echo "✓ Java already installed: $(java -version 2>&1 | head -1)"
else
    echo "Installing Java (OpenJDK 17)..."
    apt-get update -qq && apt-get install -y -qq openjdk-17-jre-headless > /dev/null 2>&1
    if command -v java &>/dev/null; then
        echo "✓ Java installed: $(java -version 2>&1 | head -1)"
    else
        echo "✗ Java installation failed — p2rank will not work"
    fi
fi

# ── 6. Install p2rank ────────────────────────────────────────────────────
P2RANK_DIR="/workspace/p2rank_2.5"
P2RANK_BIN="$P2RANK_DIR/prank"

if [ -f "$P2RANK_BIN" ]; then
    echo "✓ p2rank already installed at $P2RANK_DIR"
else
    echo "Installing p2rank 2.5..."
    P2RANK_URL="https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz"
    wget -q "$P2RANK_URL" -O /tmp/p2rank.tar.gz
    mkdir -p /workspace
    tar -xzf /tmp/p2rank.tar.gz -C /workspace
    rm /tmp/p2rank.tar.gz
    chmod +x "$P2RANK_BIN"

    if [ -f "$P2RANK_BIN" ]; then
        echo "✓ p2rank installed"
    else
        echo "✗ p2rank installation failed"
    fi
fi

# ── 7. Create output directories ─────────────────────────────────────────
echo "Creating output directories..."
mkdir -p /workspace/storage/pocket_results/fpocket_results
mkdir -p /workspace/storage/pocket_results/p2rank_results
echo "✓ Output directories ready"

# ── 8. Register Jupyter kernel ────────────────────────────────────────────
echo "Installing Jupyter kernel..."
python -m ipykernel install --user --name pockets --display-name "Python (pockets)"
echo "✓ Jupyter kernel 'pockets' registered"

# ── 9. Verify ────────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Verification"
echo "================================================================"

echo -n "fpocket:  "
if [ -f "$FPOCKET_BIN" ]; then
    "$FPOCKET_BIN" 2>&1 | head -1 || echo "installed at $FPOCKET_BIN"
else
    echo "NOT FOUND"
fi

echo -n "p2rank:   "
if [ -f "$P2RANK_BIN" ]; then
    echo "installed at $P2RANK_DIR"
else
    echo "NOT FOUND"
fi

echo -n "Java:     "
java -version 2>&1 | head -1 || echo "NOT FOUND"

conda activate "$ENV_NAME"
python -c "
import numpy, scipy, pandas, matplotlib
print(f'NumPy:      {numpy.__version__}')
print(f'SciPy:      {scipy.__version__}')
print(f'Pandas:     {pandas.__version__}')
print(f'Matplotlib: {matplotlib.__version__}')
"

echo ""
echo "================================================================"
echo " Setup Complete"
echo "================================================================"
echo ""
echo "Next steps:"
echo "  1. In Jupyter, select the 'Python (pockets)' kernel"
echo "  2. Run the Binding_Site_Dedection.ipynb notebook"
echo ""
echo "Tools installed:"
echo "  fpocket:  $FPOCKET_BIN"
echo "  p2rank:   $P2RANK_BIN"
echo ""
echo "Output directories:"
echo "  fpocket results: /workspace/storage/pocket_results/fpocket_results/"
echo "  p2rank results:  /workspace/storage/pocket_results/p2rank_results/"
echo ""
echo "Protein PDBs: $(ls /home/master_dev/Orai/*.pdb 2>/dev/null | wc -l) files in Orai/"
echo ""
