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

# ── 2. Create 'pockets' conda environment ────────────────────────────────
ENV_NAME="pockets"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Conda environment '$ENV_NAME' already exists"
else
    echo "Creating conda environment '$ENV_NAME' (Python 3.11)..."
    conda create -n "$ENV_NAME" --override-channels -c conda-forge python=3.11 -y -q
    echo "✓ Environment created"
fi

conda activate "$ENV_NAME"
echo "Python: $(python --version)"

# ── 3. Install Python dependencies ───────────────────────────────────────
echo "Installing Python packages..."
conda install -n "$ENV_NAME" --override-channels -c conda-forge \
    numpy scipy pandas matplotlib ipykernel -y -q

echo "✓ Python packages installed"

# ── Base directory for tools ──────────────────────────────────────────────
TOOLS_DIR="$HOME/tools"
mkdir -p "$TOOLS_DIR"

# ── 4. Install fpocket ───────────────────────────────────────────────────
FPOCKET_DIR="$TOOLS_DIR/fpocket"
FPOCKET_BIN="$FPOCKET_DIR/bin/fpocket"

if [ -f "$FPOCKET_BIN" ]; then
    echo "✓ fpocket already installed at $FPOCKET_BIN"
else
    echo "Installing fpocket..."

    # Install build dependencies (requires sudo)
    sudo apt-get update -qq && sudo apt-get install -y -qq \
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
            sudo apt-get install -y -qq fpocket 2>/dev/null || true
        fi
    fi
    cd "$HOME/master_dev"
fi

# ── 5. Install Java (required by p2rank) ─────────────────────────────────
if command -v java &>/dev/null; then
    echo "✓ Java already installed: $(java -version 2>&1 | head -1)"
else
    echo "Installing Java (OpenJDK 17)..."
    sudo apt-get update -qq && sudo apt-get install -y -qq openjdk-17-jre-headless > /dev/null 2>&1
    if command -v java &>/dev/null; then
        echo "✓ Java installed: $(java -version 2>&1 | head -1)"
    else
        echo "✗ Java installation failed — p2rank will not work"
    fi
fi

# ── 6. Install p2rank ────────────────────────────────────────────────────
P2RANK_DIR="$TOOLS_DIR/p2rank_2.5"
P2RANK_BIN="$P2RANK_DIR/prank"

if [ -f "$P2RANK_BIN" ]; then
    echo "✓ p2rank already installed at $P2RANK_DIR"
else
    echo "Installing p2rank 2.5..."
    P2RANK_URL="https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz"
    wget -q "$P2RANK_URL" -O /tmp/p2rank.tar.gz
    tar -xzf /tmp/p2rank.tar.gz -C "$TOOLS_DIR"
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
mkdir -p "$HOME/master_dev/storage/pocket_results/fpocket_results"
mkdir -p "$HOME/master_dev/storage/pocket_results/p2rank_results"
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
echo "  fpocket results: $HOME/master_dev/storage/pocket_results/fpocket_results/"
echo "  p2rank results:  $HOME/master_dev/storage/pocket_results/p2rank_results/"
echo ""
echo "Protein PDBs: $(ls $HOME/master_dev/Orai/*.pdb 2>/dev/null | wc -l) files in Orai/"
echo ""
