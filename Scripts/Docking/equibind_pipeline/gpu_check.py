"""GPU support validation."""
from __future__ import annotations

import subprocess

import torch

from .config import CFG


def validate_gpu() -> None:
    """Print a summary of GPU availability and the recommended device."""
    print("=" * 80)
    print("GPU SUPPORT VALIDATION")
    print("=" * 80)
    print("\nPyTorch CUDA Status:")
    print(f"  PyTorch version: {torch.__version__}")
    print(f"  CUDA available: {torch.cuda.is_available()}")

    if torch.cuda.is_available():
        print(f"  CUDA version: {torch.version.cuda}")
        print(f"  Number of GPUs: {torch.cuda.device_count()}")
        print(f"  GPU Name: {torch.cuda.get_device_name(0)}")
        print(f"  GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
        try:
            x = torch.randn(100, 100).cuda()
            y = torch.randn(100, 100).cuda()
            torch.matmul(x, y)
            print("  \u2713 GPU tensor operations working")
        except Exception as e:
            print(f"  \u2717 GPU tensor test failed: {e}")
    else:
        print("  \u26a0 No CUDA-capable GPU detected; EquiBind will use CPU (slower)")

    print("\nNVIDIA GPU Information (nvidia-smi):")
    try:
        result = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            print(result.stdout)
        else:
            print("  \u26a0 nvidia-smi not available or failed")
    except FileNotFoundError:
        print("  \u26a0 nvidia-smi command not found")
    except Exception as e:
        print(f"  \u26a0 Error running nvidia-smi: {e}")

    print("\n" + "=" * 80)
    if torch.cuda.is_available():
        print("\u2713 GPU ACCELERATION AVAILABLE")
        print(f"  Recommended setting: equibind_device = 'cuda'")
    else:
        print("\u26a0 GPU NOT AVAILABLE - Using CPU")
        print(f"  Recommended setting: equibind_device = 'cpu'")
    print(f"  Current setting: equibind_device = '{CFG.equibind_device}'")
    print("=" * 80)
