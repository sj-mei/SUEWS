# Environment Setup Guide

This guide covers Python environment setup for SUEWS development using uv, the ultra-fast Python package manager.

## Quick Start with uv (Recommended)

### Installation

```bash
# macOS/Linux
brew install uv
# or
curl -LsSf https://astral.sh/uv/install.sh | sh

# Verify installation
uv --version
```

### Setup and Usage

```bash
# Create virtual environment
uv venv

# Install packages (10-100x faster than pip!)
uv pip install pandas scipy matplotlib

# Run without activation - the modern way!
uv run python script.py
uv run pytest
uv run make test

# Or activate traditionally if preferred
source .venv/bin/activate
python script.py  # Works as normal
```

### Installing SUEWS Dependencies

```bash
# Core packages for SUEWS development
uv pip install pandas scipy matplotlib-base matplotlib-inline scikit-learn scikit-image \
    geopandas rtree openpyxl pytables psutil salem==0.3.8 floweaver==2.0.0 \
    f90nml click pydantic ipykernel jupyter_client jupyter_core \
    pytest pytest-cov ruff f90wrap==0.2.16 atmosp meson-python>=0.17.0

# Then build SUEWS
make dev
```

## Alternative: Standard Python venv

If uv is not available in your environment:

```bash
# Create and activate
python -m venv .venv
source .venv/bin/activate  # Linux/macOS

# Install packages (slower than uv)
pip install --upgrade pip
pip install pandas scipy matplotlib  # etc.

# Deactivate when done
deactivate
```

## Compiler Requirements

SUEWS requires a Fortran compiler:

```bash
# macOS
brew install gcc  # Provides gfortran

# Linux
sudo apt-get install gfortran  # Ubuntu/Debian
sudo dnf install gcc-gfortran  # Fedora/RHEL

# Verify
which gfortran
```

## Performance Comparison

| Tool | Create Env | Install 50 packages |
|------|------------|-------------------|
| uv | 1s | 5-10s |
| pip | 3s | 60-120s |

## Best Practices

1. **One environment per worktree** - Avoid conflicts
2. **Always use uv when available** - It's dramatically faster
3. **Run tests after setup** - `make test` to verify

## Quick Reference

```bash
# uv commands
uv venv                    # Create environment
uv pip install package     # Install package
uv run python script.py    # Run without activation

# venv commands (fallback)
python -m venv .venv       # Create environment
source .venv/bin/activate  # Activate
pip install package        # Install package
deactivate                 # Deactivate
```

For mamba/conda setup (complex scientific environments only), see the legacy documentation or SUEWS development guide.