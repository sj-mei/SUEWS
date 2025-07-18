# Refined UV Worktree Workflow

This guide provides an optimised workflow for using UV in git worktrees for SUEWS development, addressing specific caveats discovered during testing.

## Key Findings

1. **UV works perfectly** for creating and managing virtual environments in worktrees
2. **Python 3.13 compatibility issue**: Some packages (scipy 1.13.1) don't have wheels for Python 3.13 yet
3. **`uv run` caveat**: Attempts to sync dependencies from pyproject.toml, which can fail with newer Python versions
4. **Solution**: Use activated environments rather than `uv run` for now

## Recommended Workflow

### Step 1: Create Worktree and Environment

```bash
# Set your feature name
FEATURE="my-feature"

# Create the worktree
git worktree add worktrees/$FEATURE feature/$FEATURE || exit 1

# Navigate to the worktree
cd worktrees/$FEATURE

# Create virtual environment with UV (ultra-fast!)
uv venv

# IMPORTANT: Activate the environment (required due to Python 3.13 compatibility)
source .venv/bin/activate

# Create a marker file
echo "$FEATURE" > .worktree-name
```

### Step 2: Install Dependencies

**Critical Package Name Corrections for pip/UV:**
- Use `matplotlib` instead of `matplotlib-base` (conda name)
- Use `tables` instead of `pytables` (conda name)

```bash
# Install all requirements with correct pip names
uv pip install pandas scipy matplotlib matplotlib-inline scikit-learn scikit-image \
    geopandas rtree openpyxl tables psutil salem==0.3.8 floweaver==2.0.0 \
    f90nml click pydantic ipykernel jupyter_client jupyter_core \
    pytest pytest-cov ruff f90wrap==0.2.16 atmosp "meson-python>=0.17.0"
```

### Step 3: Build SuPy

```bash
# Build SuPy in development mode
make dev  # Auto-detects UV and uses it for installation

# Verify installation
python -c "import supy; print(f'✓ SuPy {supy.__version__} ready')"
```

### Step 4: Run Tests

```bash
# Run full test suite
make test  # All tests should pass
```

## Important Caveats

### 1. Python Version Compatibility

**Issue**: Python 3.13 is very new, and some packages don't have prebuilt wheels yet.

**Impact**: `scipy 1.13.1` attempts to build from source, requiring OpenBLAS, which fails.

**Solution**: UV automatically installs the latest compatible versions (e.g., scipy 1.16.0) which do have wheels.

### 2. UV Run Limitations

**Issue**: `uv run` attempts to sync dependencies from pyproject.toml, which can fail with version conflicts.

**Current Status**: Don't use `uv run` commands without activation. Instead:

```bash
# Always activate first
source .venv/bin/activate

# Then run commands normally
python script.py
pytest
make test
```

### 3. Package Name Differences

When converting from conda/mamba environments, note these package name changes:

| Conda/Mamba Name | pip/UV Name |
|------------------|-------------|
| matplotlib-base  | matplotlib  |
| pytables        | tables      |

## Performance Benefits

Despite the caveats, UV provides significant performance improvements:

1. **Environment creation**: ~2 seconds (vs minutes with mamba)
2. **Package installation**: 10-100x faster than pip
3. **Dependency resolution**: Modern, fast resolver
4. **Caching**: Aggressive caching reduces repeated downloads

## Testing Results

With the refined workflow:
- ✅ All 148 tests pass
- ✅ SuPy builds successfully
- ✅ Development environment fully functional
- ✅ No dependency conflicts
- ✅ Fast setup and teardown

## Quick Reference Commands

```bash
# Setup
FEATURE="my-feature"
git worktree add worktrees/$FEATURE feature/$FEATURE
cd worktrees/$FEATURE
uv venv && source .venv/bin/activate

# Install (with correct package names)
uv pip install pandas scipy matplotlib matplotlib-inline scikit-learn scikit-image \
    geopandas rtree openpyxl tables psutil salem==0.3.8 floweaver==2.0.0 \
    f90nml click pydantic ipykernel jupyter_client jupyter_core \
    pytest pytest-cov ruff f90wrap==0.2.16 atmosp "meson-python>=0.17.0"

# Build & Test
make dev
make test

# Cleanup
cd ../..
git worktree remove worktrees/$FEATURE --force
```

## Future Improvements

1. **Update pyproject.toml** to pin scipy>=1.14 when available with Python 3.13 wheels
2. **Create UV lockfile** for reproducible environments
3. **Enable `uv run`** once package compatibility improves
4. **Consider Python 3.12** for maximum compatibility if issues persist

## Conclusion

UV provides excellent performance for worktree development. The current limitation with `uv run` is temporary and related to Python 3.13 being cutting-edge. The activated environment workflow is stable and significantly faster than mamba alternatives.