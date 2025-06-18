# Docker Testing Workflow for SUEWS Environment

This workflow will help you test the fixed Dockerfile to ensure gfortran and Python work correctly together.

## Prerequisites

1. Docker installed and running
2. Access to this SUEWS repository

## Quick Test Commands

### 1. Build the Docker Image

```bash
# Navigate to the claude-dev directory
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/claude-dev/

# Build the Docker image
docker build -f Dockerfile.claude-dev -t suews-dev:test .
```

### 2. Basic Environment Test

```bash
# Run container and test basic functionality
docker run --rm -it suews-dev:test bash -c "
source ~/activate_suews.sh && 
echo 'Testing gfortran...' && 
gfortran --version && 
echo 'Testing Python...' && 
python --version && 
echo 'Testing Python packages...' && 
python -c 'import numpy, scipy, pandas; print(\"Core packages OK\")' &&
echo 'Environment test PASSED'
"
```

### 3. Interactive Testing Session

```bash
# Run container interactively for manual testing
docker run --rm -it -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test
```

Once inside the container, test:
```bash
# The environment should auto-activate, but if not:
source ~/activate_suews.sh

# Test gfortran compilation
echo 'program test; print *, "Hello Fortran!"; end program' > test.f90
gfortran test.f90 -o test_fortran
./test_fortran

# Test Python scientific stack
python -c "
import numpy as np
import scipy
import pandas as pd
import matplotlib
print('NumPy version:', np.__version__)
print('SciPy version:', scipy.__version__)
print('Pandas version:', pd.__version__)
print('Matplotlib version:', matplotlib.__version__)
print('All packages imported successfully!')
"

# Test meson-python (critical for SUEWS)
python -c "import mesonpy; print('meson-python OK')"

# Test f90wrap (critical for SUEWS)
python -c "import f90wrap; print('f90wrap OK')"
```

## SUEWS-Specific Testing

### 4. Test SUEWS Build Process

```bash
# Run with mounted SUEWS source
docker run --rm -it -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
echo 'Testing git submodule update...' &&
git submodule update --init --recursive &&
echo 'Testing SUEWS development build...' &&
make dev
"
```

### 5. Test SUEWS Full Build and Tests

```bash
# Full test including unit tests
docker run --rm -it -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
git submodule update --init --recursive &&
make clean &&
make &&
make test
"
```

## Comprehensive Test Script

Create a test script to run all tests automatically:

```bash
# Create test script
cat > test_docker_complete.sh << 'EOF'
#!/bin/bash
set -e

echo "=== SUEWS Docker Environment Test ==="
echo "1. Building Docker image..."
docker build -f Dockerfile.claude-dev -t suews-dev:test . --no-cache

echo "2. Testing basic environment..."
docker run --rm suews-dev:test bash -c "
source ~/activate_suews.sh && 
gfortran --version && 
python --version && 
python -c 'import numpy, scipy, pandas, matplotlib; print(\"All packages OK\")'
"

echo "3. Testing development tools..."
docker run --rm suews-dev:test bash -c "
source ~/activate_suews.sh && 
python -c 'import mesonpy, f90wrap; print(\"Development tools OK\")'
"

echo "4. Testing with SUEWS source..."
docker run --rm -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
git submodule update --init --recursive &&
make dev &&
echo 'SUEWS build test PASSED'
"

echo "=== All tests PASSED! ==="
EOF

# Make executable and run
chmod +x test_docker_complete.sh
./test_docker_complete.sh
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Build Fails During Micromamba Installation
```bash
# Check if you have internet access and try rebuilding with verbose output
docker build -f Dockerfile.claude-dev -t suews-dev:test . --no-cache --progress=plain
```

#### 2. gfortran Not Found
```bash
# Test if gfortran is properly installed
docker run --rm suews-dev:test which gfortran
docker run --rm suews-dev:test gfortran --version
```

#### 3. Python Environment Issues
```bash
# Check if micromamba environment exists
docker run --rm suews-dev:test bash -c "
eval \"\$(micromamba shell hook --shell=bash)\" &&
micromamba env list
"
```

#### 4. SUEWS Build Fails
```bash
# Debug build issues
docker run --rm -it -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
make clean &&
make dev -v  # Verbose output
"
```

## Performance Testing

### Memory and CPU Usage
```bash
# Monitor resource usage during build
docker stats --no-stream suews-dev-container &
docker run --name suews-dev-container --rm -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
make clean && make
"
```

### Build Time Benchmarking
```bash
# Time the build process
time docker run --rm -v "$(pwd)/../:/workspace/SUEWS" suews-dev:test bash -c "
source ~/activate_suews.sh &&
cd /workspace/SUEWS &&
make clean && make dev
"
```

## Success Criteria

The Docker environment is working correctly if:

1. ✅ Docker image builds without errors
2. ✅ gfortran compiles and runs Fortran code
3. ✅ Python 3.12 loads with all scientific packages
4. ✅ meson-python and f90wrap are available
5. ✅ SUEWS `make dev` completes successfully
6. ✅ SUEWS tests pass with `make test`

## Next Steps

Once all tests pass, you can:
1. Use the container for SUEWS development
2. Mount your local SUEWS directory for persistent changes
3. Set up VS Code dev containers for integrated development
4. Create production images with pre-built SUEWS

## Container Usage Examples

```bash
# Development session with source mounted
docker run --rm -it \
  -v "$(pwd)/../:/workspace/SUEWS" \
  -p 8080:8080 \
  suews-dev:test

# Run specific SUEWS commands
docker run --rm \
  -v "$(pwd)/../:/workspace/SUEWS" \
  suews-dev:test \
  bash -c "source ~/activate_suews.sh && cd /workspace/SUEWS && make docs"
```