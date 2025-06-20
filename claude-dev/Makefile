# SUEWS Claude Code Development Environment Testing
# Simple testing for gfortran and mamba/python readiness

.PHONY: help build test clean

# Docker image settings
IMAGE_NAME := suews-claude-dev
IMAGE_TAG := latest
FULL_IMAGE := $(IMAGE_NAME):$(IMAGE_TAG)

# Paths
SUEWS_ROOT := $(shell dirname "$(PWD)")
DOCKERFILE := Dockerfile.claude-dev

help:
	@echo "SUEWS Docker Environment - Simple Testing"
	@echo ""
	@echo "Available targets:"
	@echo "  build     - Build Docker image"
	@echo "  test      - Test gfortran and mamba/python for SUEWS development"
	@echo "  clean     - Remove Docker image"
	@echo ""
	@echo "Quick Start:"
	@echo "  make build && make test"

# Build Docker image
build:
	@echo "🏗️  Building SUEWS Docker image..."
	docker build -f $(DOCKERFILE) -t $(FULL_IMAGE) "$(SUEWS_ROOT)" --progress=plain
	@echo "✅ Image built: $(FULL_IMAGE)"

# Test SUEWS development readiness
test:
	@echo "🧪 Testing SUEWS development environment..."
	@docker run --rm $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		echo '🔧 gfortran:' && gfortran --version | head -1 && \
		echo '🐍 python:  ' && python --version && \
		echo '📦 mamba:   ' && mamba --version | head -1 && \
		echo '⚙️  meson:   ' && meson --version && \
		echo '' && \
		echo '📚 Testing key Python packages:' && \
		python -c 'import numpy, pandas, matplotlib, xarray; print(\"✅ Core packages ready\")' && \
		python -c 'import mesonpy, f90wrap; print(\"✅ Build tools ready\")' && \
		echo '' && \
		echo '🔨 Testing Fortran compilation:' && \
		echo 'program test; print *, \"Hello SUEWS!\"; end program' > test.f90 && \
		gfortran test.f90 -o test_fortran && \
		./test_fortran && \
		echo '' && \
		echo '✅ Environment ready for SUEWS development!'"

# Clean up Docker image
clean:
	@echo "🧹 Cleaning Docker image..."
	-docker rmi $(FULL_IMAGE) 2>/dev/null || echo "Image not found"
	@echo "✅ Cleanup complete"