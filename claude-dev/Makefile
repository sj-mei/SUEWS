# SUEWS Claude Code Development Environment Testing
# Quick testing and verification of Docker environment

.PHONY: help build test-env test-dev test-suews test-all clean interactive benchmark

# Docker image settings
IMAGE_NAME := suews-claude-dev
IMAGE_TAG := latest
FULL_IMAGE := $(IMAGE_NAME):$(IMAGE_TAG)

# Paths
SUEWS_ROOT := $(shell dirname "$(PWD)")
DOCKERFILE := Dockerfile.claude-dev

help:
	@echo "SUEWS Claude Code Docker Testing - Available targets:"
	@echo ""
	@echo "Build and Test:"
	@echo "  build           - Build Docker image"
	@echo "  test-env        - Test basic environment (gfortran, Python, packages)"
	@echo "  test-dev        - Test development tools (meson-python, f90wrap)"
	@echo "  test-suews      - Test SUEWS build process"
	@echo "  test-all        - Run comprehensive test suite"
	@echo ""
	@echo "Development:"
	@echo "  interactive     - Start interactive container session"
	@echo "  benchmark       - Run build performance benchmarks"
	@echo ""
	@echo "Maintenance:"
	@echo "  clean           - Remove Docker image and containers"
	@echo "  clean-all       - Clean all Docker artifacts"
	@echo ""
	@echo "Quick Start:"
	@echo "  make build && make test-all"

# Build Docker image
build:
	@echo "🏗️  Building SUEWS Claude Code Docker image..."
	@echo "📁 Build context: $(SUEWS_ROOT)"
	@echo "📄 Dockerfile: $(DOCKERFILE)"
	@echo "🏷️  Image tag: $(FULL_IMAGE)"
	docker build -f $(DOCKERFILE) -t $(FULL_IMAGE) "$(SUEWS_ROOT)" --progress=plain
	@echo "✅ Docker image built: $(FULL_IMAGE)"
	@echo ""
	@echo "To test the image interactively, run:"
	@echo "  make interactive"

# Build with no cache (clean build)
build-clean:
	@echo "🏗️  Building SUEWS Claude Code Docker image (clean build)..."
	docker build -f $(DOCKERFILE) -t $(FULL_IMAGE) "$(SUEWS_ROOT)" --no-cache --progress=plain
	@echo "✅ Docker image built: $(FULL_IMAGE)"

# Test basic environment
test-env:
	@echo "🧪 Testing basic environment..."
	docker run --rm $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		echo '=== Testing gfortran ===' && \
		gfortran --version && \
		echo '=== Testing Python ===' && \
		python --version && \
		echo '=== Testing Python packages ===' && \
		python -c 'import numpy, scipy, pandas, matplotlib; print(\"✅ Core packages OK\")' && \
		echo '✅ Basic environment test PASSED'"

# Test development tools
test-dev:
	@echo "🧪 Testing development tools..."
	docker run --rm $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		echo '=== Testing meson-python ===' && \
		python -c 'import mesonpy; print(\"✅ meson-python OK\")' && \
		echo '=== Testing f90wrap ===' && \
		python -c 'import f90wrap; print(\"✅ f90wrap OK\")' && \
		echo '=== Testing Fortran compilation ===' && \
		echo 'program test; print *, \"Hello Fortran!\"; end program' > test.f90 && \
		gfortran test.f90 -o test_fortran && \
		./test_fortran && \
		echo '✅ Development tools test PASSED'"

# Test SUEWS build process
test-suews:
	@echo "🧪 Testing SUEWS build process..."
	@echo "🐳 Using image: $(FULL_IMAGE)"
	@echo "📁 SUEWS source: $(SUEWS_ROOT)"
	@echo ""
	docker run --rm -v "$(SUEWS_ROOT):/workspace/SUEWS" -w /workspace/SUEWS $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		echo '=== Verifying environment ===' && \
		which gfortran && \
		which meson && \
		which ninja && \
		echo '' && \
		echo '=== Updating submodules ===' && \
		git submodule update --init --recursive && \
		echo '=== Cleaning previous build artifacts ===' && \
		make clean && \
		echo '=== Testing SUEWS development build ===' && \
		make dev && \
		echo '✅ SUEWS build test PASSED'"

# Test SUEWS with full build and tests
test-suews-full:
	@echo "🧪 Testing SUEWS full build and test suite..."
	docker run --rm -v "$(SUEWS_ROOT):/workspace/SUEWS" $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		cd /workspace/SUEWS && \
		git submodule update --init --recursive && \
		make clean && \
		make && \
		make test && \
		echo '✅ SUEWS full test PASSED'"

# Run all tests
test-all: test-env test-dev test-suews
	@echo ""
	@echo "🎉 All tests PASSED!"
	@echo "✅ Docker environment is ready for SUEWS development"

# Run comprehensive test including full SUEWS test suite
test-comprehensive: test-env test-dev test-suews-full
	@echo ""
	@echo "🎉 Comprehensive tests PASSED!"
	@echo "✅ Docker environment fully validated"


# Interactive container session
interactive:
	@echo "🐳 Starting interactive SUEWS development session..."
	@echo "📁 SUEWS source mounted at: /workspace/SUEWS"
	@echo "🔧 Environment auto-activated"
	@echo ""
	docker run --rm -it \
		-v "$(SUEWS_ROOT):/workspace/SUEWS" \
		-p 8080:8080 \
		$(FULL_IMAGE)

# Performance benchmark
benchmark:
	@echo "📊 Running SUEWS build performance benchmark..."
	@echo "⏱️  Timing development build process..."
	time docker run --rm -v "$(SUEWS_ROOT):/workspace/SUEWS" $(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		cd /workspace/SUEWS && \
		git submodule update --init --recursive && \
		make clean && \
		make dev"

# Monitor resource usage during build
monitor:
	@echo "📊 Monitoring resource usage during SUEWS build..."
	docker run --name suews-build-monitor --rm -d \
		-v "$(SUEWS_ROOT):/workspace/SUEWS" \
		$(FULL_IMAGE) bash -c "\
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		cd /workspace/SUEWS && \
		git submodule update --init --recursive && \
		make clean && make dev && \
		sleep 30" && \
	docker stats --no-stream suews-build-monitor && \
	echo "✅ Resource monitoring complete"

# Clean up Docker image
clean:
	@echo "🧹 Cleaning SUEWS Claude Code Docker artifacts..."
	-docker rmi $(FULL_IMAGE) 2>/dev/null || echo "Image $(FULL_IMAGE) not found"
	-docker system prune -f
	@echo "✅ Cleanup complete"

# Clean all Docker artifacts
clean-all:
	@echo "🧹 Performing comprehensive Docker cleanup..."
	-docker rmi $(FULL_IMAGE) 2>/dev/null || echo "Image $(FULL_IMAGE) not found"
	-docker image prune -a -f
	-docker container prune -f
	-docker system prune -a -f
	@echo "✅ Comprehensive cleanup complete"

# Quick verification (fast test)
verify:
	@echo "⚡ Quick verification of Docker environment..."
	@echo "🐳 Testing image: $(FULL_IMAGE)"
	@echo ""
	@docker run --rm $(FULL_IMAGE) bash -c "\
		echo '=== System Information ===' && \
		cat /etc/os-release | grep PRETTY_NAME && \
		echo '' && \
		echo '=== Conda Environment ===' && \
		source /opt/conda/etc/profile.d/conda.sh && \
		conda activate suews-dev && \
		conda info --envs | grep suews-dev && \
		echo '' && \
		echo '=== Critical Tools ===' && \
		echo -n 'gfortran: ' && gfortran --version | head -1 && \
		echo -n 'gcc:      ' && gcc --version | head -1 && \
		echo -n 'python:   ' && python --version && \
		echo -n 'meson:    ' && meson --version && \
		echo -n 'ninja:    ' && ninja --version && \
		echo '' && \
		echo '=== Python Packages ===' && \
		python -c 'import numpy; print(f\"numpy:    {numpy.__version__}\")' && \
		python -c 'import pandas; print(f\"pandas:   {pandas.__version__}\")' && \
		python -c 'import matplotlib; print(f\"matplotlib: {matplotlib.__version__}\")' && \
		python -c 'import mesonpy; print(\"meson-python: ✅\")' && \
		python -c 'import f90wrap; print(\"f90wrap:      ✅\")' && \
		echo '' && \
		echo '=== Additional Tools ===' && \
		which gh > /dev/null && echo 'GitHub CLI: ✅' || echo 'GitHub CLI: ❌' && \
		which jq > /dev/null && echo 'jq:         ✅' || echo 'jq:         ❌' && \
		which doxygen > /dev/null && echo 'doxygen:    ✅' || echo 'doxygen:    ❌' && \
		echo '' && \
		echo '✅ Quick verification PASSED'"

# Show Docker image information
info:
	@echo "🐳 Docker image information:"
	@echo "Image: $(FULL_IMAGE)"
	@echo "Dockerfile: $(DOCKERFILE)"
	@echo "SUEWS root: $(SUEWS_ROOT)"
	@echo ""
	@if docker image inspect $(FULL_IMAGE) >/dev/null 2>&1; then \
		echo "📦 Image exists:"; \
		docker image inspect $(FULL_IMAGE) --format "  Size: {{.Size}} bytes"; \
		docker image inspect $(FULL_IMAGE) --format "  Created: {{.Created}}"; \
		docker image inspect $(FULL_IMAGE) --format "  Architecture: {{.Architecture}}"; \
	else \
		echo "❌ Image not found. Run 'make build' first."; \
	fi

# Check if image exists
check:
	@if docker image inspect $(FULL_IMAGE) >/dev/null 2>&1; then \
		echo "✅ Docker image $(FULL_IMAGE) exists"; \
	else \
		echo "❌ Docker image $(FULL_IMAGE) not found"; \
		echo "Run 'make build' to create it"; \
		exit 1; \
	fi