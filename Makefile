# SUEWS Makefile - read the README file before editing

.PHONY: main clean test pip supy docs dev livehtml schema proc-csv config-ui check-dev-install mamba-dev help deactivate claude-dev claude-clean claude-start claude-stop

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON_exe = python.exe
	# F2PY_PY= /c/Users/sunt05/Anaconda2/Scripts/f2py.py
	# F2PY_EXE = $(PYTHON) $(F2PY_PY)
	TARGET=$(MODULE).pyd
else
	UNAME_S := $(shell uname -s)
	TARGET=$(MODULE).so

	ifeq ($(UNAME_S),Linux) # Linux
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

	ifeq ($(UNAME_S),Darwin) # macOS
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

endif

MODULE=SUEWS_driver

suews_dir = src/suews

docs_dir = docs

test_dir= test

release_dir = Release

makefile = Makefile.gfortran

supy_dir = src/supy

PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)

all: help

# Legacy targets for backwards compatibility
install-and-test: install test

# Display help information
help:
	@echo "SUEWS Makefile - Available targets:"
	@echo ""
	@echo "Quick Start:"
	@echo "  1. mamba/conda:  mamba activate suews-dev && make dev"
	@echo "  2. Deactivate:   make deactivate (shows command to run)"
	@echo ""
	@echo "Build and Development:"
	@echo "  dev             - Build and install SUEWS in editable mode"
	@echo "  install         - Install SUEWS to current Python environment (not editable)"
	@echo "  wheel           - Build distribution wheels"
	@echo "  claude-dev      - Set up Claude Code environment (use LOCATION=/path to specify workspace)"
	@echo "  claude-start    - Start the Claude Code sandbox (forwards arguments)"
	@echo "  claude-stop     - Stop the Claude Code sandbox"
	@echo "  claude-clean    - Remove Claude Code workspace directory"
	@echo "  clean           - Clean all build artifacts"
	@echo ""
	@echo "Legacy/Manual Commands:"
	@echo "  mamba-dev       - Build SUEWS with mamba environment check (legacy - use 'make dev')"
	@echo "  install-and-test - Install SUEWS to current Python environment + run tests (old 'make all')"
	@echo ""
	@echo "Environment Management:"
	@echo "  deactivate      - Show command to deactivate current environment"
	@echo ""
	@echo "Testing and Quality:"
	@echo "  test            - Run test suite"
	@echo ""
	@echo "Documentation:"
	@echo "  docs            - Build HTML documentation"
	@echo "  livehtml        - Start live documentation server with auto-rebuild"
	@echo "  livehtml-fast   - Start live docs server (skip supy build check)"
	@echo "  schema          - Generate JSON schema from Pydantic models"
	@echo "  proc-csv        - Process CSV files for documentation"
	@echo "  config-ui       - Start SUEWS Configuration UI server (http://localhost:8080)"
	@echo ""
	@echo "Notes:"
	@echo "  * Designed for conda/mamba environments"
	@echo "  * Automatically uses Homebrew gfortran on macOS when available"
	@echo "  * Activation commands:"
	@echo "    - mamba/conda: mamba activate suews-dev"
	@echo "  * To deactivate: use 'make deactivate' for environment-specific commands"
	@echo "  * Use 'make help' to see this help again"

# make suews driver library
suews:
	$(MAKE) -C $(suews_dir) libdriver; # make SUEWS library
	# -rm -rf *.o *.mod *.f95 *.a *.dSYM

# make supy and install locally
# NOTE: `--no-build-isolation` is used to avoid dependency issues with the editable install.
# ref: https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#editable-installs
dev:
	@echo "Building supy with development install..."
	@if [ -x "/opt/homebrew/bin/gfortran" ]; then \
		echo "Using Homebrew gfortran for better macOS compatibility"; \
		FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
	else \
		$(PYTHON) -m pip install --no-build-isolation --editable .; \
	fi

# install supy locally
install:
	$(PYTHON) -m pip install .

# make supy dist and test
test:
	$(PYTHON) -m pytest test -v --tb=short --cov=supy --cov-report=term-missing

# make supy wheels using cibuild
wheel:
	$(PYTHON) -m pip wheel --no-deps . -w wheelhouse

# Helper target to check for local supy installation
check-dev-install:
	@if ! $(PYTHON) -c "import supy" >/dev/null 2>&1; then \
		echo "ERROR: 'supy' not found. Please install it in editable mode first by running:"; \
		echo "  make dev"; \
		exit 1; \
	fi

# documentation (requires built package)
docs: check-dev-install
	@echo "Building documentation..."
	$(MAKE) -B -C $(docs_dir) html

# live html documentation (requires built package)
livehtml: check-dev-install
	@echo "Starting live documentation server..."
	$(MAKE) -B -C $(docs_dir) livehtml

# live html documentation without build check (faster for doc-only changes)
livehtml-fast:
	@echo "Starting live documentation server (skipping supy build check)..."
	$(MAKE) -B -C $(docs_dir) livehtml

# Generate JSON schema from SUEWSConfig Pydantic model
schema: check-dev-install
	@echo "Generating JSON schema from SUEWSConfig Pydantic model..."
	cd $(docs_dir) && $(PYTHON) gen_schema.py

# Process CSV files for documentation
proc-csv: check-dev-install
	@echo "Processing CSV files for documentation..."
	cd $(docs_dir) && $(PYTHON) source/related-softwares/supy/proc_var_info/gen_rst.py

# Start SUEWS Configuration UI server
config-ui:
	@echo "Starting SUEWS Configuration UI server..."
	@echo "This will serve the config UI at http://localhost:8080"
	@echo "Press Ctrl+C to stop the server"
	@echo ""
	cd $(docs_dir)/source/_static && $(PYTHON) run_server.py

# Clean all build artifacts
clean:
	@echo "Cleaning SUEWS build artifacts..."
	$(MAKE) -C $(suews_dir) clean || true
	@echo "Cleaning documentation build..."
	-cd $(docs_dir) && rm -rf build/ || true
	@echo "Cleaning Python build artifacts..."
	rm -rf build/ dist/ *.egg-info/ wheelhouse/
	@echo "Cleaning temporary files..."
	find . -name "*.pyc" -delete || true
	find . -name "__pycache__" -type d -exec rm -rf {} + || true
	find . -name "*.log" -delete || true
	@echo "Clean complete."

# this is to test cibuildwheel locally
cibw:
	CIBW_BUILD=cp312-macosx* \
	CIBW_ARCH=arm64 \
	CIBW_TEST_REQUIRES=pytest \
	CIBW_TEST_COMMAND="python -m pytest '{project}/test'" \
	pipx run cibuildwheel==2.16.5 --platform macos


# Manual mamba development (legacy target - use 'make dev' instead)
mamba-dev:
	@echo "Building supy with mamba environment..."
	@if [ -z "$$CONDA_DEFAULT_ENV" ] && [ -z "$$MAMBA_DEFAULT_ENV" ]; then \
		echo "ERROR: No mamba/conda environment detected."; \
		echo "Please activate mamba environment first:"; \
		echo "  mamba activate suews-dev"; \
		echo "Then run 'make mamba-dev' or just use 'make dev'"; \
		exit 1; \
	fi
	@if [ -x "/opt/homebrew/bin/gfortran" ]; then \
		echo "Using Homebrew gfortran for better macOS compatibility"; \
		FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
	else \
		$(PYTHON) -m pip install --no-build-isolation --editable .; \
	fi

# Deactivate current environment
deactivate:
	@echo "Attempting to deactivate current environment..."
	@if [ -n "$$CONDA_DEFAULT_ENV" ] || [ -n "$$MAMBA_DEFAULT_ENV" ]; then \
		echo "Detected conda/mamba environment: $$CONDA_DEFAULT_ENV$$MAMBA_DEFAULT_ENV"; \
		echo "Run: conda deactivate"; \
		echo "Note: You need to run 'conda deactivate' manually in your shell"; \
	elif [ -n "$$VIRTUAL_ENV" ]; then \
		echo "Detected virtual environment: $$VIRTUAL_ENV"; \
		echo "Run: deactivate"; \
		echo "Note: You need to run 'deactivate' manually in your shell"; \
	else \
		echo "No conda/mamba/virtual environment detected"; \
	fi; \
	echo ""; \
	echo "Environment deactivation commands cannot be executed from Makefiles."; \
	echo "Please run the suggested command in your shell."

# Claude Code development environment setup
# Usage: make claude-dev [LOCATION=/path/to/workspace]
# Default location: ~/claude-suews-workspace
LOCATION ?= $(HOME)/claude-suews-workspace

claude-dev:
	@echo "🚀 Setting up SUEWS Claude Code development environment..."
	@echo "📍 Target location: $(LOCATION)"
	@echo ""

	# Create workspace directory if it doesn't exist
	@mkdir -p "$(LOCATION)"

	# Clone or update SUEWS repository
	@if [ -d "$(LOCATION)/SUEWS/.git" ]; then \
		echo "📦 Updating existing SUEWS repository..."; \
		cd "$(LOCATION)/SUEWS" && git fetch --all && git pull; \
	else \
		echo "📦 Cloning SUEWS repository..."; \
		cd "$(LOCATION)" && git clone https://github.com/UMEP-dev/SUEWS.git; \
		echo "📦 Initialising submodules..."; \
		cd "$(LOCATION)/SUEWS" && git submodule init && git submodule update; \
	fi

	# The workspace will be on the default branch after cloning.
	# You can specify a different branch when running run-claude-dev.sh, e.g.:
	# ./run-claude-dev.sh --branch main
	@echo ""
	@echo "ℹ️  SUEWS workspace is on branch: $$(cd $(LOCATION)/SUEWS && git rev-parse --abbrev-ref HEAD)"

	# Copy claude-dev folder from current location if not present in target
	@if [ ! -d "$(LOCATION)/SUEWS/claude-dev" ]; then \
		if [ -d "claude-dev" ]; then \
			echo "📋 Copying claude-dev files from current location..."; \
			cp -r claude-dev "$(LOCATION)/SUEWS/"; \
		else \
			echo "❌ Error: claude-dev folder not found in current location or target"; \
			exit 1; \
		fi \
	fi

	# Run setup from the cloned location
	@echo ""
	@echo "🔧 Running Claude Code setup..."
	@cd "$(LOCATION)/SUEWS" && chmod +x claude-dev/setup-claude-dev.sh
	@cd "$(LOCATION)/SUEWS" && ./claude-dev/setup-claude-dev.sh

	@echo ""
	@echo "✅ Setup complete!"
	@echo "📂 SUEWS workspace: $(LOCATION)/SUEWS"
	@echo "🚀 To start development:"
	@echo "   cd $(LOCATION)/SUEWS && ./start-claude-dev.sh"
	@echo ""
	@echo "💡 Tip: Add to your shell profile for quick access:"
	@echo "   alias suews-claude='cd $(LOCATION)/SUEWS && ./start-claude-dev.sh'"

# Start Claude Code sandbox, passing any additional arguments
# Usage: make claude-start ARGS="--rebuild --branch main"
claude-start:
	@echo "▶️  Starting Claude Code sandbox..."
	@cd "$(LOCATION)/SUEWS" && ./start-claude-dev.sh $(ARGS)

# Stop Claude Code sandbox
claude-stop:
	@echo "⏹️  Stopping Claude Code sandbox..."
	@cd "$(LOCATION)/SUEWS" && ./stop-claude-dev.sh

# Clean up Claude Code workspace
claude-clean:
	@echo "🧹 Cleaning Claude Code workspace..."
	@if [ -d "$(LOCATION)" ]; then \
		echo "📂 Removing workspace at: $(LOCATION)"; \
		rm -rf "$(LOCATION)"; \
		echo "✅ Workspace removed"; \
	else \
		echo "ℹ️  No workspace found at: $(LOCATION)"; \
	fi
	@echo ""
	@echo "💡 Run 'make claude-dev' to create a fresh workspace"

