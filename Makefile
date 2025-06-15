# SUEWS Makefile - read the README file before editing

.PHONY: main clean test pip supy docs dev livehtml schema proc-csv config-ui check-dev-install uv-dev uv-sync uv-clean mamba-dev help deactivate uv-activate uv

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

all: install test

# Display help information
help:
	@echo "SUEWS Makefile - Available targets:"
	@echo ""
	@echo "Quick Start - Choose your environment:"
	@echo "  1. mamba/conda:  mamba activate suews-dev && make dev"
	@echo "  2. uv:           make uv-dev"
	@echo "  3. deactivate:   make deactivate (shows command to run)"
	@echo ""
	@echo "Build and Development:"
	@echo "  dev             - Smart build: auto-detects environment and uses appropriate settings"
	@echo "  uv-dev          - Create uv environment + build SUEWS (complete setup)"
	@echo "  install         - Install SUEWS (not editable)"
	@echo "  wheel           - Build distribution wheels"  
	@echo "  clean           - Clean all build artifacts"
	@echo ""
	@echo "Legacy/Manual Commands:"
	@echo "  mamba-dev       - Build SUEWS with mamba environment (manual isolated build)"
	@echo ""
	@echo "Environment Management:"
	@echo "  uv-sync         - Create and sync uv development environment"
	@echo "  uv-activate     - Show command to activate uv environment"
	@echo "  uv-clean        - Clean uv environment and build directory"
	@echo "  deactivate      - Show command to deactivate current environment"
	@echo ""
	@echo "Testing and Quality:"
	@echo "  test            - Run test suite"
	@echo ""
	@echo "Documentation:"
	@echo "  docs            - Build HTML documentation"
	@echo "  livehtml        - Start live documentation server with auto-rebuild"
	@echo "  schema          - Generate JSON schema from Pydantic models"
	@echo "  proc-csv        - Process CSV files for documentation"
	@echo "  config-ui       - Start SUEWS Configuration UI server (http://localhost:8080)"
	@echo ""
	@echo "Environment Activation:"
	@echo "  mamba/conda:  mamba activate suews-dev"
	@echo "  uv:           source /tmp/suews-uv-env/bin/activate (after make uv-dev)"
	@echo ""
	@echo "Notes:"
	@echo "  * mamba: Traditional conda environment with all dependencies"
	@echo "  * uv: Lightweight environment with faster dependency resolution"  
	@echo "  * Both environments use isolated build directories to avoid conflicts"
	@echo "  * Use 'make help' to see this help again"

# make suews driver library
suews:
	$(MAKE) -C $(suews_dir) libdriver; # make SUEWS library
	# -rm -rf *.o *.mod *.f95 *.a *.dSYM

# make supy and install locally
# NOTE: `--no-build-isolation` is used to avoid dependency issues with the editable install.
# ref: https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#editable-installs
# Automatically detects environment type and uses appropriate build directory
dev:
	@echo "Building supy with development install..."
	@if [ -n "$$CONDA_DEFAULT_ENV" ] || [ -n "$$MAMBA_DEFAULT_ENV" ]; then \
		echo "Detected conda/mamba environment: $$CONDA_DEFAULT_ENV$$MAMBA_DEFAULT_ENV"; \
		echo "Using isolated build directory for mamba environment"; \
		if [ -x "/opt/homebrew/bin/gfortran" ]; then \
			echo "Using Homebrew gfortran for better macOS compatibility"; \
			mkdir -p /tmp/suews-builds/mamba-build; \
			MESON_BUILD_ROOT="/tmp/suews-builds/mamba-build" FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
		else \
			mkdir -p /tmp/suews-builds/mamba-build; \
			MESON_BUILD_ROOT="/tmp/suews-builds/mamba-build" $(PYTHON) -m pip install --no-build-isolation --editable .; \
		fi \
	elif [ -n "$$VIRTUAL_ENV" ] && [[ "$$VIRTUAL_ENV" == *"/tmp/suews-uv-env"* ]]; then \
		echo "Detected uv environment: $$VIRTUAL_ENV"; \
		echo "Using default build approach for uv environment"; \
		if [ -x "/opt/homebrew/bin/gfortran" ]; then \
			echo "Using Homebrew gfortran for better macOS compatibility"; \
			FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
		else \
			$(PYTHON) -m pip install --no-build-isolation --editable .; \
		fi \
	else \
		echo "No conda/mamba/uv environment detected, using default build"; \
		if [ -x "/opt/homebrew/bin/gfortran" ]; then \
			echo "Using Homebrew gfortran for better macOS compatibility"; \
			FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
		else \
			$(PYTHON) -m pip install --no-build-isolation --editable .; \
		fi \
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

# UV-based development environment targets
# Cross-platform uv environment setup with isolated builds

uv-clean:
	@echo "Cleaning uv environment and build directory..."
ifeq ($(OS),Windows_NT)
	rm -rf "%TEMP%\suews-uv-env"
	rm -rf "%TEMP%\suews-builds\uv-build"
else
	rm -rf /tmp/suews-uv-env
	rm -rf /tmp/suews-builds/uv-build
endif

uv-sync: uv-clean
	@echo "Creating uv development environment outside project directory..."
ifeq ($(OS),Windows_NT)
	@echo "Environment: %TEMP%\suews-uv-env"
	@echo "Build directory: %TEMP%\suews-builds\uv-build"
	@mkdir "%TEMP%\suews-builds" 2>nul || true
	uv venv --python 3.12 "%TEMP%\suews-uv-env"
else
	@echo "Environment: /tmp/suews-uv-env"
	@echo "Build directory: /tmp/suews-builds/uv-build"
	@mkdir -p /tmp/suews-builds
	uv venv --python 3.12 /tmp/suews-uv-env
endif
	@echo "Installing build requirements..."
ifeq ($(OS),Windows_NT)
	"%TEMP%\suews-uv-env\Scripts\activate" && uv pip install pip setuptools meson-python wheel pytest f90wrap==0.2.16 numpy
else
	. /tmp/suews-uv-env/bin/activate && uv pip install pip setuptools meson-python wheel pytest f90wrap==0.2.16 numpy
endif
	@echo "Installing supy in editable mode with dev dependencies..."
ifeq ($(OS),Windows_NT)
	"%TEMP%\suews-uv-env\Scripts\activate" && python -m pip install --no-build-isolation --editable ".[dev]"
else
	@if [ -x "/opt/homebrew/bin/gfortran" ]; then \
		echo "Using Homebrew gfortran for better macOS compatibility"; \
		. /tmp/suews-uv-env/bin/activate && FC=/opt/homebrew/bin/gfortran python -m pip install --no-build-isolation --editable ".[dev]"; \
	else \
		. /tmp/suews-uv-env/bin/activate && python -m pip install --no-build-isolation --editable ".[dev]"; \
	fi
endif

# Alternative development install using uv
uv-dev: uv-sync
	@echo "SUEWS development environment ready with uv!"
	@echo "Activate with:"
ifeq ($(OS),Windows_NT)
	@echo "  %TEMP%\\suews-uv-env\\Scripts\\activate"
else
	@echo "  source /tmp/suews-uv-env/bin/activate"
endif

# Enhanced mamba development with isolated builds
mamba-dev:
	@echo "Building supy with mamba environment (isolated build)..."
ifeq ($(OS),Windows_NT)
	@echo "Build directory: %TEMP%\suews-builds\mamba-build"
	@mkdir "%TEMP%\suews-builds" 2>nul || true
	set "MESON_BUILD_ROOT=%TEMP%\suews-builds\mamba-build" && $(PYTHON) -m pip install --no-build-isolation --editable .
else
	@echo "Build directory: /tmp/suews-builds/mamba-build"
	@mkdir -p /tmp/suews-builds
	@if [ -x "/opt/homebrew/bin/gfortran" ]; then \
		echo "Using Homebrew gfortran for better macOS compatibility"; \
		MESON_BUILD_ROOT="/tmp/suews-builds/mamba-build" FC=/opt/homebrew/bin/gfortran $(PYTHON) -m pip install --no-build-isolation --editable .; \
	else \
		MESON_BUILD_ROOT="/tmp/suews-builds/mamba-build" $(PYTHON) -m pip install --no-build-isolation --editable .; \
	fi
endif

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

# Show command to activate uv environment
uv-activate:
	@if [ -d "/tmp/suews-uv-env" ]; then \
		echo "uv environment found. Activate with either:"; \
		echo "  source /tmp/suews-uv-env/bin/activate  (traditional)"; \
		echo "  cd /tmp && uv run python [script]      (uv direct)"; \
		echo "Note: uv doesn't have a direct 'activate' command, use source"; \
	else \
		echo "uv environment not found at /tmp/suews-uv-env"; \
		echo "Run 'make uv-dev' to create it first"; \
	fi

# Spawn new shell with uv environment activated
uv:
	@if [ -d "/tmp/suews-uv-env" ]; then \
		echo "Activating uv environment in new zsh shell..."; \
		echo "Type 'exit' to return to your original environment"; \
		echo 'source ~/.zshrc 2>/dev/null || true; source /tmp/suews-uv-env/bin/activate; echo "uv environment activated!"; exec zsh' > /tmp/uv-activate-tmp.zsh; \
		zsh /tmp/uv-activate-tmp.zsh; \
		rm -f /tmp/uv-activate-tmp.zsh; \
	else \
		echo "uv environment not found. Run 'make uv-dev' first"; \
	fi