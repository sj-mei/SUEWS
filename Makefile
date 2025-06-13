# SUEWS Makefile - read the README file before editing

.PHONY: main clean test pip supy docs dev livehtml schema proc-csv config-ui

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

# documentation (requires built package)
docs: 
	@echo "Building documentation (this will install supy if needed)..."
	$(PYTHON) -m pip install --no-build-isolation --editable . || echo "Build may have failed, trying docs anyway..."
	$(MAKE) -B -C $(docs_dir) html

# live html documentation (requires built package)
livehtml:
	@echo "Starting live documentation server (this will install supy if needed)..."
	$(PYTHON) -m pip install --no-build-isolation --editable . || echo "Build may have failed, trying docs anyway..."
	$(MAKE) -B -C $(docs_dir) livehtml

# Generate JSON schema from SUEWSConfig Pydantic model
schema: dev
	@echo "Generating JSON schema from SUEWSConfig Pydantic model..."
	cd $(docs_dir) && $(PYTHON) gen_schema.py

# Process CSV files for documentation
proc-csv: dev
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