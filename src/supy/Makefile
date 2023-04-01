# -*- makefile -*-
.PHONY: main clean test pip supy suews dev

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON_exe = python.exe

else
	UNAME_S := $(shell uname -s)



	ifeq ($(UNAME_S),Linux) # Linux
		PYTHON_exe=python
	endif

	ifeq ($(UNAME_S),Darwin) # macOS
		PYTHON_exe=python

	endif

endif


PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)
# All the files which include modules used by other modules (these therefore
# need to be compiled first)

MODULE = supy
SUEWS_dir = ../suews
supy_driver_dir = ../supy_driver


# install package in dev mode and do pytest
main: dev test

driver: suews
	rm -rf supy_driver
	$(MAKE) -C $(supy_driver_dir) clean
	$(MAKE) -C $(supy_driver_dir) test
	mv $(supy_driver_dir)/_supy_driver*.so $(MODULE)/
	mv $(supy_driver_dir)/supy_driver.py $(MODULE)/

# run pytest
test:
	pytest -s ${MODULE}/test

# install package in dev mode
dev: driver
	pip install -e .

# build wheel
wheel:
	python -m build . --wheel --outdir wheelhouse

# If wanted, clean all *.o files after build
clean:
	rm -rf ${MODULE}/*.so ${MODULE}/supy_driver.py ${MODULE}/*.pyc ${MODULE}/__pycache__ ${MODULE}/*.dSYM
	rm -rf ${MODULE}/test/*.pyc ${MODULE}/test/__pycache__
	rm -rf ${MODULE}/util/*.pyc ${MODULE}/util/__pycache__
	rm -rf ${MODULE}/cmd/*.pyc ${MODULE}/cmd/__pycache__
	rm -rf SuPy.log SuPy.log.*
	rm -rf build dist
	rm -rf .eggs .pytest_cache

# upload wheels to pypi using twine
upload:
	twine upload --skip-existing dist/*whl

# make fortran exe
suews:
	$(MAKE) -C $(SUEWS_dir) libdriver; # make SUEWS with the `suews` recipe

# only local testing on mac
local: export CIBUILDWHEEL = 1
local: export SETUPTOOLS_USE_DISTUTILS = stdlib
local: clean
	echo ${SETUPTOOLS_USE_DISTUTILS};
	echo ${CIBUILDWHEEL};
	pipx run cibuildwheel==2.9.0 --platform macos --config-file pyproject.toml

