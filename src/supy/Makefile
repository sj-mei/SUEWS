# -*- makefile -*-
.PHONY: main clean test pip supy

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


# install package in dev mode and do pytest
test:
	pip install -e .
	pytest -s ${MODULE}/test


# If wanted, clean all *.o files after build
clean:
	rm -rf ${MODULE}/*.so ${MODULE}/*.pyc ${MODULE}/__pycache__ ${MODULE}/*.dSYM
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
local:
	pipx run cibuildwheel==2.9.0 --platform macos --config-file pyproject-local-test.toml