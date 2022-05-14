# -*- makefile -*-
.PHONY: main clean test pip supy

SUEWS_build_base = ./build/

# system specific settings
ifeq ($(OS),Windows_NT)
	SUEWS_build_dir = $(addprefix $(SUEWS_build_base), win64)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux) # Linux
		SUEWS_build_dir = $(addprefix $(SUEWS_build_base), Linux)
	endif
	ifeq ($(UNAME_S),Darwin) # macOS
		UNAME_m := $(shell uname -m)
		ifeq ($(UNAME_m),arm64) # arm64
			SUEWS_build_dir = $(addprefix $(DIR_RELEASE_BASE_SUEWS), macOS_arm64)
		else # x64
			SUEWS_build_dir = $(addprefix $(DIR_RELEASE_BASE_SUEWS), macOS_x64)
		endif
	endif
endif



SUEWS_dir = SUEWS-SourceCode

makefile = Makefile.gfortran

SuPy_dir = supy-driver

PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)



# make fortran exe
pack:
	python pack_release.py

clean:
	rm -rf *.zip
