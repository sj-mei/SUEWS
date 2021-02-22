# SUEWS Makefile - read the README file before editing
# part of this is inspired by SPARTACUS-Surface: https://github.com/ecmwf/spartacus-surface
##################################
###     -- USEFUL PATHS --     ###
##################################

# where SUEWS source code is stored
SUEWS_DIR_SRC = ./src/

# where compiled libraries are stored
SUEWS_DIR_LIB = ./lib/

# where intermediate `.mod` are stored
SUEWS_DIR_MOD = ./mod/

# where to test built version
SUEWS_DIR_BIN = ./bin/

# where to store built binaries for release
SUEWS_DIR_RELEASE_BASE = ../Release/bin/

# system specific settings
ifeq ($(OS),Windows_NT)
	SUEWS_DIR_RELEASE = $(addprefix $(SUEWS_DIR_RELEASE_BASE), win64)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux) # Linux
		SUEWS_DIR_RELEASE = $(addprefix $(SUEWS_DIR_RELEASE_BASE), Linux)
	endif
	ifeq ($(UNAME_S),Darwin) # macOS
		SUEWS_DIR_RELEASE = $(addprefix $(SUEWS_DIR_RELEASE_BASE), macOS)
	endif
endif



LIBUTILS = $(SUEWS_DIR_LIB)/libsuewsutil.a
LIBPHYS = $(SUEWS_DIR_LIB)/libsuewsphys.a
LIBDRIVER = $(SUEWS_DIR_LIB)/libsuewsdriver.a



# version specific target name:
# determined from source code
path       := $(SUEWS_DIR_SRC)suews_ctrl_const.f95 # get program name form source code
file :=  $(if $(filter $(OS),Windows_NT),$(strip $(subst /,\,${path})),${path})
# $(info file is $(file))
grep := $(if $(filter $(OS),Windows_NT),findstr,grep)
variable   := $(shell $(grep) progname ${file})
TARGET_str := $(lastword $(subst =, ,${variable})) # program name with quotes
TARGET     := $(strip $(subst $\', ,${TARGET_str})) # program name
VERSION    := $(lastword $(subst _V, ,${TARGET}))

timestamp := $(shell date +"%Y%m%dT%H%M")
EXECUTABLE := $(strip $(TARGET))_$(strip $(timestamp))

# where to test built version
SUEWS_DIR_TEST := ../Test/BaseRun/$(VERSION)


#############################
### --- CONFIGURATION --- ###
#############################

# Use the nf-config utility, if available, to set the NETCDF_INCLUDE
# and NETCDF_LIB flags
HAVE_NFCONFIG := $(shell nf-config --version)
ifdef HAVE_NFCONFIG
$(info *** Using nf-config to obtain NetCDF flags)
NETCDF_INCLUDE=-I$(shell nf-config --fflags)
NETCDF_LIB=$(shell nf-config --flibs) -lnetcdff
ifeq ($(shell nf-config --has-nc4),yes)
NETCDF4        = 1
endif
else
$(info *** nf-config not found)
endif

# make can be invoked using "make PROFILE=<prof>" in which case your
# local configuration parameters will be obtained from
# Makefile_include.<prof>
ifndef PROFILE
$(info *** No "PROFILE" variable provided, assuming "gfortran")
PROFILE = gfortran
endif

# Include a platform-specific makefile that defines FC, FCFLAGS and
# LIBS
include	Makefile.$(PROFILE)

# Check for presence of the NETCDF_INCLUDE and NETCDF_LIB flags
ifndef NETCDF_INCLUDE
$(info *** You may need to set NETCDF_INCLUDE manually)
endif
ifndef NETCDF_LIB
$(info *** You may need to set NETCDF_LIB manually)
endif

# Add single-precision flag if SINGLE_PRECISION=1 was given on the
# "make" command line
ifdef SINGLE_PRECISION
CPPFLAGS += -DSINGLE_PRECISION
endif

# Allow the capability to write NetCDF4/HDF5 files, provided the code
# is compiled against the NetCDF4 library
ifdef NETCDF4
$(info *** Building with NetCDF4/HDF5 support)
CPPFLAGS += -DNC_NETCDF4
endif

# Consolidate flags
export FC
export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) \
	$(STATIC)	\
	$(OPTFLAGS) $(DEBUGFLAGS)
	# $(NETCDF_INCLUDE) $(OMPFLAG)
export LIBS    = $(DEBUGFLAGS) $(LDFLAGS) \
	-L$(SUEWS_DIR_LIB) -lsuewsdriver -lsuewsphys -lsuewsutil \
	$(STATICLIB)
	# $(FCLIBS) $(NETCDF_LIB) $(OMPFLAG)






####################################
###  -- OTHER CONFIGURATIONS --  ###
####################################

# options for `clean`
clean: suffix :=
clean: TARGET := $(strip $(TARGET))$(strip $(suffix))

# WRF coupling (SuMin) related options
sumin: FFLAGS += -Dwrf=1


###################################
###  -- SOURCE CODE GROUPS  --  ###
###################################


# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
# NB: DO KEEP ALL FILE NAMES in lowercase for compatability across platforms
# utilities
UTILS = suews_ctrl_const.o \
	suews_util_stringmod.o \
	suews_util_qsort.o\
	suews_util_time.o \
	suews_util_meteo.o \
	suews_util_datetime.o \
	suews_util_minpack.o

# SUEWS physics modules
PHYS =	suews_phys_narp.o \
	suews_phys_atmmoiststab.o \
	suews_phys_resist.o \
	suews_ctrl_input.o  \
	suews_phys_bluews.o \
	suews_phys_waterdist.o \
	suews_phys_evap.o \
	suews_phys_snow.o \
	suews_phys_dailystate.o \
	suews_phys_estm.o \
	suews_ctrl_output.o \
	suews_phys_anohm.o \
	suews_phys_lumps.o \
	suews_phys_anemsn.o \
    suews_phys_rslprof.o \
	suews_phys_biogenco2.o \
	suews_phys_ohm.o \
	suews_phys_solweig.o \
	suews_phys_beers.o

# SUEWS driver module (wrapper of all main calculations)
DRIVER = suews_ctrl_driver.o

# Rest of the files including modules and functions which are independent
OTHERS =  suews_ctrl_translate.o \
	suews_ctrl_error.o \
	suews_ctrl_init.o \
	suews_ctrl_calculations.o

# modules under rapid development
TEST =

# modules for WRF coupling
WRF =  suews_ctrl_sumin.o

# external libraries
# SPARTACUS
EXT_LIB = libraddriver.a \
	libradsurf.a \
	libradtool.a \
	libutilities.a \


#############################
### --- BUILD TARGETS --- ###
#############################

all: suews

help:
	@echo ${HAVE_NFCONFIG} ${NETCDF_INCLUDE} ${NETCDF_LIB}
	@echo "Usage:"
	@echo "  make PROFILE=<prof>"
	@echo "where <prof> can be either 'gfortran' or 'ifort'."
	@echo $(C_MODULE)

suews: libutil libphys libdriver exe

libutil: $(addprefix $(SUEWS_DIR_SRC), $(UTILS))
	mkdir -p $(SUEWS_DIR_LIB)
	ar r $(LIBUTILS) $^

libphys:  $(addprefix $(SUEWS_DIR_SRC), $(UTILS) $(PHYS) )
	ar r $(LIBPHYS) $^

libdriver: $(addprefix $(SUEWS_DIR_SRC), $(DRIVER))
	ar r $(LIBDRIVER) $^

exe: $(SUEWS_DIR_SRC)/suews_program$(suffix_code) \
	$(addprefix $(SUEWS_DIR_SRC), $(OTHERS)) \
	$(C_MODULE)
	@echo  $(LIBUTILS) $(LIBPHYS) $(LIBDRIVER)
	@echo $(TARGET) $(EXECUTABLE) $(OTHERS)
	@echo Building SUEWS ...
	$(FC) $(FCFLAGS) $^ -o $(SUEWS_DIR_BIN)/$(EXECUTABLE) $(LIBS)
	@echo $(TARGET) has been saved in $(SUEWS_DIR_BIN)
	cd $(SUEWS_DIR_BIN) && ln -sf $(EXECUTABLE) $(TARGET)
	mkdir -p $(SUEWS_DIR_RELEASE) && cp -f $(SUEWS_DIR_BIN)/$(EXECUTABLE) $(SUEWS_DIR_RELEASE)/$(TARGET)


# rules for f95-->f90
%.o: %$(suffix_code)
	@echo $^
	$(FC) $(FCFLAGS) -c $^ -o $@

test: test_simple test_rami4pilps test_urban

test_simple:
	cd test/simple && $(MAKE) test

test_rami4pilps:
	cd test/rami4pilps && $(MAKE) test

test_urban:
	cd test/urban && $(MAKE) test

test_urban_single:
	cd test/urban && $(MAKE) test_single

test_code:
	cd driver && $(MAKE) test_code

clean: clean-obj clean-mods clean-lib clean-bin

# dist-clean: clean-tests clean-toplevel clean-utilities clean-mods dist-clean-doc

clean-obj:
	rm -rf $(SUEWS_DIR_SRC)/*.o $(SUEWS_DIR_SRC)/*.mod

clean-lib:
	rm -rf $(SUEWS_DIR_LIB)/lib*.a

clean-mods:
	rm -rf $(SUEWS_DIR_MOD)/*.mod

clean-bin:
	rm -rf $(SUEWS_DIR_BIN)/SUEWS_*



$(C_MODULE): $(subst .o,.c, $(C_MODULE))
	$(CC) -c $^ -o $@





# If OTHERS have changed, compile them again
$(OTHERS): $(UTILS) $(PHYS) $(subst .o,$(suffix_code), $(OTHERS))

# If TEST have changed, compile them again
$(TEST): $(subst .o,$(suffix_code), $(TEST))

# WRF coupling required files
$(WRF): $(subst .o,$(suffix_code), $(WRF))



.PHONY: libutil libphys libdriver exe clean-toplevel test