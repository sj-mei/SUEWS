###################################
###  -- SOURCE CODE GROUPS  --  ###
###################################
# FC = gfortran

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
	suews_ctrl_init.o

# modules under rapid development
TEST = suews_ctrl_calculations.o

# modules for WRF coupling
WRF =  suews_ctrl_sumin.o

# external libraries
# SPARTACUS
EXT_LIB = libraddriver.a \
	libradsurf.a \
	libradtool.a \
	libutilities.a \

suffix_code = .f95
#############################
### --- BUILD TARGETS --- ###
#############################

all: build

help:
	@echo ${HAVE_NFCONFIG} ${NETCDF_INCLUDE} ${NETCDF_LIB}
	@echo "Usage:"
	@echo "  make PROFILE=<prof>"
	@echo "where <prof> can be either 'gfortran' or 'ifort'."

# build: libutil libphys libdriver

libutil: $(UTILS)
	echo $(suffix_code)
	ar r $(LIBUTIL) $(UTILS)


# rules for f95-->f90
%.o: %$(suffix_code)
	$(FC) $(FCFLAGS) -c $<

.PHONY: all help libutil %.o