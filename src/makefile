#================================================================#
#  VORONOI Makefile
# -------------------------------------------------------------- #
#
#  
#  Targets:
#    - voronoi : compile the VORONOI binary
#    - lglibs : compiles the LaGriT libraries used in VORONOI
#    - test : runs the test suite
#    - clean : destroys all build artifacts
#    - help : shows this help screen
#
#  Variables:
#    FC (default: gfortran) : Fortran MPI source compiler
#    FFLAGS : User-defined Fortran flags
#
#  Building:
#    To build VORONOI, you must first set four
#    environment variables:
#      - HDF5_DIR : /path/to/hdf5/
#      - LAGRIT_DIR : /path/to/LaGriT/
#      - PETSC_DIR : /path/to/petsc/
#      - PETSC_ARCH : dependant on sys. architecture
#
#    Further, `make lglibs` must be done before `make voronoi`.
#
#  Example:
#      $ export PETSC_DIR=/path/to/petsc 
#      $ export PETSC_ARCH=arch-darwin-c-opt
#      $ export HDF5_DIR=/path/to/hdf5/
#
#      $ make lglibs LAGRIT_DIR=../LaGriT/
#      $ make voronoi
#================================================================#

voronoi_src = ./
common_src   = ./

LG_BUILD_DIR = ./build/

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

MYFLAGS = -I. -fcray-pointer -L${HDF5_DIR}/lib -I${HDF5_DIR}/include

ifdef debug 
	MYFLAGS += ${FC_DEFINE_FLAG}DEBUG
endif

CFLAGS   =  
FFLAGS   = 
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = voronoi 

# Begin Source Block
grid_obj = \
	${common_src}grid_aux.o\
	${common_src}grid.o\
	${common_src}voronoi_constants.o\
	${common_src}user_sub.o

# End Source Block
voronoi_obj = $(grid_obj) ${voronoi_src}voronoi.o

define VORONOI_HELP
----------------------- VORONOI MAKEFILE -----------------------

This is an automatic build tool for the VORONOI tessellation
software.

USAGE:

	make [options] [target]

TARGETS:

	- voronoi : compile the VORONOI binary
	- lglibs : compiles the LaGriT libraries used in VORONOI
	- test : runs the test suite
	- clean : destroys all build artifacts
	- help : shows this help screen

VARIABLES:

	FC (default: gfortran) : Fortran MPI source compiler
	FFLAGS : User-defined Fortran flags
	HDF5_DIR : /path/to/hdf5/
	LAGRIT_DIR : /path/to/LaGriT/
	PETSC_DIR : /path/to/petsc/
	PETSC_ARCH : dependant on sys. architecture

BUILDING:

	To build VORONOI, you must first set four
	environment variables:
		- HDF5_DIR : /path/to/hdf5/
		- LAGRIT_DIR : /path/to/LaGriT/
		- PETSC_DIR : /path/to/petsc/
		- PETSC_ARCH : dependant on sys. architecture

	Further, `make lglibs` must be done before `make voronoi`.

EXAMPLE:

		> export PETSC_DIR=/path/to/petsc 
		> export PETSC_ARCH=arch-darwin-c-opt
		> export HDF5_DIR=/path/to/hdf5/

		> make lglibs LAGRIT_DIR=../LaGriT/
		> make voronoi

endef
export VORONOI_HELP

.PHONY: validate voronoi test help lglibs docs

voronoi :: validate $(voronoi_obj)
	${FLINKER} -o voronoi ${LG_BUILD_DIR}/lagrit_fdate.o $(voronoi_obj) \
	${LG_BUILD_DIR}/lg_main_lib.a ${LG_BUILD_DIR}/lg_util.a \
	-L${HDF5_DIR}/lib -I${HDF5_DIR}/include \
	${PETSC_LIB} -lhdf5_fortran ${LIBS} 

validate :
	@CHECK=0
	@[ -z "${HDF5_DIR}" ] && echo "HDF5_DIR not defined" && exit 1 || CHECK=1
	@[ -z "${PETSC_DIR}" ] && echo "PETSC_DIR not defined" && exit 1 || CHECK=1
	@[ -z "${PETSC_ARCH}" ] && echo "PETSC_ARCH not defined" && exit 1 || CHECK=1
	@[ -z "${LAGRIT_DIR}" ] && echo "LAGRIT_DIR not defined" && exit 1 || CHECK=1
	@echo "All paths are set."

test : 
	cd ../test/sanity_check && python run_tests.py

clean :: 
	rm -f ${LG_BUILD_DIR}/{*.a,*.o,*.mod}

%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

## Build LaGriT libraries and relevant object files for use in VORONOI
lglibs :
	mkdir -p ${LG_BUILD_DIR}

	make -C ${LAGRIT_DIR}/lg_util/src lib
	make -C ${LAGRIT_DIR}/src WITHEXODUS=0 lib

	cp ${LAGRIT_DIR}/src/objects/lagrit_fdate.o ${LG_BUILD_DIR}/lagrit_fdate.o
	cp ${LAGRIT_DIR}/src/lg_main_lib.a ${LG_BUILD_DIR}/lg_main_lib.a
	cp ${LAGRIT_DIR}/lg_util/src/lg_util.a ${LG_BUILD_DIR}/lg_util.a

help : 
	@echo "$$VORONOI_HELP"

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
grid.o : voronoi_constants.o user_sub.o
voronoi.o : grid_aux.o grid.o 
