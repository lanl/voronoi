cmake_minimum_required(VERSION 3.1)
project(VORONOI LANGUAGES Fortran)
set(VERSION 1.0.0)

# Set environment variables
set(CMAKE_PREFIX_PATH ${PETSC_DIR}/${PETSC_ARCH})
set(CMAKE_Fortran_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90)
set(CMAKE_Fortran_FLAGS "-fcray-pointer")

# Set source and target paths
set(EXE voronoi)
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(BIN ${CMAKE_SOURCE_DIR}/bin)
set(LIB ${CMAKE_SOURCE_DIR}/lib)

# Find required 3rd party modules
find_package(PkgConfig)
pkg_search_module(PETSC PETSc>3.1.0 IMPORTED_TARGET REQUIRED)
message(STATUS "Found PETSc ${PETSC_VERSION}")

find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")

add_definitions(${PETSC_DEFINITIONS})

set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

add_executable(
    ${EXE}
    ${SRC}/main.f90
    ${SRC}/core/constants.f90
    ${SRC}/core/Grid.f90
    ${SRC}/core/convert_accuracy.f90
    ${SRC}/compute/grid.F90
)

set_target_properties(
    ${EXE}
    PROPERTIES
    LINKER_LANGUAGE Fortran
)

set_target_properties(${EXE} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
set_target_properties(${EXE} PROPERTIES Fortran_FORMAT "FREE")

target_link_libraries(
    ${EXE}
    ${PETSC_LIBRARIES}
    ${HDF5_Fortran_LIBRARIES}
)

target_include_directories(
    ${EXE}
    PUBLIC
    ${PETSC_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

target_compile_options(
    ${EXE}
    PUBLIC
    ${PETSC_CFLAGS_OTHER}
)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

