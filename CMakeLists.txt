cmake_minimum_required(VERSION 3.1)
project(VORONOI)
enable_language(Fortran)
set(VERSION 1.0.0)

#set(CMAKE_CXX_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx)
#set(CMAKE_C_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc)
set(CMAKE_Fortran_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90)
set(CMAKE_Fortran_FLAGS "-fcray-pointer")

find_package(PkgConfig)
pkg_search_module(PETSC PETSc>3.5.0 IMPORTED_TARGET REQUIRED)

message(STATUS "Found PETSc ${PETSC_VERSION}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")

add_definitions(${PETSC_DEFINITIONS})

# -D CMAKE_C_COMPILER=${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc

set(EXE voronoi)
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(BIN ${CMAKE_SOURCE_DIR}/bin)
set(LIB ${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

add_executable(
    ${EXE}
    ${SRC}/voronoi.F90
    ${SRC}/voronoi_constants.F90
    ${SRC}/grid.F90
    ${SRC}/grid_aux.F90
    ${SRC}/convert_accuracy.f95
)

target_link_libraries(${EXE} ${PETSC_LIBRARIES})
target_include_directories(${EXE} PUBLIC ${PETSC_INCLUDE_DIRS})
target_compile_options(${EXE} PUBLIC ${PETSC_CFLAGS_OTHER})

#INSTALL(${EXE} ${FOOEXE} RUNTIME DESTINATION bin)

#include_directories(${PETSC_INCLUDES})
#set(CMAKE_C_FLAGS_INIT "")
#set(CMAKE_C_FLAGS_DEBUG_INIT "-g")
