function build_petsc()
{
    cd ./imported/petsc/
    python ./configure           \
        --download-mpich         \
        --download-hdf5=1
    #   --with-cxx=0 <- cxx compiler likely not needed
    #
}

function build_voronoi()
{
    PETSC_DIR="$(pwd)/imported/petsc"
    PETSC_ARCH="arch-darwin-c-debug"
    HDF5_DIR="${PETSC_DIR}/${PETSC_ARCH}/lib"
    LAGRIT_DIR="$(pwd)/imported/lagrit"

    make -f src/Makefile lglibs voronoi
}

build_petsc
build_voronoi

