+++
title = "Building"
description = "Compiling Voronoi and Dependencies"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "Steps to compile VORONOI on a POSIX-compliant OS or shell"
toc = true
+++


### Required software ###

* Fortran 95 compiler (`gfortran > 4.0; ifort >= 14`)
* C/C++ compiler (`gcc; icc`)
* git
* make
* Cmake >= 3.10.0

### Build steps ###

**1.0.** Ensure that all required dependencies are fulfilled.

### Building PETSc ###

**2.0.** Clone and build PETsc, on the `xsdk-0.2.0` tag:

```sh
git clone https://gitlab.com/petsc/petsc.git petsc
cd petsc
git checkout xsdk-0.2.0
```

**2.1.** Run `./configure` to configure the build with environment-specific settings. As an example,

```
$ python2 './configure' '--with-cc=gcc' '--with-cxx=g++' '--with-fc=gfortran' '--download-mpich' '--download-fblaslapack=1' '--download-hdf5=1'
```

to configure your build with requests to download MPICH, HDF5 and BLAS/LAPACK.

More information can be found [on the official PETSc website](https://www.mcs.anl.gov/petsc/documentation/installation.html).

**2.2.** After `configure` is complete, a message similar to the one below
will be printed to the console:

```
xxx=========================================================================xxx
 Configure stage complete. Now build PETSc libraries with (gnumake build):
   make PETSC_DIR=/opt/atlassian/pipelines/agent/build/petsc PETSC_ARCH=arch-linux2-c-opt all
xxx=========================================================================xxx
```

Export `PETSC_DIR` and `PETSC_ARCH` variables to your environment or
`.bash_profile` (the VORONOI build will need these) and run the `make` command
given at the bottom of the output:

```sh
export PETSC_DIR=/path/to/petsc/
export PETSC_ARCH=petsc-arch-string
export HDF5_DIR=${PETSC_DIR}/${PETSC_ARCH}/lib
make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
```

### Building LaGriT ###

**3.0** Configuring LaGriT

The VORONOI makefile will automatically generate the needed LaGriT libraries
and object files - all it needs to know is the LaGriT source directory.

If you have an existing clone of LaGriT on your system, then export that path
to the environment variable `LAGRIT_DIR`:

```sh
export LAGRIT_DIR=/path/to/LaGriT/
```

If you do not already have an existing LaGriT clone on your system,
then clone it from GitHub:

```sh
git clone https://github.com/lanl/LaGriT.git LaGriT
export LAGRIT_DIR=$(pwd)/LaGriT/
```

### Building VORONOI ###

**4.0.** First, LaGriT libraries must be compiled. Navigate into the `src/` directory 
of the VORONOI repo and run

```sh
make lglibs
```

After this completes successfully, run

```sh
make voronoi
```

The executable will be built in the current working directory and will be named `voronoi`.

## Important Notes

### Using the proper `mpiexec`

When PETSc is compiled in step 2 above, the `mpiexec` binary that is built 
is what should be used to call VORONOI.

This binary is found in `${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec`.

Calling VORONOI is then done with:

```sh
${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -np 4 voronoi [...args]
```

You can also override the system `mpiexec` by doing one of the below:

```sh
# Option 1: set an alias.
alias mpiexec=${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec

# Option 2: prepend the PATH variable.
# *May have side-effects: not recommended!*
export PATH=${PETSC_DIR}/${PETSC_ARCH}/bin:${PATH}
```

------------------------------------
