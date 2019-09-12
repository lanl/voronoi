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

```
$ git clone https://gitlab.com/petsc/petsc.git petsc
$ cd petsc
$ git checkout xsdk-0.2.0
```

**2.1.** Run `./configure` to configure the build with environment-specific settings. As an example,

```
$ ./configure --download-mpich=yes --download-metis=yes --download-parmetis=yes
```

to configure your build with requests to download MPICH, METIS, and PARMETIS.

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

```
$ export PETSC_DIR=/path/to/petsc/
$ export PETSC_ARCH=petsc-arch-string
$ make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
```

### Building HDF5 ###

Download the HDF5 source code from the [official website](https://www.hdfgroup.org/downloads/hdf5/source-code/).

You can also download via wget or cURL:

```
$ wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
```

Unpackage using:

```
$ gunzip < hdf5-1.10.5.tar.gz | tar xf -
$ cd hdf5-1.10.5
$ ./configure --prefix=/usr/local/hdf5 --enable-fortran
$ make
$ make check                # run test suite.
$ make install
$ make check-install        # verify installation.
```

where `--prefix=/usr/local/hdf5` is the build location. You may set this to any arbitrary location.

### Building LaGriT ###

**3.0** Configuring LaGriT

The VORONOI makefile will automatically generate the needed LaGriT libraries
and object files - all it needs to know is the LaGriT source directory.

If you have an existing clone of LaGriT on your system, then export that path
to the environment variable `LAGRIT_DIR`:

```
$ export LAGRIT_DIR=/path/to/LaGriT/
```

If you do not already have an existing LaGriT clone on your system,
then clone it from GitHub:

    $ git clone https://github.com/lanl/LaGriT.git LaGriT
    $ export LAGRIT_DIR=$(pwd)/LaGriT/

**3.0** Compiling LaGriT

To build the LaGriT libraries, `cd` into the `voronoi/src/` directory and run

```
make lglibs
```

### Building VORONOI ###

**4.0.** Build VORONOI by simply running

```
$ make voronoi
```

Your executable will be built in the repo `src/` dir; or, in other words, `$VORONOI_SRC_DIR/src/voronoi`.

------------------------------------
