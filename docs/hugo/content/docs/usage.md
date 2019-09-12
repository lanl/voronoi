+++
title = "Usage"
description = "Command line usage for operating VORONOI"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "Command line usage for operating VORONOI"
toc = true
+++

### Command Line Usage ###

| Command                           | Functionality                               |
|-----------------------------------|---------------------------------------------|
| -avs [INFILE.inp]                 | AVS-UCD mesh to read                        |
| -lg [INFILE.lgi]                  | LaGriT infile to run                        |
| -type [fehm,pflotran,tough2,hdf5] | filetype to write to                        |
| -cv [voronoi,median]              | control volume type                         |
| -d                                | write mesh statistics to stdout             |
| -o                                | filepath to save geometric coefficient data |
| -compress                         | coefficient compression (FEHM only)         |
| -dedud                            | epsilon coefficient removal (FEHM only)     |

More detailed information can be run by calling

    voronoi -h

### Examples ###

**1. Converting a Delaunay AVS mesh to an FEHM sparse matrix**

As VORONOI has a native AVS reader, use the flag `-avs [file]` to input any AVS mesh.
Choose the file output type using the `-type [type]` flag. In this case, we are using FEHM.
Finally, the output filename can be chosen using the `-o [output]` argument.

    mpiexec -np 4 voronoi -avs my_mesh.inp -type fehm -o my_mesh.stor

**2. Converting a LaGriT input file CMO to a median mesh**

Given a LaGriT input file, VORONOI can perform a voronoi/median/hybrid calculation
on the current mesh object (CMO) active at the close of file.

Consider the following example, where a tetrahedral Delaunay cube is generated:

`lagrit_test.lgi`:

    * TEST connect (3d) (lagrit_input_connect)
    * LaGriT input file to generate an orthogonal grid on
    * a unit cube. Just change nx,ny,nz to change the resolution.
    *
    define/nx/20
    define/ny/20
    define/nz/20
     
    cmo / create / cmo / / / tet
    rz / xyz / nx, ny, nz / 0. 0. 0. / 1. 1. 1. / 1 1 1
    cmo / setatt / cmo / imt / 1 0 0 / 1
    connect / noadd
    resetpts / itp
     
    * begin compare here
    cmo/status
    cmo/printatt//-all-/minmax
    quality
    finish

The control volume on the final mesh can then calculated with the `-lg [infile]` flag: 

    mpiexec -np 4 voronoi -lg lagrit_test.lgi -type tough2 -o MESH


