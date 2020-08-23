+++
title = "Getting Started"
description = "Quick tutorial on getting up and running fast"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "Quick tutorial on getting up and running fast"
toc = true
+++

### Execution Overview ###

Follow the documentation on [building VORONOI](../building/) and run the test
suite to verify build integrity.

General usage for executing VORONOI follows this pattern:

    mpiexec -np N voronoi [commands...]
    
where `N` is the number of MPI nodes desired.

A list of available commands can be found by running:

    voronoi -h

### Examples ###

To generate a VORONOI tessellation in the FEHM format from an AVS mesh, run

    mpiexec -np 4 voronoi -avs my_mesh.inp -type fehm -o my_mesh.stor
    
To build a median tessellation in a PFLOTRAN format with a LaGriT infile, run

    mpiexec -np 4 voronoi -lg infile.lgi -type pflotran -o geom.uge -cv median

------------------------------------