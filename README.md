![alt text](documentation/logo.png "Voronoi")

[![pipeline status](https://gitlab.com/sft-meshing/voronoi/badges/master/pipeline.svg)](https://gitlab.com/sft-meshing/voronoi/commits/master)

VORONOI is a parallel, scalable control volume tessellation generator, built for use in subsurface flow and transport solvers.

The program accepts a Delaunay mesh composed of one of four element types: triangles, tetrahedrals, quads, or hexes.

It then finds the Voronoi dual or median of that mesh, and writes it in a solver-specific format:

1. FEHM (`.stor`)
2. PFLOTRAN (`.uge`)
3. TOUGH (`MESH`)
4. HDF5 (for general use; `.h5`)

VORONOI works as a standalone MPI application or as a serial command in the
[LaGriT mesh generator](lagrit.lanl.gov).

### Features ###

* Accurate calculation of Voronoi/Median cell volume and cell face area in 2D and 3D
* Reads AVS and LaGriT input files (LGI)
* Outputs to FEHM, PFLOTRAN, TOUGH, and HDF5
* Outputs visualization of control volume cells
* Assigns ROCK coefficients based on node attributes (TOUGH)
* Automatic element assignment of isotropic / anisotropic permeability coefficients (TOUGH)
* Removes duplicate and negligable coefficients (FEHM)
* Built with PETSc for parallel execution and sparse matrix data types
Cross-platform & open source

### Building ###

For building instructions, see the [Build page](http://lanl.github.io/voronoi/build.html).

### License ###

VORONOI is open-source software licensed under the 3-Clause BSD License.

*This is open source software; you can redistribute it and/or modify it under
the terms of the BSD-3 License. If software is modified to produce derivative
works, such modified software should be clearly marked, so as not to confuse
it with the version available from LANL. Full text of the BSD-3 License can be
found in the LICENSE.md file of the repository.*

### Contact ###

For help and support, contact Daniel Livingston (livingston@lanl.gov), Carl Gable (gable@lanl.gov), or Satish Karra (satkarra@lanl.gov).

To report a bug or feature request, open an issue on the GitHub Issue Tracker.



