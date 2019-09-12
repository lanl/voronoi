+++
title = "About"
description = "Overview of Voronoi"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "VORONOI was built from the ground up for subsurface flow and transport simulators; with fast, accurate conversion of unstructured finite element meshes into Voronoi and median tessellations"
toc = true
+++

### Overview ###

VORONOI is a parallel, scalable control volume tessellation generator, built for use in subsurface flow and transport solvers.

The program accepts a Delaunay mesh composed of one of four element types: triangles, tetrahedrals, quads, or hexes.

It then finds the Voronoi dual or median of that mesh, and writes it in a solver-specific format:

1. FEHM (`.stor`)
2. PFLOTRAN (`.uge`)
3. TOUGH (`MESH`)
4. HDF5 (for general use; `.h5`)

### Features ###

* Accurate calculation of Voronoi/Median cell volume and cell face area in 2D and 3D
* Reads AVS and LaGriT input files (.lgi)
* Outputs to FEHM, PFLOTRAN, TOUGH, and HDF5
* Outputs visualization of control volume cells
* Assigns ROCK coefficients based on node attributes (TOUGH)
* Automatic element assignment of isotropic / anisotropic permeability coefficients (TOUGH)
* Removes duplicate and negligable coefficients (FEHM)
* Built with PETSc for parallel execution and sparse matrix data types
* Cross-platform & open source