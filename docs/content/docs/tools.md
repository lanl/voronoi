+++
title = "Tools"
description = "Supplementary utilities to help visualize and process results"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "Supplementary utilities to help visualize and process results"
toc = true
+++

### Visualizing the sparse matrix ###

*Found in:* `viewSparseMatrix.py`

This script visualizes a sparse matrix for a given mesh -
this is a handy way to quickly view the i,j connections between nodes.

Voronoi cell volumes are represented on the matrix diagonal, while the
off-diagonals represent: (Voronoi area_ij/interface_ij length).

In the current implementation, an NxN sparse matrix is converted to an
NxN image - be mindful of very large meshes as to not run into memory
issues.

References:

[1] https://ieeexplore.ieee.org/document/112361/ <br/>
[2] https://en.wikipedia.org/wiki/Sparse_matrix

### Mesh conversion to AVS ###

*Found in:* `createMeshFromData.py`

This script converts a file of node (x,y,z) values and element (i,j,k,...) indices into
an AVS-UCD mesh.

Node input file should be in the form:

    x0 y0 z0
    x1 y1 z1
    ...
    xN yN zN

Element input file should be in the form:

    i1 j1 k1 ...
    i2 j2 k2 ...
    ...
    iN jN kN ...

A delimiter between entries may be specified with the `-d` argument. Defaults to space (`' '`).

It is recommended that the element file reference the node list with 1-based indexing; that is,
the element input file should reference the node `x0 y0 z0` as `1`.
If you use zero- or N-based indexing, use the `--index` flag to indicate this.

This script will automatically assume that an element file with 3 integers on a line refers to
a triangle element, and 4 integers refers to tetrahedrons. If you wish to manually specify what the
element type is (i.e., quad) then use `--type ['tet', 'tri', 'quad', 'hex']`. Note that only one element
type per file is allowed - mixing of element types in a single mesh is not supported in Voronoi.

If you only have a list of nodes, this script will still write out a file - but no element connectivity
will be defined. You can import the mesh into LaGriT and triangulate it, or use the SciPy Delaunay function.

If you wish to view the created mesh, it is recommended that you use ParaView: https://www.paraview.org

### Testing strong / weak MPI scaling ###

*Found in:* `viewMPIscaling.py`

This script tests against strong or weak scaling, and formats the results into:

1. A set of .PNG plots, generated with Matplotlib, and/or
2. A table of the form

```
    =================================
       n         t(n)       E_fehm(n) ...
    ---------------------------------
       1        9.763s     1.0000
       2        4.901s     0.7540
       4        2.731s     0.7540
       8        1.220s     0.7540
       16       0.813s     0.7540
       32       0.763s     0.7540
       64       0.763s     0.7540
```

