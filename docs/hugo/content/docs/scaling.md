+++
title = "Parallel Performance"
description = "Scaling capabilities of VORONOI"
date = 2018-10-22T12:17:12-06:00
weight = 20
draft = false
bref = "Steps to compile VORONOI on a POSIX-compliant OS or shell"
toc = true
+++

### Performance on high-element meshes ###

**Example 1**

Example name: TSA250_9
Nodes: 40,127,977
Elements: 82,304,132
Element type: triangle
Input file size: 16583177171
Time spent on reading: 4m1s
Time spent on coefficient calculations: 30.01 s
Time spent on writing: 10m6s
Total time: 14m35s
Maximum memory usage: 1.134 GB <-- check this
MPI cores used: 64

*All examples were run on a __ processor with 88 cores and 252 GB RAM; x86_64 GNU/Linux Ubuntu 16.04*

------------------------------------
