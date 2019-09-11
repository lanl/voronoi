Introduction
============

VORONOI is a parallel 

Citing VORONOI 
---------------
`Satish Karra and Zhuolin Qu. VORONOI: A parallel Delauney to Voronoi mesh convertor for subsurface codes. Technical Report No. LA-UR-XX-XXXXX, Los Alamos National Laboratory, NM, USA.`

Where can one get VORONOI?
-----------------------------
Developer version of VORONOI can be downloaded from ..... 

v1.0 can be download from .....

Installation
------------
VORONOI uses the PETSc toolkit for parallel datastructures, linear solvers, timestepping as well as I/O. 

PETSc
^^^^^
PETSc_ can be installed by following the instructions below:

.. _PETSc: https://www.mcs.anl.gov/petsc/

1. Clone PETSc:

    .. code-block:: bash
        
        git clone https://bitbucket.org/petsc/petsc petsc
        cd petsc
        git checkout xsdk-0.2.0


2. Set PETSC_DIR and PETSC_ARCH environmental variables (see `here <http://www.mcs.anl.gov/petsc/documentation/installation.html#vars>`_). Then type:

    .. code-block:: bash

       ./configure --download-mpich=yes --download-metis=yes --download-parmetis=yes


3. Compile PETSc using the instructions you see on the screen after configuring. You should see something similar to:

    .. code-block:: bash

        xxx=========================================================================xxx
         Configure stage complete. Now build PETSc libraries with (gnumake build):
            make PETSC_DIR=/Users/satkarra/src/petsc-3.7-release PETSC_ARCH=macx-nodebug all
        xxx=========================================================================xxx

__ PETSc_

VORONOI
^^^^^^^^^^
VORONOI can then be installed using make:

.. code-block:: bash
    
    make voronoi


About this  manual
------------------

This manual comprises of information on: algorithms used; installing VORONOI; setting up input files 

Contributors
-------------
- Satish Karra
- Zhuolin Qu

Copyright information
----------------------

