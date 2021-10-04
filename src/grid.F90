module Grid_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

   use Grid_Aux_module
   use Voronoi_Constants_module

   implicit none

   private

   public :: GridRead, &
             mpi_print, &
             ScatterGatherCoordinates, &
             AllocateConnectMatrix, &
             CalculateGeometricCoeff, &
             CalculateGeometricCoeff3D, &
             CreateEdgeMatrix, &
             GridWriteFEHM, &
             GridWriteTOUGH2, &
             GridWriteHDF5, &
             GridWritePFLOTRAN, &
             CreateConnMatrix, &
             ReadSTORFile, &
             PrintInitialConditions, &
             PrintMeshAttributes, &
             Diagnostics

contains

   subroutine PrintInitialConditions(grid, atts, rank, size)

#include "petsc/finclude/petscvec.h"
      use petscvec

      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      character(len=1) :: h1 = "="
      character(len=150) :: input_string
      character(len=50) :: l
      character(len=18) :: ftype
      character(len=7) :: cvtype
      integer :: w1 = 65
      PetscInt :: rank, size, io_rank = 0

      if (rank == io_rank) then

45       format(4X, A41, 1x, A20)
46       format(4X, A51, 1x, i10)
49       format(4X, A60, 1x, L1)

         l = repeat(' ', 50)

         select case (grid%outtype)
         case (1)
            ftype = '    FEHM STOR file'
         case (2)
            ftype = '   TOUGH MESH file'
         case (3)
            ftype = 'PFLOTRAN grid file'
         case (4)
            ftype = '         HDF5 file'
         case (5)
            ftype = '  Stanford Polygon'
         case default
            ftype = '   UNKNOWN (ERR 1)'
         end select

         select case (grid%cv_type)
         case (1)
            cvtype = 'Voronoi'
         case (2)
            cvtype = ' Median'
         case (3)
            cvtype = ' Hybrid'
         case default
            cvtype = 'UNKNOWN'
         end select

         if (grid%avs_flag .EQV. PETSC_TRUE) then
            input_string = grid%avs_str
         else if (grid%lg_flag .EQV. PETSC_TRUE) then
            input_string = grid%lg_str
         else
            input_string = 'ERROR: Undefined input'
         endif

         print *, char(10), 'PARAMETERS'
         print *, repeat(h1, w1)
         print 46, 'MPI Nodes:'//l, size
         print 45, 'Input:'//l, trim(input_string)
         print 45, 'Output:'//l, trim(grid%dump_str)
         print 45, 'Filetype:'//l, ftype
         print 45, 'Control volume:'//l, cvtype

         if (grid%outtype == 1) then
            print *, char(10), 'FEHM SETTINGS'
            print *, repeat(h1, w1)
            print 49, 'Coefficient compression:'//l, atts%coefficients_are_compressed
            print 49, 'Coefficient dedudding:'//l, atts%coefficients_are_dedudded
         else if (grid%outtype == 2) then
            print *, char(10), 'TOUGH SETTINGS'
            print *, repeat(h1, w1)
            print '(a)', 'Materials: nil'
         endif

      endif

   end subroutine PrintInitialConditions

   subroutine PrintMeshAttributes(grid, rank, size)

      implicit none

      type(grid_type) :: grid
      character(len=1) :: h1 = "="
      character(len=50) :: l
      character(len=11) :: etype
      integer :: w1 = 65
      PetscInt :: rank, size, io_rank = 0

      if (rank == io_rank) then

45       format(4X, A41, 1x, A20)
46       format(4X, A51, 1x, i10)

         l = repeat(' ', 50)

         select case (grid%ndim)
         case (2)
            etype = '   triangle'
         case (3)
            etype = 'tetrahedron'
         case default
            etype = '    UNKNOWN'
         end select

         print *, char(10), 'MESH ATTRIBUTES'
         print *, repeat(h1, w1)
         print 46, 'Nodes:'//l, grid%num_pts_global
         print 46, 'Elements:'//l, grid%num_elems_global
         print 45, 'Element type:'//l, etype
         !print 45,'Attributes (node):'//l,'imt1 ipopt askem'
         !print 45,'Attributes (cell):'//l,'imt1 ipopt'

      endif

   end subroutine PrintMeshAttributes

   subroutine mpi_print(message, rank, io_rank)
#include "petsc/finclude/petscvec.h"
      use petscvec

      implicit none

      ! implement io level as optional
      PetscInt :: rank, io_rank, zero, verbosity
      character(len=*) :: message

      common verbosity

      zero = 0

      if (verbosity == 0) return
      if (rank == io_rank) print *, message

   end subroutine mpi_print

   integer function ParseDimension(infile, length)
      !
      ! Scans an AVS file for element types,
      ! and returns dimensionality as an integer
      !
      ! Author: Daniel Livingston
      ! Date: 10/18/2017
      !

      implicit none

      character(len=150) :: infile
      character(len=20), dimension(3) :: string
      integer :: length, fileid
      integer :: nnodes, n, dim
      integer :: data(5)

      ! Open file and read header
      fileid = 61
      open (fileid, file=infile(1:length))
      read (fileid, *) data

      ! read point count
      nnodes = data(1)

      ! skip past nodes
      do n = 1, nnodes
         read (fileid, *)
      enddo

      ! This is the first line with cells -
      ! the third column here tells the cell type
      read (fileid, *) string

      ! Detect cell type
      if (trim(string(3)) == 'tet') then
         dim = 3
         !print*,'   Element type: tetrahedron'
      else if (trim(string(3)) == 'tri') then
         dim = 2
         !print*,'   Element type: triangle'
      else
         dim = 2
         !print*,'   Element type: unknown'
      end if

      close (fileid)

      ParseDimension = dim

   end function ParseDimension

! ************************************************************************** !

   subroutine ReadSTORFile(grid, rank, size, infile)
      !
      ! Reads a STOR file into a grid object
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 08/27/2018
      !

#include "petsc/finclude/petscvec.h"
      use petscvec

      implicit none

      type(grid_type) :: grid
      character(len=100) :: infile, trash
      PetscInt :: i, j, k, delta, rank, size, fileid = 432, io_rank = 0
      PetscReal, allocatable :: voronoi_volumes(:), real_line(:), geom_coeffs(:)
      PetscInt, allocatable :: row_count(:), integer_line(:), row_entries(:), coeff_pointers(:)
      PetscInt :: ncoefs, neq, ncoef_sum, num_area_coef, ncon_max, coeff_count, junk_lines, nwc

      if (rank == io_rank) then

         ! ==================================
         ! 1. Begin reading file
         ! ==================================

         ! FEHM STOR has a five-column format. Hence the step size.
         allocate (integer_line(5))
         allocate (real_line(5))

         ! Open the input file ... duh
         open (unit=fileid, file=trim(infile))

         ! Read and ignore comments
         read (fileid, *) trash  ! fehmstor ascir8i4 LaGriT Sparse Matrix Voronoi Coefficients
         read (fileid, *) trash  ! Tue Jul 31 15:18:23 2018

         ! Grab the five-integer header
         read (fileid, *) nwc, neq, ncoef_sum, num_area_coef, ncon_max
         ncoefs = ncoef_sum - (neq + 1)

         if (num_area_coef /= 1) then
            print *, 'ERROR PARSING STOR FILE: non-scalar areas are not supported'
            call EXIT(1) ! standardize error codes plz
         endif

         ! Read in the Voronoi volumes
         allocate (voronoi_volumes(neq))
         do i = 1, neq, 5
            j = i + 4
            if (j > neq) j = neq
            read (fileid, *) voronoi_volumes(i:j)
         enddo

         ! Read in the row count - i.e. what columns are represented in each row
         allocate (row_count(neq + 1))
         do i = 1, neq + 1, 5
            j = i + 4
            k = 5

            ! Row count may bleed into the next block - compensate for this
            if (j > neq + 1) then
               j = neq + 1
               k = (neq + 1) - i + 1
            endif

            read (fileid, *) integer_line
            row_count(i:j) = integer_line(1:k)
         enddo

         ! Read in the explicitly represented row-column pairs in the sparse matrix
         allocate (row_entries(ncoefs))

         ! If there are any row entries that bled into the last captured line,
         ! then we want to parse them into the array before reading again
         delta = 5 - k
         if (delta /= 0) row_entries(1:delta) = integer_line(k + 1:5)

         do i = 1 + delta, ncoefs, 5
            j = i + 4
            if (j > ncoefs) j = ncoefs
            read (fileid, *) row_entries(i:j)
         enddo

         ! This array contains 'pointers' to the geometric coeffs. array
         coeff_count = ncoef_sum + neq
         allocate (coeff_pointers(ncoefs))

         do i = 1, ncoefs, 5
            j = i + 4
            k = 5
            if (j > ncoefs) then
               j = ncoefs
               k = modulo(j, 5)
            endif
            read (fileid, *) integer_line
            coeff_pointers(i:j) = integer_line(1:k)
         enddo

         ! Calculate how many lines of 0 will be present next
         ! There are NEQ+1 0's
         ! We also want to ignore the NEQ pointers to Vor. volumes
         junk_lines = ceiling((neq + 1 - (5 - k))/5.0) + ceiling(neq/5.0)
         do i = 1, junk_lines
            read (fileid, *) trash
         enddo

         ! Read the geometric coefficients into an array
         allocate (geom_coeffs(nwc))

         do i = 1, nwc, 5
            j = i + 4
            k = 5
            if (j > nwc) then
               j = nwc
               k = modulo(j, 5)
            endif
            read (fileid, *) real_line(1:k)
            geom_coeffs(i:j) = real_line(1:k)
         enddo

         print *, geom_coeffs

         ! ==================================
         ! 2. Pack data into a grid object
         ! ==================================

         ! TODO
         ! Configure vectors and Mats that GridWrite depends on
         ! grid%degree_tot,grid%edgematrix,grid%adjmatrix_full,grid%vertex_ids
         ! grid%degree, grid%num_pts_local, grid%num_pts_total
         ! It may (probably will) just be easier to refactor the write subroutine

         !grid%num_pts_local = 0
         !grid%num_pts_global = 0

         !d_nnz1 = (/(i,i=grid%num_pts_local,1,-1)/)

         !call VecGetArrayF90(grid%connections,vec_ptr,ierr);CHKERRQ(ierr)
         !d_nnz = min(d_nnz1,int(vec_ptr))
         !o_nnz = min(int(vec_ptr),grid%num_pts_global-grid%num_pts_local)
         !d_nnz2 = min(int(vec_ptr),grid%num_pts_local)
         !o_nnz2 = min(int(vec_ptr),grid%num_pts_global-grid%num_pts_local)
         !call VecRestoreArrayF90(grid%connections,vec_ptr,ierr);CHKERRQ(ierr)

         !d_nnz = 0
         !o_nnz = 0
         !d_nnz2 = 0
         !o_nnz2 = 0

         ! Create the adjacency matrix
         !call MatCreateAIJ(PETSC_COMM_WORLD,grid%num_pts_local,grid%num_pts_local, &
         !                                      PETSC_DETERMINE,PETSC_DETERMINE, &
         !                                      PETSC_NULL_INTEGER,d_nnz, &
         !                                      PETSC_NULL_INTEGER,o_nnz, &
         !                                      grid%adjmatrix,ierr);CHKERRQ(ierr)

         ! Clear allocations
         deallocate (integer_line)
         deallocate (real_line)
         deallocate (voronoi_volumes)
         deallocate (row_count)
         deallocate (row_entries)
         deallocate (geom_coeffs)
         deallocate (coeff_pointers)

      endif

   end subroutine ReadSTORFile

! ************************************************************************** !
   subroutine GridRead(grid, rank, size)
      !
      ! Reads a grid in avs format
      !
      ! Author: Satish Karra and Zhuolin Qu, LANL
      ! Date: 07/23/2015
      !

#include "petsc/finclude/petscvec.h"
      use petscvec

      implicit none

      integer :: imt_index
      PetscReal :: tmp_att(5)
      PetscReal, dimension(:), allocatable :: imt_vector

      type(grid_type) :: grid
      character(len=3) :: tcell
      character(len=MAXSTRINGLENGTH) :: cmo
      character(len=5) :: tmp_str
      PetscInt :: fileid, ivertex, iatt, irank, remainder
      PetscInt :: temp_int1, temp_int2, temp_int3, num_to_read, num_to_read_save, ielem
      PetscInt :: data(5)
      PetscInt :: ifrac
      PetscInt :: num_elems, num_pts, num_atts
      PetscInt :: num_pts_local, num_pts_local_save
      PetscInt :: num_elems_local, num_elems_local_save
      PetscErrorCode :: ierr
      PetscReal, allocatable :: temp_real_array(:), temp_real(:)
      PetscInt, allocatable :: temp_int_array(:), temp_int(:)
      PetscInt :: rank, size
      PetscInt :: io_rank = 0
      PetscReal, pointer :: vec_ptr(:)
      PetscInt, allocatable :: pos(:)
      PetscInt, allocatable :: sendcounts(:)
      PetscInt :: length
      PetscInt :: dim_skip
      PetscInt :: ierror, ilen, itype
      integer*8 :: nsdtopo_8, nsdgeom_8, npoints_8, ntets_8 ! For getting LaGrit 64-bit integer data
      PetscInt :: nsdtopo, nsdgeom, npoints, ntets

      pointer(ixic, xic)
      pointer(iyic, yic)
      pointer(izic, zic)
      pointer(ipitet, itet_8)
      pointer(ipitettyp, itettyp_8)
      pointer(ipimt1, imt1_8)

      PetscReal xic(*), yic(*), zic(*)
      integer*8 itet_8(*), itettyp_8(*), imt1_8(*)
      PetscInt, allocatable :: itet(:), itettyp(:)
      PetscReal, allocatable :: imt1(:)

      if (rank == io_rank) then

         ! Check if LaGriT infile was passed
         if (grid%lg_flag .eqv. PETSC_TRUE) then

            if (grid%avs_flag .eqv. PETSC_TRUE) then
               print *, '\nWARNING: Both AVS and LaGriT infiles passed to Voronoi.'
               print *, '         Defaulting to LaGriT.\n'
            endif

            call initlagrit('silent', ' ', ' ')
            !call dotask('read/avs/' // trim(grid%avs_str) // /mo1'; finish',ierror)
            call dotask('infile '//trim(grid%lg_str)//'; finish', ierror)

            call cmo_get_name(cmo, ierror)
            call cmo_get_info('ndimensions_topo', cmo, nsdtopo_8, ilen, itype, ierror)
            call cmo_get_info('ndimensions_geom', cmo, nsdgeom_8, ilen, itype, ierror)
            call cmo_get_info('nnodes', cmo, npoints_8, ilen, itype, ierror)
            call cmo_get_info('xic', cmo, ixic, ilen, itype, ierror)
            call cmo_get_info('yic', cmo, iyic, ilen, itype, ierror)
            call cmo_get_info('zic', cmo, izic, ilen, itype, ierror)
            call cmo_get_info('nelements', cmo, ntets_8, ilen, itype, ierror)
            call cmo_get_info('itet', cmo, ipitet, ilen, itype, ierror)
            call cmo_get_info('itettyp', cmo, ipitettyp, ilen, itype, ierror)
            call cmo_get_info('imt1', cmo, ipimt1, ilen, itype, ierr)

            ! Re-cast LaGrit integer*8 to PetscInt
            nsdtopo = nsdtopo_8
            nsdgeom = nsdgeom_8
            npoints = npoints_8
            ntets = ntets_8

            ! TODO: Dimension should be a parsing of itettyp
            grid%ndim = nsdtopo

            allocate (itet(ntets*(grid%ndim + 1)))
            allocate (itettyp(ntets*(grid%ndim + 1)))
            allocate (imt1(npoints))

            ! Re-cast integer*8 to PestcInt
            itet(1:ntets*(grid%ndim + 1)) = itet_8(1:ntets*(grid%ndim + 1))
            itettyp(1:ntets) = itettyp_8(1:ntets)
            imt1(1:npoints) = imt1_8(1:npoints)

            temp_int1 = npoints
            temp_int2 = ntets

#if DEBUG
            call dotask('mmprint; finish', ierror) ! LG memory info
#endif
         else
            length = len(trim(grid%avs_str))
            grid%ndim = ParseDimension(grid%avs_str, length)

            fileid = 86
            open (fileid, file=trim(grid%avs_str))
            read (fileid, *) data

            temp_int1 = data(1) ! read point count
            temp_int2 = data(2) ! read simplex count
            temp_int3 = data(3) ! read node attribute count
         endif

#if DEBUG
         print *, 'Rank:', rank
         print *, 'Number of points:', temp_int1
         print *, 'Number of elements:', temp_int2
         print *, 'Dimension:', grid%ndim
#endif
      endif

      ! Broadcast the number of dimensions to the rest of the nodes
      call MPI_Bcast(grid%ndim, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! The number of indexes each simplex requires
      ! E.g., the 2D case: [1][n1][n2][n3][2][n1][n2][n3]
      dim_skip = grid%ndim + 2
      call MPI_Bcast(temp_int1, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! Assign grid->num_pts_global to point count
      num_pts = temp_int1
      grid%num_pts_global = num_pts

      call MPI_Bcast(temp_int2, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! Assign grid->num_elems_global to simplex count
      num_elems = temp_int2
      grid%num_elems_global = num_elems

      ! Broadcast node attribute count
      call MPI_Bcast(temp_int3, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      num_atts = temp_int3

#if DEBUG
      print *, 'Rank:', rank, '       ', 'Total number of points:', num_pts
      print *, 'Rank:', rank, '       ', 'Total number of elements:', num_elems
#endif

      ! divide vertices across ranks
      num_pts_local = num_pts/size
      num_pts_local_save = num_pts_local
      remainder = num_pts - num_pts_local*size
      if (rank < remainder) num_pts_local = num_pts_local + 1

      grid%num_pts_local = num_pts_local

#if DEBUG
      print *, 'Rank:', rank, '       ', 'Local number of points:', num_pts_local
#endif

      ! create a parallel vector with a stride of 3, storing the coordinates
      call VecCreate(PETSC_COMM_WORLD, grid%coordinates_vec, ierr); CHKERRQ(ierr)
      call VecSetType(grid%coordinates_vec, VECMPI, ierr); CHKERRQ(ierr)
      call VecSetSizes(grid%coordinates_vec, grid%num_pts_local*3, PETSC_DECIDE, ierr); CHKERRQ(ierr)
      call VecSetBlockSize(grid%coordinates_vec, 3, ierr); CHKERRQ(ierr)

      ! create a PETSc vector to store Voronoi cell centers
      call VecCreate(PETSC_COMM_WORLD, grid%lcc, ierr); CHKERRQ(ierr)
      call VecSetType(grid%lcc, VECMPI, ierr); CHKERRQ(ierr)
      call VecSetSizes(grid%lcc, grid%num_pts_local*3, PETSC_DECIDE, ierr); CHKERRQ(ierr)
      call VecSetBlockSize(grid%lcc, 3, ierr); CHKERRQ(ierr)

!=====================================================================
!  Read the vertex coordinates
!  -- record the coordinates in a vector to prepare for scatter/gather
!=====================================================================

      allocate (grid%vertex_ids(grid%num_pts_local))
      allocate (sendcounts(size))
      allocate (pos(size))
      grid%vertex_ids = 0

      allocate (temp_real_array(grid%num_pts_global*4))

      ! Read coordinates through io_rank and pass to other ranks
      if (rank == io_rank) then
         !allocate(temp_real_array(grid%num_pts_global*4))
         allocate (temp_real(4))

         ! Read coordinates and store
         do ivertex = 1, grid%num_pts_global

            if (grid%lg_flag .eqv. PETSC_TRUE) then
               temp_real = (/float(ivertex), real(xic(ivertex)), real(yic(ivertex)), real(zic(ivertex))/)
            else
               read (fileid, *) temp_real
            endif

            temp_real_array((ivertex - 1)*4 + 1:(ivertex - 1)*4 + 4) = temp_real

         enddo
         deallocate (temp_real)
      endif

      call MPI_Bcast(temp_real_array, grid%num_pts_global*4, MPI_DOUBLE, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      num_to_read_save = grid%num_pts_global/size
      remainder = grid%num_pts_global - num_to_read_save*size

      pos = 0
      do irank = 0, size - 2
         num_to_read = num_pts_local_save
         if (irank < remainder) num_to_read = num_to_read_save + 1
         pos(irank + 2) = pos(irank + 1) + num_to_read*4
      enddo

      do irank = 0, size - 1
         num_to_read = num_pts_local_save
         if (irank < remainder) num_to_read = num_to_read_save + 1
         sendcounts(irank + 1) = num_to_read*4
      enddo

      allocate (temp_real(grid%num_pts_local*4))
      !if (rank == io_rank) &
      !  call MPI_Scatterv(temp_real_array,sendcounts,pos,MPI_DOUBLE_PRECISION,temp_real,&
      !       grid%num_pts_global*4,MPI_DOUBLE_PRECISION,io_rank,MPI_COMM_WORLD,ierr);CHKERRQ(ierr)

      ! TODO: This is the 'poor man's MPI_Scatterv.
      ! This should be replaced with Scatterv soon.
      temp_real = temp_real_array(pos(rank + 1) + 1:pos(rank + 1) + sendcounts(rank + 1))

      if (rank == io_rank) deallocate (temp_real_array)
      call VecGetArrayF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)

      ! Store vertex ids and coordinates
      do ivertex = 1, grid%num_pts_local
         grid%vertex_ids(ivertex) = int(temp_real((ivertex - 1)*4 + 1))
         vec_ptr((ivertex - 1)*3 + 1:(ivertex - 1)*3 + 3) = temp_real((ivertex - 1)*4 + 2:(ivertex - 1)*4 + 4)
      enddo

      call VecRestoreArrayF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)
      deallocate (temp_real)
      deallocate (pos)
      deallocate (sendcounts)

#if DEBUG
      write (string, *) grid%num_pts_local
      write (word, *) rank
      string = trim(adjustl(string))//' vertices stored on p'//trim(adjustl(word))
      print *, trim(string)
#endif

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'vertices.out', viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%coordinates_vec, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

!===========================================
!  Read the element connectivity
!  -- recored in a vector and scatter/gather
!===========================================

      !divide elements across ranks
      num_elems_local = num_elems/size
      num_elems_local_save = num_elems_local
      remainder = num_elems - num_elems_local*size
      if (rank < remainder) num_elems_local = num_elems_local + 1

      grid%num_elems_local = num_elems_local

#if DEBUG
      print *, 'Rank:', rank, '       ', 'Local number of elements:', num_elems_local
#endif

      allocate (grid%elem_ids(grid%num_elems_local))
      allocate (grid%elem_connectivity(grid%num_elems_local, grid%ndim + 1))
      allocate (sendcounts(size))
      allocate (pos(size))

      ! Initialize vectors for PFLOTRAN cell centers and volumes
      if (grid%outtype == 3) allocate (grid%cell_cc(grid%num_elems_local, 3))
      if (grid%outtype == 3) allocate (grid%cell_vol(grid%num_elems_local))

      grid%elem_ids = 0

      ! read the connectivity through io_rank and pass to other ranks
      if (rank == io_rank) then

         ! allocate arrays to hold all tet/tri's, and the current file-read tet/tri
         allocate (temp_int_array(grid%num_elems_global*dim_skip))
         allocate (temp_int(dim_skip))

         do ielem = 1, grid%num_elems_global

            ! Construct grid matrix from LaGriT pointers...
            if (grid%lg_flag .eqv. PETSC_TRUE) then
               if (grid%ndim == 2) then
                  temp_int = (/ielem, &
                               int(itet(3*(ielem - 1) + 1)), &
                               int(itet(3*(ielem - 1) + 2)), &
                               int(itet(3*(ielem - 1) + 3))/)
               else
                  temp_int = (/ielem, &
                               int(itet(4*(ielem - 1) + 1)), &
                               int(itet(4*(ielem - 1) + 2)), &
                               int(itet(4*(ielem - 1) + 4)), &
                               int(itet(4*(ielem - 1) + 3))/)
               endif
               ! Or, construct from AVS infile
            else
               read (fileid, *) temp_int(1), ifrac, tcell, temp_int(2:grid%ndim + 2)
            endif

            temp_int_array((ielem - 1)*dim_skip + 1:(ielem - 1)*dim_skip + dim_skip) = temp_int

         enddo

         deallocate (temp_int)
      endif

      num_to_read_save = grid%num_elems_global/size
      remainder = grid%num_elems_global - num_to_read_save*size

      pos = 0
      do irank = 0, size - 2
         num_to_read = num_to_read_save
         if (irank < remainder) num_to_read = num_to_read_save + 1
         pos(irank + 2) = pos(irank + 1) + num_to_read*dim_skip!4
      enddo

      do irank = 0, size - 1
         num_to_read = num_to_read_save
         if (irank < remainder) num_to_read = num_to_read_save + 1
         sendcounts(irank + 1) = num_to_read*dim_skip
      enddo

      allocate (temp_int(grid%num_elems_local*dim_skip))
      call MPI_Scatterv(temp_int_array, sendcounts, pos, MPI_INT, &
                        temp_int, sendcounts, MPI_INT, 0, &
                        MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) deallocate (temp_int_array)

      ! Assign simplex # to elem_ids; simplex nodes to elem_connectivity
      do ielem = 1, grid%num_elems_local
         grid%elem_ids(ielem) = temp_int((ielem - 1)*dim_skip + 1)
         grid%elem_connectivity(ielem, :) = temp_int((ielem - 1)*dim_skip + 2:(ielem - 1)*dim_skip + dim_skip)
      enddo

      deallocate (temp_int)
      deallocate (pos)
      deallocate (sendcounts)

      if (allocated(itet)) deallocate (itet)
      if (allocated(itettyp)) deallocate (itettyp)

#if DEBUG
      do ielem = 1, grid%num_elems_local
         print *, rank, grid%elem_ids(ielem), grid%elem_connectivity(ielem, :)
      enddo
#endif

      !===========================================
      !  Read the imt values for a mesh
      !  -- convert to a character array and store
      !===========================================

      allocate (grid%imt_values(num_pts))

      if (rank == io_rank) then
         ! if we're reading an AVS file...
         if (grid%lg_flag .eqv. PETSC_TRUE) then
            call getMatFromIMT(imt1, num_pts, grid%imt_values)
         else
            ! if there is no zone flag, AND we're in TOUGH2 mode...
            if ((grid%zone_flag .EQV. PETSC_FALSE) .AND. (grid%is_tough .EQV. PETSC_TRUE)) then

               ! if the AVS file actually HAS attributes...
               if (num_atts .gt. 0) then
                  allocate (imt_vector(num_pts))
                  imt_index = 0

                  ! figure out which column is imt
                  do iatt = 1, num_atts + 1
                     read (fileid, *) tmp_str
                     if ((tmp_str(1:3) == 'imt') .AND. (imt_index .lt. 1)) then
                        imt_index = iatt
                     endif
                  enddo

                  ! verify that there actually *is* an imt field
                  if (imt_index .gt. 0) then
                     ! iterate over rows & assign the correct column to imt_vector
                     do iatt = 1, num_pts
                        read (fileid, *) tmp_att(1), tmp_att(2), tmp_att(3), tmp_att(4), tmp_att(5)
                        imt_vector(iatt) = tmp_att(imt_index)
                     enddo

                     ! convert to a materials string
                     call getMatFromIMT(imt_vector, num_pts, grid%imt_values)
                  else
                     grid%imt_values = '    1'
                  endif

                  deallocate (imt_vector)
               ! if num_atts <= 0
               else
                  grid%imt_values = '    1'
               endif
            endif
            close (fileid)
         endif
      endif

      if (allocated(imt1)) deallocate (imt1)

   end subroutine GridRead

! ************************************************************************** !
   subroutine ScatterGatherCoordinates(grid, rank, size)
      !
      ! Counts the local vertices that are needed for calculations,
      ! and scatter/gather from the global parallel vector of coordinates
      ! to local needed coordinates
      ! & counting the connections meanwhile
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/21/2015
      !

#include "petsc/finclude/petscvec.h"
      use petscvec

      implicit none

      type(grid_type) :: grid
      PetscInt :: rank, size, i
      PetscInt :: temp(grid%ndim + 1)
      PetscScalar :: one = 1, two = 2
      PetscInt, pointer :: int_array_pointer(:)
      PetscInt, allocatable :: int_array(:), int_array2(:)
      PetscInt, allocatable :: int_array3(:), int_array4(:)
      PetscInt :: vertex_count, count, vertex_id
      PetscErrorCode :: ierr
      VecScatter :: vec_scatter
      IS :: is_scatter, is_gather
      PetscInt :: nd1 ! number_of_dimensions + 1

      nd1 = grid%ndim + 1

      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, grid%num_pts_global, grid%connections, ierr); CHKERRQ(ierr)
      call VecSet(grid%connections, one, ierr); CHKERRQ(ierr)

      vertex_count = grid%num_elems_local*nd1
      allocate (int_array_pointer(vertex_count))

      ! put all the connectivity into a single vector
      ! and calcuate a rough estimate of the connections
      do i = 1, grid%num_elems_local
         temp = grid%elem_connectivity(i, :)
         int_array_pointer((i - 1)*nd1 + 1:(i - 1)*nd1 + nd1) = temp

         call VecSetValue(grid%connections, temp(1) - 1, two, ADD_VALUES, ierr); CHKERRQ(ierr)
         call VecSetValue(grid%connections, temp(2) - 1, two, ADD_VALUES, ierr); CHKERRQ(ierr)
         call VecSetValue(grid%connections, temp(3) - 1, two, ADD_VALUES, ierr); CHKERRQ(ierr)
         if (grid%ndim == 3) call VecSetValue(grid%connections, temp(4) - 1, two, ADD_VALUES, ierr); CHKERRQ(ierr)
      enddo

      call VecAssemblyBegin(grid%connections, ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(grid%connections, ierr); CHKERRQ(ierr)

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'connections.out', viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%connections, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

      ! sort the vertex ids
      allocate (int_array(vertex_count))
      int_array = int_array_pointer
      allocate (int_array2(vertex_count))
      int_array2 = (/(i, i=1, vertex_count)/)
      deallocate (int_array_pointer)
      nullify (int_array_pointer)
      int_array2 = int_array2 - 1
      call PetscSortIntWithPermutation(vertex_count, int_array, int_array2, ierr); CHKERRQ(ierr)
      int_array2 = int_array2 + 1

      ! remove the duplicates
      allocate (int_array3(vertex_count))
      allocate (int_array4(vertex_count))
      allocate (grid%vertex_ids_map(vertex_count))

      int_array3 = 0
      int_array4 = 0
      int_array3(1) = int_array(int_array2(1))
      count = 1
      int_array4(int_array2(1)) = count

      do i = 2, vertex_count
         vertex_id = int_array(int_array2(i))

         if (vertex_id > int_array3(count)) then
            count = count + 1
            int_array3(count) = vertex_id
         endif

         int_array4(int_array2(i)) = count
      enddo

      vertex_count = count
      deallocate (int_array)
      deallocate (int_array2)

      ! update the local vertices set
      grid%num_pts_need = vertex_count

      ! record the map from global to local
      grid%vertex_ids_map = int_array4

      deallocate (int_array4)
      allocate (grid%vertex_ids_need(grid%num_pts_need))
      grid%vertex_ids_need = int_array3(1:grid%num_pts_need)

      ! get needed coordinates from global parallel coordinate vector
      ! IS for gather operatioin -- need local numbering
      allocate (int_array(vertex_count))
      int_array = (/(i, i=0, vertex_count - 1)/)
      call ISCreateBlock(PETSC_COMM_SELF, 3, vertex_count, int_array, PETSC_COPY_VALUES, is_gather, ierr); CHKERRQ(ierr)
      deallocate (int_array)
      ! IS for scatter operation -- need global numbering
      allocate (int_array(vertex_count))
      int_array = int_array3(1:vertex_count) - 1
      deallocate (int_array3)

      call ISCreateBlock(PETSC_COMM_SELF, 3, vertex_count, int_array, PETSC_COPY_VALUES, is_scatter, ierr); CHKERRQ(ierr)
      deallocate (int_array)

      ! allocate the space for local vertices needed on current processor
      call VecCreate(PETSC_COMM_SELF, grid%coordinates_local, ierr); CHKERRQ(ierr)
      call VecSetType(grid%coordinates_local, VECSEQ, ierr); CHKERRQ(ierr)
      call VecSetSizes(grid%coordinates_local, vertex_count*3, PETSC_DECIDE, ierr); CHKERRQ(ierr)
      call VecSetBlockSize(grid%coordinates_local, 3, ierr); CHKERRQ(ierr)

      ! scatter and gather
      call VecScatterCreate(grid%coordinates_vec, is_scatter, grid%coordinates_local, is_gather, vec_scatter, ierr); CHKERRQ(ierr)
      call ISDestroy(is_scatter, ierr); CHKERRQ(ierr)
      call ISDestroy(is_gather, ierr); CHKERRQ(ierr)

      call VecScatterBegin(vec_scatter, grid%coordinates_vec, grid%coordinates_local, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
      call VecScatterEnd(vec_scatter, grid%coordinates_vec, grid%coordinates_local, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

      call VecScatterDestroy(vec_scatter, ierr); CHKERRQ(ierr)

      ! TOUGH2 ---------------------------------------
      !if (grid%is_tough .EQV. PETSC_FALSE) call VecDestroy(grid%coordinates_vec,ierr);CHKERRQ(ierr)
      ! TOUGH2 END -----------------------------------

#if DEBUG
      write (filename, '("vertices_local_",I1,".out")') rank
      call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%coordinates_local, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      print *, 'Rank: ', rank, '       ', 'After scatter/gather, local number of points:', grid%num_pts_need
#endif

   end subroutine ScatterGatherCoordinates

! ************************************************************************** !

   subroutine AllocateConnectMatrix(grid, atts)
      !
      ! Allocates the adjacent matrix
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/22/2015
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat
      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscReal, pointer :: vec_ptr(:)
      PetscErrorCode :: ierr
      PetscInt :: i
      PetscInt :: d_nnz1(grid%num_pts_local)
      PetscInt :: d_nnz(grid%num_pts_local), o_nnz(grid%num_pts_local)
      PetscInt :: d_nnz2(grid%num_pts_local), o_nnz2(grid%num_pts_local)

      d_nnz1 = (/(i, i=grid%num_pts_local, 1, -1)/)

      call VecGetArrayF90(grid%connections, vec_ptr, ierr); CHKERRQ(ierr)
      d_nnz = min(d_nnz1, int(vec_ptr))
      o_nnz = min(int(vec_ptr), grid%num_pts_global - grid%num_pts_local)
      d_nnz2 = min(int(vec_ptr), grid%num_pts_local)
      o_nnz2 = min(int(vec_ptr), grid%num_pts_global - grid%num_pts_local)
      call VecRestoreArrayF90(grid%connections, vec_ptr, ierr); CHKERRQ(ierr)

      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz, &
                        PETSC_DECIDE, o_nnz, grid%adjmatrix, ierr); CHKERRQ(ierr)
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz2, &
                        PETSC_DECIDE, o_nnz2, grid%edgematrix, ierr); CHKERRQ(ierr)
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz2, &
                        PETSC_DECIDE, o_nnz2, grid%adjmatrix_full, ierr); CHKERRQ(ierr)
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz2, &
                        PETSC_DECIDE, o_nnz2, grid%connectivity, ierr); CHKERRQ(ierr)
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz2, &
                        PETSC_DECIDE, o_nnz2, grid%connect_area, ierr); CHKERRQ(ierr)

      ! TOUGH2
      !if ((grid%is_tough .EQV. PETSC_TRUE) .OR. (atts%are_on)) then
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz, &
                        PETSC_DECIDE, o_nnz, grid%adjmatrix_area, ierr); CHKERRQ(ierr)
      call MatCreateAIJ(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_local, &
                        PETSC_DETERMINE, PETSC_DETERMINE, &
                        PETSC_DECIDE, d_nnz, &
                        PETSC_DECIDE, o_nnz, grid%adjmatrix_len, ierr); CHKERRQ(ierr)

      ! TODO: FIX!!!
      call MatSetOption(grid%adjmatrix_area, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      call MatSetOption(grid%adjmatrix_len, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      !endif
      ! END TOUGH2

      ! TODO: The below three lines are NOT a permanent fix -
      ! MEMORY MISMANAGEMENT MAY OCCUR IF LEFT IN PLACE
      ! TEMPORARY fix to get around Petsc errors for different debugging
      call MatSetOption(grid%adjmatrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      call MatSetOption(grid%adjmatrix_full, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      call MatSetOption(grid%edgematrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      call MatSetOption(grid%connect_area, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
      call MatSetOption(grid%connectivity, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)

   end subroutine AllocateConnectMatrix

! ************************************************************************** !
   subroutine CalculateGeometricCoeff(grid, rank, atts)
      !
      ! Calculates the geometric coefficients from grid and put the values in the adjacent matrix
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/22/2015
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat

      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscInt  :: i, temp(3), rank
      PetscReal :: cell_area(3), cell_vol(3), cell_len(3)
      PetscErrorCode :: ierr

      temp = 0

      do i = 1, grid%num_elems_local

         temp = grid%elem_connectivity(i, :)

         ! calculate the areas,volumes and edge lengths associated with each element stored on processors
         call PhaseAreaVolumeLength(grid, rank, cell_area, cell_vol, cell_len, temp(1), temp(2), temp(3), i, atts)

         ! dump the results into the adjacent matrix grid%adjmatrix
         ! dump areas/lengths
         call MatSetValue(grid%adjmatrix, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_area(1)/cell_len(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_area(2)/cell_len(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_area(3)/cell_len(3), ADD_VALUES, ierr); CHKERRQ(ierr)

         !if ((grid%is_tough .EQV. PETSC_TRUE) .OR. (atts%are_on)) then
         call MatSetValue(grid%adjmatrix_area, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_area(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_area(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_area(3), ADD_VALUES, ierr); CHKERRQ(ierr)

         call MatSetValue(grid%adjmatrix_area, temp(1) - 1, temp(1) - 1, cell_vol(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, temp(2) - 1, temp(2) - 1, cell_vol(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, temp(3) - 1, temp(3) - 1, cell_vol(3), ADD_VALUES, ierr); CHKERRQ(ierr)

         ! dump lengths of connections, needed to be kept for final dumpping to TOUGH2 file
         call MatSetValue(grid%adjmatrix_len, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_len(1)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_len(2)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_len(3)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         !endif

         ! dump volumes
         call MatSetValue(grid%adjmatrix, temp(1) - 1, temp(1) - 1, cell_vol(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, temp(2) - 1, temp(2) - 1, cell_vol(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, temp(3) - 1, temp(3) - 1, cell_vol(3), ADD_VALUES, ierr); CHKERRQ(ierr)
      enddo

      call MatAssemblyBegin(grid%adjmatrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      ! TOUGH2
      !if ((grid%is_tough .EQV. PETSC_TRUE) .OR. (atts%are_on)) then
      call MatAssemblyBegin(grid%adjmatrix_len, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix_len, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%adjmatrix_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      !endif
      ! TOUGH2 END

      call VecDestroy(grid%coordinates_local, ierr); CHKERRQ(ierr)

      ! elem%connectivity needs to stay alive for PFLOTRAN dump
      if (grid%outtype /= 3) deallocate (grid%elem_connectivity)
      deallocate (grid%vertex_ids_need)
      deallocate (grid%vertex_ids_map)

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'adjmatrix.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%adjmatrix, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

   end subroutine CalculateGeometricCoeff

! ************************************************************************** !
   subroutine CalculateGeometricCoeff3D(grid, rank, atts)
      !
      ! Calculates the geometric coefficients from grid and put the values in the adjacent matrix
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 07/17/2017
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat

      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscInt  :: i, rank
      PetscInt, dimension(:), allocatable :: temp
      PetscReal, dimension(:), allocatable :: v1, v2, v3
      PetscErrorCode :: ierr

      PetscReal :: cell_area(6), cell_vol(4), cell_len(6)

      allocate (v1(3), v2(3), v3(3))
      allocate (temp(grid%ndim + 1))

      temp = 0

      do i = 1, grid%num_elems_local

         ! temp pulls a single element (tet/tri) from the array
         temp = grid%elem_connectivity(i, :)

         ! calculate the areas,volumes and edge lengths associated with each element stored on processors
         call PhaseAreaVolumeLength3D(grid, rank, cell_area, cell_vol, cell_len, temp, i, atts)

         ! write out area / length coefficients
         call MatSetValue(grid%adjmatrix, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_area(1)/cell_len(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_area(2)/cell_len(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_area(3)/cell_len(3), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(1), temp(4)) - 1, &
                          max0(temp(1), temp(4)) - 1, cell_area(4)/cell_len(4), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(2), temp(4)) - 1, &
                          max0(temp(2), temp(4)) - 1, cell_area(5)/cell_len(5), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, min0(temp(3), temp(4)) - 1, &
                          max0(temp(3), temp(4)) - 1, cell_area(6)/cell_len(6), ADD_VALUES, ierr); CHKERRQ(ierr)

         ! TOUGH2 ------------------------------------------------------------------
         !if ((grid%is_tough .EQV. PETSC_TRUE) .OR. (atts%are_on .EQV. PETSC_TRUE)) then
         call MatSetValue(grid%adjmatrix_area, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_area(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_area(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_area(3), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(1), temp(4)) - 1, &
                          max0(temp(1), temp(4)) - 1, cell_area(4), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(2), temp(4)) - 1, &
                          max0(temp(2), temp(4)) - 1, cell_area(5), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, min0(temp(3), temp(4)) - 1, &
                          max0(temp(3), temp(4)) - 1, cell_area(6), ADD_VALUES, ierr); CHKERRQ(ierr)

         call MatSetValue(grid%adjmatrix_area, temp(1) - 1, temp(1) - 1, cell_vol(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, temp(2) - 1, temp(2) - 1, cell_vol(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, temp(3) - 1, temp(3) - 1, cell_vol(3), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_area, temp(4) - 1, temp(4) - 1, cell_vol(4), ADD_VALUES, ierr); CHKERRQ(ierr)

         ! dump lengths of connections, needed to be kept for final dumpping to TOUGH2 file
         call MatSetValue(grid%adjmatrix_len, min0(temp(1), temp(2)) - 1, &
                          max0(temp(1), temp(2)) - 1, cell_len(1)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(2), temp(3)) - 1, &
                          max0(temp(2), temp(3)) - 1, cell_len(2)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(1), temp(3)) - 1, &
                          max0(temp(1), temp(3)) - 1, cell_len(3)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(1), temp(4)) - 1, &
                          max0(temp(1), temp(4)) - 1, cell_len(4)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(2), temp(4)) - 1, &
                          max0(temp(2), temp(4)) - 1, cell_len(5)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix_len, min0(temp(3), temp(4)) - 1, &
                          max0(temp(3), temp(4)) - 1, cell_len(6)/2.d0, INSERT_VALUES, ierr); CHKERRQ(ierr)
         !endif

         ! dump volumes ()
         call MatSetValue(grid%adjmatrix, temp(1) - 1, temp(1) - 1, cell_vol(1), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, temp(2) - 1, temp(2) - 1, cell_vol(2), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, temp(3) - 1, temp(3) - 1, cell_vol(3), ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValue(grid%adjmatrix, temp(4) - 1, temp(4) - 1, cell_vol(4), ADD_VALUES, ierr); CHKERRQ(ierr)
      enddo

      call MatAssemblyBegin(grid%adjmatrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      ! TOUGH2
      !if ((grid%is_tough .EQV. PETSC_TRUE) .OR. (atts%are_on)) then
      call MatAssemblyBegin(grid%adjmatrix_len, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix_len, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%adjmatrix_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      !endif
      ! TOUGH2 END

      call VecDestroy(grid%coordinates_local, ierr); CHKERRQ(ierr)

! elem%connectivity needs to stay alive for PFLOTRAN dump
      if (grid%outtype /= 3) deallocate (grid%elem_connectivity)
      deallocate (grid%vertex_ids_need)
      deallocate (grid%vertex_ids_map)

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'adjmatrix.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%adjmatrix, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

   end subroutine CalculateGeometricCoeff3D

!**************************************************************************

   subroutine PhaseAreaVolumeLength3D(grid, rank, cell_area, cell_vol, cell_len, idx, ielem, atts)
      !
      ! Calculates the volumes and areas for readin elements and nodes
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 12/7/2017
      !

#include "petsc/finclude/petscvec.h"

      use petscvec
      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscInt :: rank
      PetscInt :: idx(:)
      PetscInt :: idx1, idx2, idx3, idx4 ! local index of element vertices
      PetscInt :: i, ielem
      PetscScalar, dimension(3) :: v1, v2, v3, v4, nn
      PetscScalar, dimension(3) :: e12, e13, e14, e23, e24, e34
      PetscScalar, dimension(3) :: C123, C124, C134, C234, C1234
      PetscReal :: a(6)
      PetscReal :: cell_area(6), cell_vol(4), cell_len(6)
      PetscReal, pointer :: vec_ptr(:)
      PetscErrorCode :: ierr

      idx1 = idx(1)
      idx2 = idx(2)
      idx3 = idx(3)
      idx4 = idx(4)

!  Fetch the coordinates
      call VecGetArrayF90(grid%coordinates_local, vec_ptr, ierr); CHKERRQ(ierr)
      v1 = vec_ptr((grid%vertex_ids_map((ielem - 1)*4 + 1) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*4 + 1) - 1)*3 + 3)
      v2 = vec_ptr((grid%vertex_ids_map((ielem - 1)*4 + 2) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*4 + 2) - 1)*3 + 3)
      v3 = vec_ptr((grid%vertex_ids_map((ielem - 1)*4 + 3) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*4 + 3) - 1)*3 + 3)
      v4 = vec_ptr((grid%vertex_ids_map((ielem - 1)*4 + 4) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*4 + 4) - 1)*3 + 3)
      call VecRestoreArrayF90(grid%coordinates_local, vec_ptr, ierr); CHKERRQ(ierr)

!  calculate the positive normal direction of the element if no longer on the same fracture
      call Cross(nn, v2 - v1, v3 - v1) ! possibly depreciated
      nn = nn/Norm(nn) ! possibly depreciated

      ! Get the bisectors for each of the cell edges
      e12 = 0.5d0*(v1 + v2)
      e13 = 0.5d0*(v1 + v3)
      e14 = 0.5d0*(v1 + v4)
      e23 = 0.5d0*(v2 + v3)
      e24 = 0.5d0*(v2 + v4)
      e34 = 0.5d0*(v3 + v4)

      ! Find the circumcenter for each of the external faces
      call CalcCircumcenter(grid, C123, v1, v2, v3)
      call CalcCircumcenter(grid, C124, v1, v2, v4)
      call CalcCircumcenter(grid, C134, v1, v3, v4)
      call CalcCircumcenter(grid, C234, v2, v3, v4)

      ! Get the centroid of the tetrahedron
      call TetrahedronCentroid(grid, C1234, v1, v2, v3, v4)

      ! -- Calculate Voronoi face areas --------------------------------------
      ! After dividing a single tetrahedron into its 'sub-tets', there are
      ! 12 faces generated, internal and facing each other relative to the
      ! larger tet. Six of these faces are degenerate. This leaves 6 faces,
      ! which we calculate the area of, as a 2D triangle in 3D space.

      cell_area(1) = TriangleFaceArea(v2, v1, e12, C124, C123, C1234) ! e12 C123 C124 CC
      cell_area(2) = TriangleFaceArea(v3, v2, e23, C234, C123, C1234) ! e23 C234 C123 CC
      cell_area(3) = TriangleFaceArea(v3, v1, e13, C123, C134, C1234) ! e13 C134 C123 CC
      cell_area(4) = TriangleFaceArea(v4, v1, e14, C134, C124, C1234) ! e14 C134 C124 CC
      cell_area(5) = TriangleFaceArea(v4, v2, e24, C124, C234, C1234) ! e24 C234 C124 CC
      cell_area(6) = TriangleFaceArea(v4, v3, e34, C234, C134, C1234) ! e34 C134 C234 CC

      ! -- Calculate the Voronoi volume ---------------------------------------
      ! At each vertex of a tetrahedon, there are six sub-tets that are formed.
      ! This calculates the volume of each sub-tet, then adds the volume of each
      ! to find the volume a single node contributes.

      a(1) = Volume(v1, e12, C124, C1234)
      a(2) = Volume(v1, C124, e14, C1234)
      a(3) = Volume(v1, C134, e13, C1234)
      a(4) = Volume(v1, e14, C134, C1234)
      a(5) = Volume(v1, C123, e12, C1234)
      a(6) = Volume(v1, e13, C123, C1234)

      cell_vol(1) = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)

      a(1) = Volume(v2, C124, e12, C1234)
      a(2) = Volume(v2, C234, e24, C1234)
      a(3) = Volume(v2, e23, C234, C1234)
      a(4) = Volume(v2, e12, C123, C1234)
      a(5) = Volume(v2, C123, e23, C1234)
      a(6) = Volume(v2, e24, C124, C1234)

      cell_vol(2) = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)

      a(1) = Volume(v3, C134, e34, C1234)
      a(2) = Volume(v3, e34, C234, C1234)
      a(3) = Volume(v3, C234, e23, C1234)
      a(4) = Volume(v3, e13, C134, C1234)
      a(5) = Volume(v3, e23, C123, C1234)
      a(6) = Volume(v3, C123, e13, C1234)

      cell_vol(3) = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)

      a(1) = Volume(v4, C134, e14, C1234)
      a(2) = Volume(v4, e14, C124, C1234)
      a(3) = Volume(v4, C124, e24, C1234)
      a(4) = Volume(v4, e24, C234, C1234)
      a(5) = Volume(v4, C234, e34, C1234)
      a(6) = Volume(v4, e34, C134, C1234)

      cell_vol(4) = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)

      ! -- Calculate the edge lengths -----------------------------
      ! Find the edge length between vertex v_i and vertex v_j

      cell_len(1) = Norm((/v1(1) - v2(1), v1(2) - v2(2), v1(3) - v2(3)/)) ! | v1 - v2 |
      cell_len(2) = Norm((/v2(1) - v3(1), v2(2) - v3(2), v2(3) - v3(3)/)) ! | v2 - v3 |
      cell_len(3) = Norm((/v1(1) - v3(1), v1(2) - v3(2), v1(3) - v3(3)/)) ! | v1 - v3 |
      cell_len(4) = Norm((/v1(1) - v4(1), v1(2) - v4(2), v1(3) - v4(3)/)) ! | v1 - v4 |
      cell_len(5) = Norm((/v2(1) - v4(1), v2(2) - v4(2), v2(3) - v4(3)/)) ! | v2 - v4 |
      cell_len(6) = Norm((/v3(1) - v4(1), v3(2) - v4(2), v3(3) - v4(3)/)) ! | v3 - v4 |

      ! Store tetrahedron center
      if (grid%outtype == 3) grid%cell_cc(ielem, :) = C1234
      if (grid%outtype == 3) grid%cell_vol(ielem) = Volume(v1, v2, v3, v4)

      ! -- Diagnostics collection -----------------------------
      ! If diagnostics are turned on, start collecting attributes

      if (atts%are_on .EQV. PETSC_TRUE) then
         atts%face_area = atts%face_area - ( &
                          Area(v1, v2, v3, (v4 - v1)/Norm(v4 - v1)) + &
                          Area(v2, v3, v4, (v4 - v1)/Norm(v4 - v1)) + &
                          Area(v1, v3, v4, (v4 - v1)/Norm(v4 - v1)) + &
                          Area(v1, v2, v4, (v4 - v1)/Norm(v4 - v1)))

         atts%vol = atts%vol - Volume(v1, v2, v3, v4)

         ! Count negative areas and volumes for Diagnostics() allocation
         ! TODO: this is NOT true
         do i = 1, 6
            if (cell_area(i)/cell_len(i) > 0) then
               atts%neg_coeff = atts%neg_coeff + 1
            endif
            if (i <= 4) then
               if (cell_vol(i) < 0) then
                  atts%neg_vol = atts%neg_vol + 1
               endif
            endif
         enddo
      endif

   end subroutine PhaseAreaVolumeLength3D

!**************************************************************************

   PetscReal function TriangleFaceArea(v_j, v_i, e_ij, C_1, C_2, C1234)

#include "petsc/finclude/petscvec.h"
      use petscvec
      implicit none

      PetscReal, dimension(3) :: v_i, v_j, e_ij, C_1, C_2, C1234
      PetscReal, dimension(3) :: a_1, a_2, nn

      nn = (v_j - v_i)/Norm((v_j - v_i))

      call Cross(a_1, C_1 - e_ij, C1234 - e_ij)
      call Cross(a_2, C1234 - e_ij, C_2 - e_ij)

      TriangleFaceArea = dot_product(0.5d0*(a_1 + a_2), nn)

   end function TriangleFaceArea

!**************************************************************************

   subroutine PhaseAreaVolumeLength(grid, rank, cell_area, cell_vol, cell_len, idx1, idx2, idx3, ielem, atts)
      !
      ! Calculates the volumes and areas for each readin elements/triangle
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/22/2015
      !

#include "petsc/finclude/petscvec.h"

      use petscvec
      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscInt :: rank
      PetscInt :: idx1, idx2, idx3 ! local index of element vertices
      PetscInt :: ielem
      PetscScalar, dimension(3) :: v1, v2, v3, nn, lcc, e12, e13, e23
      PetscReal :: a(6)
      PetscReal :: cell_area(3), cell_vol(3), cell_len(3), tmp(3)
      PetscReal, pointer :: vec_ptr(:)
      PetscErrorCode :: ierr

!  Fetch the coordinates
      call VecGetArrayF90(grid%coordinates_local, vec_ptr, ierr); CHKERRQ(ierr)
      v1 = vec_ptr((grid%vertex_ids_map((ielem - 1)*3 + 1) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*3 + 1) - 1)*3 + 3)
      v2 = vec_ptr((grid%vertex_ids_map((ielem - 1)*3 + 2) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*3 + 2) - 1)*3 + 3)
      v3 = vec_ptr((grid%vertex_ids_map((ielem - 1)*3 + 3) - 1)*3 + 1:(grid%vertex_ids_map((ielem - 1)*3 + 3) - 1)*3 + 3)
      call VecRestoreArrayF90(grid%coordinates_local, vec_ptr, ierr); CHKERRQ(ierr)

!  calculate the positive normal direction of the element if no longer on the same fracture
      call Cross(nn, v2 - v1, v3 - v1)
      nn = nn/Norm(nn)

      call CalcCircumcenter(grid, lcc, v1, v2, v3)
      !if (grid%outtype == 3) grid%cell_cc(ielem,:) = lcc

      e12 = 0.5d0*(v1 + v2) ! midpoint of cell edges
      e23 = 0.5d0*(v2 + v3)
      e13 = 0.5d0*(v1 + v3)

      ! calculate the Voronoi phase area--vec
      cell_area(1) = VoronoiPhaseArea(lcc, e12, v1, v2, nn) !idx1,idx2
      cell_area(2) = VoronoiPhaseArea(lcc, e23, v2, v3, nn) !idx2,idx3
      cell_area(3) = VoronoiPhaseArea(lcc, e13, v3, v1, nn) !idx1,idx3

!----------------------------------------------------------
! calculate the Voronoi volume pieces--area
!  the order of the vertex area(v1,v2,v3) matters! define the counterclockwise (righthand-rule) as
!  the positive direction, i.e. area can be negative when the order is in clockwise.
!  Permutation is allowed, but keep the relative position

      a(2) = Area(e12, lcc, v1, nn)
      a(3) = a(2) !Area(e12,v2,lcc,nn)
      a(4) = Area(e23, lcc, v2, nn)
      a(5) = a(4) !Area(e23,v3,lcc,nn)
      a(6) = Area(e13, lcc, v3, nn)
      a(1) = a(6) !Area(e13,v1,lcc,nn)

      cell_vol(1) = a(1) + a(2) !idx1,idx1
      cell_vol(2) = a(3) + a(4) !idx2,idx2
      cell_vol(3) = a(5) + a(6) !idx3,idx3

      cell_len(1) = Norm((/v1(1) - v2(1), v1(2) - v2(2), v1(3) - v2(3)/))
      cell_len(2) = Norm((/v2(1) - v3(1), v2(2) - v3(2), v2(3) - v3(3)/))
      cell_len(3) = Norm((/v1(1) - v3(1), v1(2) - v3(2), v1(3) - v3(3)/))

      ! -- Diagnostics collection -----------------------------
      ! If diagnostics are turned on, start collecting attributes

      if (atts%are_on .EQV. PETSC_TRUE) then
         atts%face_area = atts%face_area + Area(v1, v2, v3, nn)
         atts%vol = atts%vol + Area(v1, v2, v3, nn)
      endif

      ! Calculate the unsigned area of a face
      if (grid%outtype == 3) then
         call Cross(tmp, v2 - v1, v3 - v1)
         grid%cell_vol(ielem) = 0.5d0*Norm(tmp)
      endif

   end subroutine PhaseAreaVolumeLength

!**************************************************************************
   subroutine CalcCircumcenter(grid, circ, v1, v2, v3)
      !
      ! Calculates the circumcenter for the Delauney element
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/08/2015
      !

#include "petsc/finclude/petscvec.h"

      use petscvec
      implicit none

      type(grid_type) :: grid
      PetscReal, dimension(3) :: v1, v2, v3, circ, mm
      PetscReal :: coef(3), denom

      if (grid%cv_type == 1) then
         ! Voronoi cell center
         call Cross(mm, v1 - v2, v2 - v3)
         denom = 2.d0*Norm(mm)**2

         coef(1) = Norm(v2 - v3)**2*dot_product(v1 - v2, v1 - v3)/denom
         coef(2) = Norm(v1 - v3)**2*dot_product(v2 - v1, v2 - v3)/denom
         coef(3) = Norm(v1 - v2)**2*dot_product(v3 - v1, v3 - v2)/denom

         circ = coef(1)*v1 + coef(2)*v2 + coef(3)*v3
      else
         ! Median cell center
         circ = (v1 + v2 + v3)/3.0
      endif

   end subroutine CalcCircumcenter

!**************************************************************************
   subroutine TetrahedronCentroid(grid, cent, v1, v2, v3, v4)
      !
      ! Calculates the circumsphere of a tetrahedron
      !
      ! The circumsphere for a tetrahedron is given by a set
      ! of determinates of matrices:
      !
      !     x0 = det(Dx) / 2*det(a), y0 = det(Dy) / 2*det(a), z0 = det(Dz) / 2*det(a)
      !
      ! where Di and a are 4x4 matrices constructed from vectors v1,..,v4.
      ! The set of equations can be found at Wolfram MathWorld,
      ! under 'Circumsphere'.
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 09/05/2017
      !

#include "petsc/finclude/petscvec.h"

      use petscvec
      implicit none

      type(grid_type) :: grid
      PetscReal, dimension(4, 4) :: MatA, MatDx, MatDy, MatDz
      PetscReal, dimension(3) :: v1, v2, v3, v4, cent
      PetscReal, dimension(4) :: vsum
      PetscReal :: a, Dx, Dy, Dz

      ! If (voronoi)...
      if (grid%cv_type == 1) then
         ! Populate the 'a' matrix
         MatA(1:3, 1) = v1
         MatA(1:3, 2) = v2
         MatA(1:3, 3) = v3
         MatA(1:3, 4) = v4
         MatA(4, :) = 1

         ! square and sum x,y,z in each vector
         vsum(1) = v1(1)**2 + v1(2)**2 + v1(3)**2
         vsum(2) = v2(1)**2 + v2(2)**2 + v2(3)**2
         vsum(3) = v3(1)**2 + v3(2)**2 + v3(3)**2
         vsum(4) = v4(1)**2 + v4(2)**2 + v4(3)**2

         ! Fill the Dx matrix
         MatDx(1, :) = vsum
         MatDx(2, :) = MatA(2, :)
         MatDx(3, :) = MatA(3, :)
         MatDx(4, :) = 1

         ! Fill the Dy matrix
         MatDy(1, :) = vsum
         MatDy(2, :) = MatA(1, :)
         MatDy(3, :) = MatA(3, :)
         MatDy(4, :) = 1

         ! Fill the Dz matrix
         MatDz(1, :) = vsum
         MatDz(2, :) = MatA(1, :)
         MatDz(3, :) = MatA(2, :)
         MatDz(4, :) = 1

         ! Calculate the determinates
         a = Det4x4(MatA(:, 1), MatA(:, 2), MatA(:, 3), MatA(:, 4))
         Dx = Det4x4(MatDx(:, 1), MatDx(:, 2), MatDx(:, 3), MatDx(:, 4))
         Dy = Det4x4(MatDy(:, 1), MatDy(:, 2), MatDy(:, 3), MatDy(:, 4))*(-1)
         Dz = Det4x4(MatDz(:, 1), MatDz(:, 2), MatDz(:, 3), MatDz(:, 4))

         ! Find the centroid
         cent(1) = Dx/(2*a)
         cent(2) = Dy/(2*a)
         cent(3) = Dz/(2*a)

         ! else (median)...
      else
         cent = (v1 + v2 + v3 + v4)/4.0
      endif

   end subroutine TetrahedronCentroid

!**************************************************************************
   PetscReal function Det4x4(v1, v2, v3, v4)
      ! Calculates the determinant of a 4x4 matrix
      !
      ! For sequential rows in the matrix, v1,v2,v3,v4,
      ! return the determinate
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 10/05/2017
      !

      implicit none

      PetscReal, dimension(4) :: v1, v2, v3, v4
      PetscReal :: c1, c2, c3, c4

      c1 = v1(1)*Det3x3(v2(2:4), v3(2:4), v4(2:4))
      c2 = v1(2)*Det3x3(v2((/1, 3, 4/)), v3((/1, 3, 4/)), v4((/1, 3, 4/)))
      c3 = v1(3)*Det3x3(v2((/1, 2, 4/)), v3((/1, 2, 4/)), v4((/1, 2, 4/)))
      c4 = v1(4)*Det3x3(v2(1:3), v3(1:3), v4(1:3))

      Det4x4 = c1 - c2 + c3 - c4

   end function Det4x4

!**************************************************************************
   PetscReal function Det3x3(v1, v2, v3)
      ! Calculates the determinant of a 3x3 matrix
      !
      ! For sequential rows in the matrix, v1,v2,v3,
      ! return the determinate
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 10/05/2017
      !

      implicit none

      PetscReal, dimension(3) :: v1, v2, v3
      PetscReal :: c1, c2, c3

      c1 = v1(1)*(v2(2)*v3(3) - v3(2)*v2(3))
      c2 = v1(2)*(v2(1)*v3(3) - v3(1)*v2(3))
      c3 = v1(3)*(v2(1)*v3(2) - v3(1)*v2(2))

      Det3x3 = c1 - c2 + c3

   end function Det3x3

!**************************************************************************
   subroutine Cross(v12, v1, v2)
      !
      ! Calculates the cross product of v1 and v2
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/08/2015
      !

      implicit none

      PetscReal, dimension(3) :: v1, v2, v12

      v12(1) = v1(2)*v2(3) - v2(2)*v1(3)
      v12(2) = v2(1)*v1(3) - v1(1)*v2(3)
      v12(3) = v1(1)*v2(2) - v2(1)*v1(2)

   end subroutine Cross

!**************************************************************************
   PetscReal function VoronoiPhaseArea(center, mid, v1, v2, nn)
      !
      ! Calculates the Voronoi Phase Area
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/08/2015
      !

      implicit none

      PetscReal, dimension(3) :: center, mid, v1, v2, temp, enorm, nn
      enorm = (v2 - v1)/Norm(v2 - v1)
      call Cross(temp, nn, enorm)
      VoronoiPhaseArea = dot_product(temp, center - mid)

   end function VoronoiPhaseArea

!**************************************************************************
   PetscReal function Area(v1, v2, v3, nn)
      !
      ! Calculates the Area
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/08/2015
      !

      implicit none

      PetscReal, dimension(3) :: v1, v2, v3, temp, nn

      call Cross(temp, v2 - v1, v3 - v1)

      Area = 0.5d0*dot_product(temp, nn)

   end function Area

!**************************************************************************

   PetscReal function Volume(v1, v2, v3, v4)
      !
      ! Calculates the volume of a tetrahedon
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 09/05/2017
      !

      implicit none

      PetscReal, dimension(3) :: v1, v2, v3, v4
      PetscReal, dimension(3) :: v21, v31, v41, temp

      v21 = v2 - v1
      v31 = v3 - v1
      v41 = v4 - v1

      call Cross(temp, v31, v41)
      Volume = dot_product(v21, temp)/6.d0

   end function Volume
!**************************************************************************

   PetscReal function Norm(v)
      !
      ! Calculates the 2-norm of a vector
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/08/2015
      !

      implicit none

      PetscReal :: v(3)
      Norm = dsqrt(sum(v**2))

   end function Norm

!**************************************************************************
   subroutine CreateEdgeMatrix(grid, rank)
      !
      ! Creates the exact degree vector and degree_tot vector
      ! and edge matrix from adjacentmatrix
      ! to prepare for the output
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/22/2015
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat
      implicit none

      type(grid_type) :: grid
      PetscInt :: i, j, ncols
      PetscInt :: idx, rank, istart, iend, istart_mat1, istart_mat2, iend_mat
      PetscScalar :: zero = 0, one = 1
      PetscInt :: zero_t = 0, one_t = 1
      PetscInt, allocatable :: cols(:)
      PetscScalar, allocatable :: vals(:), temp_cols(:), temp_vals(:)
      PetscErrorCode :: ierr
      PetscScalar :: temp
      PetscReal :: max_degree

      call VecCreateMPI(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_global, grid%degree, ierr); CHKERRQ(ierr)
      call VecSet(grid%degree, zero, ierr); CHKERRQ(ierr)
      call VecMax(grid%connections, PETSC_NULL_INTEGER, max_degree, ierr); CHKERRQ(ierr)
      call VecDestroy(grid%connections, ierr); CHKERRQ(ierr)

      ! Create the 'mapping' vector: each element at index i represents the count of connections for element i
      call VecCreateMPI(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_global, grid%connect_map, ierr); CHKERRQ(ierr)
      call VecSet(grid%connect_map, zero, ierr); CHKERRQ(ierr)
      call VecGetOwnershipRange(grid%connect_map, istart, iend, ierr); CHKERRQ(ierr)
      call MatGetOwnershipRange(grid%connectivity, istart_mat1, iend_mat, ierr); CHKERRQ(ierr)
      call MatGetOwnershipRange(grid%connectivity, istart_mat2, iend_mat, ierr); CHKERRQ(ierr)

      if (grid%ndim == 3) then
         allocate (cols(int(max_degree) + 1))
         allocate (vals(int(max_degree) + 1))
      else
         allocate (cols(int(max_degree)))
         allocate (vals(int(max_degree)))
      end if

      idx = 0

! form the upper triangular sparse matrix
      do i = 1, grid%num_pts_local
         cols = 0

         ! TOUGH2
         if (grid%is_tough .EQV. PETSC_TRUE) then
            call MatGetRow(               &
               grid%adjmatrix_area,       &
               grid%vertex_ids(i) - 1,    &
               ncols,                     &
               cols,                      &
               vals,                      &
               ierr                       &
            ); CHKERRQ(ierr)
         else
            call MatGetRow(               &
               grid%adjmatrix,            &
               grid%vertex_ids(i) - 1,    &
               ncols,                     &
               cols,                      &
               vals,                      &
               ierr                       &
            ); CHKERRQ(ierr)
         endif
         ! TOUGH2 END

         temp = ncols

         allocate (temp_cols(ncols))
         allocate (temp_vals(ncols))

         temp_cols = cols(1:ncols) + 1
         temp_vals = vals(1:ncols)

         ! Store the connectivity graph (i->(j_1,...,j_N))
         ! Each row i corresponds to Voronoi cell i, and each column j is the connection between cell i and cell j
         !call MatSetValues(grid%connectivity,one,(/i-1/),ncols,(/(j,j=1,ncols)/),temp_cols,INSERT_VALUES,ierr);CHKERRQ(ierr)
         call MatSetValues(grid%connectivity, one_t, (/i - 1/) + istart_mat1, ncols, &
                           (/(j, j=1, ncols)/), temp_cols, INSERT_VALUES, ierr); CHKERRQ(ierr)

         ! Store the area of the connectivity (i->(j_1,...,j_N))
         ! Each row i corresponds to Voronoi cell i, and each column j is the Voronoi area orthogonal to connection i->j
         !call MatSetValues(grid%connect_area,one,(/i-1/),ncols,(/(j,j=1,ncols)/),temp_vals,INSERT_VALUES,ierr);CHKERRQ(ierr)
         call MatSetValues(grid%connect_area, one_t, (/i - 1/) + istart_mat2, ncols, &
                           (/(j, j=1, ncols)/), temp_vals, INSERT_VALUES, ierr); CHKERRQ(ierr)

         ! Store the number of connections for element i, in vector element i
         call VecSetValue(grid%connect_map, (i - 1) + istart, temp, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call VecSetValue(grid%degree, grid%vertex_ids(i) - 1, temp, ADD_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValues(grid%edgematrix, one_t, grid%vertex_ids(i) - 1, ncols, cols(1:ncols), temp_cols, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValues(grid%adjmatrix_full, one_t, grid%vertex_ids(i) - 1, ncols, cols(1:ncols), temp_vals, INSERT_VALUES, ierr); CHKERRQ(ierr)

         temp_cols = int(grid%vertex_ids(i))

         ! make a symmetric matrix with zero diagonal
         call MatSetValues(grid%edgematrix, ncols - 1, cols(2:ncols), one_t, &
                           grid%vertex_ids(i) - 1, temp_cols(2:ncols), INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatSetValues(grid%adjmatrix_full, ncols - 1, cols(2:ncols), one_t, &
                           grid%vertex_ids(i) - 1, temp_vals(2:ncols), INSERT_VALUES, ierr); CHKERRQ(ierr)

         call MatSetValue(             &
            grid%adjmatrix_full,       &
            grid%vertex_ids(i) - 1,    &
            grid%vertex_ids(i) - 1,    &
            zero,                      &
            INSERT_VALUES,             &
            ierr                       &
         ); CHKERRA(ierr)

         ! TOUGH2
         if (grid%is_tough .EQV. PETSC_TRUE) then
            call MatRestoreRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, cols, PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
         else
            call MatRestoreRow(grid%adjmatrix, grid%vertex_ids(i) - 1, ncols, cols, PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
         endif
         ! END TOUGH2

         deallocate (temp_cols)
         deallocate (temp_vals)
      enddo

      ! =============================================================================================== !
      !   Assemble vectors and matrices
      ! =============================================================================================== !

      call VecAssemblyBegin(grid%degree, ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(grid%degree, ierr); CHKERRQ(ierr)

      call VecAssemblyBegin(grid%connect_map, ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(grid%connect_map, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%connectivity, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%connectivity, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%connect_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%connect_area, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%edgematrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%edgematrix, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      call MatAssemblyBegin(grid%adjmatrix_full, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%adjmatrix_full, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

      deallocate (vals)
      deallocate (cols)

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'connect_map.out', viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%connect_map, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'connectivity.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%connectivity, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'connect_area.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%connect_area, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

! calculate the grid%degree_tot based on grid%edgematrix

      call VecCreateMPI(PETSC_COMM_WORLD, grid%num_pts_local, grid%num_pts_global, grid%degree_tot, ierr); CHKERRQ(ierr)
      call VecSet(grid%degree_tot, zero, ierr); CHKERRQ(ierr)

      do i = 1, grid%num_pts_local
         call MatGetRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
         temp = ncols
         call VecSetValue(grid%degree_tot, grid%vertex_ids(i) - 1, temp, INSERT_VALUES, ierr); CHKERRQ(ierr)
         call MatRestoreRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, cols, PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
      enddo

      call VecAssemblyBegin(grid%degree_tot, ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(grid%degree_tot, ierr); CHKERRQ(ierr)

#if DEBUG
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'degree.out', viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%degree, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'degree_tot.out', viewer, ierr); CHKERRQ(ierr)
      call VecView(grid%degree_tot, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'edge_list.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%edgematrix, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'adjmatrix_full.out', viewer, ierr); CHKERRQ(ierr)
      call MatView(grid%adjmatrix_full, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
#endif

   end subroutine CreateEdgeMatrix

!**************************************************************************
   subroutine GridWriteFEHM(grid, atts, rank, size)
      !
      ! Dumps all the information into FEHM .stor file
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 07/26/2015
      !

      ! row_count -> 'Count for Each Row' (S4 in STOR file descriptor)
      !   Deprecates `degree_local`.

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat
      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      character(len=30) :: fdate
      PetscInt :: fileid, i, rank, irank, size
      PetscInt :: rec_tot, rec_local, NEQ, coeff_count
      PetscInt :: io_rank = 0, num_area_coef = 1
      PetscInt :: temp_int, pos1, ncon_max, ncon_max_local, mask_delta, mask_delta_local
      PetscReal :: degree_tot_max
      PetscReal, pointer :: vec_ptr(:), vec_ptr1(:)
      PetscInt, allocatable :: sendcounts_pts(:), pos_pts(:), tmp_map(:)
      PetscInt, allocatable :: sendcounts_edge(:), pos_edge(:)
      PetscReal, allocatable :: vol_local(:), geom_coeffs(:), voronoi_volumes(:)
      PetscScalar, allocatable :: edge_temp(:)
      PetscReal, allocatable :: edge_local(:)
      PetscReal, allocatable :: adj_temp(:), adj_local(:)
      PetscInt, allocatable :: degree_1(:), degree_2(:), degree_local(:), degree_local2(:)
      PetscReal, allocatable :: data_local_real(:), data_global_real(:)
      PetscInt, allocatable :: num_pts_0(:), rec_local_0(:)
      PetscInt, allocatable :: indx(:), indx_save(:), coeff_pointers(:)
      PetscInt, allocatable :: sendcounts(:), pos(:), off_diagonal_map(:), on_diagonal_map(:)
      PetscInt, allocatable :: row_count(:), row_entries(:), pointer_stride(:)
      PetscBool, allocatable :: nil_mask(:), full_mask(:)
      PetscInt :: ncols
      Vec :: mat_diag
      PetscErrorCode :: ierr

      NEQ = grid%num_pts_local

      call VecMax(grid%degree_tot, PETSC_NULL_INTEGER, degree_tot_max, ierr); CHKERRQ(ierr)

      ! get the vol_local
      allocate (vol_local(grid%num_pts_local))
      call VecCreate(PETSC_COMM_WORLD, mat_diag, ierr); CHKERRQ(ierr)
      call VecSetSizes(mat_diag, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
      call VecSetType(mat_diag, VECMPI, ierr); CHKERRQ(ierr)
      call MatGetDiagonal(grid%adjmatrix, mat_diag, ierr); CHKERRQ(ierr)
      call MatDestroy(grid%adjmatrix, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      vol_local = vec_ptr
      call VecRestoreArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      call VecDestroy(mat_diag, ierr); CHKERRQ(ierr)

      ! get the edge_local & adj_local
      call VecGetArrayReadF90(grid%degree_tot, vec_ptr, ierr); CHKERRQ(ierr)
      rec_local = int(sum(vec_ptr))
      allocate (edge_local(rec_local))
      allocate (adj_local(rec_local))
      allocate (full_mask(rec_local))

      pos1 = 0
      temp_int = int(maxval(vec_ptr))
      allocate (edge_temp(temp_int))
      allocate (adj_temp(temp_int))
      allocate (nil_mask(temp_int)) ! mask

      allocate (row_count(NEQ + 1))

      ! -------------------------------------------------------------
      ! CAPTURE CONNECTIVITY && AREA COEFFICIENTS
      !
      ! Iterate over each row in the sparse matrix
      ! Here, we will pull connectivity (i.e. where entries exist), and
      ! the corresponding off-diagonal coefficients

      ncon_max_local = 0
      mask_delta = 0
      mask_delta_local = 0

      do i = 1, grid%num_pts_local
         temp_int = int(vec_ptr(i))

         ! Capture connectivity
         call MatGetRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)
         edge_local(pos1 + 1:pos1 + temp_int) = int(edge_temp(1:temp_int))
         call MatRestoreRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)

         ! Capture coefficients
         call MatGetRow(grid%adjmatrix_full, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         adj_local(pos1 + 1:pos1 + temp_int) = adj_temp(1:temp_int)
         call MatRestoreRow(grid%adjmatrix_full, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)

         ! Store number of connections
         row_count(i + 1) = temp_int

         ! Generate a mask for 'removing' entries below a threshold.
         ! Namely, (i) all non-zeros should be preserved, EXCEPT FOR
         ! (ii) non-zeros on the diagonal (i.e. Voronoi volumes)
         if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) then
            nil_mask = PETSC_TRUE
            where (abs(adj_temp(1:temp_int)) < 1e-6) nil_mask(1:temp_int) = PETSC_FALSE
            where (int(edge_temp(1:temp_int)) == grid%vertex_ids(i)) nil_mask(1:temp_int) = PETSC_TRUE

            ! Fill full masking vector with captured subset
            full_mask(pos1 + 1:pos1 + temp_int) = nil_mask(1:temp_int)
            ncon_max_local = max(ncon_max_local, count(nil_mask(1:temp_int)))

            !row_count(i+1) = count(nil_mask(1:temp_int))
            vec_ptr(i) = count(nil_mask(1:temp_int))
            mask_delta_local = mask_delta_local + (temp_int - count(nil_mask(1:temp_int)))
         endif

         pos1 = pos1 + temp_int
      enddo

      call MPI_Reduce(rec_local, rec_tot, 1, MPI_INT, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(mask_delta_local, mask_delta, 1, MPI_INT, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(ncon_max_local, ncon_max, 1, MPI_INT, MPI_MAX, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! Below is probably not necessary as only io_rank needs access to rec_tot
      call MPI_Bcast(rec_tot, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Bcast(ncon_max, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Bcast(mask_delta, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (sendcounts_pts(size))
         allocate (sendcounts_edge(size))
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, sendcounts_pts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(rec_local, 1, MPI_INT, sendcounts_edge, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (pos_pts(size))
         allocate (pos_edge(size))
         pos_pts(1) = 0
         do irank = 1, size - 1
            pos_pts(irank + 1) = pos_pts(irank) + sendcounts_pts(irank)
         enddo
         pos_edge(1) = 0
         do irank = 1, size - 1
            pos_edge(irank + 1) = pos_edge(irank) + sendcounts_edge(irank)
         enddo
         if (atts%coefficients_are_dedudded .eqv. PETSC_FALSE) rec_tot = pos_edge(size) + sendcounts_edge(size) !!!!!!!!!!

         ! Update NCOEF based on mask
         !if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) rec_tot = 33!count(full_mask)
         !if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) rec_local = count(full_mask) !!! WARNING: IS NOT TRUE FOR PARALLEL!!!
         !if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) degree_tot_max = ncon_max
      endif

      if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) rec_local = count(full_mask)
      if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) degree_tot_max = ncon_max

      ! calculate degree_local
      allocate (degree_local(grid%num_pts_local))
      call MPI_Scatter(pos_edge, 1, MPI_INT, degree_local(1), 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! Calculate and populate the `Count for Each Row` vector
      row_count(1) = NEQ + 1
      do i = 2, grid%num_pts_local
         degree_local(i) = degree_local(i - 1) + int(vec_ptr(i - 1))
         row_count(i) = row_count(i) + row_count(i - 1)
      enddo

      degree_local = degree_local + 1 + grid%num_pts_global

      ! calculate degree_local2 = degree_diff_local+degree_local+1
      allocate (degree_1(grid%num_pts_local))
      allocate (degree_2(grid%num_pts_local))
      allocate (degree_local2(grid%num_pts_local))
      degree_1 = int(vec_ptr)
      call VecGetArrayReadF90(grid%degree, vec_ptr1, ierr); CHKERRQ(ierr)
      degree_2 = int(vec_ptr1)
      call VecRestoreArrayReadF90(grid%degree, vec_ptr1, ierr); CHKERRQ(ierr)
      degree_local2 = degree_local + degree_1 - degree_2 + 1
      deallocate (degree_1)
      deallocate (degree_2)

      call VecRestoreArrayReadF90(grid%degree_tot, vec_ptr, ierr); CHKERRQ(ierr)

#if DEBUG
      print *, rank, 'vol_local=', vol_local
      print *, rank, 'degree_local', degree_local
      print *, rank, 'edge_lcoal=', edge_local
      print *, rank, 'adj_local=', adj_local
      print *, rank, 'degree_local2=', degree_local2
#endif

      allocate (data_local_real(grid%num_pts_local*3 + rec_local*2 + 2))
      data_local_real(1:grid%num_pts_local) = vol_local
      data_local_real(grid%num_pts_local + 1:grid%num_pts_local*2) = degree_local
      data_local_real(grid%num_pts_local*2 + rec_local + 1:grid%num_pts_local*3 + rec_local) = degree_local2
      data_local_real(grid%num_pts_local*3 + rec_local*2 + 1) = grid%num_pts_local
      data_local_real(grid%num_pts_local*3 + rec_local*2 + 2) = rec_local

      ! Apply mask if compression
      if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) then
         data_local_real(grid%num_pts_local*3 + rec_local + 1:grid%num_pts_local*3 + rec_local*2) = pack(adj_local, full_mask)
         data_local_real(grid%num_pts_local*2 + 1:grid%num_pts_local*2 + rec_local) = pack(edge_local, full_mask)
      else
         data_local_real(grid%num_pts_local*3 + rec_local + 1:grid%num_pts_local*3 + rec_local*2) = adj_local
         data_local_real(grid%num_pts_local*2 + 1:grid%num_pts_local*2 + rec_local) = edge_local
      endif

      ! deallocate(degree_diff_local)
      deallocate (vol_local)
      deallocate (degree_local)
      deallocate (edge_local)
      deallocate (degree_local2)
      deallocate (adj_local)

      call MPI_Reduce(rec_local, rec_tot, 1, MPI_INT, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! prepare for gathering the data
      if (rank == io_rank) then
         allocate (sendcounts(size))
         allocate (pos(size))

         sendcounts = sendcounts_pts*3 + sendcounts_edge*2 + 2
         pos = pos_pts*3 + pos_edge*2 + (/(2*(i - 1), i=1, size)/)

         allocate (data_global_real(grid%num_pts_global*3 + rec_tot*2 + 2*size))
         deallocate (sendcounts_pts)
         deallocate (pos_pts)
         deallocate (sendcounts_edge)
         deallocate (pos_edge)
      endif

      if (atts%coefficients_are_dedudded .eqv. PETSC_TRUE) then
         call MPI_Gather((grid%num_pts_local*3 + rec_local*2 + 2), 1, MPI_INT, &
                         sendcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

         if (rank == io_rank) then
            pos = 0
            if (size > 1) pos(2:size) = (/(sum(sendcounts(1:i)), i=1, size)/)
         endif
      endif

      ! gather the data to io_rank
      call MPI_Gatherv(data_local_real, grid%num_pts_local*3 + rec_local*2 + 2, &
                       MPI_DOUBLE_PRECISION, data_global_real, sendcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! print to file
      if (rank == io_rank) then

         fileid = 49
         allocate (num_pts_0(size))
         allocate (rec_local_0(size))
         allocate (indx(size))
         allocate (voronoi_volumes(grid%num_pts_global))
         allocate (pointer_stride(size))
         allocate (row_entries(size))
         allocate (geom_coeffs(rec_tot))
         allocate (coeff_pointers(rec_tot))
         allocate (off_diagonal_map(rec_tot))
         allocate (on_diagonal_map(size))

55       format(5(1PE20.12))
56       format(5i15)

         ! -------------------------------------------------------------
         ! Break data_global_real into discrete arrays
         ! This is for later manipulation, but also benefits readability
         ! This needs to be significantly refactored. Ugh.

         indx = (/pos(2:size) - 1, pos(size) + sendcounts(size) - 1/)
         indx = indx - 2*mask_delta
         num_pts_0 = int(data_global_real(indx))
         rec_local_0 = int(data_global_real(indx + 1))

         indx = pos + 1
         voronoi_volumes = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/) ! <- what?

         allocate (indx_save(size))
         indx = indx + num_pts_0
         indx_save = indx + num_pts_0
         pointer_stride = (/(int(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)
         row_entries = (/(int(data_global_real(indx_save(i):indx_save(i) - 1 + rec_local_0(i))), i=1, size)/)
         indx = indx_save
         deallocate (indx_save)

         indx = indx + rec_local_0
         off_diagonal_map = (/(i, i=1, rec_tot)/)
         on_diagonal_map = (/(int(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)
         indx = indx + num_pts_0
         geom_coeffs = (/(-(data_global_real(indx(i):indx(i) - 1 + rec_local_0(i))), i=1, size)/)

         ! --------------------------------------------------------------------------------

         ! Open outfile for writing
         if (grid%dump_flag .EQV. PETSC_TRUE) then
            open (fileid, file=trim(grid%dump_str))
         else
            open (fileid, file='voronoi.stor')
         endif

         ! -------------------------------------------------------------
         ! Write out the header
         write (fileid, '(a)') 'fehmstor ascir8i4 Sparse Matrix Voronoi Coefficients'
         write (fileid, *) fdate()
         write (fileid, 56) rec_tot, grid%num_pts_global, rec_tot + grid%num_pts_global + 1, num_area_coef, int(degree_tot_max)

         ! -------------------------------------------------------------
         ! Write out the Voronoi volumes
         write (fileid, 55) voronoi_volumes

         ! -------------------------------------------------------------
         ! Write out the count for each row, row entries

         write (fileid, 56) pointer_stride, rec_tot + 1 + grid%num_pts_global, row_entries

         ! -------------------------------------------------------------
         ! Write out the 'pointers' and N_vertices zeros

         allocate (tmp_map(rec_tot))

         ! rec tot may not be stable under epsilon removal
         if (atts%coefficients_are_compressed .eqv. PETSC_FALSE) then
            off_diagonal_map = (/(i, i=1, rec_tot)/)
            coeff_count = rec_tot
         else
            call CompressCoefficients(geom_coeffs, tmp_map, off_diagonal_map, rec_tot, coeff_count)
         endif

         write (fileid, 56) off_diagonal_map, (0, i=1, grid%num_pts_global + 1)

         ! -------------------------------------------------------------
         ! Write out the pointers indicating which represent volumes
         write (fileid, 56) on_diagonal_map

         ! -------------------------------------------------------------
         ! Write out the geometric coefficients

         if (atts%coefficients_are_compressed .eqv. PETSC_FALSE) then
            write (fileid, 55) geom_coeffs
         else
            write (fileid, 55) (geom_coeffs(tmp_map(i)), i=1, coeff_count)
         endif

         deallocate (indx)
         deallocate (data_global_real)
         deallocate (tmp_map)
         deallocate (pointer_stride)
         deallocate (voronoi_volumes)
         deallocate (row_entries)

         close (fileid)
      endif

   end subroutine GridWriteFEHM

!**************************************************************************

   subroutine GridWriteHDF5(grid, rank, size)

      ! Writes out mesh to a compressed HDF5 file format.
      ! The structure of the file format is:
      !   |_  material
      !   |   geom
      !    \_ x
      !    \_ y
      !    \_ z
      !   |   cell
      !    \_ volume
      ! AREA
      ! LEN
      ! CONNECTIVITY

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat
      use HDF5

      implicit none

      type(grid_type) :: grid
      PetscInt :: rank, size, ncols, istart, iend, i0, i1
      PetscInt :: temp_int
      PetscInt :: io_rank = 0
      INTEGER :: hdferr
      INTEGER(HID_T) :: file, space2, dset3 ! handles
      INTEGER(HSIZE_T), DIMENSION(1:2) :: dims ! size read/write buffer
      INTEGER :: i
      integer :: ierr ! error handler

      PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
      PetscReal, allocatable :: vol_local(:), x(:), y(:), z(:), x0(:), y0(:), z0(:)
      PetscReal, allocatable :: length_matrix(:, :)
      PetscReal, allocatable :: area_local(:), temp_connections(:), area_tmp(:), area_matrix(:, :)
      PetscReal, allocatable :: length_local(:), length_tmp(:), volumes_full(:)
      PetscInt, allocatable :: recvcounts(:), pos(:)
      Vec :: mat_diag

      ! TODO: add connectivity

      ! =====================================================
      ! GET VORONOI CELL VOLUMES
      ! =====================================================

      ! Pull the Voronoi cell volumes from PETSc matrix 'adjmatrix'
      allocate (vol_local(grid%num_pts_local))
      call VecCreate(PETSC_COMM_WORLD, mat_diag, ierr); CHKERRQ(ierr)
      call VecSetSizes(mat_diag, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
      call VecSetType(mat_diag, VECMPI, ierr); CHKERRQ(ierr)
      call MatGetDiagonal(grid%adjmatrix, mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)

      vol_local = vec_ptr

      call VecRestoreArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      call VecDestroy(mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(grid%degree, vec_ptr, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (volumes_full(grid%num_pts_global))
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = 0
         if (size > 1) pos(2:size) = (/(sum(recvcounts(1:i)), i=1, size)/)
      endif

      call MPI_Gatherv(vol_local, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, volumes_full, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         deallocate (recvcounts)
         deallocate (pos)
      endif

      ! =====================================================
      ! GET VORONOI X,Y,Z CENTROIDS
      ! =====================================================

      allocate (x0(grid%num_pts_local))
      allocate (y0(grid%num_pts_local))
      allocate (z0(grid%num_pts_local))

      !  Get coordinate vector (arranged as: (/xyzxyzxyzxyz.../))
      call VecGetArrayF90(grid%coordinates_vec, vec_ptr2, ierr); CHKERRQ(ierr)

      x0(:) = vec_ptr2(1::3) ! Get x1, skip 3, get x2, ...
      y0(:) = vec_ptr2(2::3) ! Get y1, skip 3, get x2, ...
      z0(:) = vec_ptr2(3::3) ! Get z1, skip 3, get z2, ...

      call VecRestoreArrayF90(grid%coordinates_vec, vec_ptr2, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (x(grid%num_pts_global))
         allocate (y(grid%num_pts_global))
         allocate (z(grid%num_pts_global))
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = 0
         if (size > 1) pos(2:size) = (/(sum(recvcounts(1:i)), i=1, size)/)
      endif

      call MPI_Gatherv(x0, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, x, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      call MPI_Gatherv(y0, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, y, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      call MPI_Gatherv(z0, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, z, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      deallocate (x0)
      deallocate (y0)
      deallocate (z0)

      if (rank == io_rank) then
         deallocate (recvcounts)
         deallocate (pos)
      endif

      ! ====================================================
      !  Parse out area data
      !  (upper triangular represented as full matrix)
      ! ====================================================
      call MatGetOwnershipRange(grid%adjmatrix_area, istart, iend, ierr); CHKERRQ(ierr)
      allocate (area_local((iend - istart)*grid%num_pts_global))
      allocate (temp_connections(grid%num_pts_global))

      do i = 1, iend - istart

         temp_connections = 0.
         temp_int = int(vec_ptr(i))
         i0 = ((i - 1)*grid%num_pts_global) + 1
         i1 = ((i - 1)*grid%num_pts_global) + grid%num_pts_global

         call MatGetRow(grid%adjmatrix_area, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr); CHKERRQ(ierr)
         area_local(i0:i0 + (i + istart) - 1) = 0.
         area_local(i0 + (i + istart) - 1:i1) = temp_connections(1:temp_int)
         call MatRestoreRow(grid%adjmatrix_area, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr)
      enddo

      deallocate (temp_connections)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = pos*grid%num_pts_global
         recvcounts = recvcounts*grid%num_pts_global
         allocate (area_tmp(grid%num_pts_global*grid%num_pts_global))
      endif

      call MPI_Gatherv(area_local, (iend - istart)*grid%num_pts_global, &
                       MPI_DOUBLE_PRECISION, area_tmp, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (area_matrix(grid%num_pts_global, grid%num_pts_global))
         area_matrix = reshape(area_tmp, (/grid%num_pts_global, grid%num_pts_global/))

         deallocate (area_tmp)
         deallocate (pos)
         deallocate (recvcounts)
      endif

      deallocate (area_local)

      ! ====================================================
      !  Parse out edge length data
      !  (upper triangular represented as full matrix)
      ! ====================================================
      call MatGetOwnershipRange(grid%adjmatrix_len, istart, iend, ierr); CHKERRQ(ierr)
      allocate (length_local((iend - istart)*grid%num_pts_global))
      allocate (temp_connections(grid%num_pts_global))

      do i = 1, iend - istart

         temp_connections = 0.
         temp_int = int(vec_ptr(i))
         i0 = ((i - 1)*grid%num_pts_global) + 1
         i1 = ((i - 1)*grid%num_pts_global) + grid%num_pts_global

         call MatGetRow(grid%adjmatrix_len, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr); CHKERRQ(ierr)
         length_local(i0:i0 + (i + istart) - 1) = 0. ! zero out unfilled section of vector
         length_local(i0 + (i + istart) - 1:i1) = temp_connections(1:temp_int)
         call MatRestoreRow(grid%adjmatrix_len, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr)
      enddo

      deallocate (temp_connections)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = pos*grid%num_pts_global
         recvcounts = recvcounts*grid%num_pts_global
         allocate (length_tmp(grid%num_pts_global*grid%num_pts_global))
      endif

      call MPI_Gatherv(length_local, (iend - istart)*grid%num_pts_global, &
                       MPI_DOUBLE_PRECISION, length_tmp, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (length_matrix(grid%num_pts_global, grid%num_pts_global))
         length_matrix = reshape(length_tmp, (/grid%num_pts_global, grid%num_pts_global/))

         deallocate (length_tmp)
         deallocate (pos)
         deallocate (recvcounts)
      endif

      deallocate (length_local)

      ! =====================================================
      ! BEGIN WRITING TO HDF5 FILE
      ! =====================================================

      if (rank == io_rank) then
         ! Initialize FORTRAN interface.
         call h5open_f(hdferr)

         ! Initialize and open the file
         call h5fcreate_f(trim(grid%dump_str), H5F_ACC_TRUNC_F, file, hdferr)

         ! Write out volumes
         dims = (/grid%num_pts_global, 1/)
         ! Create a dataspace for volumes and write out vector
         call h5screate_simple_f(1, dims, space2, hdferr)
         call h5dcreate_f(file, 'VOLUMES', H5T_NATIVE_REAL, space2, dset3, hdferr)
         call h5dwrite_f(dset3, H5T_NATIVE_REAL, real(volumes_full), dims, hdferr)
         call h5dclose_f(dset3, hdferr)
         call h5sclose_f(space2, hdferr)

         ! Create a dataspace for X and write out vector
         call h5screate_simple_f(1, dims, space2, hdferr)
         call h5dcreate_f(file, 'X', H5T_NATIVE_REAL, space2, dset3, hdferr)
         call h5dwrite_f(dset3, H5T_NATIVE_REAL, real(x), dims, hdferr)
         call h5dclose_f(dset3, hdferr)
         call h5sclose_f(space2, hdferr)

         ! Create a dataspace for Y and write out vector
         call h5screate_simple_f(1, dims, space2, hdferr)
         call h5dcreate_f(file, 'Y', H5T_NATIVE_REAL, space2, dset3, hdferr)
         call h5dwrite_f(dset3, H5T_NATIVE_REAL, real(y), dims, hdferr)
         call h5dclose_f(dset3, hdferr)
         call h5sclose_f(space2, hdferr)

         ! Create a dataspace for Z and write out vector
         call h5screate_simple_f(1, dims, space2, hdferr)
         call h5dcreate_f(file, 'Z', H5T_NATIVE_REAL, space2, dset3, hdferr)
         call h5dwrite_f(dset3, H5T_NATIVE_REAL, real(z), dims, hdferr)
         call h5dclose_f(dset3, hdferr)
         call h5sclose_f(space2, hdferr)

         ! Create a dataspace for connection area and write out matrix
         dims = (/grid%num_pts_global, grid%num_pts_global/)
         CALL h5screate_simple_f(2, dims, space2, hdferr)
         CALL h5dcreate_f(file, 'AREAS', H5T_NATIVE_REAL, space2, dset3, hdferr)
         CALL h5dwrite_f(dset3, H5T_NATIVE_REAL, real(area_matrix), dims, hdferr)
         CALL h5dclose_f(dset3, hdferr)
         CALL h5sclose_f(space2, hdferr)

         ! Create a dataspace for edge lengths and write out matrix
         call h5screate_simple_f(2, dims, space2, hdferr)
         call h5dcreate_f(file, 'LENGTHS', H5T_NATIVE_REAL, space2, dset3, hdferr)
         call h5dwrite_f(dset3, H5T_NATIVE_REAL, real(length_matrix), dims, hdferr)
         call h5dclose_f(dset3, hdferr)
         call h5sclose_f(space2, hdferr)

         ! Close and release resources.
         call h5fclose_f(file, hdferr)
         deallocate (vol_local)
         deallocate (x)
         deallocate (y)
         deallocate (z)
         deallocate (area_matrix)
         deallocate (length_matrix)
      endif

   end subroutine GridWriteHDF5

!**************************************************************************

   subroutine GridWriteTOUGH2(grid, rank, size, ELNE_len)
      !
      ! Writes out mesh to TOUGH2 file format
      !
      ! Authors: Daniel Livingston, Mikey Hannon, Manual Sentis, and Zhuolin Qu
      ! Date: 04/11/20018
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat

      implicit none

      type(grid_type) :: grid
      PetscInt :: ncols, istart, iend
      PetscInt :: fileid, fileid_zone, rank, irank, size
      PetscInt :: rec_tot_tri, rec_local_tri
      PetscInt :: io_rank = 0
      PetscInt :: i, j, k, ELNE_len
      character(len=ELNE_len), allocatable :: ELE(:)
      character(len=5), allocatable        :: MAT(:)
      character(len=4)   :: nnum
      PetscInt :: NN, countmat, isot
      PetscInt, allocatable :: MM(:), recvcounts(:)
      PetscReal :: AHTX, temp_real, cosine, distance
      PetscReal, allocatable :: x(:), y(:), z(:)
      PetscInt :: temp_int, temp_int2, pos1, pos2, vsize
      PetscReal, pointer :: vec_ptr(:)
      PetscInt, allocatable :: sendcounts_pts(:), pos_pts(:)
      PetscInt, allocatable :: sendcounts_edge(:), pos_edge(:)
      PetscReal, allocatable :: vol_local(:), vol_global(:)
      PetscReal, allocatable :: edge_temp(:), edge_local(:)
      PetscInt, allocatable :: edge_global(:), degree_global(:)
      PetscReal, allocatable :: adj_temp_len(:), adj_local_len(:), adj_global_len(:)
      PetscReal, allocatable :: adj_temp_area(:), adj_local_area(:), adj_global_area(:)
      PetscInt, allocatable :: degree_local(:)
      PetscReal, allocatable :: data_local_real(:), data_global_real(:)
      PetscInt, allocatable :: num_pts_0(:), rec_local_tri_0(:)
      PetscInt, allocatable :: indx(:)
      PetscInt, allocatable :: sendcounts(:), pos(:)
      PetscReal, allocatable :: connectivity(:, :), connect_area(:, :), connect_map(:), connect_map_local(:), temp_connections(:)
      PetscReal, allocatable :: areax_tmp(:), areax_local(:), connectivity_local(:), connectivity_tmp(:)
      logical :: print_isot_warning = .true.
      Vec :: mat_diag
      PetscErrorCode :: ierr

      ! get the vol_local
      allocate (vol_local(grid%num_pts_local))
      call VecCreate(PETSC_COMM_WORLD, mat_diag, ierr); CHKERRQ(ierr)
      call VecSetSizes(mat_diag, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
      call VecSetType(mat_diag, VECMPI, ierr); CHKERRQ(ierr)

      call MatGetDiagonal(grid%adjmatrix_area, mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)

      vol_local = vec_ptr

      call VecRestoreArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      call VecDestroy(mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(grid%degree, vec_ptr, ierr); CHKERRQ(ierr)

      rec_local_tri = int(sum(vec_ptr - 1))
      allocate (edge_local(rec_local_tri))
      allocate (adj_local_len(rec_local_tri))
      allocate (adj_local_area(rec_local_tri))

      temp_int = int(maxval(vec_ptr))
      allocate (edge_temp(temp_int + 10))
      allocate (adj_temp_area(temp_int))
      allocate (adj_temp_len(temp_int))

      pos1 = 0
      pos2 = 0

      ! ==================================================
      !  Prepare connectivity matrices
      ! ==================================================

      call VecGetOwnershipRange(grid%connect_map, istart, iend, ierr); CHKERRQ(ierr)
      allocate (connect_map_local(iend - istart))
      call VecGetValues(grid%connect_map, iend - istart, (/(i, i=istart, iend)/), connect_map_local, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) allocate (connect_map(grid%num_pts_global))

      call MPI_Gatherv(connect_map_local, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, connect_map, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         deallocate (pos)
         deallocate (recvcounts)
      endif

      if (rank == io_rank) vsize = int(maxval(connect_map))
      call MPI_Bcast(vsize, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      do i = 1, grid%num_pts_local
         if (grid%vertex_ids(i) /= grid%num_pts_global) then

            temp_int = int(vec_ptr(i))

            call MatGetRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)
            edge_local(pos1 + 1:pos1 + temp_int - 1) = int(edge_temp(2:temp_int))
            call MatRestoreRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)

            pos1 = pos1 + temp_int - 1
            temp_int2 = int(vec_ptr(i) - 1)

            call MatGetRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_len, ierr); CHKERRQ(ierr)
            adj_local_len(pos2 + 1:pos2 + temp_int2) = adj_temp_len(1:temp_int2)
            call MatRestoreRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_len, ierr); CHKERRQ(ierr)

            call MatGetRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_area, ierr); CHKERRQ(ierr)
            adj_local_area(pos2 + 1:pos2 + temp_int2) = adj_temp_area(2:temp_int)
            call MatRestoreRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_area, ierr); CHKERRQ(ierr)

            pos2 = pos2 + temp_int2
         endif
      enddo

      ! ====================================================
      !  Parse out connectivity data
      ! ====================================================
      call MatGetOwnershipRange(grid%connectivity, istart, iend, ierr); CHKERRQ(ierr)
      allocate (connectivity_local((iend - istart)*vsize))
      allocate (temp_connections(vsize))

      do i = 1, iend - istart
         temp_connections = 0.
         call MatGetRow(grid%connectivity, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr); CHKERRQ(ierr)
         connectivity_local(((i - 1)*vsize) + 1:((i - 1)*vsize) + vsize) = temp_connections
         call MatRestoreRow(grid%connectivity, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr)
      enddo

      deallocate (temp_connections)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = pos*vsize
         recvcounts = recvcounts*vsize
         allocate (connectivity_tmp(grid%num_pts_global*vsize))
      endif

      call MPI_Gatherv(connectivity_local, (iend - istart)*vsize, &
                       MPI_DOUBLE_PRECISION, connectivity_tmp, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (connectivity(grid%num_pts_global, vsize))
         connectivity = reshape(connectivity_tmp, (/grid%num_pts_global, vsize/), order=(/2, 1/))

         deallocate (connectivity_tmp)
         deallocate (pos)
         deallocate (recvcounts)
      endif

      deallocate (connectivity_local)

      ! ====================================================
      !  Parse out AREAX data
      ! ====================================================
      call MatGetOwnershipRange(grid%connect_area, istart, iend, ierr); CHKERRQ(ierr)
      allocate (areax_local((iend - istart)*vsize))
      allocate (temp_connections(vsize))

      do i = 1, iend - istart
         temp_connections = 0.
         call MatGetRow(grid%connect_area, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr); CHKERRQ(ierr)
         areax_local(((i - 1)*vsize) + 1:((i - 1)*vsize) + vsize) = temp_connections
         call MatRestoreRow(grid%connect_area, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr)
      enddo

      deallocate (temp_connections)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = pos*vsize
         recvcounts = recvcounts*vsize
         allocate (areax_tmp(grid%num_pts_global*vsize))
      endif

      call MPI_Gatherv(areax_local, (iend - istart)*vsize, &
                       MPI_DOUBLE_PRECISION, areax_tmp, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (connect_area(grid%num_pts_global, vsize))
         connect_area = reshape(areax_tmp, (/grid%num_pts_global, vsize/), order=(/2, 1/))

         deallocate (areax_tmp)
         deallocate (pos)
         deallocate (recvcounts)
      endif

      deallocate (areax_local)

      ! ====================================================

      if (rank == io_rank) then
         allocate (sendcounts_pts(size))
         allocate (sendcounts_edge(size)) ! for the upper triangular case without diagonal element
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, sendcounts_pts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(rec_local_tri, 1, MPI_INT, sendcounts_edge, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (pos_pts(size))
         allocate (pos_edge(size))
         pos_pts(1) = 0
         do irank = 1, size - 1
            pos_pts(irank + 1) = pos_pts(irank) + sendcounts_pts(irank)
         enddo
         pos_edge(1) = 0
         do irank = 1, size - 1
            pos_edge(irank + 1) = pos_edge(irank) + sendcounts_edge(irank)
         enddo
         rec_tot_tri = pos_edge(size) + sendcounts_edge(size)
      endif

      ! calculate degree_local
      allocate (degree_local(grid%num_pts_local))
      call MPI_Scatter(pos_edge, 1, MPI_INT, degree_local(1), 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      do i = 2, grid%num_pts_local
         degree_local(i) = degree_local(i - 1) + int(vec_ptr(i - 1) - 1)
      enddo

      degree_local = degree_local + 1

      call VecRestoreArrayReadF90(grid%degree, vec_ptr, ierr); CHKERRQ(ierr)

#if DEBUG
      print *, rank, 'vol_local=', vol_local
      print *, rank, 'degree_local', degree_local
      print *, rank, 'edge_lcoal=', edge_local
      print *, rank, 'adj_local_len=', adj_local_len
      print *, rank, 'adj_local_area=', adj_local_area
      print *, rank, 'rec_local_tri=', rec_local_tri
#endif

      allocate (data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 2))
      data_local_real(1:grid%num_pts_local) = vol_local
      data_local_real(grid%num_pts_local + 1:grid%num_pts_local*2) = degree_local
      data_local_real(grid%num_pts_local*2 + 1:grid%num_pts_local*2 + rec_local_tri) = edge_local

      call VecGetArrayReadF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)

      data_local_real(grid%num_pts_local*2 + rec_local_tri + 1:grid%num_pts_local*3 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 1), i=1, grid%num_pts_local)/)
      data_local_real(grid%num_pts_local*3 + rec_local_tri + 1:grid%num_pts_local*4 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 2), i=1, grid%num_pts_local)/)
      data_local_real(grid%num_pts_local*4 + rec_local_tri + 1:grid%num_pts_local*5 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 3), i=1, grid%num_pts_local)/)

      call VecRestoreArrayReadF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)

      data_local_real(grid%num_pts_local*5 + rec_local_tri + 1:grid%num_pts_local*5 + rec_local_tri*2) = adj_local_len
      data_local_real(grid%num_pts_local*5 + rec_local_tri*2 + 1:grid%num_pts_local*5 + rec_local_tri*3) = adj_local_area
      data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 1) = grid%num_pts_local
      data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 2) = rec_local_tri

      deallocate (vol_local)
      deallocate (degree_local)
      deallocate (edge_local)
      deallocate (adj_local_len)
      deallocate (adj_local_area)

      call VecDestroy(grid%coordinates_vec, ierr); CHKERRQ(ierr)

      ! prepare for gathering the data
      if (rank == io_rank) then
         allocate (sendcounts(size))
         allocate (pos(size))
         sendcounts = sendcounts_pts*5 + sendcounts_edge*3 + 2
         pos = pos_pts*5 + pos_edge*3 + (/(2*(i - 1), i=1, size)/)
         allocate (data_global_real(grid%num_pts_global*5 + rec_tot_tri*3 + 2*size))
         deallocate (sendcounts_pts)
         deallocate (pos_pts)
         deallocate (sendcounts_edge)
         deallocate (pos_edge)
      endif

      ! GAther the data to io_rank
      call MPI_Gatherv(data_local_real, grid%num_pts_local*5 + rec_local_tri*3 + 2, &
                       MPI_DOUBLE_PRECISION, data_global_real, sendcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      deallocate (data_local_real)

9110  format(A80)
9140  format(i10, 9i11)
9150  format(A4)
9160  format(10X, I11)

      ! print to file
      if (rank == io_rank) then
         fileid_zone = 51

         ! ---------------------------------------------
         ! Build the ELEME materials char array
         ! ---------------------------------------------

         countmat = 0
         allocate (MAT(grid%num_pts_global))

         ! If a zone file is provided, then...
         if (grid%zone_flag .EQV. PETSC_TRUE) then
            open (fileid_zone, file=grid%zone_str)

            do i = 1, 100000 ! TOUGH allows for a max n of nodes
               read (fileid_zone, 9150) nnum
               if (nnum .eq. 'nnum') then
                  countmat = countmat + 1
                  read (fileid_zone, 9160) NN
                  allocate (MM(NN))
                  read (fileid_zone, 9140) (MM(j), j=1, NN)
                  MAT((/(MM(j), j=1, NN)/)) = Material(countmat)
                  deallocate (MM)
               endif
               if (nnum .eq. 'stop') exit
            enddo
            close (fileid_zone)

            ! else, read from the IMT values array
         else

            ! Does this need a do loop? MAT = grid%imt_values
            do i = 1, grid%num_pts_global
               MAT(i) = grid%imt_values(i)
            enddo

         endif

         AHTX = 0.d0
         allocate (ELE(grid%num_pts_global))
         allocate (num_pts_0(size))
         allocate (rec_local_tri_0(size))
         allocate (indx(size))

         indx = (/pos(2:size) - 1, pos(size) + sendcounts(size) - 1/)
         num_pts_0 = int(data_global_real(indx))
         rec_local_tri_0 = int(data_global_real(indx + 1))

         ! global volumes
         indx = pos + 1
         deallocate (sendcounts)
         deallocate (pos)

         vol_global = (/(dabs(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)

         allocate (x(grid%num_pts_global))
         allocate (y(grid%num_pts_global))
         allocate (z(grid%num_pts_global))
         allocate (adj_global_len(rec_tot_tri))
         allocate (adj_global_area(rec_tot_tri))
         allocate (degree_global(grid%num_pts_global))

         indx = indx + num_pts_0
         degree_global = (/(int(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)
         indx = indx + num_pts_0
         edge_global = (/(int(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)

         indx = indx + rec_local_tri_0
         x = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         y = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         z = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         adj_global_len = (/(dabs(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)
         indx = indx + rec_local_tri_0
         adj_global_area = (/(dabs(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)

         deallocate (data_global_real)
         deallocate (num_pts_0)
         deallocate (rec_local_tri_0)
         deallocate (indx)

         ! ---------------------------------------------
         ! Begin the formatting and writing process
         ! ---------------------------------------------

         fileid = 49

         if (grid%dump_flag .EQV. PETSC_TRUE) then
            open (fileid, file=grid%dump_str)
         else
            open (fileid, file='MESH_VORONOI')
         endif

         ! write out the ELEME block
         write (fileid, 9110) 'ELEME----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8'

         do i = 1, grid%num_pts_global
            ! Get the element name
            call idxToELNE(ELE(i), i, ELNE_len)

            ! Assign AHTX
            ! TODO: Is this right?
            AHTX = connect_area(i, 1)

            ! and write line to file
            write (fileid, 9110) elemeTOUGH2(ELE(i), MAT(i), x(i), y(i), z(i), vol_global(i), AHTX)
         enddo

         ! write out the CONNE block
         write (fileid, *)
         write (fileid, 9110) 'CONNE----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8'

         pos1 = 1

         ! BETAX is the cosine between u, the connection vector, and Fg, the gravitational vector
         ! It is calculated using the following formula:
         ! cos(theta) = u * Fg / (||u|| * ||Fg||)
         !            = (9.81 * u_z) / (9.81 * sqrt(u_x^2 + u_y^2 + u_z^2))
         !            = (u_z) / (sqrt(u_x^2 + u_y^2 + u_z^2))

         !do i=1,grid%num_pts_global-1
         !do k=degree_global(i),degree_global(i+1)-1

         !  j = edge_global(k)

         ! Are we on the same node?
         !  if (j .eq. i) CYCLE

         ! Calculate the cosine of Fg and the connection vector
         !   temp_real = (x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
         !   cosine = (z(i)-z(j))/dsqrt(temp_real)

         !   if ((adj_global_len(pos1) .eq. 0.d0) .or. (adj_global_area(pos1) .lt. 1.d-12)) then
         !     pos1 = pos1 + 1
         !     cycle
         !   endif

         !    isot = getISOT((/ x(i), y(i), z(i) /), (/ x(j), y(j), z(j) /))

         !    write(fileid,9110) conneTOUGH2(ELE(i),ELE(j),isot,adj_global_len(pos1),adj_global_len(pos1),adj_global_area(pos1),cosine)
         !    pos1 = pos1 + 1
         !  enddo
         !enddo

         ! Iterate over all points...
         do i = 1, grid%num_pts_global

            ! Iterate over each connection...
            do k = 2, int(connect_map(i))

               j = int(connectivity(i, k))

               ! Are we on the same node?
               if (j .eq. i) CYCLE

               ! Calculate the cosine of Fg and the connection vector
               temp_real = (x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2
               cosine = (z(i) - z(j))/dsqrt(temp_real)

               ! TODO: This should be taken from len. NOT computed on the fly!
               distance = 0.5*((z(j) - z(i))**2.+(y(j) - y(i))**2.+(x(j) - x(i))**2.)**0.5

               ! TODO: Doesn't work for small values!!
               if (connect_area(i, k) .lt. 0.0001) CYCLE

               ! TODO: This should be parallel processed. NOT computed on the fly (single core!)
               isot = getISOT((/x(i), y(i), z(i)/), (/x(j), y(j), z(j)/))

               if ((print_isot_warning .eqv. .true.) .and. (isot == 0)) then
                  print *, 'WARNING: At least one non-conforming connection vector. Setting ISOT to 0.'
                  print_isot_warning = .false.
               endif

               write (fileid, 9110) conneTOUGH2(ELE(i), ELE(j), isot, distance, distance, connect_area(i, k), cosine)
            enddo
         enddo

         ! Write empty line at EOF & close
         write (fileid, *) ''
         close (fileid)

         deallocate (MAT)
         deallocate (ELE)
         deallocate (vol_global)
         !  deallocate(degree_global)
         deallocate (edge_global)
         deallocate (x)
         deallocate (y)
         deallocate (z)
         deallocate (adj_global_len)
         deallocate (adj_global_area)

      endif

   end subroutine GridWriteTOUGH2

   subroutine CompressCoefficients(values, map, indx, N, max_index)
      !
      ! Subroutine for compressing FEHM STOR file
      ! geometric coefficients
      !
      ! Usage:
      !
      !  values(indx1(i))
      !
      ! Parameters:
      !
      !    values (Nx1 array) - contains geometric coeff. values
      !    N (integer) - size of values
      !    indx1 (Nx1 array) (output) - remapped pointers to coeff. indices
      !

      implicit none

      PetscInt  :: N, max_index, indx(N), map(N)
      PetscReal :: values(N)
      real*8 :: x0, x1, epsilon = 0.00001, ascend = -1.0
      real*8, allocatable :: values_8(:)
      integer*8, allocatable :: map_8(:), indx_8(:)
      integer*8 :: N_8, i, j

      allocate (map_8(N))
      allocate (values_8(N))
      allocate (indx_8(N))

      map_8 = (/(i, i=1, N)/)
      indx_8 = 0
      j = 1

      ! Type promotion must be done to accomodate the <T>*8 LaGriT uses
      N_8 = N
      values_8 = values

      ! Call a sorting function on the geometric coefficients - uses `map` as a key
      ! hpsort1 : in namespace `LaGriT::lg_util`
      call hpsort1(N_8, values_8, ascend, map_8)

      ! Iterate over the sorted vector, capturing 'unique' values.
      ! Unique values are defined as elements in the sorted vector
      ! that have a .gt. `epsilon` difference than the previous vector element.
      ! `j` is the iterator counting the unique values.
      x0 = values_8(map_8(1))
      do i = 1, N_8
         x1 = values_8(map_8(i))
         if (abs(x1 - x0) > epsilon) then
            j = j + 1
            values_8(map_8(j)) = x1
         endif
         indx_8(map_8(i)) = j
         x0 = x1
      enddo

      ! Promote <T>*8 variables back to PETSc precision
      max_index = j
      values = values_8
      indx = indx_8
      map = map_8

      deallocate (map_8)
      deallocate (values_8)
      deallocate (indx_8)

   end subroutine CompressCoefficients

!**************************************************************************

   character(len=5) function Material(countmat)
      ! Returns a 5-char string indicating the material type of
      ! the fracture.
      !
      ! Author: Zhuolin Qu, LANL
      ! Revised: Daniel Livingston, LANL
      ! Revision: 03/01/2018
      !
      implicit none

      PetscInt :: countmat

      SELECT CASE (countmat)
      CASE (1)
         Material = 'clayr'
      CASE (2)
         Material = 'waste'
      CASE (3)
         Material = 'bento'
      CASE (4)
         Material = 'inwas'
      CASE (5)
         Material = 'edzal'
      CASE (6)
         Material = 'backf'
      CASE (7)
         Material = 'inben'
      CASE (8)
         Material = 'inbac'
      CASE (9)
         Material = 'boun1'
      CASE (10)
         Material = 'bounf'
      CASE (11)
         Material = 'sourc'
      CASE (12)
         Material = 'mbent'
      CASE (13)
         Material = 'medzw'
      CASE (14)
         Material = 'mlbkf'
      CASE (15)
         Material = 'msbkf'
      CASE (16)
         Material = 'edzal'
      CASE (17)
         Material = 'mbcez'
      CASE (18)
         Material = 'edzat'
      CASE (19)
         Material = 'mbeca'
      CASE (20)
         Material = 'mbkad'
      CASE DEFAULT
         Material = '?????'
      END SELECT
   end function

!**************************************************************************

   logical function epsEqual(value1, value2)
      !
      ! Returns 'true' if values are equal
      ! within some epsilon.
      !

      implicit none

      PetscReal :: value1, value2
      PetscReal :: eps = 2E-3
      epsEqual = (abs(abs(value2) - abs(value1)) .le. eps)
      return
   end function epsEqual

!**************************************************************************

   PetscInt function getISOT(v1, v2)
      !
      ! Given two points composing a vector, this function
      ! determines which axis or plane the vector lies on.
      ! It then returns an ISOT value accordingly.
      !
      ! TODO: Implement (theta,phi) for off-axis rotation
      !

      implicit none

      PetscReal :: v1(3)
      PetscReal :: v2(3)
      PetscReal :: inorm(3)
      PetscReal :: c0, c1

      c0 = 0.0
      c1 = 1.0

      inorm = v2 - v1
      inorm = inorm/norm2(inorm)

      if (epsEqual(inorm(1), c1) .eqv. .true.) then
         ! parallel to cardinal X
         getISOT = 1
      elseif (epsEqual(inorm(2), c1) .eqv. .true.) then
         ! parallel to cardinal Y
         getISOT = 2
      elseif (epsEqual(inorm(3), c1) .eqv. .true.) then
         ! parallel to cardinal Z
         getISOT = 3
      elseif (epsEqual(inorm(1), c0) .eqv. .true.) then
         ! within YZ plane
         getISOT = 2
      elseif (epsEqual(inorm(2), c0) .eqv. .true.) then
         ! within XZ plane
         getISOT = 1
      elseif (epsEqual(inorm(3), c0) .eqv. .true.) then
         ! within XY plane
         getISOT = 1
      else
         getISOT = 0
      endif

      return

   end function getISOT

!**************************************************************************

   subroutine getMatFromIMT(val, val_size, final_chr)

      !
      ! Given a vector of imt floating point values, presum-
      ! ably sourced from an AVS file, 'return' a vector of
      ! the same length with char(5) material names instead.
      ! Used for TOUGH2 ELEME materials when ZONE file isn't
      ! present.
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 03/20/2018
      !

      implicit none

      integer :: i, j, val_size, final_size
      PetscReal :: min_val, max_val
      PetscReal :: val(val_size), unique(val_size)
      PetscReal, dimension(:), allocatable :: final
      !character(len=2) :: tmp
      character(len=5) :: tmp
      character(len=5), dimension(val_size) :: final_chr

878   format(I0.1)

      ! ---------------------------------------------
      ! Get the unique values in imt vector
      ! ---------------------------------------------

      i = 0

      min_val = minval(val) - 1
      max_val = maxval(val)

      ! Sort imt values
      do while (min_val < max_val)
         i = i + 1
         min_val = minval(val, mask=val > min_val)
         unique(i) = min_val
      enddo

      ! Effectively trim unique to generate a
      ! 'de-duped' and sorted vector
      allocate (final(i), source=unique(1:i))

      ! ---------------------------------------------
      ! Map each attribute to its position in the
      ! final vector, then to a materials string
      ! ---------------------------------------------

      print *, 'WARNING: Not properly configured for IMT > E+01'
      print *, '         Also not configured for atts other than IMT.'

      final_size = i

      ! Iterate over each node attribute...
      do i = 1, val_size
         ! iterate over each de-duped attribute...
         do j = 1, final_size
            ! if we've found the index, then...
            if (val(i) == final(j)) then
               ! store in the character array & break
               write (tmp, 878) j
               final_chr(i) = adjustr(tmp)
               EXIT
            endif
         enddo
      enddo

      deallocate (final)

   end subroutine getMatFromIMT

!**************************************************************************

   character(len=80) function conneTOUGH2(elne1, elne2, isot, d1, d2, ax, bx)

      !
      ! Returns an CONNE line in the TOUGH2 filespec
      ! See the TOUGH2 manual for more information, pp 172:
      !   http://esd1.lbl.gov/FILES/research/projects/tough/documentation/TOUGH2_V2_Users_Guide.pdf
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 03/20/2018
      !

      implicit none

      character(len=10) :: tmp

      ! Begin variable definitions
      character(len=5) :: elne1, elne2 ! 1st & 2nd element names
      integer :: isot                 ! abs. permeability for elem. materials (1,2, or 3)
      PetscReal :: d1, d2              ! dst from elems to their common interface
      PetscReal :: ax, bx              ! interface area ; cos. of angle bt. grav vec & elem bndry
      PetscReal :: zero = 0.0

85    format(6ES10.4)
86    format(I1)
87    format(6E10.3)

      !TODO: ONLY 87 outputs values .lt. 0. Fix formatting!

      ! Set element names
      conneTOUGH2 = repeat(' ', 80)
      conneTOUGH2(1:5) = elne1
      conneTOUGH2(6:10) = elne2

      ! Write out permeability ID
      write (tmp, 86) isot
      conneTOUGH2(30:31) = tmp(1:2)

      ! Write out floating point values
      write (tmp, 85) d1
      conneTOUGH2(31:40) = tmp
      write (tmp, 85) d2
      conneTOUGH2(41:50) = tmp
      write (tmp, 85) ax
      conneTOUGH2(51:60) = tmp

      ! BETAX may be skipped if it is zero
      ! TOUGH assumes that blank columns for BETAX are 0.
      if (epsEqual(bx, zero) .eqv. .false.) then
         write (tmp, 87) bx
         conneTOUGH2(61:70) = tmp
      endif

      return

   end function conneTOUGH2

!**************************************************************************

   character(len=80) function elemeTOUGH2(elne, ma, x, y, z, volx, ahtx)

      !
      ! Returns an ELEME line in the TOUGH2 filespec
      ! See the TOUGH2 manual for more information, pp 171:
      !   http://esd1.lbl.gov/FILES/research/projects/tough/documentation/TOUGH2_V2_Users_Guide.pdf
      !
      ! Author: Daniel Livingston, LANL
      ! Date: 03/20/2018
      !

      implicit none

      character(len=10) :: tmp

      ! Mutable element variables
      character(len=5) :: elne   ! element name
      character(len=5) :: ma     ! material ID corresponding to a ROCKS entry
      PetscReal :: volx, ahtx    ! element volume & interface area
      PetscReal :: x, y, z       ! grid block center

75    format(6E10.4)

      ! Write out element & material names
      elemeTOUGH2 = repeat(' ', 80)
      elemeTOUGH2(1:5) = elne
      elemeTOUGH2(16:20) = ma

      ! Don't write out junk values
      if (abs(volx) < 1e-15) volx = 0.0
      if (abs(ahtx) < 1e-15) ahtx = 0.0
      if (abs(x) < 1e-15) x = 0.0
      if (abs(y) < 1e-15) y = 0.0
      if (abs(z) < 1e-15) z = 0.0

      ! Write out floating point values
      write (tmp, 75) volx
      elemeTOUGH2(21:30) = tmp
      write (tmp, 75) ahtx
      elemeTOUGH2(31:40) = tmp
      write (tmp, 75) x
      elemeTOUGH2(51:60) = tmp
      write (tmp, 75) y
      elemeTOUGH2(61:70) = tmp
      write (tmp, 75) z
      elemeTOUGH2(71:80) = tmp

      return

   end function elemeTOUGH2

!**************************************************************************

   subroutine idxToELNE(result, idx, length)

      !
      ! Given an integer, returns a five-letter code where:
      !  - The first three chars are letters or numbers (A-Z;a-z;0-9)
      !  - The last two chars are numbers (0-9)
      !
      ! This subroutine 'chops off' the last two digits of idx to store
      ! as the last two chars of 'result'.
      ! For the remainder of idx, we take the modulus of idx to the sum of
      ! all allowable characters (26+26+10), and assign an ASCII character
      ! to it based on where it falls.
      !

      implicit none

      integer :: num, cchr, remainder, idx, length, nn
      integer :: prefix, postfix, i, permutations
      character :: result(length)
      character :: digit(1)
      character(len=1) :: FMT
      character(len=length - 3) :: tmp

      if ((length .lt. 5) .OR. (length .gt. 9)) then
         print *, 'ERROR: length must be in range [5,9]'
         print *, 'Defaulting to length = 5'
         length = 5
      endif

      ! [A-Z; a-z; 0-9]
      permutations = 26 + 26 + 10

      ! Count of integers in final ELNE; len(NE)
      nn = length - 3

      ! prefix captures all numbers except the last two
      prefix = int(floor(idx/real(10**nn)))

      ! postfix captures the last two numbers
      postfix = int(mod(idx, 10**nn))

      ! Add 1 so we still get 'A' in the ASCII conversion
      num = prefix + 1

      i = 0
      result = 'A'

      ! If idx is less than 100, then
      ! do not bother with the loop
      if (prefix .eq. 0) then
         num = -1
      endif

      ! The amount of times this loop interates reflects
      ! the number of characters in result it replaces,
      ! with an upper bound of 3 ('AAA')
      do while (num > 0)
         num = num - 1
         remainder = mod(num, int(permutations))
         cchr = remainder

         ! A - Z
         if (cchr .lt. 26) then
            cchr = cchr + 65

            ! When this is true, we are moving into
            ! the second character space.
            ! Make sure that it doesn't repeat 'A'
            if (i .gt. 0) then
               cchr = cchr + 1
            endif

         ! a - z
         elseif ((cchr .ge. 26) .AND. (cchr .lt. 52)) then
            cchr = cchr - 27 + 97 + 1

         ! 0 - 9
         elseif ((cchr .ge. 52)) then
            cchr = cchr - 52 + 48
         endif

         ! Convert to ASCII char
         digit = char(cchr)

         ! Which index of 'A' should we replace?
         result(3 - i) = digit(1)

         ! Do we have any iterations left?
         num = int((num - remainder)/real(permutations))
         i = i + 1
      enddo

      write (FMT, "(I1)") nn
      write (tmp, "(I0."//ADJUSTL(FMT)//")") postfix

      ! Why does it need to have 1 added?
      ! 'Access at index 0' error if i = 1, nn
      result(4:length) = (/tmp(1:1), tmp(2:2)/)!(/ (tmp(i+1:i+1), i = 0, nn-1) /)

   end subroutine idxToELNE

!**************************************************************************

   subroutine CreateConnMatrix(grid, size, rank)
      !
      ! Create two connectivity matrices for elements
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 08/21/2015
      !

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
      use petscis

      implicit none

      type(grid_type) :: grid
      PetscErrorCode :: ierr
      PetscInt :: i, size, a1, b1, a2, b2, a3, b3, num_elems_local_save, rank
      Mat :: submat1, submat2, submat3
      PetscInt :: conn(grid%ndim + 1)
      PetscInt :: io_rank = 0
      PetscScalar :: v1(1), v2(1), v3(1)
      IS :: irow1, icol1, irow2, icol2, irow3, icol3
      PetscScalar :: temp_int

      call mpi_print('   Building connectivity matrices', rank, io_rank)

      call MatDuplicate(grid%edgematrix, MAT_SHARE_NONZERO_PATTERN, grid%conn1, ierr); CHKERRQ(ierr)
      call MatDuplicate(grid%edgematrix, MAT_SHARE_NONZERO_PATTERN, grid%conn2, ierr); CHKERRQ(ierr)

      ! MatGetSubMatrices <-> MatCreateSubMatrices

      num_elems_local_save = grid%num_elems_global/size

      do i = 1, num_elems_local_save + 1
         if (i .le. grid%num_elems_local) then
            conn = grid%elem_connectivity(i, :)
            temp_int = int(grid%elem_ids(i))

            a1 = min0(conn(1), conn(2)); b1 = max0(conn(1), conn(2))
            a2 = min0(conn(1), conn(3)); b2 = max0(conn(1), conn(3))
            a3 = min0(conn(2), conn(3)); b3 = max0(conn(2), conn(3))

            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/a1 - 1/), PETSC_COPY_VALUES, irow1, ierr); CHKERRQ(ierr)
            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/b1 - 1/), PETSC_COPY_VALUES, icol1, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow1, icol1, MAT_INITIAL_MATRIX, submat1, ierr); CHKERRQ(ierr)
            call MatGetValues(submat1, 1, (/0/), 1, (/0/), v1, ierr); CHKERRQ(ierr)
            call ISDestroy(irow1, ierr); CHKERRQ(ierr)
            call ISDestroy(icol1, ierr); CHKERRQ(ierr)

            if (int(v1(1)) == 0) then
               call MatSetValues(grid%conn1, 1, a1 - 1, 1, b1 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            else
               call MatSetValues(grid%conn2, 1, a1 - 1, 1, b1 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            endif

            call MatAssemblyBegin(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
            call MatAssemblyEnd(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/a2 - 1/), PETSC_COPY_VALUES, irow2, ierr); CHKERRQ(ierr)
            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/b2 - 1/), PETSC_COPY_VALUES, icol2, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow2, icol2, MAT_INITIAL_MATRIX, submat2, ierr); CHKERRQ(ierr)
            call MatGetValues(submat2, 1, (/0/), 1, (/0/), v2, ierr); CHKERRQ(ierr)

            if (int(v2(1)) == 0) then
               call MatSetValues(grid%conn1, 1, a2 - 1, 1, b2 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            else
               call MatSetValues(grid%conn2, 1, a2 - 1, 1, b2 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            endif

            call MatAssemblyBegin(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
            call MatAssemblyEnd(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
            call ISDestroy(irow2, ierr); CHKERRQ(ierr)
            call ISDestroy(icol2, ierr); CHKERRQ(ierr)

            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/a3 - 1/), PETSC_COPY_VALUES, irow3, ierr); CHKERRQ(ierr)
            call ISCreateGeneral(PETSC_COMM_SELF, 1, (/b3 - 1/), PETSC_COPY_VALUES, icol3, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow3, icol3, MAT_INITIAL_MATRIX, submat3, ierr); CHKERRQ(ierr)
            call MatGetValues(submat3, 1, (/0/), 1, (/0/), v3, ierr); CHKERRQ(ierr)

            if (int(v3(1)) == 0) then
               call MatSetValues(grid%conn1, 1, a3 - 1, 1, b3 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            else
               call MatSetValues(grid%conn2, 1, a3 - 1, 1, b3 - 1, temp_int, INSERT_VALUES, ierr); CHKERRQ(ierr)
            endif

            call MatAssemblyBegin(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
            call MatAssemblyEnd(grid%conn1, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
            call ISDestroy(irow3, ierr); CHKERRQ(ierr)
            call ISDestroy(icol3, ierr); CHKERRQ(ierr)

         else

            call ISCreateGeneral(PETSC_COMM_SELF, 0, (/1/), PETSC_COPY_VALUES, irow1, ierr); CHKERRQ(ierr)
            call ISCreateGeneral(PETSC_COMM_SELF, 0, (/1/), PETSC_COPY_VALUES, icol1, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow1, icol1, MAT_INITIAL_MATRIX, submat1, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow1, icol1, MAT_INITIAL_MATRIX, submat2, ierr); CHKERRQ(ierr)
            call MatGetSubMatrices(grid%conn1, 1, irow1, icol1, MAT_INITIAL_MATRIX, submat3, ierr); CHKERRQ(ierr)
            call ISDestroy(irow1, ierr); CHKERRQ(ierr)
            call ISDestroy(icol1, ierr); CHKERRQ(ierr)

         endif
      enddo

      call MatDestroy(submat1, ierr); CHKERRQ(ierr)
      call MatDestroy(submat2, ierr); CHKERRQ(ierr)
      call MatDestroy(submat3, ierr); CHKERRQ(ierr)
      call MatAssemblyBegin(grid%conn2, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(grid%conn2, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
      deallocate (grid%elem_connectivity)

   end subroutine CreateConnMatrix

!**************************************************************************

!**************************************************************************
   subroutine GridWritePFLOTRAN(grid, rank, size)
      !
      ! Dumps all the information into PFLOTRAN input file
      !
      ! Author: Zhuolin Qu, LANL
      ! Date: 08/10/2015
      !

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat

      implicit none

      type(grid_type) :: grid
      PetscBool :: keep_i_to_i_connections = .false.
      PetscInt :: ncols, connections_count, istart, iend
      PetscInt :: fileid, rank, irank, size
      PetscInt :: rec_tot_tri, rec_local_tri
      PetscInt :: io_rank = 0
      PetscInt :: i, j, k
      PetscInt, allocatable :: recvcounts(:)
      PetscReal, allocatable :: x(:), y(:), z(:)
      PetscInt :: temp_int, temp_int2, pos1, pos2, vsize
      PetscReal, pointer :: vec_ptr(:)
      PetscInt, allocatable :: sendcounts_pts(:), pos_pts(:)
      PetscInt, allocatable :: sendcounts_edge(:), pos_edge(:)
      PetscReal, allocatable :: vol_local(:), vol_global(:)
      PetscReal, allocatable :: edge_temp(:), edge_local(:)
      PetscInt, allocatable :: edge_global(:), degree_global(:)
      PetscReal, allocatable :: adj_temp_len(:), adj_local_len(:), adj_global_len(:)
      PetscReal, allocatable :: adj_temp_area(:), adj_local_area(:), adj_global_area(:)
      PetscInt, allocatable :: degree_local(:)
      PetscReal, allocatable :: data_local_real(:), data_global_real(:)
      PetscInt, allocatable :: num_pts_0(:), rec_local_tri_0(:)
      PetscInt, allocatable :: indx(:)
      PetscInt, allocatable :: sendcounts(:), pos(:)
      PetscReal, allocatable :: connectivity(:, :), connectivity_local(:), connectivity_tmp(:)
      PetscReal, allocatable :: connect_area(:, :), connect_map(:), connect_map_local(:), temp_connections(:)
      Vec :: mat_diag
      PetscErrorCode :: ierr

      ! get the vol_local
      allocate (vol_local(grid%num_pts_local))
      call VecCreate(PETSC_COMM_WORLD, mat_diag, ierr); CHKERRQ(ierr)
      call VecSetSizes(mat_diag, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
      call VecSetType(mat_diag, VECMPI, ierr); CHKERRQ(ierr)

      call MatGetDiagonal(grid%adjmatrix_area, mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)

      vol_local = vec_ptr

      call VecRestoreArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      call VecDestroy(mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(grid%degree, vec_ptr, ierr); CHKERRQ(ierr)

      rec_local_tri = int(sum(vec_ptr - 1))
      allocate (edge_local(rec_local_tri))
      allocate (adj_local_len(rec_local_tri))
      allocate (adj_local_area(rec_local_tri))

      temp_int = int(maxval(vec_ptr))
      allocate (edge_temp(temp_int + 10))
      allocate (adj_temp_area(temp_int))
      allocate (adj_temp_len(temp_int))

      pos1 = 0
      pos2 = 0

      ! ==================================================
      !  Prepare connectivity matrices
      ! ==================================================

      call VecGetOwnershipRange(grid%connect_map, istart, iend, ierr); CHKERRQ(ierr)
      allocate (connect_map_local(iend - istart))
      call VecGetValues(grid%connect_map, iend - istart, (/(i, i=istart, iend)/), connect_map_local, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) allocate (connect_map(grid%num_pts_global))

      call MPI_Gatherv(connect_map_local, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, connect_map, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         deallocate (pos)
         deallocate (recvcounts)
      endif

      if (rank == io_rank) vsize = int(maxval(connect_map))
      call MPI_Bcast(vsize, ONE_INTEGER_MPI, MPI_INTEGER, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      do i = 1, grid%num_pts_local
         if (grid%vertex_ids(i) /= grid%num_pts_global) then

            temp_int = int(vec_ptr(i))

            call MatGetRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)
            edge_local(pos1 + 1:pos1 + temp_int - 1) = int(edge_temp(2:temp_int))
            call MatRestoreRow(grid%edgematrix, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, edge_temp, ierr); CHKERRQ(ierr)

            pos1 = pos1 + temp_int - 1
            temp_int2 = int(vec_ptr(i) - 1)

            call MatGetRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_len, ierr); CHKERRQ(ierr)
            adj_local_len(pos2 + 1:pos2 + temp_int2) = adj_temp_len(1:temp_int2)
            call MatRestoreRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_len, ierr); CHKERRQ(ierr)

            call MatGetRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_area, ierr); CHKERRQ(ierr)
            adj_local_area(pos2 + 1:pos2 + temp_int2) = adj_temp_area(2:temp_int)
            call MatRestoreRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp_area, ierr); CHKERRQ(ierr)

            pos2 = pos2 + temp_int2
         endif
      enddo

      ! ====================================================
      !  Parse out connectivity data and area
      ! ====================================================
      call MatGetOwnershipRange(grid%connectivity, istart, iend, ierr); CHKERRQ(ierr)
      allocate (connectivity_local((iend - istart)*vsize))
      allocate (connect_area(grid%num_pts_global, vsize))
      allocate (temp_connections(vsize))

      do i = 1, iend - istart
         temp_connections = 0.
         call MatGetRow(grid%connectivity, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr); CHKERRQ(ierr)
         connectivity_local(((i - 1)*vsize) + 1:((i - 1)*vsize) + vsize) = temp_connections
         call MatRestoreRow(grid%connectivity, (i - 1) + istart, ncols, PETSC_NULL_INTEGER, temp_connections, ierr)
      enddo

      deallocate (temp_connections)

      if (rank == io_rank) then
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(iend - istart, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(istart, 1, MPI_INT, pos, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = pos*vsize
         recvcounts = recvcounts*vsize
         allocate (connectivity_tmp(grid%num_pts_global*vsize))
      endif

      call MPI_Gatherv(connectivity_local, (iend - istart)*vsize, &
                       MPI_DOUBLE_PRECISION, connectivity_tmp, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (connectivity(grid%num_pts_global, vsize))
         connectivity = reshape(connectivity_tmp, (/grid%num_pts_global, vsize/), order=(/2, 1/))

         deallocate (connectivity_tmp)
         deallocate (pos)
         deallocate (recvcounts)
      endif

      if (rank == io_rank) then
         allocate (sendcounts_pts(size))
         allocate (sendcounts_edge(size)) ! for the upper triangular case without diagonal element
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, sendcounts_pts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Gather(rec_local_tri, 1, MPI_INT, sendcounts_edge, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         allocate (pos_pts(size))
         allocate (pos_edge(size))
         pos_pts(1) = 0
         do irank = 1, size - 1
            pos_pts(irank + 1) = pos_pts(irank) + sendcounts_pts(irank)
         enddo
         pos_edge(1) = 0
         do irank = 1, size - 1
            pos_edge(irank + 1) = pos_edge(irank) + sendcounts_edge(irank)
         enddo
         rec_tot_tri = pos_edge(size) + sendcounts_edge(size)
      endif

      ! calculate degree_local
      allocate (degree_local(grid%num_pts_local))
      call MPI_Scatter(pos_edge, 1, MPI_INT, degree_local(1), 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      do i = 2, grid%num_pts_local
         degree_local(i) = degree_local(i - 1) + int(vec_ptr(i - 1) - 1)
      enddo

      degree_local = degree_local + 1

      call VecRestoreArrayReadF90(grid%degree, vec_ptr, ierr); CHKERRQ(ierr)

#if DEBUG
      print *, rank, 'vol_local=', vol_local
      print *, rank, 'degree_local', degree_local
      print *, rank, 'edge_lcoal=', edge_local
      print *, rank, 'adj_local_len=', adj_local_len
      print *, rank, 'adj_local_area=', adj_local_area
      print *, rank, 'rec_local_tri=', rec_local_tri
#endif

      allocate (data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 2))
      data_local_real(1:grid%num_pts_local) = vol_local
      data_local_real(grid%num_pts_local + 1:grid%num_pts_local*2) = degree_local
      data_local_real(grid%num_pts_local*2 + 1:grid%num_pts_local*2 + rec_local_tri) = edge_local

      call VecGetArrayReadF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)

      data_local_real(grid%num_pts_local*2 + rec_local_tri + 1:grid%num_pts_local*3 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 1), i=1, grid%num_pts_local)/)
      data_local_real(grid%num_pts_local*3 + rec_local_tri + 1:grid%num_pts_local*4 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 2), i=1, grid%num_pts_local)/)
      data_local_real(grid%num_pts_local*4 + rec_local_tri + 1:grid%num_pts_local*5 + rec_local_tri) = &
         (/(vec_ptr((i - 1)*3 + 3), i=1, grid%num_pts_local)/)

      call VecRestoreArrayReadF90(grid%coordinates_vec, vec_ptr, ierr); CHKERRQ(ierr)

      data_local_real(grid%num_pts_local*5 + rec_local_tri*1 + 1:grid%num_pts_local*5 + rec_local_tri*2) = adj_local_len
      data_local_real(grid%num_pts_local*5 + rec_local_tri*2 + 1:grid%num_pts_local*5 + rec_local_tri*3) = adj_local_area
      data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 1) = grid%num_pts_local
      data_local_real(grid%num_pts_local*5 + rec_local_tri*3 + 2) = rec_local_tri

      deallocate (vol_local)
      deallocate (degree_local)
      deallocate (edge_local)
      deallocate (adj_local_len)
      deallocate (adj_local_area)

      call VecDestroy(grid%coordinates_vec, ierr); CHKERRQ(ierr)

      ! prepare for gathering the data
      if (rank == io_rank) then
         allocate (sendcounts(size))
         allocate (pos(size))
         sendcounts = sendcounts_pts*5 + sendcounts_edge*3 + 2
         pos = pos_pts*5 + pos_edge*3 + (/(2*(i - 1), i=1, size)/)
         allocate (data_global_real(grid%num_pts_global*5 + rec_tot_tri*3 + 2*size))
         deallocate (sendcounts_pts)
         deallocate (pos_pts)
         deallocate (sendcounts_edge)
         deallocate (pos_edge)
      endif

      ! GAther the data to io_rank
      call MPI_Gatherv(data_local_real, grid%num_pts_local*5 + rec_local_tri*3 + 2, &
                       MPI_DOUBLE_PRECISION, data_global_real, sendcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      deallocate (data_local_real)

9140  format(I11, 2X, E20.12, 2X, E20.12, 2X, E20.12, 2X, E20.12)
9150  format(I11, 2X, I11, 2X, E20.12, 2X, E20.12, 2X, E20.12, 2X, E20.12)

      ! print to file
      if (rank == io_rank) then
         allocate (num_pts_0(size))
         allocate (rec_local_tri_0(size))
         allocate (indx(size))
         indx = (/pos(2:size) - 1, pos(size) + sendcounts(size) - 1/)
         num_pts_0 = int(data_global_real(indx))
         rec_local_tri_0 = int(data_global_real(indx + 1))
         ! global volumes
         indx = pos + 1
         deallocate (sendcounts)
         deallocate (pos)
         vol_global = (/(dabs(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)
         allocate (x(grid%num_pts_global))
         allocate (y(grid%num_pts_global))
         allocate (z(grid%num_pts_global))
         allocate (adj_global_len(rec_tot_tri))
         allocate (adj_global_area(rec_tot_tri))
         allocate (degree_global(grid%num_pts_global))
         indx = indx + num_pts_0
         degree_global = (/(int(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i))), i=1, size)/)
         indx = indx + num_pts_0
         edge_global = (/(int(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)
         indx = indx + rec_local_tri_0
         x = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         y = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         z = (/(data_global_real(indx(i):indx(i) - 1 + num_pts_0(i)), i=1, size)/)
         indx = indx + num_pts_0
         adj_global_len = (/(dabs(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)
         indx = indx + rec_local_tri_0
         adj_global_area = (/(dabs(data_global_real(indx(i):indx(i) - 1 + rec_local_tri_0(i))), i=1, size)/)

         deallocate (data_global_real)
         deallocate (num_pts_0)
         deallocate (rec_local_tri_0)
         deallocate (indx)
         fileid = 49
         open (fileid, file=trim(grid%dump_str))

         ! write id_N,x_N,y_N,z_N,volume_N
         write (fileid, '(A11,2X,I11)') 'CELLS', grid%num_pts_global
         do i = 1, grid%num_pts_global
            write (fileid, 9140) i, x(i), y(i), z(i), vol_global(i)
         enddo

         ! Determine number of connections
         ! Inelegant and should be deprecated for a faster solution
         connections_count = int(sum(connect_map))
         if (keep_i_to_i_connections .eqv. .false.) connections_count = connections_count - grid%num_pts_global

         ! write id_up_N,id_dn_N,x_M,y_M,z_M,area_M
         !write(fileid,*)
         write (fileid, '(A11,2X,I11)') 'CONNECTIONS', connections_count
         pos1 = 1
         do i = 1, grid%num_pts_global!-1
            do k = 2, int(connect_map(i))
               !do k=degree_global(i),degree_global(i+1)-1
               !j=edge_global(k)
               j = int(connectivity(i, k))
               if ((keep_i_to_i_connections .eqv. .false.) .AND. (j .eq. i)) cycle
               write (fileid, 9150) i, j, (x(i) + x(j))*0.5d0, (y(i) + y(j))*0.5d0, (z(i) + z(j))*0.5d0, adj_global_area(pos1)
               pos1 = pos1 + 1
            enddo
         enddo

         deallocate (vol_global)
         deallocate (edge_global)
         deallocate (x)
         deallocate (y)
         deallocate (z)
         deallocate (adj_global_len)
         deallocate (adj_global_area)
         close (fileid)
      endif

   end subroutine GridWritePFLOTRAN

!===========================================================

!   subroutine orderPointsInPolygon(c, v1, v2, v3, v4)
!      !
!      ! Given four points, sort them in a clockwise-ordering
!      ! according to the RHR.
!      !
!      implicit none
!      logical :: mask(4)
!      PetscInt :: c(4), l(1), i
!      PetscReal :: v1(3), v2(3), v3(3), v4(3)
!      PetscReal :: angles(4), centroid(3)

!      mask = .true.
!      centroid = 0.25*(v1 + v2 + v3 + v4)

!      angles = (/ATAN2(v1(2) - centroid(2), v1(1) - centroid(1)), &
!                 ATAN2(v2(2) - centroid(2), v2(1) - centroid(1)), &
!                 ATAN2(v3(2) - centroid(2), v3(1) - centroid(1)), &
!                 ATAN2(v4(2) - centroid(2), v4(1) - centroid(1))/)

!      do i = 1, 4
!         l = MINLOC(angles, mask)
!         c(i) = l(1)
!         mask(l(1)) = .false.
!      enddo

!   end subroutine orderPointsInPolygon

!===========================================================

   subroutine Diagnostics(grid, atts, rank, size)
      ! minimal hides top banner

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

      use petscvec
      use petscmat
      implicit none

      type(grid_type) :: grid
      type(diag_atts) :: atts
      PetscInt :: rank, size

      PetscReal :: tmp, area_local, face_area, vol
      PetscReal :: area_min, area_max, len_min, len_max, area_len_min, area_len_max
      PetscInt :: mint, temp_int, rec_local, negative_coeffs_local, negative_coeffs
      PetscBool :: hide_banner = PETSC_TRUE

      character(len=50) :: l
      character(len=1) :: h1, h2
      character(len=20) :: element_type
      PetscInt :: w1, w2, i, j, k, pos1, io_rank = 0
      PetscErrorCode :: ierr

      PetscInt :: row_idx, ncols, max_connections, min_connections, max_connections_local, min_connections_local
      PetscInt, allocatable :: recvcounts(:), pos(:)
      PetscReal, pointer :: vec_ptr(:)
      PetscReal, allocatable :: vor_vol(:)
      PetscReal, allocatable :: adj_temp(:), adj_local(:), vol_local(:)
      Vec :: mat_diag

      !allocate(adj_temp(grid%num_pts_global))
      !allocate(adj_local(grid%num_pts_global))
      !allocate(vor_vol(grid%num_pts_local))

      ! Get the node volumes
      ! get the vol_local
      allocate (vol_local(grid%num_pts_local))
      call VecCreate(PETSC_COMM_WORLD, mat_diag, ierr); CHKERRQ(ierr)
      call VecSetSizes(mat_diag, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
      call VecSetType(mat_diag, VECMPI, ierr); CHKERRQ(ierr)
      call MatGetDiagonal(grid%adjmatrix, mat_diag, ierr); CHKERRQ(ierr)
      call VecGetArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)

      vol_local = vec_ptr

      call VecRestoreArrayReadF90(mat_diag, vec_ptr, ierr); CHKERRQ(ierr)
      call VecDestroy(mat_diag, ierr); CHKERRQ(ierr)

      ! Get edge and adj local
      call VecGetArrayReadF90(grid%degree_tot, vec_ptr, ierr); CHKERRQ(ierr)

      rec_local = int(sum(vec_ptr))
      allocate (adj_local(rec_local))

      temp_int = int(maxval(vec_ptr))
      allocate (adj_temp(temp_int))

      !call VecRestoreArrayReadF90(grid%degree_tot,vec_ptr,ierr);CHKERRQ(ierr)

      ! Get the area/len coefficients
      row_idx = 1 ! depreciated
      j = 1 ! neg_vols counter
      k = 1 ! neg_coeff counter

      area_max = -100000.0
      len_max = area_max
      area_len_max = area_max

      area_min = 100000.0
      len_min = area_min
      area_len_min = area_min

      ! ==================================
      ! CAPTURE VOLUMES
      ! ==================================

      if (rank == io_rank) then
         allocate (vor_vol(grid%num_pts_global))
         allocate (recvcounts(size))
         allocate (pos(size))
      endif

      call MPI_Gather(grid%num_pts_local, 1, MPI_INT, recvcounts, 1, MPI_INT, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         pos = 0
         if (size > 1) pos(2:size) = (/(sum(recvcounts(1:i)), i=1, size)/)
      endif

      call MPI_Gatherv(vol_local, grid%num_pts_local, &
                       MPI_DOUBLE_PRECISION, vor_vol, recvcounts, &
                       pos, MPI_DOUBLE_PRECISION, io_rank, &
                       MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      if (rank == io_rank) then
         deallocate (recvcounts)
         deallocate (pos)
      endif

      call MPI_Reduce(atts%vol, vol, 1, MPI_DOUBLE, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! ==================================
      ! AREA / LENGTH STATISTICS
      ! ==================================
      max_connections_local = -1
      min_connections_local = 1e4
      pos1 = 0
      do i = 1, grid%num_pts_local
         temp_int = int(vec_ptr(i))

         ! Capture coefficients
         call MatGetRow(grid%adjmatrix_full, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         adj_local(pos1 + 1:pos1 + temp_int) = adj_temp(1:temp_int)
         max_connections_local = max(max_connections_local, count(abs(adj_temp(1:temp_int)) .gt. 1e-12))
         min_connections_local = min(min_connections_local, count(abs(adj_temp(1:temp_int)) .gt. 1e-12))

         negative_coeffs_local = negative_coeffs_local + count(adj_temp(1:temp_int) .lt. 0.0)

         call MatRestoreRow(grid%adjmatrix_full, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         pos1 = pos1 + temp_int
      enddo

      call MPI_Reduce(minval(adj_local), area_len_min, 1, MPI_DOUBLE, MPI_MIN, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(maxval(adj_local), area_len_max, 1, MPI_DOUBLE, MPI_MAX, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(max_connections_local, max_connections, 1, MPI_INT, MPI_MAX, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(min_connections_local, min_connections, 1, MPI_INT, MPI_MIN, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(negative_coeffs_local, negative_coeffs, 1, MPI_INT, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! ==================================
      ! AREA STATISTICS
      ! ==================================
      pos1 = 0
      adj_local = 0.
      do i = 1, grid%num_pts_local
         temp_int = int(vec_ptr(i))

         ! Capture coefficients
         call MatGetRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         adj_local(pos1 + 1:pos1 + temp_int) = adj_temp(1:temp_int)
         call MatRestoreRow(grid%adjmatrix_area, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         pos1 = pos1 + temp_int
         area_local = area_local + sum(adj_temp(1:temp_int))

      enddo

      call MPI_Reduce(area_local, atts%vor_area, 1, MPI_DOUBLE, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(atts%face_area, face_area, 1, MPI_DOUBLE, MPI_SUM, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(minval(adj_local), area_min, 1, MPI_DOUBLE, MPI_MIN, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(maxval(adj_local), area_max, 1, MPI_DOUBLE, MPI_MAX, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      ! ==================================
      ! LENGTH STATISTICS
      ! ==================================
      pos1 = 0
      adj_local = 0.
      do i = 1, grid%num_pts_local
         temp_int = int(vec_ptr(i))

         ! Capture coefficients
         call MatGetRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         adj_local(pos1 + 1:pos1 + temp_int) = adj_temp(1:temp_int)
         call MatRestoreRow(grid%adjmatrix_len, grid%vertex_ids(i) - 1, ncols, PETSC_NULL_INTEGER, adj_temp, ierr); CHKERRQ(ierr)
         pos1 = pos1 + temp_int
      enddo

      call MPI_Reduce(minval(adj_local), len_min, 1, MPI_DOUBLE, MPI_MIN, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
      call MPI_Reduce(maxval(adj_local), len_max, 1, MPI_DOUBLE, MPI_MAX, io_rank, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)

      !call EXIT(0)

      ! Iterate over every sparse matrix row
      ! do i=1,grid%num_pts_local
      !
      !   ! ==================================
      !   ! AREA / LENGTH STATISTICS
      !   ! ==================================

      !   call MatGetRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)
      !   adj_local(1:grid%num_pts_global) = adj_temp(1:grid%num_pts_global)
      !   call MatRestoreRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)

      !   do ii=2,grid%num_pts_global
      !     temp_real = adj_local(ii)

      !     area_len_min = min(area_len_min,temp_real)
      !     area_len_max = max(area_len_max,temp_real)

      !     if ((temp_real < -1e-8) .AND. (k <= 10)) then
      !       neg_coeff_val(k) = -temp_real     ! Value
      !       neg_coeff_idx(k,1) = i            ! Column
      !       neg_coeff_idx(k,2) = ii           ! Row
      !       k = k + 1
      !     endif
      !   enddo

      !   ! ==================================
      !   ! AREA STATISTICS
      !   ! ==================================

      !   call MatGetRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)
      !   adj_local(1:grid%num_pts_global) = adj_temp(1:grid%num_pts_global)
      !   call MatRestoreRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)

      !   area_min = min(area_min,minval(adj_local(i+1:)))
      !   area_max = max(area_max,maxval(adj_local(i+1:)))

      !   atts%vor_area = atts%vor_area + sum(adj_local(i+1:))

      !   ! ==================================
      !   ! LENGTH STATISTICS
      !   ! ==================================

      !   call MatGetRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)
      !   adj_local(1:grid%num_pts_global) = adj_temp(1:grid%num_pts_global)
      !   call MatRestoreRow(grid%adjmatrix_full,i-1,ncols,PETSC_NULL_INTEGER,adj_temp,ierr);CHKERRQ(ierr)

      !   len_min = min(len_min,minval(adj_local(i+1:)))
      !   len_max = max(len_max,maxval(adj_local(i+1:)))

      !   ! ==================================
      !   ! OTHER STATISTICS
      !   ! ==================================

      !   max_node_connections = max(max_node_connections,count(abs(adj_local(i+1:)) .gt. 1e-9))

      !   ! Capture negative Voronoi volumes
      !   if ((vor_vol(i) < -1e-8) .AND. (j <= 10)) then
      !     neg_vols_val(j) = vor_vol(i)
      !     neg_vols_idx(j) = i
      !     j = j + 1
      !   endif

      ! enddo

      if (rank == io_rank) then
         if (grid%ndim == 3) then
            element_type = 'tetrahedral'
         else
            element_type = 'triangle'
         endif

         ! Necessary to 'left-flush' lines
         l = repeat(' ', 50)

         h1 = '='    ! Header 1 linebreak
         h2 = '-'    ! Header 2 linebreak

         w1 = 65     ! Width of header 1
         w2 = w1 - 3 ! Width of header 2

         tmp = -0.0000456278
         mint = 12

         ! Text: float
44       format(4X, A46, 1x, E15.7)
         ! Text: integer
46       format(4X, A51, 1x, i10)
         ! Text
47       format(4X, A)

         ! Print out header
         if (hide_banner .eqv. PETSC_FALSE) then
            print *, char(10), ' ', repeat('-', (w1 - 15)/2), '[ DIAGNOSTICS ]', repeat('-', (w1 - 15)/2)
         endif

         ! Begin outputting diagnostics
         print *, char(10), ' CONNECTIONS'
         print *, repeat(h1, w1)
         print 46, 'Minimum Cell Connections:'//l, min_connections
         print 46, 'Maximum Cell Connections:'//l, max_connections

         print *, char(10), ' VOLUME'
         print *, repeat(h1, w1)
         print 44, 'Min. Voronoi volume:'//l, minval(vor_vol)
         print 44, 'Max. Voronoi volume:'//l, maxval(vor_vol)
         print 47, repeat(h2, w2)
         print 44, 'Total Voronoi volume:'//l, sum(vor_vol)
         print 44, 'Total '//trim(element_type)//' volume:'//l, vol!atts%vol

         print *, char(10), ' COEFFICIENTS'
         print *, repeat(h1, w1)
         print 44, 'Min. Voronoi area:'//l, area_min
         print 44, 'Max. Voronoi area:'//l, area_max
         print 47, repeat(h2, w2)
         print 44, 'Min. Voronoi edge length:'//l, len_min
         print 44, 'Max. Voronoi edge length:'//l, len_max
         print 47, repeat(h2, w2)
         print 44, 'Min. area/len coefficient:'//l, area_len_min
         print 44, 'Max. area/len coefficient:'//l, area_len_max
         print 47, repeat(h2, w2)
         print 44, 'Total Voronoi area:'//l, atts%vor_area
         print 44, 'Total triangle face area:'//l, face_area

         print *, char(10), ' WARNINGS'
         print *, repeat(h1, w1)
         print 46, 'Negative Voronoi volumes:'//l, count(vor_vol .lt. 0.0)!j-1

         ! Write out where negative volumes are
         !if (j > 1) then
         !  print*,''
         !  ! k = ceil(log10(y+1))
         !  do i=1,atts%neg_vol
!
         !    if (i > 10) then
         !      print*,repeat(' ',(w1-1)/2),'...'
         !      exit
         !    endif
!
         !    write(tmp_str,'(i10)') neg_vols_idx(i)
         !    print 44,'    - Node '//trim(tmp_str)//':'//l,neg_vols_val(i)
         !  enddo
         !  print*,''
         !endif

         print 47, repeat(h2, w2)
         print 46, 'Negative area/len coefficients:'//l, negative_coeffs!k-1

         ! Write out where negative coefficients are
         !if (k > 1) then
         !  print*,''
         !  do i=1,atts%neg_coeff
!
         !    if (i > 10) then
         !      print*,repeat(' ',(w1-1)/2),'...'
         !      exit
         !    endif
         !
         !    write(tmp_str,'(i8)') neg_coeff_idx(i,1)
         !    write(tmp_str2,'(i8)') neg_coeff_idx(i,2)
         !    print 44,'    - Row '//trim(tmp_str)//', column '//tmp_str2//':'//l,neg_coeff_val(i)
         !  enddo
         !  print*,''
         !endif

         print 47, repeat(h2, w2)
      endif

   end subroutine Diagnostics

end module Grid_module
