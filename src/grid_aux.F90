module Grid_Aux_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
  use petscvec
  use petscmat
  use petscis

  implicit none

  private 
  
  ! Contains attributes about mesh for diagnostics
  type, public :: diag_atts
    PetscBool :: are_on
    PetscBool :: coefficients_are_compressed
    PetscBool :: coefficients_are_dedudded

    PetscReal :: min_area
    PetscReal :: max_area
    
    PetscReal :: min_edge
    PetscReal :: max_edge

    PetscReal :: vor_area
    PetscReal :: face_area

    PetscInt :: neg_vol
    PetscInt :: neg_coeff

    PetscReal :: vol
  end type diag_atts

  type, public :: point3d_type
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point3d_type
  
  type, public :: grid_type
    character(len=150) :: avs_str, lg_str, out_str, zone_str, dump_str ! argument passing strings
    character(len=150) :: aperture_file
    character(len=150) :: outfile                                      ! stor file output
    character*132 :: aperture_attribute
    character(len=5), dimension(:), allocatable :: imt_values

    PetscBool :: lg_flag, avs_flag, out_flag, zone_flag, dump_flag, verbose_flag
    PetscBool :: adjust_aperture
    PetscBool :: is_tough
    PetscInt :: cv_type ! control volume type
    PetscInt :: verbose, outtype
    PetscInt, pointer :: vertex_ids(:)
    PetscInt, pointer :: vertex_ids_need(:)
    PetscInt, pointer :: vertex_ids_map(:)
    PetscInt, pointer :: elem_ids(:)
    PetscInt, pointer :: elem_connectivity(:,:)

    PetscReal, pointer :: apertures(:)
    PetscReal, pointer :: attribute(:)
    PetscReal, pointer :: cell_cc(:,:) ! PFLOTRAN
    PetscReal, pointer :: cell_vol(:)  ! PFLOTRAN

    PetscInt :: num_elems_global  ! Number of cells in the entire domain
    PetscInt :: num_elems_local   ! Number of elements locally
    PetscInt :: num_elems_max     ! Max of all the number of elements locally
    PetscInt :: num_pts_local
    PetscInt :: num_pts_need
    PetscInt :: num_pts_global
    PetscInt :: ndim         ! The dimensionality (2 or 3) of the mesh elements
    Vec :: lcc               ! (parallel vector)  Voronoi cell centers
    Vec :: coordinates_local ! (series vector)    local stored coordinates depending on the connectivity on current processor
    Vec :: coordinates_vec   ! (parallel vector)  store all the coordinates in a 3*num_pts long vector
    Vec :: connections       ! (parallel vector)  rough estimate of the connections 
    Vec :: degree            ! (upper triangular/parallel vector) true degree of connections for vertices
    Vec :: degree_tot        ! (full matrix/parallel vector) exact degree of connections for vertices 
    Vec :: connect_map       ! TOUGH2 (perhaps PFLOTRAN) : get connectivity matr
    Mat :: connectivity      ! TOUGH2 (perhaps PFLOTRAN) : get connectivity matr
    Mat :: connect_area      ! TOUGH2 (perhaps PFLOTRAN) : get connectivity matr
    Mat :: adjmatrix         ! (upper triangular) sparse parallel matrix: ajacent matrix that contains the geometric coefficients
    Mat :: adjmatrix_len     ! TOUGH2
    Mat :: adjmatrix_area    ! TOUGH2
    Mat :: adjmatrix_full    ! (full matrix/parallel matrix) ajacent matrix that contains the geometric coefficients
    Mat :: edgematrix        ! (full matrix/parallel matrix) records the edge list

    Mat :: conn1             ! PFLOTRAN
    Mat :: conn2             ! PFLOTRAN

  end type grid_type

  public :: GridCreate, &
            DiagCreate, &
            GridDestroy

contains

! ************************************************************************** !

function GridCreate()
  ! 
  ! Creates an grid object
  ! 
  ! Author: Satish Karra & Zhuolin Qu, LANL
  ! Date: 07/21/2015
  ! 

  implicit none
  
  type(grid_type), pointer :: GridCreate
  type(grid_type), pointer :: grid

  allocate(grid)

  nullify(grid%vertex_ids)
  nullify(grid%vertex_ids_need)
  nullify(grid%vertex_ids_map)
  nullify(grid%elem_ids)
  nullify(grid%elem_connectivity)

  nullify(grid%cell_cc)  ! PFLOTRAN
  nullify(grid%cell_vol) ! PFLOTRAN

  !grid%aperture_attribute = " "
  grid%is_tough = PETSC_FALSE
  grid%outtype = 0
  grid%num_elems_global = 0 
  grid%num_elems_local = 0 
  grid%num_elems_max = 0 
  grid%connections = PETSC_NULL_VEC
  grid%degree = PETSC_NULL_VEC
  grid%degree_tot = PETSC_NULL_VEC
  grid%connect_map = PETSC_NULL_VEC
  grid%adjmatrix = PETSC_NULL_MAT
  grid%adjmatrix_full = PETSC_NULL_MAT

  ! TOUGH2
  grid%adjmatrix_area = PETSC_NULL_MAT
  grid%adjmatrix_len = PETSC_NULL_MAT
  ! TOUGH2 END
  
  grid%edgematrix = PETSC_NULL_MAT

  GridCreate => grid
  
end function GridCreate

! ************************************************************************** !

function DiagCreate()

  implicit none

  type(diag_atts), pointer :: DiagCreate
  type(diag_atts), pointer :: atts

  allocate(atts)

  atts%min_area = 99999.0
  atts%max_area = -99999.0
  atts%min_edge = 99999.0
  atts%max_edge = -99999.0
  atts%vor_area = 0.0
  atts%face_area = 0.0
  atts%vol = 0.0
  atts%neg_vol = 0
  atts%neg_coeff = 0
  atts%are_on = PETSC_TRUE
  atts%coefficients_are_dedudded = PETSC_FALSE
  atts%coefficients_are_compressed = PETSC_FALSE

  DiagCreate => atts

end function DiagCreate

! ************************************************************************** !

subroutine TestMethod(grid)
  implicit none
  
  type(grid_type), pointer :: grid

  print*,'hello world'

end subroutine TestMethod


! ************************************************************************** !

subroutine GridDestroy(grid)
  ! 
  ! Destroys an grid object
  ! 
  ! Author: Satish Karra & Zhuolin Qu, LANL
  ! Date: 07/21/2015
  ! 

  implicit none
  
  type(grid_type), pointer :: grid
  PetscErrorCode :: ierr

  if (grid%connections /= PETSC_NULL_VEC) then
    call VecDestroy(grid%connections,ierr);CHKERRQ(ierr)
  endif
  if (grid%adjmatrix /= PETSC_NULL_MAT) then
    call MatDestroy(grid%adjmatrix,ierr);CHKERRQ(ierr)
  endif

  ! BEGIN TOUGH2
  if (grid%adjmatrix_len /= PETSC_NULL_MAT) then
    call MatDestroy(grid%adjmatrix_len,ierr);CHKERRQ(ierr)
  endif
  if (grid%adjmatrix_area /= PETSC_NULL_MAT) then
    call MatDestroy(grid%adjmatrix_area,ierr);CHKERRQ(ierr)
  endif
  ! END TOUGH2

  if (grid%adjmatrix_full /= PETSC_NULL_MAT) then
    call MatDestroy(grid%adjmatrix_full,ierr);CHKERRQ(ierr)
  endif
  if (grid%edgematrix /= PETSC_NULL_MAT) then
    call MatDestroy(grid%edgematrix,ierr);CHKERRQ(ierr)
  endif
  if (grid%coordinates_local /= PETSC_NULL_VEC) then
    call VecDestroy(grid%coordinates_local,ierr);CHKERRQ(ierr)
  endif
  if (grid%coordinates_vec /=PETSC_NULL_VEC) then
    call VecDestroy(grid%coordinates_vec,ierr);CHKERRQ(ierr)
  endif
  if (grid%degree /= PETSC_NULL_VEC) then
    call VecDestroy(grid%degree,ierr);CHKERRQ(ierr)
  endif
  if (grid%degree_tot /= PETSC_NULL_VEC) then
    call VecDestroy(grid%degree_tot,ierr);CHKERRQ(ierr)
  endif 

  if (.not.associated(grid)) return
  deallocate(grid%vertex_ids)
  deallocate(grid%elem_ids)

end subroutine GridDestroy

end module Grid_Aux_module
