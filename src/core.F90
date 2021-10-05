module VoronoiGridCore

   use Grid_Aux_module
   use Voronoi_Constants_module
   implicit none

private

public :: &
    GetVolumes

contains

    function GetVolumes(grid, rank, size) result(volumes)
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

        use petscvec
        use petscmat
        implicit none

        type(grid_type) :: grid
        Vec :: vec_local, vec_global
        VecScatter :: ctx
        PetscInt :: i, rank, size
        PetscInt :: io_rank = 0
        PetscReal, allocatable :: volumes(:)
        PetscReal, pointer :: ptr(:)
        PetscErrorCode :: ierr

        ! Create the PetscVec that will store volumes
        ! Recall that volumes are stored on the matrix diagonal
        call VecCreate(PETSC_COMM_WORLD, vec_local, ierr); CHKERRQ(ierr)
        call VecSetSizes(vec_local, grid%num_pts_local, grid%num_pts_global, ierr); CHKERRQ(ierr)
        call VecSetType(vec_local, VECMPI, ierr); CHKERRQ(ierr)
        call MatGetDiagonal(grid%adjmatrix_area, vec_local, ierr); CHKERRQ(ierr)

        ! Gather all data onto a single rank (rank = 0)
        call VecScatterCreateToZero(vec_local, ctx, vec_global, ierr); CHKERRQ(ierr)
        call VecScatterBegin(vec_local, vec_global, INSERT_VALUES, SCATTER_FORWARD, ctx, ierr); CHKERRQ(ierr)
        call VecScatterEnd(vec_local, vec_global, INSERT_VALUES, SCATTER_FORWARD, ctx, ierr); CHKERRQ(ierr)
        call VecGetArrayF90(vec_global, ptr, ierr); CHKERRQ(ierr)

        ! Move all data to a single rank
        if (rank == 0) then
            allocate(volumes(grid%num_pts_global))

            ! TODO: is copying the most effective/efficient way
            ! of working around the destructors?
            do i = 1, grid%num_pts_global
                volumes(i) = ptr(i)
            enddo
        endif

        ! Invoke destructors for all created contexts
        ! Only the standard Fortran array (volumes) should remain
        call VecRestoreArrayF90(vec_global, ptr, ierr); CHKERRQ(ierr)
        call VecScatterDestroy(ctx, ierr); CHKERRQ(ierr)
        call VecDestroy(vec_global, ierr); CHKERRQ(ierr)
        call VecRestoreArrayReadF90(vec_local, ptr, ierr); CHKERRQ(ierr)
        call VecDestroy(vec_local, ierr); CHKERRQ(ierr)
    end function

end module VoronoiGridCore