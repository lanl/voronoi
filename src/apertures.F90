module VoronoiApertures

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

   use Grid_Aux_module
   use Voronoi_Constants_module
   use VoronoiGridCore, only: GetVolumes

    implicit none

    private

    public :: &
        VoronoiApplyApertures, &
        VoronoiReadApertures

contains

    subroutine VoronoiApplyApertures(grid, apertures, attribute, rank, size)

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

        use petscvec
        use petscmat
        implicit none

        type(grid_type) :: grid
        character(len=*) :: attribute
        Vec :: apertures
        PetscInt :: rank, size

        PetscReal, allocatable :: volumes(:)

        volumes = GetVolumes(grid, rank, size)
        ! xyz = GetCellCenters()
        ! edges = GetEdgeList()
        ! areas = GetAreas()
        ! lengths = GetLengths()

        ! 0. take grid and...
        ! 1. take attribute vector and...
        ! 1. apply to areas
        ! 2. apply to volumes

        if (allocated(volumes)) then
            deallocate(volumes)
        endif

    end subroutine

    subroutine VoronoiReadApertures(filename, vec_apertures)

#include "petsc/finclude/petscvec.h"
        use petscvec
        implicit none

        character(len=*) :: filename
        character(len=100) :: tmp_string, att_name
        Vec :: vec_apertures
        PetscInt :: rank, io_rank
        PetscInt :: i, file_unit, index, ierr
        PetscInt :: num_apertures
        PetscBool :: file_exists
        PetscReal :: value
        PetscReal, allocatable :: apertures(:)

        call VecCreate(PETSC_COMM_WORLD, vec_apertures, ierr); CHKERRQ(ierr)
        call VecSetType(vec_apertures, VECMPI, ierr); CHKERRQ(ierr)

        if (rank == io_rank) then
            file_unit = 15
            inquire(file=filename, exist=file_exists)

            if (file_exists .eqv. .false.) then
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, 'File ' // trim(filename) // ' does not exist')
                call PETSCABORT(PETSC_COMM_WORLD, ierr); CHKERRA(ierr)
            endif

            ! Read over open file until EOF
            ! Fill in zdepths using the first column as
            ! an index
            open(unit=file_unit, file=filename)
            do
                read(file_unit,'(A)',IOSTAT=ierr) tmp_string
                if (ierr > 0) then
                    ! Malformed file format / broken stream / etc.
                    close(file_unit)
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, 'Something went wrong during file read')
                    call PETSCABORT(PETSC_COMM_WORLD, ierr); CHKERRA(ierr)
                else if (ierr < 0) then
                    ! EOF
                    exit
                endif

                ! Ignore comments and empty lines
                if ((tmp_string(1:1) == '#') .OR. (trim(tmp_string(1:1)) == '')) then
                    cycle
                endif

                ! Read in the attribute name from control file of form:
                ! attribute: imt1
                if (tmp_string(1:10) == 'attribute:') then
                    att_name = trim(tmp_string(11:len(tmp_string)))
                    cycle
                endif

                ! Read in the attribute name from control file of form:
                ! apertures: 4
                if (tmp_string(1:10) == 'apertures:') then
                    tmp_string = trim(tmp_string(11:))
                    read(tmp_string, *) index
                    allocate(apertures(index))
                    cycle
                endif

                ! Read in attribute_id, aperture_val
                read(tmp_string,*) index, value
                apertures(index) = value
            enddo
            close(file_unit)

            call VecAssemblyBegin(vec_apertures, ierr); CHKERRQ(ierr)

            call VecSetSizes( &
                vec_apertures, &
                num_apertures, &
                num_apertures, &
                ierr &
            ); CHKERRQ(ierr)

            call VecSetValues( &
                vec_apertures, &
                num_apertures, &
                (/(i, i=0, num_apertures-1)/), &
                apertures, &
                INSERT_VALUES, &
                ierr &
            ); CHKERRQ(ierr)

            deallocate(apertures)

            call VecAssemblyEnd(vec_apertures, ierr); CHKERRQ(ierr)
        endif

        call PetscBarrier(vec_apertures, ierr); CHKERRQ(ierr)
   end subroutine

end module VoronoiApertures