!=====================================
! Main driver for voronoi calculations
! from given delauney mesh
! All done in parallel
! Authors: Zhuolin Qu and Satish Karra
! Addl. coding by Daniel Livingston
! Los Alamos National Laboratory
! Started Date: July 6, 2015
! Update: October 2017
!=====================================

program voronoi

#include "petsc/finclude/petscsys.h"

   use petscsys
   use Grid_Aux_module
   use Grid_module

   implicit none
   character(len=*), parameter :: version = '1.0.0'
   character(len=*), parameter :: date    = 'Oct. 2021'
   character(len=10) :: cv_str
   character(len=150) :: tmp_str

   ! grid_type defined in grid_aux.F90
   type(grid_type), pointer :: grid
   type(diag_atts), pointer :: atts
   PetscErrorCode :: ierr
   PetscBool :: flg, run_tests, verbosity_bool, help_bool
   PetscBool :: cv_flag
   PetscMPIInt :: rank, size
   PetscInt :: io_rank = 0
   PetscReal :: time_start, time_finish
   PetscInt :: out_type, verbosity

   common verbosity, verbosity_bool

   verbosity = 2 ! Set default verbosity level

   !---------------- Begin main block --------------------------
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
   call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
   call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr); CHKERRQ(ierr)

   grid => GridCreate() ! Creates a grid object (grid_aux.F90)
   atts => DiagCreate() ! Creates a diagnostic attributes object

   !----------! Parse command line for flags and input !----------!
   ! Get input type - avs or LaGriT?
   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-avs', &
                              grid%avs_str, grid%avs_flag, ierr); CHKERRQ(ierr)

   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-lg', &
                              grid%lg_str, grid%lg_flag, ierr); CHKERRQ(ierr)

   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-zone', &
                              grid%zone_str, grid%zone_flag, ierr); CHKERRQ(ierr)

   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-o', &
                              grid%dump_str, grid%dump_flag, ierr); CHKERRQ(ierr)

   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-cv', &
                              cv_str, cv_flag, ierr); CHKERRQ(ierr)

   ! Get verbosity level
   call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-v', &
                           verbosity, verbosity_bool, ierr); CHKERRQ(ierr)

   ! Get diagnostics bool
   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-d', flg, &
                            atts%are_on, ierr); CHKERRQ(ierr)

   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-compress', flg, &
                            atts%coefficients_are_compressed, ierr); CHKERRQ(ierr)

   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-dedud', flg, &
                            atts%coefficients_are_dedudded, ierr); CHKERRQ(ierr)

   ! Run tests?
   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-t', flg, &
                            run_tests, ierr); CHKERRQ(ierr)

   ! Show 'man pages'?
   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-help', flg, &
                            help_bool, ierr); CHKERRQ(ierr)

   ! Capture desired output filetype
   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-type', &
                              grid%out_str, grid%out_flag, ierr); CHKERRQ(ierr)

   if (rank == io_rank) call banner()

   ! Launch the test suite and quit this instance of the program
   if ((run_tests .EQV. PETSC_TRUE) .AND. (rank == io_rank)) then
      call execute_command_line("python test/sparse_matrix_compare.py", wait=.false.)
      stop ": PLEASE WAIT: Loading test suite..."
   endif

   if (rank == io_rank) then
      if (atts%coefficients_are_dedudded .EQV. PETSC_TRUE) then
         print *, 'ERROR: Coefficient removal is not yet MPI-compatible.'
         call EXIT(1)
      endif
   endif

   ! Configure the control volume type
   ! Setting this as a int prevents a more
   ! expensive string bitwise comparison in loops
   if (cv_flag .EQV. PETSC_TRUE) then
      select case (cv_str(:7))
      case ('voronoi')
         grid%cv_type = 1
      case ('median ')
         grid%cv_type = 2
      case ('hybrid ')
         grid%cv_type = 3
      case default
         grid%cv_type = 1
         print *, 'Unable to parse ', cv_str, '\n', &
            'Defaulting to Voronoi control volume.'
      end select
   else
      grid%cv_type = 1
   endif

   ! The output type is stored in out_flag.
   ! Parse this flag to determine the output type.
   if (grid%out_flag .EQV. PETSC_TRUE) then
      select case (grid%out_str(:4))
      case ('pflo')
         out_type = 3
         grid%outtype = 3
         grid%is_tough = PETSC_TRUE
         tmp_str = 'pflotran.uge'
      case ('fehm')
         out_type = 1
         grid%outtype = 1
         tmp_str = 'voronoi.stor'
      case ('toug')
         out_type = 2
         grid%is_tough = PETSC_TRUE
         grid%outtype = 2
         tmp_str = 'VORONOI_MESH'
      case ('hdf5')
         out_type = 4
         grid%outtype = 4
         tmp_str = 'voronoi.h5'
      case default
         if (rank == io_rank) print *, 'UNKNOWN OUTPUT FLAG: DEFAULTING TO FEHM'
         out_type = 1
         grid%outtype = 1
         tmp_str = 'voronoi.stor'
      end select
   else
      out_type = 1
      grid%outtype = 1
      tmp_str = 'voronoi.stor'
   endif

   if (help_bool .EQV. PETSC_TRUE) then
      if (rank == io_rank) call print_help()
   endif

   ! Check if no flags were passed
   if ((grid%lg_flag .eqv. PETSC_FALSE) .AND. (grid%avs_flag .eqv. PETSC_FALSE)) then
      if (rank == io_rank) call short_help()
      return
   endif

   if (grid%dump_flag .eqv. PETSC_FALSE) then
      grid%dump_flag = PETSC_TRUE
      grid%dump_str = tmp_str
   endif

   call PrintInitialConditions(grid, atts, rank, size)

   !----------! Begin the main program !----------!
   ! Read the input file and begin internal structuring
   !call mpi_print('> Reading coordinates/elements...',rank,io_rank)
   call GridRead(grid, rank, size) ! Reads a grid in avs format

   if (rank == io_rank) call cpu_time(time_start) ! Begin sys clock timer

   call PrintMeshAttributes(grid, rank, size)

   call mpi_print('', rank, io_rank)
   call mpi_print('=================================================================', rank, io_rank)
   call mpi_print('RUNNING PROGRAM...', rank, io_rank)
   call mpi_print('-----------------------------------------------------------------', rank, io_rank)

   call mpi_print('> Scatter/gather the coordinates...', rank, io_rank)
   call ScatterGatherCoordinates(grid, rank, size)
   call mpi_print('> Building connection matrix...', rank, io_rank)
   call AllocateConnectMatrix(grid, atts)

   ! Calculate geometric coefficients
   call mpi_print('> Calculating Voronoi tesselation...', rank, io_rank)

   ! Choose function based on dimensionality
   if (grid%ndim == 3) then
      call CalculateGeometricCoeff3D(grid, rank, atts)
   else
      call CalculateGeometricCoeff(grid, rank, atts)
   endif

   call mpi_print('> Reconstructing the full matrices...', rank, io_rank)
   call CreateEdgeMatrix(grid, rank)

   call mpi_print('> Writing to file...', rank, io_rank)
   call mpi_print('> Computing mesh statistics...', rank, io_rank)
   call mpi_print('=================================================================', rank, io_rank)

   if (atts%are_on) call Diagnostics(grid, atts, rank, size)

   if (rank == io_rank) call cpu_time(time_finish)
   if (rank == io_rank) call save_time(time_start, time_finish, 2)

   ! Determine what the desired output format is, then write to it
   if (out_type == 1) then
      call GridWriteFEHM(grid, atts, rank, size)
      !call TEMPGridWriteFEHM(grid,atts,rank,size)
   elseif (out_type == 2) then
      call GridWriteTOUGH2(grid, rank, size, 5)
   elseif (out_type == 3) then
      !call CreateConnMatrix(grid,size,rank)
      call GridWritePFLOTRAN(grid, rank, size)
   elseif (out_type == 4) then
      call GridWriteHDF5(grid, rank, size)
   endif

   call GridDestroy(grid)

   call PetscFinalize(ierr); CHKERRQ(ierr)

contains

   subroutine short_help()
      print '(a)', 'ERROR: No mesh has been passed as input.'
      print '(a)', ''
      print '(a)', 'USAGE'
      print '(a)', '     voronoi [-avs | -lg FILE] [-o FILE] [-zone FILE] [-type OUTTYPE]'
      print '(a)', '             [-help] [-d] [-t] [-compress] [-cv CONTROL_VOLUME]'
      print '(a)', ''
      print '(a)', '      For more information, run with the help flag (-help):'
      print '(a)', ''
      print '(a)', '             ./voronoi -help'
      call EXIT(0)
   end subroutine short_help

   subroutine print_help()
      print '(a)', 'SYNOPSIS'
      print '(a)', '     voronoi [-avs | -lg FILE] [-o FILE] [-type outtype] [-v VAL]'
      print '(a)', '             [-help] [-d] [-t] [-nc]'
      print '(a)', ''
      print '(a)', 'DESCRIPTION'
      print '(a)', '     -avs FILE'
      print '(a)', '             Configure Voronoi to receive FILE as input. FILE must be'
      print '(a)', '             an AVS-UCD mesh with either triangle or tetrahedral elements.'
      print '(a)', ''
      print '(a)', '     -lg FILE'
      print '(a)', '             Runs LaGriT infile FILE to completion, then performs a Voronoi'
      print '(a)', '             tesselation on the current mesh object.'
      print '(a)', '             For more information on LaGriT, visit: http://lagrit.lanl.gov'
      print '(a)', ''
      print '(a)', '     -o FILE'
      print '(a)', '             Define the output file to be written to. Defaults to "voronoi.stor".'
      print '(a)', ''
      print '(a)', '     -cv CONTROL_VOLUME'
      print '(a)', '             Performs calculations on a Voronoi or Median control volume.'
      READ *
      print '(a)', '     -type outtype'
      print '(a)', '             Supported arguments for outtype are: "fehm", "pflotran", "tough2", and "hdf5".'
      print '(a)', '             If outtype is not specified, the program will default to FEHM.'
      print '(a)', ''
      print '(a)', '     -v VAL'
      print '(a)', '             Verbosity level - change level of printed console output.'
      print '(a)', '                     -v 2 (default): Run with all console output'
      print '(a)', '                     -v 1: Hide the banner; print the rest of program execution'
      print '(a)', '                     -v 0: Hide all console output'
      print '(a)', ''
      print '(a)', '     -d'
      print '(a)', '             Displays a diagnostics screen at the end of the program, generating'
      print '(a)', '             information on the tesselation and MESH quality.'
      print '(a)', ''
      print '(a)', '     -help'
      print '(a)', '             Displays this help screen. This screen will also appear if -avs'
      print '(a)', '             or -lg are not passed as arguments.'
      print '(a)', ''
      print '(a)', '     -t'
      print '(a)', '             Run test suite to verify the integrity of the build and configuration.'
      print '(a)', '             This runs the program with many different mesh geometries as input'
      print '(a)', '             and compares output against "gold standard" files.'
      print '(a)', ''
      print '(a)', '     -compress'
      print '(a)', '             Generate a STOR file with area coefficient compression. That is,'
      print '(a)', '             duplicate area coefficients will be removed.'
      print '(a)', ''
      READ *
      print '(a)', 'EXAMPLES'
      print '(a)', '     The simplest test case is to build a Voronoi tesselation from an AVS mesh:'
      print '(a)', ''
      print '(a)', '             ./voronoi -avs test/2D/quad/quad6.inp'
      print '(a)', ''
      print '(a)', '     This command will output a voronoi.stor STOR file in the directory that Voronoi lives.'
      call EXIT(0)
   end subroutine print_help

   subroutine banner()
      PetscInt :: verbosity
      common verbosity
      if (verbosity < 2) return
      print '(a)',  ''
      print '(a)',  '                                            |'
      print '(a)',  ' _      ___   ___   ___   _      ___   _    |  A massively parallel Voronoi'
      print '(a)',  '\\\  / /// \ |||_) /// \ |||\ | /// \ |||   |  tesselation generator'
      print '(a)',  ' \\\/  \\\_/ ||| \ \\\_/ ||| \| \\\_/ |||   |  '
      print '(5a)', '                                            |  v', version, ' (', date, ')'
      print '(a)',  '                                            |  (c) 2021 Los Alamos Natl. Lab'
      print '(a)',  '                                            |'
      print '(a)',  ''
   end subroutine banner

   subroutine save_time(start, finish, flag)
#include "petsc/finclude/petscvec.h"
      use petscvec
      PetscReal :: start, finish
      PetscInt :: flag

      if (flag == 1) then
         print *, ''
         print '("Finished in",f6.3," seconds.")', finish - start
      elseif (flag == 2) then
         print *, ''
         print '("Finished in",f6.3," seconds.")', finish - start
         open (826, file='._runtime')
         write (826, *) (finish - start)
         close (826)
      endif

   end subroutine save_time

end program
