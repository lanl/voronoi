module Voronoi_Constants_module 

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: ndim = 3 ! deprecated
  PetscInt, parameter, public :: MAXSTRINGLENGTH = 512
  PetscInt, parameter, public :: MAXWORDLENGTH = 32
  
  PetscInt, parameter, public :: ZERO_INTEGER = 0
  PetscInt, parameter, public :: ONE_INTEGER = 1
  PetscInt, parameter, public :: TWO_INTEGER = 2
  PetscInt, parameter, public :: THREE_INTEGER = 3
  PetscInt, parameter, public :: FOUR_INTEGER = 4
  PetscInt, parameter, public :: FIVE_INTEGER = 5
  PetscInt, parameter, public :: SIX_INTEGER = 6
  PetscInt, parameter, public :: SEVEN_INTEGER = 7
  PetscInt, parameter, public :: EIGHT_INTEGER = 8
  PetscInt, parameter, public :: NINE_INTEGER = 9
  PetscInt, parameter, public :: TEN_INTEGER = 10
  PetscInt, parameter, public :: ELEVEN_INTEGER = 11
  PetscInt, parameter, public :: TWELVE_INTEGER = 12
  PetscInt, parameter, public :: NEG_ONE_INTEGER = -1
  
  PetscMPIInt, parameter, public :: ZERO_INTEGER_MPI = ZERO_INTEGER
  PetscMPIInt, parameter, public :: ONE_INTEGER_MPI = ONE_INTEGER
  PetscMPIInt, parameter, public :: TWO_INTEGER_MPI = TWO_INTEGER
  PetscMPIInt, parameter, public :: THREE_INTEGER_MPI = THREE_INTEGER
  PetscMPIInt, parameter, public :: FOUR_INTEGER_MPI = FOUR_INTEGER
  PetscMPIInt, parameter, public :: SIX_INTEGER_MPI = SIX_INTEGER
  PetscMPIInt, parameter, public :: SEVEN_INTEGER_MPI = SEVEN_INTEGER
  PetscMPIInt, parameter, public :: TWELVE_INTEGER_MPI = TWELVE_INTEGER
  PetscMPIInt, parameter, public :: MAXSTRINGLENGTH_MPI = MAXSTRINGLENGTH

  ! uninitialized values
  PetscInt, parameter, public :: UNINITIALIZED_INTEGER = -999
  PetscReal, parameter, public :: UNINITIALIZED_DOUBLE = -999.d0

end module Voronoi_Constants_module
