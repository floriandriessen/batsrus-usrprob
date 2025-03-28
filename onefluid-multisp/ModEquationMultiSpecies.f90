!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModSingleFluid
  use ModExtraVariables

  implicit none

  save

  character(len=*), parameter :: &
       NameEquationFile = 'ModEquationMultiSpecies.f90'

  ! This equation module contains the multi-species MHD equations for a planet
  ! with an atmosphere consisting of ions (H+) and neutrals (H).
  ! NOTE: BATSRUS is only aware of ions in multi-species, to account for the
  !       neutrals extra computations are done in ModUser.f90
  character(len=*), parameter :: NameEquation='Exoplanet multi-species (M)HD'

  ! Number of variables without energy:
  integer, parameter :: nVar = 11

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_    = 1,    &
       HpRho_  = 2,    & ! Protons
       H1Rho_  = 3,    & ! Neutrals
       PaSc_   = 4,    & ! Passive scalar
       RhoUx_  = 5,    &
       RhoUy_  = 6,    &
       RhoUz_  = 7,    &
       Bx_     = 8,    &
       By_     = 9,    &
       Bz_     = 10,   &
       p_      = nVar, &
       Energy_ = nVar+1

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = [Rho_]
  integer, parameter :: iRhoUx_I(nFluid) = [RhoUx_]
  integer, parameter :: iRhoUy_I(nFluid) = [RhoUy_]
  integer, parameter :: iRhoUz_I(nFluid) = [RhoUz_]
  integer, parameter :: iP_I(nFluid)     = [p_]

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+1) = [ &
       1.0, & ! Rho_
       1.0, & ! HpRho_
       1.0, & ! H1Rho_
       0.0, & ! passive scalar PaSc_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! p_
       1.0 ]  ! Energy_

  ! The names of the variables used in i/o
  character(len=4) :: NameVar_V(nVar+1) = [ &
       'Rho ', & ! Rho_
       'Hp  ', & ! HpRho_
       'H1  ', & ! H1Rho_
       'PaSc', & ! PaSc_
       'Mx  ', & ! RhoUx_
       'My  ', & ! RhoUy_
       'Mz  ', & ! RhoUz_
       'Bx  ', & ! Bx_
       'By  ', & ! By_
       'Bz  ', & ! Bz_
       'p   ', & ! p_
       'e   ' ]  ! Energy_

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! Scalars to be advected.
  integer, parameter :: ScalarFirst_ = HpRho_, ScalarLast_ = PaSc_

end module ModVarIndexes
