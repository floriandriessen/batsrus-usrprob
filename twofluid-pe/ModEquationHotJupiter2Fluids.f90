!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, Redefine1 => Pe_

  implicit none

  save

  ! This equation module contains the multi-fluid MHD equations for a planet
  ! with an atmosphere consisting of ions (H+) and neutrals (H). Electrons are
  ! separately treated via pressure.
  character(len=*), parameter :: NameEquation='Exoplanet two-fluid + Pe (M)HD'

  ! Number of variables without energy:
  integer, parameter :: nVar = 14

  ! Adapt fluid settings from ModSingleFluid
  integer, parameter :: nFluid    = 2
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 1
  logical, parameter :: IsMhd     = .true.  ! the first fluid (ion) obeys MHD
  real               :: MassFluid_I(nFluid) = [ 1.0, 1.0 ]

  ! Tack these strings to MHD state variable names in TECplot output
  character(len=2), parameter :: NameFluid_I(nFluid) = [ 'Hp', 'H1' ]

  ! Named indexes for State_VGB and other variables.
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  !
  ! NOTE: the first fluid MHD variables need the default name for internal code
  !       calls. For user routines a placeholder name can be used instead.
  integer, parameter :: &
       Rho_      =  1, HpRho_ = 1, &
       RhoUx_    =  2, Ux_ = 2, HpRhoUx_ = 2, HpUx_ = 2, &
       RhoUy_    =  3, Uy_ = 3, HpRhoUy_ = 3, HpUy_ = 3, &
       RhoUz_    =  4, Uz_ = 4, HpRhoUz_ = 4, HpUz_ = 4, &
       Bx_       =  5, &
       By_       =  6, &
       Bz_       =  7, &
       p_        =  8, HpP_ = 8, &
       H1Rho_    =  9, &
       H1RhoUx_  = 10, H1Ux_ = 10, &
       H1RhoUy_  = 11, H1Uy_ = 11, &
       H1RhoUz_  = 12, H1Uz_ = 12, &
       H1P_      = 13, &
       Pe_       = 14, &
       Energy_   = nVar+1, HpEnergy_ = nVar+1, &
       H1Energy_ = nVar+2

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_-1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = [HpRho_,   H1Rho_]
  integer, parameter :: iRhoUx_I(nFluid) = [HpRhoUx_, H1RhoUx_]
  integer, parameter :: iRhoUy_I(nFluid) = [HpRhoUy_, H1RhoUy_]
  integer, parameter :: iRhoUz_I(nFluid) = [HpRhoUz_, H1RhoUz_]
  integer, parameter :: iP_I(nFluid)     = [HpP_,     H1P_]

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = [ &
       1.0, & ! HpRho_
       0.0, & ! HpRhoUx_
       0.0, & ! HpRhoUy_
       0.0, & ! HpRhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! HpP_
       1.0, & ! H1Rho_
       0.0, & ! H1RhoUx_
       0.0, & ! H1RhoUy_
       0.0, & ! H1RhoUz_
       1.0, & ! H1P_
       1.0, & ! Pe_
       1.0, & ! HpEnergy_
       1.0 ]  ! H1Energy_

  ! The names of the variables used in i/o
  character(len=5) :: NameVar_V(nVar+nFluid) = [ &
       'HpRho', & ! HpRho_
       'HpMx ', & ! HpRhoUx_
       'HpMy ', & ! HpRhoUy_
       'HpMz ', & ! HpRhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'HpP  ', & ! HpP_
       'H1Rho', & ! H1Rho_
       'H1Mx ', & ! H1RhoUx_
       'H1My ', & ! H1RhoUy_
       'H1Mz ', & ! H1RhoUz_
       'H1P  ', & ! H1P_
       'Pe_  ', & ! Pe_
       'HpE  ', & ! HpEnergy_
       'H1E  ' ]  ! H1Energy_

  ! Scalars to be advected (none)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

end module ModVarIndexes
