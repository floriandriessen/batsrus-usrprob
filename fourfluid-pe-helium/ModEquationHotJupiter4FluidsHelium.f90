!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, Redefine1 => Pe_

  implicit none

  save

  ! This equation module contains the multi-fluid MHD equations for a planet
  ! with an atmosphere consisting of ions (H+, He+) and neutrals (H, He).
  ! Electron fluid is treated via pressure.
  character(len=*), parameter :: NameEquation='Exoplanet four-fluid + Pe MHD'

  ! Number of variables without energy:
  integer, parameter :: nVar = 24

  ! Adapt fluid settings from ModSingleFluid
  integer, parameter :: nFluid    = 4
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 2
  logical, parameter :: IsMhd     = .false.
  real               :: MassFluid_I(nFluid) = [ 1.0, 4.0, 1.0, 4.0 ]

  ! Tack these strings to MHD state variable names in TECplot output
  character(len=3), parameter :: &
       NameFluid_I(nFluid) = [ 'Hp ', 'Hep', 'H1 ', 'He4' ]

  ! Named indexes for State_VGB and other variables.
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_      =  1, HpRho_ = 1, &
       RhoUx_    =  2, Ux_ = 2, HpRhoUx_ = 2, HpUx_ = 2, &
       RhoUy_    =  3, Uy_ = 3, HpRhoUy_ = 3, HpUy_ = 3, &
       RhoUz_    =  4, Uz_ = 4, HpRhoUz_ = 4, HpUz_ = 4, &
       Bx_       =  5, &
       By_       =  6, &
       Bz_       =  7, &
       p_        =  8, HpP_ = 8, &
       HepRho_   =  9, &
       HepRhoUx_ = 10, HepUx_ = 10, &
       HepRhoUy_ = 11, HepUy_ = 11, &
       HepRhoUz_ = 12, HepUz_ = 12, &
       HepP_     = 13, &
       H1Rho_    = 14, &
       H1RhoUx_  = 15, H1Ux_ = 15, &
       H1RhoUy_  = 16, H1Uy_ = 16, &
       H1RhoUz_  = 17, H1Uz_ = 17, &
       H1P_      = 18, &
       He4Rho_   = 19, &
       He4RhoUx_ = 20, He4Ux_ = 20, &
       He4RhoUy_ = 21, He4Uy_ = 21, &
       He4RhoUz_ = 22, He4Uz_ = 22, &
       He4P_     = 23, &
       Pe_       = 24, &
       Energy_   = nVar+1, HpEnergy_ = nVar+1, &
       HepEnergy_= nVar+2, &
       H1Energy_ = nVar+3, &
       He4Energy_= nVar+4

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_-1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)   = [HpRho_,   HepRho_,   H1Rho_,   He4Rho_],   &
       iRhoUx_I(nFluid) = [HpRhoUx_, HepRhoUx_, H1RhoUx_, He4RhoUx_], &
       iRhoUy_I(nFluid) = [HpRhoUy_, HepRhoUy_, H1RhoUy_, He4RhoUy_], &
       iRhoUz_I(nFluid) = [HpRhoUz_, HepRhoUz_, H1RhoUz_, He4RhoUz_], &
       iP_I(nFluid)     = [HpP_,     HepP_,     H1P_,     He4P_]

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
       1.0, & ! HepRho_
       0.0, & ! HepRhoUx_
       0.0, & ! HepRhoUy_
       0.0, & ! HepRhoUz_
       1.0, & ! HepP_
       1.0, & ! H1Rho_
       0.0, & ! H1RhoUx_
       0.0, & ! H1RhoUy_
       0.0, & ! H1RhoUz_
       1.0, & ! H1P_
       1.0, & ! He4Rho_
       0.0, & ! He4RhoUx_
       0.0, & ! He4RhoUy_
       0.0, & ! He4RhoUz_
       1.0, & ! He4P_
       1.0, & ! Pe_
       1.0, & ! HpEnergy_
       1.0, & ! HepEnergy_
       1.0, & ! H1Energy_
       1.0 ]  ! He4Energy_

  ! The names of the variables used in i/o
  character(len=6) :: NameVar_V(nVar+nFluid) = [ &
       'HpRho ', & ! HpRho_
       'HpMx  ', & ! HpRhoUx_
       'HpMy  ', & ! HpRhoUy_
       'HpMz  ', & ! HpRhoUz_
       'Bx    ', & ! Bx_
       'By    ', & ! By_
       'Bz    ', & ! Bz_
       'HpP   ', & ! HpP_
       'HepRho', & ! HepRho_
       'HepMx ', & ! HepRhoUx_
       'HepMy ', & ! HepRhoUy_
       'HepMz ', & ! HepRhoUz_
       'HepP  ', & ! HepP_
       'H1Rho ', & ! H1Rho_
       'H1Mx  ', & ! H1RhoUx_
       'H1My  ', & ! H1RhoUy_
       'H1Mz  ', & ! H1RhoUz_
       'H1P   ', & ! H1P_
       'He4Rho', & ! He4Rho_
       'He4Mx ', & ! He4RhoUx_
       'He4My ', & ! He4RhoUy_
       'He4Mz ', & ! He4RhoUz_
       'He4P  ', & ! He4P_
       'Pe_   ', & ! Pe_
       'HpE   ', & ! HpEnergy_
       'HepE  ', & ! HepEnergy_
       'H1E   ', & ! H1Energy_
       'He4E  ' ]  ! He4Energy_

  ! Scalars to be advected (none)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

end module ModVarIndexes
