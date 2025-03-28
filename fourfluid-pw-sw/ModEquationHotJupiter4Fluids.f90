!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables

  implicit none

  save

  ! This equation module contains the multi-fluid MHD equations for a planet
  ! with an atmosphere consisting of ions (H+) and neutrals (H) and a stellar
  ! wind consisting of hot ions (H+) and cold neutrals (H)
  character(len=*), parameter :: NameEquation='Star-exoplanet four-fluid MHD'

  ! Number of variables without energy:
  integer, parameter :: nVar = 23

  ! Adapt fluid settings from ModSingleFluid
  integer, parameter :: nFluid    = 4
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 2
  logical, parameter :: IsMhd     = .false.
  real               :: MassFluid_I(nFluid) = [ 1.0, 1.0, 1.0, 1.0 ]

  ! Tack these strings to MHD state variable names in TECplot output
  character(len=4), parameter :: NameFluid_I(nFluid) = &
       [ 'SwHp', 'Hp  ', 'SwH1', 'H1  ' ]

  ! Named indexes for State_VGB and other variables.
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_        =  1, SwHpRho_ = 1, &
       RhoUx_      =  2, Ux_ = 2, SwHpRhoUx_ = 2, SwHpUx_ = 2, &
       RhoUy_      =  3, Uy_ = 3, SwHpRhoUy_ = 3, SwHpUy_ = 3, &
       RhoUz_      =  4, Uz_ = 4, SwHpRhoUz_ = 4, SwHpUz_ = 4, &
       Bx_         =  5, &
       By_         =  6, &
       Bz_         =  7, &
       p_          =  8, SwHpP_ = 8, &
       HpRho_      =  9, &
       HpRhoUx_    = 10, HpUx_ = 10, &
       HpRhoUy_    = 11, HpUy_ = 11, &
       HpRhoUz_    = 12, HpUz_ = 12, &
       HpP_        = 13, &
       SwH1Rho_    = 14, &
       SwH1RhoUx_  = 15, SwH1Ux_ = 15, &
       SwH1RhoUy_  = 16, SwH1Uy_ = 16, &
       SwH1RhoUz_  = 17, SwH1Uz_ = 17, &
       SwH1P_      = 18, &
       H1Rho_      = 19, &
       H1RhoUx_    = 20, H1Ux_ = 20, &
       H1RhoUy_    = 21, H1Uy_ = 21, &
       H1RhoUz_    = 22, H1Uz_ = 22, &
       H1P_        = 23, &
       Energy_     = nVar+1, SwHpEnergy_ = nVar+1, &
       HpEnergy_   = nVar+2,                       &
       SwH1Energy_ = nVar+3,                       &
       H1Energy_   = nVar+4

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_-1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = [SwHpRho_, HpRho_, SwH1Rho_, H1Rho_]
  integer, parameter :: iRhoUx_I(nFluid) = &
       [SwHpRhoUx_, HpRhoUx_, SwH1RhoUx_, H1RhoUx_]
  integer, parameter :: iRhoUy_I(nFluid) = &
       [SwHpRhoUy_, HpRhoUy_, SwH1RhoUy_, H1RhoUy_]
  integer, parameter :: iRhoUz_I(nFluid) = &
       [SwHpRhoUz_, HpRhoUz_, SwH1RhoUz_, H1RhoUz_]
  integer, parameter :: iP_I(nFluid)     = [SwHpP_, HpP_, SwH1P_, H1P_]

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = [ &
       1.0, & ! SwHpRho_
       0.0, & ! SwHpRhoUx_
       0.0, & ! SwHpRhoUy_
       0.0, & ! SwHpRhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! SwHpP_
       1.0, & ! HpRho_
       0.0, & ! HpRhoUx_
       0.0, & ! HpRhoUy_
       0.0, & ! HpRhoUz_
       1.0, & ! HpP_
       1.0, & ! SwH1Rho_
       0.0, & ! SwH1RhoUx_
       0.0, & ! SwH1RhoUy_
       0.0, & ! SwH1RhoUz_
       1.0, & ! SwH1P_
       1.0, & ! H1Rho_
       0.0, & ! H1RhoUx_
       0.0, & ! H1RhoUy_
       0.0, & ! H1RhoUz_
       1.0, & ! H1P_
       1.0, & ! SwHpEnergy_
       1.0, & ! HpEnergy_
       1.0, & ! SwH1Energy_
       1.0 ]  ! H1Energy_

  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = [ &
       'SwHpRho', & ! SwHpRho_
       'SwHpMx ', & ! SwHpRhoUx_
       'SwHpMy ', & ! SwHpRhoUy_
       'SwHpMz ', & ! SwHpRhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'SwHpP  ', & ! SwHpP_
       'HpRho  ', & ! HpRho_
       'HpMx   ', & ! HpRhoUx_
       'HpMy   ', & ! HpRhoUy_
       'HpMz   ', & ! HpRhoUz_
       'HpP    ', & ! HpP_
       'SwH1Rho', & ! SwH1Rho_
       'SwH1Mx ', & ! SwH1RhoUx_
       'SwH1My ', & ! SwH1RhoUy_
       'SwH1Mz ', & ! SwH1RhoUz_
       'SwH1P  ', & ! SwH1P_
       'H1Rho  ', & ! H1Rho_
       'H1Mx   ', & ! H1RhoUx_
       'H1My   ', & ! H1RhoUy_
       'H1Mz   ', & ! H1RhoUz_
       'H1P    ', & ! H1P_
       'SwHpE  ', & ! SwHpEnergy_
       'HpE    ', & ! HpEnergy_
       'SwH1E  ', & ! SwH1Energy_
       'H1E    ' ]  ! H1Energy_

  ! Scalars to be advected (none)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

end module ModVarIndexes
