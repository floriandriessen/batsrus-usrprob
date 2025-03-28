!==============================================================================
!  Multi-species (M)HD model for atmospheric escape of extrasolar gas giants.
!==============================================================================
! The planetary atmosphere consists of the following species:
!   > singly ionised hydrogen
!   > neutral hydrogen
! with heat capacity ratio of 5/3. Global charge neutrality holds. Strictly
! speaking this is not a multi-species model for BATSRUS as multi-species =
! multiple ion species. The correct inclusion of neutrals is in this ModUser.
!
! Calculating tau by interpolating neutral density to a uniform 3D grid. Each
! processor implements a portion of the grid, then in ModBatsrusMethods.f90 I
! sum the contribution of each processor into a master grid of densities, which
! are then used to sum along x to get tau(x,y,z). This uniform grid of tau is
! then interpolated back to the AMR grid in calc_sources.
! The forces (Coriolis and tidal) and the input stellar wind as implemented by
! Stephen Carolan
!
! HISTORY (AAV, SC)
! 19 Jan 2021 : fixed issue that code is not working correctly after restart
!               (although I don't understand why...)
! 22 Feb 2021 : Optimising way to read the SW coefficients, to reduce number of
!               open/read/close wind.dat files.
! 19 Mar 2021 : implementing passive scalar. AAV
! 26 Mar 2021 : editing the sources for passive scalar.
! 22 Apr 2021 : adding collisional cooling / collisional ionisation
! 21 May 2021 : adding new variables for plotting
! 10 Jun 2021 : changing the resolution of the interpolated grid
! 03 Jun 2022 : (re)incorporating magnetic field
!
! HISTORY (FLO)
!   Dec 2022 : major formatting, display, and deletion of unused code
!   Jan 2023 : bugfix in the treatment of ion and neutral species; neutrals
!              have to be handled in ModUser for correct body values
!   Sep 2023 : standardise code following other versions for thesis student
!==============================================================================
module ModUser

  use BATL_lib, ONLY: test_start, test_stop
  use ModUserEmpty,                             &
       IMPLEMENTED1  => user_read_inputs,       &
       IMPLEMENTED2  => user_io_units,          &
       IMPLEMENTED3  => user_normalization,     &
       IMPLEMENTED4  => user_init_session,      &
       IMPLEMENTED5  => user_set_face_boundary, &
       IMPLEMENTED6  => user_set_cell_boundary, &
       IMPLEMENTED7  => user_set_ics,           &
       IMPLEMENTED8  => user_calc_sources,      &
       IMPLEMENTED9  => user_set_plot_var,      &
       IMPLEMENTED10 => user_action,            &
       IMPLEMENTED11 => user_get_log_var,       &
       IMPLEMENTED12 => user_update_states

  ! List of public methods
  include 'user_module.h'

  ! Make public in order to add to ModBatsrusMethods
  public :: create_uniform_3d_grid_tau

  real,              parameter :: VersionUserModule = 6.0
  character (len=*), parameter :: NameUserModule = &
       'Multi-species (M)HD simulation for extrasolar gas giants'

  ! Orbital properties planet
  real, public :: Rorbit_au=0.0, OmegaOrbit_dim
  real         :: Rorbit, OmegaOrbit, PeriodOrbitDays

  ! Star properties
  real :: Mstar=0.0, Rstar=0.0, B0StarCgs=0.0, MdotSwMsunYr=0.0, TempSwSi=0.0
  real :: gStar, MdotSw, B0Star

  ! User switches to add physics to the problem
  logical :: UseCoriolis=.false., UseTidal=.false., UseStellarWind=.false.

  ! User strings for problem specific settings
  character(len=30) :: NameStar='empty', TypeBfieldStar='none'

  ! Used for the IC planetary wind
  real :: Vrad0PwIo=0.0, VinfPwIo=0.0, BetaIndex=1.0, VinfPw, Vrad0Pw

  ! Used for radiation hydrodynamics
  real :: FxuvCgs=0.0, EphotCgs=0.0, SigmaNu0Cgs, EpsilonNu, EphotEv

  ! Coefficients for simplifying numerical factors in rate equations
  real :: PhotoIonisationCoef, ElectronImpactCoef, RecombinationCoef
  real :: HeatPhotoCoef, CoolElectronImpactCoef, CoolLyalphaCoef

  ! Passive scalar for planetary outflow (<0) and stellar wind (>0)
  real, parameter :: pScalarPw = -1000.0, pScalarSw = 1000.0

  ! Polynomial coefficients for stellar wind injected in the grid
  integer, parameter :: OrderCoef=9
  real               :: SwCoef_I(OrderCoef)

  ! Used for printing new variables
  real, allocatable :: PlotExp_mTau_CB(:,:,:,:), PlotHeat_CB(:,:,:,:), &
       PlotCoolLya_CB(:,:,:,:), PlotCoolColl_CB(:,:,:,:),              &
       PlotRatePhotoIon_CB(:,:,:,:), PlotRateCollIon_CB(:,:,:,:),      &
       PlotRateRec_CB(:,:,:,:)

  ! Interpolation grid: use even numbers to avoid interpolation to fall in xyz
  ! in-between 2 blocks
  integer, parameter :: size_3d_i=200, size_3d_j=200, size_3d_k=128
  real               :: x_array(size_3d_i),y_array(size_3d_j),z_array(size_3d_k)
  real               :: array_3d_dens_plot(size_3d_i,size_3d_j,size_3d_k)

  ! Optical depth related: public so can be added to modBatsrusMethods
  real, public :: array_3d_dens(size_3d_i, size_3d_j, size_3d_k)=0.0
  real, public :: array_3d_Exp_mTau(size_3d_i, size_3d_j, size_3d_k)=1.0

  ! Neutral specie properties
  real, parameter :: BodyNNeuSpeciesDim = 2.39e11, MassNeuSpecies = 1.0
  real :: BodyRhoNeuSpecies, BodyPNeuSpecies

contains

!==============================================================================
! Read the user variables given in PARAM.in with #USERINPUT.
!==============================================================================
  subroutine user_read_inputs

    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModUtilities, ONLY: lower_case
    use BATL_lib,     ONLY: iProc

    implicit none

    character(len=100) :: NameCommand

    logical :: DoTest
    character (len=*), parameter :: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do
      if (.not.read_line() ) EXIT
      if (.not.read_command(NameCommand)) CYCLE

      select case(NameCommand)

      case('#USERINPUTEND')
        if (iProc == 0) write(*,*) 'USERINPUTEND'
        EXIT

      ! Extra source terms for MHD; by default set to False
      case('#USR_SOURCES')
        call read_var('UseCoriolis', UseCoriolis)
        call read_var('UseTidal', UseTidal)

      ! Stellar mass and radius and planet orbital distance
      case('#USR_SYSTEMPROPERTIES')
        call read_var('Mstar', Mstar)
        call read_var('Rstar', Rstar)
        call read_var('Rorbit_au', Rorbit_au)

      ! Stellar flux (cgs) at orbit and interaction energy
      case('#USR_FXUV')
        call read_var('FxuvCgs', FxuvCgs)
        call read_var('EphotEv', EphotEv)

      ! Planetary atmosphere
      case('#USR_INIT_PL_ATM')
        call read_var('Vrad0Pw', Vrad0PwIo)
        call read_var('VinfPw', VinfPwIo)
        call read_var('BetaIndex', BetaIndex)

      ! Stellar wind properties
      case('#USR_INJECT_SW')
        call read_var('UseStellarWind', UseStellarWind)

        ! Read in the sw-velocity.inp, stellar mass loss, and
        ! wind temperature (K)
        if (UseStellarWind) then
          call read_var('NameStar', NameStar)
          call lower_case(NameStar)
          call read_var('MdotSw', MdotSwMsunYr)
          call read_var('TemperatureSw', TempSwSi)
          call read_stellar_wind_coeff(OrderCoef, SwCoef_I)
          call read_var('TypeBfieldStar', TypeBfieldStar)

          if ('TypeBfieldStar' /= 'none') call read_var('B0StarCgs', B0StarCgs)
        endif

      case default
        call stop_mpi('read_inputs: unrecognized command: '//NameCommand)
      end select
    enddo

    call test_stop(NameSub, DoTest)

  end subroutine user_read_inputs

!==============================================================================
! Conversion between units if not default input of SI. The following should
! always be true: No2Si_V*Si2Io_V = No2Io_V.
! This subroutine is read after each #RUN in PARAM.in.
!==============================================================================
  subroutine user_io_units

    use ModPhysics, ONLY: Io2Si_V, Si2Io_V, No2Si_V, No2Io_V, Io2No_V,       &
         NameTecUnit_V, rPlanetSi, UnitX_, UnitU_, UnitRho_, UnitT_, UnitN_, &
         UnitP_, UnitB_, UnitRhoU_, UnitEnergydens_, UnitPoynting_, UnitJ_,  &
         UnitElectric_, UnitTemperature_, UnitDivB_, UnitAngle_
    use ModConst,   ONLY: cRadToDeg

    implicit none

    logical :: DoTest
    character (len=*), parameter :: NameSub = 'user_io_units'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    !\
    ! Set I/O units in terms of SI units if they differ
    ! As a default use SI units, so below only the differences need to be set
    ! For example, input densities in g/cm^3. So density has to be multiplied
    ! by 1.0E+3 to be transfomed in SI units kg/m^3.
    !/
    Io2Si_V(UnitX_)           = rPlanetSi ! Rplanet -> m
    Io2Si_V(UnitU_)           = 1.0E+3    ! km/s -> m/s
    Io2Si_V(UnitRho_)         = 1.0E+3    ! g/cm^3 -> kg/m^3
    Io2Si_V(UnitT_)           = 1.0       ! s
    Io2Si_V(UnitN_)           = 1.0E+6    ! #/cm^3 -> #/m^3
    Io2Si_V(UnitP_)           = 1.0E-1    ! dyne/cm^2 -> Pa
    Io2Si_V(UnitB_)           = 1.0E-4    ! Gauss -> Tesla
    Io2Si_V(UnitRhoU_)        = 1.0E+1    ! g/cm^2/s -> kg/m^2/s
    Io2Si_V(UnitEnergydens_)  = 1.0E-1    ! erg/cm^3 -> J/m^3
    Io2Si_V(UnitPoynting_)    = 1.0E-3    ! erg/cm^2/s -> J/m^2/s
    Io2Si_V(UnitJ_)           = 1.0E-6    ! microA/m^2
    Io2Si_V(UnitElectric_)    = 1.0       ! V/m
    Io2Si_V(UnitTemperature_) = 1.0       ! K
    Io2Si_V(UnitDivB_)        = 1.0E-2    ! Gauss/cm
    Io2Si_V(UnitAngle_)       = cRadToDeg ! degrees -> rad

    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V * Si2Io_V
    Io2No_V = 1/No2Io_V

    !\
    ! Set strings for writing Tecplot output
    !/
    NameTecUnit_V(UnitX_)          = '(R_p)'
    NameTecUnit_V(UnitRho_)        = '(g/cm^3)'
    NameTecUnit_V(UnitU_)          = '(km/s)'
    NameTecUnit_V(UnitN_)          = '(amu/cm^3)'
    NameTecUnit_V(UnitP_)          = '(dyn/cm^2)'
    NameTecUnit_V(UnitB_)          = '(G)'
    NameTecUnit_V(UnitRhoU_)       = '(g/cm^2/s)'
    NameTecUnit_V(UnitEnergyDens_) = '(erg/cm^3)'
    NameTecUnit_V(UnitJ_)          = '(`mA/m^2)'
    NameTecUnit_V(UnitDivB_)       = '(G/cm)'
    NameTecUnit_V(UnitAngle_)      = '(deg)'
    NameTecUnit_V(UnitPoynting_)   = '(erg/cm^2/s)'

    call test_stop(NameSub, DoTest)

  end subroutine user_io_units

!==============================================================================
! Do user normalization. To activate in Param.in set #NORMALIZATION with USER.
!==============================================================================
  subroutine user_normalization

    use ModVarindexes, ONLY: MassFluid_I
    use ModPhysics,    ONLY: No2Si_V, UnitX_, UnitU_, UnitRho_, rPlanetSi, &
         Gamma, BodyNDim_I, BodyTDim_I
    use ModConst,      ONLY: cBoltzmann, cProtonMass

    implicit none

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_normalization'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Three independent normalisation units (in SI)
    No2Si_V(UnitX_)   = rPlanetSi
    No2Si_V(UnitU_)   = sqrt(Gamma*cBoltzmann*BodyTDim_I(1)/cProtonMass)
    No2Si_V(UnitRho_) = 1e6*cProtonMass * ( MassFluid_I(1)*BodyNDim_I(1) &
         + MassNeuSpecies*BodyNNeuSpeciesDim )

    call test_stop(NameSub, DoTest)

  end subroutine user_normalization

!==============================================================================
! Initialize new global variables to be used within a session. This routine is
! read at the first iteration and after a restart.
!==============================================================================
  subroutine user_init_session

    use ModMain,       ONLY: Body1_
    use ModVarIndexes, ONLY: Rho_, HpRho_, H1Rho_, P_
    use ModPhysics,    ONLY: Si2No_V, Io2No_V, No2Io_V, UnitX_, UnitU_, &
         UnitN_,  UnitB_, UnitTemperature_, UnitEnergyDens_, UnitMass_, &
         UnitT_, mSun, BodyRho_I, BodyP_I, BodyTDim_I, FaceState_VI
    use ModConst,      ONLY: cGravitation, cAU, cSecondPerYear, cSecondPerDay,&
         cEV, cErg
    use ModNumConst,   ONLY: cPi

    implicit none

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    !\
    ! Deriving other initial variables
    !/

    ! Stellar mass-loss rate; Msun/yr to kg/s to normalised
    MdotSw = MdotSwMsunYr * mSun / cSecondPerYear &
         * Si2No_V(UnitMass_) / Si2No_V(UnitT_)

    ! G*Mstar; m^3/s^2 to normalised
    gStar = cGravitation * Mstar*mSun * Si2No_V(UnitU_)**2 * Si2No_V(UnitX_)

    ! Orbital distance; au to normalised
    Rorbit = Rorbit_au * cAU * Si2No_V(UnitX_)

    ! Orbit quantities: Omega in (rad/s) and Period in (days)
    OmegaOrbit_dim   = sqrt(cGravitation * Mstar * mSun / (Rorbit_au*cAU)**3)
    OmegaOrbit       = OmegaOrbit_dim / Si2No_V(UnitT_)
    PeriodOrbitDays  = 2.0*cPi / OmegaOrbit_dim / cSecondPerDay

    ! Radiation field quantities:
    ! Photon energy; eV to erg
    EphotCgs = EphotEv * cEV / cErg

    ! Fraction photon energy to be deposited as heat
    EpsilonNu = 1.0 - (13.6/EphotEv)

    ! Photo-ionisation cross-section Hydrogen, Spitzer (1998) eq. 5-6/table 5.1
    SigmaNu0Cgs = 6e-18 * (13.6/EphotEv)**3

    ! Planetary wind velocity; km/s to normalised
    VinfPw  = VinfPwIo * Io2No_V(UnitU_)
    Vrad0Pw = Vrad0PwIo * Io2No_V(UnitU_)

    ! Stellar magnetic field; Gauss to normalised
    B0Star = B0StarCgs * Io2No_V(UnitB_)

    ! Do correct handling of MHD mass density equation and MHD pressure by
    ! including neutrals as separate from ions (and electrons)
    BodyRhoNeuSpecies = BodyNNeuSpeciesDim * Io2No_V(UnitN_) * MassNeuSpecies
    BodyPNeuSpecies   = BodyNNeuSpeciesDim * Io2No_V(UnitN_) &
         * BodyTDim_I(1) * Io2No_V(UnitTemperature_)

    ! Add neutral species on top of ion species for total fluid in body
    BodyRho_I(1) = BodyRho_I(1) + BodyRhoNeuSpecies
    BodyP_I(1)   = BodyP_I(1) + BodyPNeuSpecies

    ! Amend body faces to account for neutral species
    FaceState_VI(Rho_,Body1_)   = BodyRho_I(1)
    FaceState_VI(HpRho_,Body1_) = BodyRho_I(1) - BodyRhoNeuSpecies
    FaceState_VI(H1Rho_,Body1_) = BodyRhoNeuSpecies
    FaceState_VI(P_,Body1_)     = BodyP_I(1)

    ! Numerical factors for mass-loading and collisional rate processes with
    ! proper normalisation for later use in source term computation:
    ! Photo-ionisation coefficient (1/s)
    PhotoIonisationCoef = SigmaNu0Cgs * FxuvCgs / EphotCgs / Io2No_V(UnitT_)

    ! Heating coefficient (erg/cm3/s)
    HeatPhotoCoef = EpsilonNu * SigmaNu0Cgs * FxuvCgs &
         * No2Io_V(UnitN_) * Io2No_V(UnitEnergyDens_) / Io2No_V(UnitT_)

    ! Case-B recombination coefficient (cm3/s): Storey & Hummer (1995), table 1
    RecombinationCoef = 2.7e-13 * No2Io_V(UnitN_) / Io2No_V(UnitT_)

    ! e-impact rate (cm3/s): Voronov (1997), eq. 1 + table 1
    ElectronImpactCoef = 0.291e-7 * No2Io_V(UnitN_) / Io2No_V(UnitT_)

    ! e-impact cooling rate (erg/cm3/s)
    CoolElectronImpactCoef = 0.291e-7 * (13.6*cEV/cErg) * No2Io_V(UnitN_)**2 &
         * Io2No_V(UnitEnergyDens_) / Io2No_V(UnitT_)

    ! Lyman alpha cooling rate (erg/cm3/s): Black (1981), table 3
    CoolLyalphaCoef = 7.5d-19 * No2Io_V(UnitN_)**2 &
         * Io2No_V(UnitEnergyDens_) / Io2No_V(UnitT_)

    ! Generate the grid for the optical depth computation
    call make_uniform_rtgrid

    call test_stop(NameSub, DoTest)

  end subroutine user_init_session

!==============================================================================
! Boundary conditions for a single inner boundary face of a body. Active when
! #INNERBOUNDARY is set to USER. BATSRUS is block cartesian and values inside
! the boundary face must be passed back in cartesian coordinates.
! Values that must be set are: Rho_, Ux_, Uy_, Uz_, p_, Bx_, By_, Bz_
!==============================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModFaceBoundary, ONLY: VarsTrueFace_V
    use ModVarIndexes,   ONLY: nVar, Rho_, HpRho_, H1Rho_, PaSc_, Ux_, Uz_, P_
    use ModPhysics,      ONLY: BodyP_I, BodyRho_I

    implicit none

    real, intent(out) :: VarsGhostFace_V(nVar)

    ! By default apply floating condition to all state variables
    VarsGhostFace_V = VarsTrueFace_V

    ! Fix mass density and pressure to the respective body value
    VarsGhostFace_V(Rho_)   = BodyRho_I(1)
    VarsGhostFace_V(HpRho_) = BodyRho_I(1) - BodyRhoNeuSpecies
    VarsGhostFace_V(H1Rho_) = BodyRhoNeuSpecies
    VarsGhostFace_V(P_)     = BodyP_I(1)

    ! Passive scalar for planetary outflow
    VarsGhostFace_V(PaSc_) = pScalarPw * VarsGhostFace_V(Rho_)

    ! Change sign for velocities to force 0 velocity at the boundary
    VarsGhostFace_V(Ux_:Uz_) = -VarsTrueFace_V(Ux_:Uz_)

  end subroutine user_set_face_boundary

!==============================================================================
! Boundary conditions for the edges of the simulation domain. Active when in
! PARAM.in the #OUTERBOUNDARY is  'userfloat' or 'userInjectSW' for edges.
!==============================================================================
  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, found)

    use ModVarIndexes,   ONLY: Rho_, HpRho_, H1Rho_, PaSc_, RhoUx_, RhoUy_, &
         RhoUz_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_, MassFluid_I
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitX_, UnitU_, &
         UnitTemperature_, rSun
    use ModGeometry,     ONLY: Xyz_DGB
    use ModAdvance,      ONLY: State_VGB
    use ModCellBoundary, ONLY: iMin, iMax, jMin, jMax, kMin, kMax
    use ModNumConst,     ONLY: cPi
    use BATL_lib,        ONLY: nI, nJ, nK

    implicit none

    integer, intent(in)          :: iBlock, iSide
    logical, intent(out)         :: found
    character(len=*), intent(in) :: TypeBc

    ! Local variables
    integer :: i, j, k
    real    :: x, y, z, rFrameStar, rFrameStarRstar, logrp, SwVrad, NumDens
    real    :: BradStar, BthetaStar, CosTheta, SinTheta, CosPhi, SinPhi

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    select case(TypeBc)

    ! Below is a floating boundary condition, taken from ModCellBoundary
    ! with an if statement that checks for inflow on a specific iSide of
    ! the grid. If there is an inflow this sets the inflow momentum to 0
    ! at that boundary.
    case('userfloat')
      found = .true.

      select case(iSide)
      case(1) ! The minimum X boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,1,j,k,iBlock)

          if (State_VGB(RhoUx_,1,j,k,iBlock) > 0.0) then
            State_VGB(RhoUx_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo

      case(2) ! The maximum X boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)

          if (State_VGB(RhoUx_,nI,j,k,iBlock) < 0.0) then
            State_VGB(RhoUx_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo

      case(3) ! The minimum Y boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,1,k,iBlock)

          if (State_VGB(RhoUy_,i,1,k,iBlock) > 0.0) then
            State_VGB(RhoUy_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo

      case(4) ! The maximum Y boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,nJ,k,iBlock)

          if (State_VGB(RhoUy_,i,nJ,k,iBlock) < 0.0) then
            State_VGB(RhoUy_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo

      case(5) ! The minimum Z boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,1,iBlock)

          if (State_VGB(RhoUz_,i,j,1,iBlock) > 0.0) then
            State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo

      case(6) ! The maximum Z boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,nK,iBlock)

          if (State_VGB(RhoUz_,i,j,nK,iBlock) < 0.0) then
            State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
          endif
        enddo; enddo; enddo
      end select

    ! Only works on the negative x boundary. This BC is used to inject the SW
    case('userInjectSW')
      found = .true.

      if (iSide /= 1 .or. .not.UseStellarWind) &
           call stop_mpi('Wrong iSide in user_set_cell_boundary or &
           UseStellarWind set to false, userSW')

      select case(iSide)
      case(1) ! The minimum X boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax

          x = Xyz_DGB(x_,i,j,k,iBlock)
          y = Xyz_DGB(y_,i,j,k,iBlock)
          z = Xyz_DGB(z_,i,j,k,iBlock)
          rFrameStar = sqrt( (Rorbit + x)**2 + y**2 + z**2 )

          ! Radial coordinate in stellar frame in units of Rstar
          rFrameStarRstar = rFrameStar * No2Si_V(UnitX_) / (Rstar*rSun)

          ! Setup the stellar wind based on the sw-velocity.inp file
          logrp = log10(rFrameStarRstar)

          SwVrad = SwCoef_I(1) + SwCoef_I(2) * logrp &
               + SwCoef_I(3) * logrp**2 + SwCoef_I(4) * logrp**3 &
               + SwCoef_I(5) * logrp**4 + SwCoef_I(6) * logrp**5 &
               + SwCoef_I(7) * logrp**6 + SwCoef_I(8) * logrp**7 &
               + SwCoef_I(9) * logrp**8

          ! The fit for vradial is given in 10^[V[cm/s]]. First multiply by
          ! 1e-2 to make it in m/s (SI) and then multiply by Si2No_V to convert
          ! to normalised units
          SwVrad = 10**SwVrad * 1e-2 * Si2No_V(UnitU_)

          ! Stellar wind velocities are purely radial
          State_VGB(Ux_:Uz_,i,j,k,iBlock) = &
               SwVrad * [ (Rorbit + x), y, z ] / rFrameStar

          ! Add phi component to account for the stellar wind being radial in
          ! the inertial (stellar) frame, but it has an azimuthal component
          ! in the rotating (planet) frame because of the orbiting planet
          State_VGB(Ux_,i,j,k,iBlock) = State_VGB(Ux_,i,j,k,iBlock) &
               + y * OmegaOrbit
          State_VGB(Uy_,i,j,k,iBlock) = State_VGB(Uy_,i,j,k,iBlock) &
               - (Rorbit + x) * OmegaOrbit

          ! Convert velocities to momenta for the total fluid
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
               State_VGB(Ux_:Uz_, i,j,k,iBlock) * State_VGB(Rho_,i,j,k,iBlock)

          ! Density and pressure of incoming stellar wind (spherical symmetry)
          State_VGB(Rho_,i,j,k,iBlock) = &
               MdotSw / (4.0*cPi * rFrameStar**2 * SwVrad)

          ! Set no neutrals in stellar wind
          State_VGB(HpRho_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
      	  State_VGB(H1Rho_,i,j,k,iBlock) = 0.0

      	  ! Passive scalar for stellar wind > 0
          State_VGB(PaSc_,i,j,k,iBlock) =pScalarSw*State_VGB(Rho_,i,j,k,iBlock)

          ! Isothermal pressure
          NumDens = 2.0*State_VGB(HpRho_,i,j,k,iBlock) / MassFluid_I(1) &
               + State_VGB(H1Rho_,i,j,k,iBlock) / MassNeuSpecies

          State_VGB(P_,i,j,k,iBlock) = NumDens * TempSwSi &
               * Si2No_V(UnitTemperature_)

          ! Possible stellar magnetic field, always converted from spherical
          ! to Cartesian to fill in hydro array
          select case(TypeBfieldStar)
          case('none')
            ! Non-magnetic star
            State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0

          case('radial')
            ! Stellar magnetic field is radial
            ! Monopolar field strength decays with (Rstar/r)^2
            BradStar = B0Star / rFrameStarRstar**2

            State_VGB(Bx_:Bz_,i,j,k,iBlock) = BradStar &
                 * [ (Rorbit + x), y, z ] / rFrameStar

          case('dipole')
            ! Stellar magnetic field is dipolar, no tilt
            CosTheta = z / rFrameStar
            SinTheta = sqrt((Rorbit + x)**2 + y**2) / rFrameStar
            CosPhi   = (Rorbit + x) / sqrt((Rorbit + x)**2 + y**2)
            SinPhi   = y / sqrt((Rorbit + x)**2 + y**2)

            ! Dipolar field strength decays with (Rstar/r)^3
            BradStar = B0Star * CosTheta / rFrameStarRstar**3
            BthetaStar  = 0.5*B0Star * SinTheta / rFrameStarRstar**3

            State_VGB(Bx_,i,j,k,iBlock) = SinTheta * CosPhi * BradStar &
                 + CosTheta * CosPhi * BthetaStar

            State_VGB(By_,i,j,k,iBlock) = SinTheta * SinPhi * BradStar &
                 + CosTheta * SinPhi * BthetaStar

            State_VGB(Bz_,i,j,k,iBlock) = CosTheta * BradStar &
                 - SinTheta * BthetaStar
          end select ! TypeBfieldStar
        enddo; enddo; enddo
      end select ! Iside
    end select ! BC type

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_set_cell_boundary

!==============================================================================
! Calculates the initial conditions of the grid for the planetary wind outflow.
! After a restart, this subroutine will not be accessed.
!==============================================================================
  subroutine user_set_ics(iBlock)

    use ModVarIndexes, ONLY: Rho_, HpRho_, H1Rho_, PaSc_, RhoUx_, RhoUz_, P_
    use ModPhysics,    ONLY: rBody, BodyRho_I, BodyP_I, Gamma
    use ModGeometry,   ONLY: Xyz_DGB, true_cell, r_BLK
    use ModAdvance,    ONLY: State_VGB
    use BATL_lib,      ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    implicit none

    integer, intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k
    real    :: x0, y0, z0, rp, Vrad

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI

      x0 = Xyz_DGB(x_,i,j,k,iBlock)
      y0 = Xyz_DGB(y_,i,j,k,iBlock)
      z0 = Xyz_DGB(z_,i,j,k,iBlock)
      rp = r_BLK(i,j,k,iBlock)

      if (.not.true_cell(i,j,k,iBlock) .or. rp < rBody) CYCLE

      ! IC of the planetary outflow (beta profile fit in python code by SC)
      Vrad = VinfPw &
           * (1.0 - (1.0 - (Vrad0Pw/VinfPw)**BetaIndex) * rBody/rp)**BetaIndex

      ! Mass density for each fluid from mass conservation
      State_VGB(HpRho_,i,j,k,iBlock) = (BodyRho_I(1) - BodyRhoNeuSpecies) &
            * Vrad0Pw/Vrad * (rBody/rp)**2

      State_VGB(H1Rho_,i,j,k,iBlock) = BodyRhoNeuSpecies * Vrad0Pw/Vrad &
           * (rBody/rp)**2

      State_VGB(Rho_,i,j,k,iBlock) = State_VGB(HpRho_,i,j,k,iBlock) &
           + State_VGB(H1Rho_,i,j,k,iBlock)

      ! Passive scalar for planetary outflow < 0
      State_VGB(PaSc_,i,j,k,iBlock) = pScalarPw * State_VGB(Rho_,i,j,k,iBlock)

      ! Thermal gas pressure for fluid is a polytrope
      State_VGB(P_,i,j,k,iBlock) = &
           BodyP_I(1) * (State_VGB(Rho_,i,j,k,iBlock)/BodyRho_I(1))**Gamma

      ! Momentum of fluid (radial component, spherically symmetric)
      State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
           State_VGB(Rho_,i,j,k,iBlock) * Vrad * [ x0, y0, z0 ] / rp
    enddo; enddo; enddo

    ! Initialising the uniform grid of densities for optical depth calculation
    if (iBlock == 1) call init_uniform_rho_grid

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_set_ics

!==============================================================================
! Calculate extra source terms from chemical reactions and external forces.
!==============================================================================
  subroutine user_calc_sources(iBlock)

    use ModVarIndexes,  ONLY: Rho_, HpRho_, H1Rho_, PaSc_, RhoUx_, RhoUy_, &
         RhoUz_, P_, Energy_, MassFluid_I
    use ModPhysics,     ONLY: No2Si_V, No2Io_V, Io2No_V, Si2No_V,    &
         UnitRho_, UnitU_, UnitX_, UnitN_, UnitTemperature_, UnitT_, &
         UnitEnergyDens_, mSun, rBody
    use ModGeometry,    ONLY: Xyz_DGB, true_cell, r_BLK
    use ModAdvance,     ONLY: State_VGB, Source_VC
    use ModInterpolate, ONLY: interpolate_scalar
    use ModConst,       ONLY: cKToEV
    use BATL_lib,       ONLY: nI, nJ, nK

    implicit none

    integer, intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k
    real :: x, y, z, rFrameStar
    real :: GravityStar_D(3), ForceCentrifugal_D(3), SrhoSpecies, SE
    real :: Temperature, NumDens, NumDensHp, NumDensH1, NumDensElectron
    real :: PhotoIonisationRate, ElectronImpactRate, RecombinationRate, Ufrac
    real :: HeatPhotoIon, CoolLyalpha, CoolElectronImpact, ExpMinTau

    logical :: DoTest, DoExtrapolate=.false.
    character(len=*), parameter :: NameSub = 'user_calc_source'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    !******************
    ! Atomic processes
    !******************
    do k = 1,nK; do j = 1,nJ; do i = 1,nI
      if (.not.true_cell(i,j,k,iBlock)) CYCLE

      ! Cell centred coordinates
      x = Xyz_DGB(x_,i,j,k,iBlock)
      y = Xyz_DGB(y_,i,j,k,iBlock)
      z = Xyz_DGB(z_,i,j,k,iBlock)

      ! Only calculate sources for regions with neutral atoms (Passive<0)
      if (State_VGB(PaSc_,i,j,k,iBlock) / State_VGB(Rho_,i,j,k,iBlock) &
           < 0.98*pScalarSw) then

        NumDensHp       = State_VGB(HpRho_,i,j,k,iBlock) / MassFluid_I(1)
        NumDensH1       = State_VGB(H1Rho_,i,j,k,iBlock) / MassNeuSpecies
        NumDensElectron = NumDensHp
        NumDens         = NumDensHp + NumDensH1 + NumDensElectron

        ! To find attenuation in a particular cell, interpolate exp(-tau) from
        ! uniform grid to the BATSRUS AMR grid
        if (r_BLK(i,j,k,iBlock) <= rBody) then
          ExpMinTau = 0.0
        else
          ExpMinTau = interpolate_scalar(array_3d_Exp_mTau(:,:,:),           &
                                        3, [ 1, 1, 1 ],                      &
                                        [ size_3d_i, size_3d_j, size_3d_k ], &
                                        [ x, y, z ],                         &
                                        x_array(:), y_array(:), z_array(:),  &
                                        DoExtrapolate)
        endif

        ! Plasma temperature in Kelvin
        Temperature =  State_VGB(P_,i,j,k,iBlock) / NumDens &
             * No2Si_V(UnitTemperature_)

        ! 1. Photo-ionisation rate by radiation: H1 + hv -> H+ + e- + Eheat
        PhotoIonisationRate = PhotoIonisationCoef * ExpMinTau

        ! 1.1 Photo-ionising heat carried away by the electron
        HeatPhotoIon = HeatPhotoCoef * NumDensH1 * ExpMinTau

        ! 2. Electron-impact ionisation: H1 + e- -> H+ + 2e-
        ! Fit from Voronov (1997), eq. 1 + table 1
        Ufrac = 13.6 / (cKToEV * Temperature)
        ElectronImpactRate = ElectronImpactCoef * NumDensElectron &
             * Ufrac**0.39 / (0.232 + Ufrac) * exp(-Ufrac)

        ! 2.1 Cooling from ejected bound electron carrying away heat
        CoolElectronImpact = CoolElectronImpactCoef * NumDensElectron &
             * NumDensH1 * Ufrac**0.39 / (0.232 + Ufrac) * exp(-Ufrac)

        ! 3. Radiative (ion-electron) recombination: H+ + e- -> H1 + hv
        ! Fit from Storey & Hummer (1995), table 1
        RecombinationRate = RecombinationCoef * NumDensElectron &
             * (1.0e4/Temperature)**0.9

        ! 4. Lyman alpha cooling from Black (1981), table 3
        CoolLyalpha = CoolLyalphaCoef * NumDensElectron * NumDensH1 &
             * exp(-118348.0/Temperature)

        ! 5. Update hydro equations

        ! Rate per volume at which H+ mass is produced from H1 mass
        SrhoSpecies = (PhotoIonisationRate + ElectronImpactRate)  &
             * State_VGB(H1Rho_,i,j,k,iBlock) - RecombinationRate &
             * State_VGB(HpRho_,i,j,k,iBlock)

        Source_VC(HpRho_,i,j,k) = Source_VC(HpRho_,i,j,k) + SrhoSpecies
        Source_VC(H1Rho_,i,j,k) = Source_VC(H1Rho_,i,j,k) - SrhoSpecies

        ! Rate per volume at which total energy is added to the fluid
        SE = HeatPhotoIon - CoolLyalpha - CoolElectronImpact

        Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + SE

        ! Store for printing variables
        PlotRateRec_CB(i,j,k,iBlock)      = RecombinationRate
        PlotRateCollIon_CB(i,j,k,iBlock)  = ElectronImpactRate
        PlotRatePhotoIon_CB(i,j,k,iBlock) = PhotoIonisationRate
        PlotHeat_CB(i,j,k,iBlock)         = HeatPhotoIon
        PlotCoolLya_CB(i,j,k,iBlock)      = CoolLyalpha
        PlotCoolColl_CB(i,j,k,iBlock)     = CoolElectronImpact
        PlotExp_mTau_CB(i,j,k,iBlock)     = ExpMinTau
      else
        ! This is the stellar wind part, when passive scalar is >0 [or >980]
        PlotRateRec_CB(i,j,k,iBlock)      = 0.0
        PlotRateCollIon_CB(i,j,k,iBlock)  = 0.0
        PlotRatePhotoIon_CB(i,j,k,iBlock) = 0.0
        PlotHeat_CB(i,j,k,iBlock)         = 0.0
        PlotCoolLya_CB(i,j,k,iBlock)      = 0.0
        PlotCoolColl_CB(i,j,k,iBlock)     = 0.0
        PlotExp_mTau_CB(i,j,k,iBlock)     = 1.0
      endif
    enddo; enddo; enddo

    !**************************************************
    ! Tidal gravity (star gravity + centrifugal force)
    !**************************************************
    if (UseTidal) then
      ! NOTE: for notation convenience forces are computed per unit mass
      do k = 1,nK; do j = 1,nJ; do i = 1,nI
        if (.not.true_cell(i,j,k,iBlock)) CYCLE

        x = Xyz_DGB(x_,i,j,k,iBlock)
        y = Xyz_DGB(y_,i,j,k,iBlock)
        z = Xyz_DGB(z_,i,j,k,iBlock)
        rFrameStar = sqrt( (Rorbit + x)**2 + y**2 + z**2 )

        ! Stellar gravity: -GMstar/R^2, acts in -r direction
        GravityStar_D = -gStar * [ (Rorbit + x), y, z ] / rFrameStar**3

        ! Update momenta and total energy
        Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
             + State_VGB(Rho_,i,j,k,iBlock) * GravityStar_D

        Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
             + sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) * GravityStar_D)

        ! Centrifugal force: -rho * Omega x (Omega x r), in +r direction
        ! NOTE: has only x,y component if rotation axis is along z
        ForceCentrifugal_D = OmegaOrbit**2 * [ (Rorbit + x), y, 0.0 ]

        ! Update momenta and total energy
        Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
             + State_VGB(Rho_,i,j,k,iBlock) * ForceCentrifugal_D

        Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
             + sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) * ForceCentrifugal_D)
      enddo; enddo; enddo
    endif

    !**************************************
    ! Coriolis force: -2*rho * (Omega x v)
    !**************************************
    if (UseCoriolis) then
      ! NOTE: has only x,y component if rotation axis is along z
      !       AND does not perform fictitious work
      do k = 1,nK; do j = 1,nJ; do i = 1,nI
        if (.not.true_cell(i,j,k,iBlock)) CYCLE

        Source_VC(RhoUx_,i,j,k) = Source_VC(RhoUx_,i,j,k) &
             + 2.0*OmegaOrbit * State_VGB(RhoUy_,i,j,k,iBlock)

        Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
             - 2.0*OmegaOrbit * State_VGB(RhoUx_,i,j,k,iBlock)
      enddo; enddo; enddo
    endif

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_calc_sources

!==============================================================================
! Perform hydrodynamic update and set hydrogen mass density on radiation grid.
!==============================================================================
  subroutine user_update_states(iBlock)

    use ModUpdateState, ONLY: update_state_normal

    implicit none

    integer, intent(in) :: iBlock

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Update hydrodynamic state on each block
    call update_state_normal(iBlock)

    ! Interpolate updated hydrogen mass density onto radiation grid
    call create_uniform_3d_grid_rho(iBlock)

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_update_states

!==============================================================================
! This subroutine used for plotting new variables with #SAVEPLOT beyond the
! default MHD variables. Variables should be allocated at startup.
! For example, to add PARAM.in:
! #SAVEPLOT
!   1             	nPlotFiles
!   z=0 VAR tec   	StringPlot
!   1            	  DnSavePlot
!   -1            	DtSavePlot
!   {MHD} tau heat cool rion rrec   NameVars
!   {default}       NamePars
!==============================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, PlotVar_G, &
                               PlotVarBody, UsePlotVarBody, NameTecVar,   &
                               NameTecUnit, NameIdlUnit, IsFound)

    use ModPhysics, ONLY: No2Io_V, UnitRho_, UnitT_, UnitEnergyDens_
    use BATL_lib,   ONLY: nI, nJ, nK

    implicit none

    integer,          intent(in)    :: iBlock
    character(len=*), intent(in)    :: NameVar
    logical,          intent(in)    :: IsDimensional
    real,             intent(out)   :: PlotVar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,             intent(out)   :: PlotVarBody
    logical,          intent(out)   :: UsePlotVarBody
    character(len=*), intent(inout) :: NameTecVar
    character(len=*), intent(inout) :: NameTecUnit
    character(len=*), intent(inout) :: NameIdlUnit
    logical,          intent(out)   :: IsFound

    ! Local variables
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    IsFound = .true.

    select case(NameVar)
    case('tau')
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotExp_mTau_CB(:,:,:,iBlock)
      !NameTecUnit = ' '
      !NameIdlUnit = ' '
      NameTecVar = 'exp(-t)'

    case('heat') ! Heating rate
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotHeat_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
      NameIdlUnit = 'erg/cm^3/s'
      NameTecVar  = 'Heat'

    case('clya') ! Cooling rate Ly alpha
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotCoolLya_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
      NameIdlUnit = 'erg/cm^3/s'
      NameTecVar  = 'Cool_Lya'

    case('ccoll') ! Cooling rate collisions
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotCoolColl_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
      NameIdlUnit = 'erg/cm^3/s'
      NameTecVar  = 'Cool_coll'

    case('rion')
      PlotVar_G(1:nI,1:nJ,1:nK) = &
           PlotRatePhotoIon_CB(:,:,:,iBlock) + PlotRateCollIon_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
      NameTecVar  = 'R_ion'

    case('rphoto')
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotRatePhotoIon_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
      NameTecVar  = 'R_phion'

    case('rcoll')
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotRateCollIon_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
      NameTecVar  = 'R_collion'

    case('rrec')
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotRateRec_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
      NameTecVar  = 'R_rec'

    case default
      IsFound = .false.
      call stop_mpi(NameSub//': unimplemented variable='//NameVar)
    end select

    UsePlotVarBody = .true.
    PlotVarBody    = 0.0

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_set_plot_var

!==============================================================================
! Perform additional manipulations at specific instances of the simulation.
!==============================================================================
  subroutine user_action(NameAction)

    use ModMpi,   ONLY: MPI_REAL, MPI_SUM, mpi_reduce_real_array
    use BATL_lib, ONLY: iProc, nProc, nBlock, iComm, nI, nJ, nK, MaxBlock

    implicit none

    character(len=*), intent(in) :: NameAction

    ! Local variables
    integer :: iBlock, iError

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_action'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    select case(NameAction)
    case('initialize module')
      ! Allocate arrays for including in tecplot output
      if (.not.allocated(PlotExp_mTau_CB)) then
        allocate(PlotExp_mTau_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotHeat_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotCoolLya_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotCoolColl_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotRatePhotoIon_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotRateCollIon_CB(1:nI,1:nJ,1:nK,MaxBlock))
        allocate(PlotRateRec_CB(1:nI,1:nJ,1:nK,MaxBlock))
      endif

    case('clean module')
      ! At end give memory back to computer - Fortran 90 practise
      if (allocated(PlotExp_mTau_CB))     deallocate(PlotExp_mTau_CB)
      if (allocated(PlotHeat_CB))         deallocate(PlotHeat_CB)
      if (allocated(PlotCoolLya_CB))      deallocate(PlotCoolLya_CB)
      if (allocated(PlotCoolColl_CB))     deallocate(PlotCoolColl_CB)
      if (allocated(PlotRatePhotoIon_CB)) deallocate(PlotRatePhotoIon_CB)
      if (allocated(PlotRateCollIon_CB))  deallocate(PlotRateCollIon_CB)
      if (allocated(PlotRateRec_CB))      deallocate(PlotRateRec_CB)

    case('initial condition done')
      ! Initialising grid of densities after startup OR restart
      do iBlock = 1,nBlock
        call create_uniform_3d_grid_rho(iBlock)
      enddo

      if (nProc > 1) then
        call mpi_reduce_real_array(array_3d_dens, size(array_3d_dens), &
                                   MPI_SUM, 0, iComm, iError)

        if (iProc == 0) call create_uniform_3d_grid_tau

        call MPI_BCAST(array_3d_Exp_mTau, size(array_3d_Exp_mTau), &
                       MPI_REAL, 0, iComm, iError)

        ! Reseting for next timestep
        array_3d_dens = 0.0
      endif

    case('write progress')
      ! Output some extra model setup information
      call print_model_info
    end select

    call test_stop(NameSub, DoTest)

  end subroutine user_action

!==============================================================================
! This routine can add extra variables of interest to the .log file created
! with the #SAVELOGFILE header in the PARAM.in
!==============================================================================
  subroutine user_get_log_var(VarValue, NameVar, Radius)

    use ModMain,            ONLY: nI, nJ, nK, nBlock, Unused_B, UseB0
    use ModAdvance,         ONLY: State_VGB, tmp1_BLK
    use ModB0,              ONLY: B0_DGB
    use ModVarIndexes,      ONLY: Rho_, Bx_, By_, Bz_, Ux_, Uy_, Uz_
    use ModPhysics,         ONLY: No2Si_V, UnitX_, UnitU_, UnitMass_, UnitB_
    use ModGeometry,        ONLY: r_BLK
    use ModWriteLogSatFile, ONLY: calc_sphere
    use ModNumConst,        ONLY: cPi
    use BATL_lib,           ONLY: Xyz_DGB

    real, intent(out)            :: VarValue
    character(len=*), intent(in) :: NameVar
    real, intent(in), optional   :: Radius

    ! Local variables
    integer :: iBlock, i, j, k
    real    :: FullB_DG(3,0:nI+1,0:nJ+1,0:nK+1), FullV_DG(3,0:nI+1,0:nJ+1,0:nK+1)
    real    :: xloc, yloc, rplane, Bxloc, Byloc, Brloc, Bploc
    real    :: rholoc, vxloc, vyloc, vrloc, vploc

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_get_log_var'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    VarValue = 0.0

    select case(NameVar)
    case("ubflx") ! Unsigned magnetic flux

      do iBlock = 1,nBlock

        if (Unused_B(iBlock)) CYCLE

        FullB_DG = State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBlock)

        if (UseB0) then
          FullB_DG = FullB_DG + B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBlock)
        endif

        ! Compute Br
        do k = 0,nK+1; do j = 0,nJ+1; do i = 0,nI+1
          tmp1_BLK(i,j,k,iBlock) = &
               sum( FullB_DG(:,i,j,k) * Xyz_DGB(:,i,j,k,iBlock)) &
               / r_BLK(i,j,k,iBlock )
        enddo; enddo; enddo
      enddo

      ! To get the unsigned Br, take absolute value
      VarValue = calc_sphere('integrate', 360, Radius, abs(tmp1_BLK))

      ! Handle conversion to IO units; avoids working with UserUnit_V
      VarValue = VarValue * ( No2Si_V(UnitB_)*No2Si_V(UnitX_)**2 )
      ! print*, "Logfile check:", VarValue, NameVar, Radius

    case('pbflx') ! Magnetic torque planet

      do iBlock = 1,nBlock

        if (Unused_B(iBlock)) CYCLE

        FullB_DG = State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBlock)

        if (UseB0) then
          FullB_DG = FullB_DG + B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBlock)
        endif

        ! Compute Br and Bphi to get magnetic torque Jdot_mag
        do k = 0,nK+1; do j = 0,nJ+1; do i = 0,nI+1
          Brloc = sum( FullB_DG(:,i,j,k) * Xyz_DGB(:,i,j,k,iBlock)) &
               / r_BLK(i,j,k,iBlock )

          xloc =  Xyz_DGB(x_,i,j,k,iBlock)
          yloc =  Xyz_DGB(y_,i,j,k,iBlock)
          rplane = sqrt( xloc**2.0 + yloc**2.0 )
          Bxloc = FullB_DG(Bx_,i,j,k)
          Byloc = FullB_DG(By_,i,j,k)

          Bploc = (-Bxloc*yloc + Byloc*xloc) / rplane

          tmp1_BLK(i,j,k,iBlock) = -rplane * Brloc * Bploc/(4.0*cPi)
        enddo; enddo; enddo
      enddo

      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)

      ! Handle conversion to IO units; avoids working with UserUnit_V
      VarValue = VarValue * ( No2Si_V(UnitMass_)*No2Si_V(UnitU_)**2.0 )

    case('pwflx') ! Wind material torque planet

      do iBlock = 1,nBlock

        if (Unused_B(iBlock)) CYCLE

        FullV_DG = State_VGB(Ux_:Uz_,0:nI+1,0:nJ+1,0:nK+1,iBlock)

        ! Compute vr and vphi to get wind torque Jdot_wind
        do k = 0,nK+1; do j = 0,nJ+1; do i = 0,nI+1
          rholoc = State_VGB(Rho_,i,j,k,iBlock)

          vrloc = sum( FullV_DG(:,i,j,k) * Xyz_DGB(:,i,j,k,iBlock)) &
               / r_BLK(i,j,k,iBlock )

          xloc =  Xyz_DGB(x_,i,j,k,iBlock)
          yloc =  Xyz_DGB(y_,i,j,k,iBlock)
          rplane = sqrt( xloc**2.0 + yloc**2.0 )
          vxloc = FullV_DG(Ux_,i,j,k)
          vyloc = FullV_DG(Uy_,i,j,k)

          vploc = (-vxloc*yloc + vyloc*xloc) / rplane

          tmp1_BLK(i,j,k,iBlock) = rplane * rholoc * vrloc * vploc
        enddo; enddo; enddo
      enddo

      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)

      ! Handle conversion to IO units; avoids working with UserUnit_V
      VarValue = VarValue * ( No2Si_V(UnitMass_)*No2Si_V(UnitU_)**2.0 )
      !print*, "Logfile check:", VarValue, NameVar, Radius

    case default
       call stop_mpi('Unknown user logvar='//NameVar)
    end select

    call test_stop(NameSub, DoTest)

  end subroutine user_get_log_var

!==============================================================================
! Read stellar wind velocity for injection into the grid. The velocity profile
! is created in 1D models of polytropic winds, whose results are fitted by a
! 8th degree polynomial.
!==============================================================================
  subroutine read_stellar_wind_coeff(OrderCoef, SwCoef_I)

    use ModIoUnit, ONLY: io_unit_new

    implicit none

    integer, intent(in)    :: OrderCoef
    real,    intent(inout) :: SwCoef_I(OrderCoef)

    ! Local variables
    integer            :: iUnit, i, iError, line, pos
    character(len=100) :: buffer, label
    logical            :: exist=.false.

    ! Return an unused unit number for extended use
    iUnit = io_unit_new()

    open(iUnit, file='sw-velocity.inp', status='old')
    line = 0

    ! iError is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected. iError is zero otherwise.
    iError = 0
    do while (iError == 0)
      read(iUnit,*, iostat=iError) buffer
      line = line + 1

      ! Find the first instance of whitespace. Split label and comment.
      pos = scan(buffer, ' ')
      label = buffer(1:pos)

      !write (*,*) line, pos, buffer

      if (label == NameStar) then
        read(iUnit,*) (SwCoef_I(i), i=1,OrderCoef)
        exist  = .true.
        iError = -1
      endif
    enddo
    close(iUnit)

    if (.not.exist) then
      call stop_mpi('NameStar does not exist in sw-velocity.inp. &
           Make sure NameStar is writen in lower case!.')
    endif

  end subroutine read_stellar_wind_coeff

!==============================================================================
! Setup a uniform base grid for performing the ray tracing to compute the
! optical depth integral. To be used within the source term calculation.
!==============================================================================
  subroutine make_uniform_rtgrid

    use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2

    implicit none

    ! Local variables
    integer :: itemp, points_left, points_right
    real    :: res_mid, res_right, res_left
    real    :: length_mid, length_left, length_right

    ! Grid is in normalised units, but since Io2No(x_)=1, it does not matter

    ! === x-axis ===
    length_mid   = 2.0 - (-2.0)
    length_left  = -2.0 - x1
    length_right = x2 - 2.0

    points_right = int(length_right*(size_3d_i+2-41)/(length_left+length_right))
    points_left  = size_3d_i+2-41 - points_right

    ! 40 deltas within +/- 2, the same in all three axis
    res_mid   = length_mid/(41-1)
    res_left  = length_left/(points_left-1)
    res_right = length_right/(points_right-1)

    ! Generate x-grid
    x_array(1:points_left) = &
        (/ ( x1 + res_left*(itemp-1.0), itemp=1,points_left) /)

    x_array(points_left:points_left+41) = &
        (/ ( -2.0 + res_mid*(itemp-points_left), &
        itemp=points_left,points_left+41) /)

    x_array(points_left+41:size_3d_i) = &
        (/ ( 2.0 + res_right*(itemp -(points_left+40)), &
        itemp=points_left+41,size_3d_i) /)

    ! === y-axis ====
    length_left  = -2.0 - y1
    length_right = y2 - 2.0

    points_right = int(length_right*(size_3d_j+2-41)/(length_left+length_right))
    points_left  = size_3d_j+2-41 - points_right

    res_left  = length_left/(points_left-1)
    res_right = length_right/(points_right-1)

    y_array(1:points_left) = &
        (/ ( y1 + res_left*(itemp-1.0), itemp=1,points_left) /)

    y_array(points_left:points_left+41) = &
        (/ ( -2. + res_mid*(itemp-points_left), &
        itemp=points_left,points_left+41) /)

    y_array(points_left+41:size_3d_j) = &
        (/ ( 2. + res_right*(itemp -(points_left+40)), &
        itemp=points_left+41,size_3d_j) /)

    ! === z-axis ===
    length_left  = -2.0 - z1
    length_right = z2 - 2.0

    points_right = int(length_right*(size_3d_k+2-41)/(length_left+length_right))
    points_left  = size_3d_k+2-41 - points_right

    res_left  = length_left/(points_left-1)
    res_right = length_right/(points_right-1)

    z_array(1:points_left) = &
        (/ ( z1 + res_left*(itemp-1.0), itemp=1,points_left) /)

    z_array(points_left:points_left+41) = &
        (/ ( -2.0 + res_mid*(itemp-points_left), &
        itemp=points_left,points_left+41) /)

    z_array(points_left+41:size_3d_k) = &
        (/ ( 2.0 + res_right*(itemp -(points_left+40)), &
        itemp=points_left+41,size_3d_k) /)

  end subroutine make_uniform_rtgrid

!==============================================================================
! Initialisation of the uniform grid of density and optical depth that will be
! used in the first instance of user_calc_sources.
!==============================================================================
  subroutine init_uniform_rho_grid

    use ModPhysics, ONLY: BodyRhoSpecies_I

    implicit none

    ! Local variables
    integer :: i, j, k
    real    :: r_vel, radius

    do k = 1,size_3d_k; do j = 1,size_3d_j; do i = 1,size_3d_i
      radius = sqrt(x_array(i)**2 + y_array(j)**2 + z_array(k)**2)

      if (radius <= 1.0) then
        array_3d_dens(i,j,k) = BodyRhoNeuSpecies
        cycle
      end if

      r_vel = (Vrad0Pw / (VinfPw*(1.0 - 1.0/radius)**BetaIndex)) &
          * (1.0/radius)**2
     array_3d_dens(i,j,k) = r_vel * BodyRhoNeuSpecies
   enddo; enddo; enddo

    ! This will create my initial array of tau
    call create_uniform_3d_grid_tau()

    ! Reset for next timestep
    array_3d_dens(1:size_3d_i,1:size_3d_j,1:size_3d_k) = 0.0

  end subroutine init_uniform_rho_grid

!==============================================================================
! Interpolate the values of density from the AMR grid to the uniform grid. The
! computation is performed on the MHD grid spanned by each individual block.
! The densities for the entire grid are computed in ModBatsrusMethods.f90.
!==============================================================================
  subroutine create_uniform_3d_grid_rho(iBlock)

    use ModVarIndexes,  ONLY: H1Rho_
    use ModPhysics,     ONLY: rBody
    use ModGeometry,    ONLY: Xyz_DGB
    use ModAdvance,     ONLY: State_VGB
    use ModInterpolate, ONLY: interpolate_scalar
    use BATL_lib,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         CoordMin_DB, CoordMax_DB

    implicit none

    integer , intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k, iMinX, iMaxX, iMinY, iMaxY, iMinZ, iMaxZ
    logical :: DoExtrapolate=.true.

    ! Instead of looping through all i,j,k only loop through the indices
    ! that will be interpolated in this block

    ! Minimum index where x_array > x_min_iblock
    iMinX = minloc(x_array, 1, mask = x_array > CoordMin_DB(x_, iBlock))

    ! Maximum index where x_array <= x_max_iblock
    iMaxX = maxloc(x_array, 1, mask = x_array <= CoordMax_DB(x_, iBlock))

    ! Same procedure for y,z direction
    iMinY = minloc(y_array, 1, mask = y_array > CoordMin_DB(y_, iBlock))
    iMaxY = maxloc(y_array, 1, mask = y_array <= CoordMax_DB(y_, iBlock))
    iMinZ = minloc(z_array, 1, mask = z_array > CoordMin_DB(z_, iBlock))
    iMaxZ = maxloc(z_array, 1, mask = z_array <= CoordMax_DB(z_, iBlock))

    do k = iMinZ,iMaxZ; do j = iMinY,iMaxY; do i = iMinX,iMaxX

      ! No interpolation inside planet
      if (sqrt(x_array(i)**2 + y_array(j)**2 + z_array(k)**2) <= rBody) then
        array_3d_dens(i,j,k) = BodyRhoNeuSpecies
        CYCLE
      endif

      ! There are two ways of interpolating: considering the ghost cells or not
      array_3d_dens(i,j,k) =                                                 &
        interpolate_scalar(State_VGB(H1Rho_,:,:,:,iBlock),                   &
      	                   3, (/MinI, MinJ, MinK/), (/MaxI, MaxJ, MaxK/),    &
                           (/x_array(i), y_array(j), z_array(k)/),           &
                           Xyz_DGB(x_,MinI:MaxI,1        ,1        ,iBlock), &
                           Xyz_DGB(y_,1        ,MinJ:MaxJ,1        ,iBlock), &
                           Xyz_DGB(z_,1        ,1        ,MinK:MaxK,iBlock), &
                           DoExtrapolate)
    enddo; enddo; enddo ! loop in ijk

  end subroutine create_uniform_3d_grid_rho

!==============================================================================
! Calculate the sum over n(Hydrogen) to get optical depth. This subroutine is
! called from ModBatsrusMethods.f90 after message passing to iProc=0. Thus,
! the value of array_3d_Exp_mTau is only complete for iProc=0.
!==============================================================================
  subroutine create_uniform_3d_grid_tau

    use ModVarIndexes, ONLY: H1Rho_
    use ModPhysics,    ONLY: No2Si_V, No2Io_V, UnitX_, UnitN_

    implicit none

    ! Local variables
    integer :: idum
    real    :: cell_size, array_3d_tau(size_3d_i, size_3d_j, size_3d_k), cTau

    ! Constant optical depth (factor 100 is to go from m to cm)
    cTau = SigmaNu0Cgs * No2Io_V(UnitN_) * 1e2 * No2Si_V(UnitX_)

    ! Summing density in all cells along x to get optical depth
    ! This is zero in the first point, at the left edge of the grid
    array_3d_tau(1,:,:) = array_3d_dens(1,:,:) * (x_array(2) - x_array(1))

    do idum = 2,size_3d_i
      cell_size = x_array(idum) - x_array(idum-1)
      array_3d_tau(idum,:,:) = array_3d_dens(idum,:,:) * cell_size &
           + array_3d_tau(idum-1,:,:)
    enddo

    ! Normalisation of optical depth
    array_3d_tau      = array_3d_tau * cTau/MassNeuSpecies
    array_3d_Exp_mTau = exp(-array_3d_tau)

  end subroutine create_uniform_3d_grid_tau

!==============================================================================
! Print some extra information of interest to the screen when called for.
!==============================================================================
  subroutine print_model_info

    use ModMain,        ONLY: IsStandAlone, UseRotatingFrame, UseRotatingBC, &
         UseRayTrace
    use ModVarIndexes,  ONLY: HpRho_, H1Rho_, ScalarFirst_, ScalarLast_, &
         PaSc_, MassFluid_I
    use ModPhysics,     ONLY: No2Si_V, No2Io_V, UnitX_, UnitRho_, UnitU_, &
         UnitN_, UnitP_, UnitB_, UnitMass_, UnitT_, UnitRhoU_, UnitJ_,    &
         UnitDivB_, UnitAngle_, UnitTemperature_, UnitEnergyDens_,        &
         UnitElectric_, UnitPoynting_, BodyRhoSpecies_I, BodyRho_I,       &
         BodyNDim_I, BodyP_I, rPlanetSi, RotPeriodSi, OmegaBody
    use ModGeometry,    ONLY: TypeGeometry
    use ModAdvance,     ONLY: UseNonConservative
    use ModPlanetConst, ONLY: mPlanet_I, rPlanet_I, Jupiter_
    use CON_planet,     ONLY: MassPlanet
    use ModConst,       ONLY: cSecondPerDay
    use ModIO,          ONLY: restart

    if (UseNonConservative) &
         call stop_mpi('Only runs with conservative energy formulation!')

    if (Mstar == 0.0) &
         call stop_mpi('Missing #USR_SYSTEMPROPERTIES in PARAM.in. &
         Maybe due to a restart of the code?')

    if (FxuvCgs == 0.0) call stop_mpi('Missing #USR_FXUV in PARAM.in. &
         Maybe due to a restart of the code?')

    if (UseStellarWind .and. NameStar == 'empty') &
         call stop_mpi('Missing #USR_INJECT_SW in PARAM.in. &
         Maybe due to a restart of the code?')

    write(*,'(2X, A)') ' **** Planetary properties:'
    write(*,'(5X,A23,ES10.3,A10)') 'Rplanet = ', rPlanetSi, ' m'
    write(*,'(5X,A23,F10.3,A10)')  &
         'Rplanet = ', rPlanetSi / rPlanet_I(Jupiter_), ' Rjup'
    write(*,'(5X,A23,ES10.3,A10)') 'Mplanet = ', MassPlanet, ' kg'
    write(*,'(5X,A23,F10.3,A10)')  &
         'Mplanet = ', MassPlanet / mPlanet_I(Jupiter_), ' Mjup'
    write(*,'(5X,A23,ES10.3,A10)') &
         'Planet Rotation Period = ', RotPeriodSi/cSecondPerDay, ' days'
    write(*,'(5X,A23,ES10.3)') &
         'Planet rotation rate = ', OmegaBody*No2Si_V(UnitT_)
    write(*,'(5X,A23,f10.3,A10)')  'Rorbit = ', Rorbit_au, ' au'
    write(*,'(5X,A23,es10.3,A10)') 'OmegaOrbit = ', OmegaOrbit_dim,' 1/s'
    write(*,'(5X,A23,es10.3,A10)') 'PeriodOrbit = ', PeriodOrbitDays,' days'
    write(*,'(5X,A23,es10.3,A10)') 'OmegaOrbit = ', OmegaOrbit,'normalised'
    write(*,'(5X,A23,es10.3,A10)') 'Rorbit = ', Rorbit,' normalised'
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Stellar properties:'
    write(*,*) ' Star: ', NameStar
    write(*,'(5X,A23,f10.3,A10)') 'Mstar = ', Mstar, ' Msun'
    write(*,'(5X,A23,f10.3,A10)') 'Rstar = ', Rstar, ' Rsun'
    write(*,'(5X,A23,es10.3,A10)') 'Twind = ', TempSwSi, ' K'
    write(*,'(5X,A23,es10.3,A10)') 'Mdot = ', MdotSwMsunYr, ' Msun/yr'
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Radiation properties:'
    write(*,'(5X,A23,es10.3,A10)') 'FxuvCgs = ', FxuvCgs,' erg/cm2/s'
    write(*,'(5X,A23,es10.3,A10)') 'SigmaNu0 = ', SigmaNu0Cgs,' cm^2'
    write(*,'(5X,A23,es10.3,A10)') 'EpsilonNu = ', EpsilonNu,'  '
    write(*,'(5X,A23,es10.3,A10)') 'EphotEv = ', EphotEv,' eV'
    write(*,'(5X,A23,es10.3,A10)') 'EphotCgs = ', EphotCgs,' erg'
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Other fluid properties:'
    write(*,'(5X,A23,es10.3,A10)') 'BodyRho_I(fluid=1)  = ', BodyRho_I(1),' '
    write(*,'(5X,A23,es10.3,A10)') 'BodyP_I(fluid=1) = ', BodyP_I(1),' '

    write(*,'(5X,A23,i10)') 'size_3d_i = ', size_3d_i
    write(*,'(5X,A23,i10)') 'size_3d_j = ', size_3d_j
    write(*,'(5X,A23,i10)') 'size_3d_k = ', size_3d_k
    write(*,'(2X, A)') ' '
    ! NOTE: rPlanetSi, MassBodySi, and RotPeriodSi are read in #PLANET

    write(*,'(2X, A)') ' **** Normalisation variables:'
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitX_) = '  , No2Si_V(UnitX_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitU_) = '  , No2Si_V(UnitU_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitRho_) = ', No2Si_V(UnitRho_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitN_) = '  , No2Si_V(UnitN_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitP_) = '  , No2Si_V(UnitP_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitB_) = '  , No2Si_V(UnitB_)
    write(*,'(5X,A23,ES10.3)') &
         'No2Si_V(UnitTemperature_) = ', No2Si_V(UnitTemperature_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitMass_) = ', No2Si_V(UnitMass_)
    write(*,'(5X,A23,ES10.3)') &
         'No2Si_V(UnitEnergyDens_) = ', No2Si_V(UnitEnergyDens_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitT_) = ', No2Si_V(UnitT_)
    write(*,'(5X,A23,ES10.3)') &
         'No2Si_V(UnitPoynting_) = ', No2Si_V(UnitPoynting_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitRhoU_) = ', No2Si_V(UnitRhoU_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitJ_) = ', No2Si_V(UnitJ_)
    write(*,'(5X,A23,ES10.3)') &
         'No2Si_V(UnitElectric_) = ', No2Si_V(UnitElectric_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitDivB_) = ', No2Si_V(UnitDivB_)
    write(*,'(5X,A23,ES10.3)') 'No2Si_V(UnitAngle_) = ', No2Si_V(UnitAngle_)
    write(*,'(5X,A23,ES10.3)') &
         'No2Si_V(energy) = ', No2Si_V(UnitEnergydens_)*No2Si_V(UnitX_)**3
    write(*,*)
    write(*,'(5X,A23,ES10.3)') 'No2Io_V(UnitX_) = '  , No2Io_V(UnitX_)
    write(*,'(5X,A23,ES10.3)') 'No2Io_V(UnitU_) = '  , No2Io_V(UnitU_)
    write(*,'(5X,A23,ES10.3)') 'No2Io_V(UnitRho_) = ', No2Io_V(UnitRho_)
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Logical variables:'
    write(*,*) 'IsStandAlone       = ', IsStandAlone
    write(*,*) 'UseRotatingFrame   = ', UseRotatingFrame
    write(*,*) 'UseRotatingBC      = ', UseRotatingBC
    write(*,*) 'UseNonConservative = ', UseNonConservative
    write(*,*) 'Restart            = ', restart
    write(*,*) 'TypeGeometry       = ', TypeGeometry
    write(*,*) 'UseRaytrace        = ', UseRaytrace
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Beta profile fit for planetary wind IC:'
    write(*,'(5X,A23,f10.3,A10)') &
         'Vrad0Pw = ', Vrad0Pw*No2Io_V(UnitU_),' km/s'
    write(*,'(5X,A23,F10.3,A10)') &
         'VinfPw = ', VinfPw*No2Io_V(UnitU_),' km/s'
    write(*,'(5X,A23,f10.3)') 'BetaIndex = ', BetaIndex

    if (UseStellarWind) then
      write(*,'(2X, A)') ' **** Polynomial coefficients injected stellar wind'
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(1)= ', SwCoef_I(1)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(2)= ', SwCoef_I(2)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(3)= ', SwCoef_I(3)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(4)= ', SwCoef_I(4)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(5)= ', SwCoef_I(5)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(6)= ', SwCoef_I(6)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(7)= ', SwCoef_I(7)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(8)= ', SwCoef_I(8)
      write(*,'(5X,A23,f10.3)') 'SwCoef_I(9)= ', SwCoef_I(9)
    endif

  end subroutine print_model_info

end module ModUser
