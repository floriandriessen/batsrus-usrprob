!\
! Multi-fluid MHD model for atmospheric escape in (magnetic) Hot Jupiters.
! The fluid consists of singly ionised hydrogen and neutral hydrogen with heat
! capacity ratio of 5/3. The free electrons are assumed to be in thermal
! equilibrium with the ions and global charge neutrality holds.
!
! Fluids interact amongst each other by means of hard-sphere collisions. Other
! physical processes included are species production and loss from ionisation
! (photo and collision) and radiative recombination, and charge exchange.
! The optical depth for the ionisation is computed via ray tracing on a uniform
! grid which is interpolated back on the AMR grid--see ModBatsrusMethods.
!
! HISTORY
!   Apr 2023 : init from two-fluid MHD code
!/
module ModUser

  use BATL_lib, ONLY: test_start, test_stop
  use ModUserEmpty,                               &
       IMPLEMENTED1  => user_read_inputs,         &
       IMPLEMENTED2  => user_io_units,            &
       IMPLEMENTED3  => user_normalization,       &
       IMPLEMENTED4  => user_init_session,        &
       IMPLEMENTED5  => user_set_face_boundary,   &
       IMPLEMENTED6  => user_set_cell_boundary,   &
       IMPLEMENTED7  => user_set_ics,             &
       IMPLEMENTED8  => user_calc_sources,        &
       IMPLEMENTED9  => user_set_plot_var,        &
       IMPLEMENTED10 => user_action,              &
       IMPLEMENTED11 => user_init_point_implicit, &
       IMPLEMENTED12 => user_update_states

  ! List of public methods
  include 'user_module.h'

  ! Make public in order to add to ModBatsrusMethods
  public :: create_uniform_3d_grid_tau

  real,             parameter :: VersionUserModule = 1.0
  character(len=*), parameter :: NameUserModule = &
       'Four-fluid MHD simulation for Hot Jupiter with star'

  ! Orbital properties planet
  real, public :: RorbitAu=0.0, OmegaOrbitSi
  real         :: Rorbit, OmegaOrbit, PeriodOrbitDays

  ! Star properties
  real :: Mstar=0.0, Rstar=0.0, gStar

  ! User switches to add physics to the problem
  logical :: UseCoriolis=.false., UseTidal=.false., UsePwChargeExchange=.false.

  ! Used for the IC planetary wind
  real :: Vrad0PwIo=0.0, VinfPwIo=0.0, BetaIndex=1.0, VinfPw, Vrad0Pw

  ! Used for radiation hydrodynamics
  real :: FxuvCgs=0.0, EphotCgs=0.0, SigmaNu0Cgs, EpsilonNu, EphotEv

  ! Used for printing new variables
  real, allocatable :: &
       PlotExp_mTau_CB(:,:,:,:), PlotHeat_CB(:,:,:,:),            &
       PlotCoolLya_CB(:,:,:,:), PlotCoolColl_CB(:,:,:,:),         &
       PlotRatePhotoIon_CB(:,:,:,:), PlotRateCollIon_CB(:,:,:,:), &
       PlotRateRec_CB(:,:,:,:), PlotChExcHpH1_CB(:,:,:,:)

  ! Radiation grid uses even numbers of points to avoid interpolation to fall
  ! in-between two blocks
  integer, parameter :: nRadGridX=200, nRadGridY=200, nRadGridZ=128
  real :: RadGridX_C(nRadGridX), RadGridY_C(nRadGridY), RadGridZ_C(nRadGridZ)

  ! Optical depth computation: public so can be added to modBatsrusMethods
  real, public :: RhoRadGrid_C(nRadGridX, nRadGridY, nRadGridZ)=0.0
  real, public :: ExpTauRadGrid_C(nRadGridX, nRadGridY, nRadGridZ)=1.0

  ! NEW Flo
  real :: SwVrad0, SwVinf, SwMdot

contains

!==============================================================================
! Read the user variables given in PARAM.in with #USERINPUT.
!==============================================================================
  subroutine user_read_inputs

    use ModReadParam, ONLY: read_line, read_command, read_var
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
        call read_var('UsePwChargeExchange', UsePwChargeExchange)

      ! Stellar mass and radius and planet orbital distance
      case('#USR_SYSTEMPROPERTIES')
        call read_var('Mstar', Mstar)
        call read_var('Rstar', Rstar)
        call read_var('RorbitAu', RorbitAu)

      ! Stellar flux (cgs) at orbit and interaction energy
      case('#USR_FXUV')
        call read_var('FxuvCgs', FxuvCgs)
        call read_var('EphotEv', EphotEv)

      ! Planetary atmosphere
      case('#USR_INIT_PL_ATM')
        call read_var('Vrad0Pw', Vrad0PwIo)
        call read_var('VinfPw', VinfPwIo)
        call read_var('BetaIndex', BetaIndex)

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
         Gamma_I, BodyTDim_I, BodyNDim_I
    use ModConst,      ONLY: cBoltzmann, cProtonMass

    implicit none

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_normalization'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Three independent normalisation units (in SI) - using planet quantities
    No2Si_V(UnitX_)   = rPlanetSi
    No2Si_V(UnitU_)   = sqrt(Gamma_I(2)*cBoltzmann*BodyTDim_I(2)/cProtonMass)
    No2Si_V(UnitRho_) = 1.0e6*cProtonMass &
         * (MassFluid_I(2)*BodyNDim_I(2) + MassFluid_I(4)*BodyNDim_I(4))

    call test_stop(NameSub, DoTest)

  end subroutine user_normalization

!==============================================================================
! Initialize new global variables to be used within a session. This routine is
! read at the first iteration and after a restart.
!==============================================================================
  subroutine user_init_session

    use ModPhysics,  ONLY: No2Si_V, Si2No_V, Io2No_V, No2Io_V, UnitX_, &
         UnitU_, UnitN_, UnitB_, UnitT_, UnitMass_, mSun, rSun
    use ModConst,    ONLY: cGravitation, cAU, cSecondPerDay, cSecondPerYear
    use ModNumConst, ONLY: cPi

    implicit none

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    !\
    ! Deriving other initial variables
    !/

    ! G*Mstar; m^3/s^2 -> normalised
    gStar = cGravitation * Mstar*mSun * Si2No_V(UnitU_)**2 * Si2No_V(UnitX_)

    ! Orbital distance; au -> normalised
    Rorbit = RorbitAu * cAU * Si2No_V(UnitX_)

    ! Orbit quantities: Omega in (rad/s) and Period in (days)
    OmegaOrbitSi    = sqrt(cGravitation * Mstar * mSun / (RorbitAu*cAU)**3)
    OmegaOrbit      = OmegaOrbitSi / Si2No_V(UnitT_)
    PeriodOrbitDays = 2.0*cPi / OmegaOrbitSi / cSecondPerDay

    ! Radiation field quantities:
    ! Photon energy; eV to erg
    EphotCgs = EphotEv / 624146012218.78162

    ! Fraction photon energy to be deposited as heat
    EpsilonNu = 1.0 - (13.6/EphotEv)

    ! Photo-ionisation cross-section Hydrogen, Spitzer (1998) eq. 5-6/table 5.1
    SigmaNu0Cgs = 6e-18 * (13.6/EphotEv)**3

    ! Planetary wind velocity; km/s to normalised
    VinfPw = VinfPwIo * Io2No_V(UnitU_)
    Vrad0Pw = Vrad0PwIo * Io2No_V(UnitU_)

    ! NEW Flo
    SwVinf = 12.0 * VinfPw
    SwMdot = 2e-13 * mSun/cSecondPerYear * Si2No_V(UnitMass_) / Si2No_V(UnitT_)

    ! Generate the grid for the optical depth computation
    call make_uniform_radiation_grid

    call test_stop(NameSub, DoTest)

  end subroutine user_init_session

!==============================================================================
! Boundary conditions for a single inner boundary face of a body. Active when
! #INNERBOUNDARY is set to USER. BATSRUS is block cartesian and values inside
! the boundary face must be passed back in cartesian coordinates.
! Values that must be set are: Rho_, Ux_, Uy_, Uz_, p_, Bx_, By_, Bz_
!==============================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModFaceBoundary, ONLY: VarsTrueFace_V, iSide
    use ModVarIndexes,   ONLY: nVar, nFluid
    use ModMultiFluid,   ONLY: iRho_I, iP_I, iUx_I, iUy_I, iUz_I, iFluid, &
         select_fluid, iUx, iUy, iUz
    use ModPhysics,      ONLY: BodyP_I, BodyRho_I

    implicit none

    real, intent(out) :: VarsGhostFace_V(nVar)

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_face_boundary'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! By default apply floating condition to all state variables of each fluid
    VarsGhostFace_V = VarsTrueFace_V

    !****************
    ! Planetary wind
    !****************
    ! Fix mass density and pressure of planet fluid to the body value
    VarsGhostFace_V(iRho_I(2)) = BodyRho_I(2)
    VarsGhostFace_V(iP_I(2))   = 2.0*BodyP_I(2)
    VarsGhostFace_V(iRho_I(4)) = BodyRho_I(4)
    VarsGhostFace_V(iP_I(4))   = BodyP_I(4)

    ! Change sign for velocities to force 0 velocity at the boundary
    VarsGhostFace_V(iUx_I(2)) = -VarsTrueFace_V(iUx_I(2))
    VarsGhostFace_V(iUy_I(2)) = -VarsTrueFace_V(iUy_I(2))
    VarsGhostFace_V(iUz_I(2)) = -VarsTrueFace_V(iUz_I(2))
    VarsGhostFace_V(iUx_I(4)) = -VarsTrueFace_V(iUx_I(4))
    VarsGhostFace_V(iUy_I(4)) = -VarsTrueFace_V(iUy_I(4))
    VarsGhostFace_V(iUz_I(4)) = -VarsTrueFace_V(iUz_I(4))

    !**************
    ! Stellar wind
    !**************
    ! Small density and pressure to avoid planet to be a source of stellar wind
    VarsGhostFace_V(iRho_I(1)) = BodyRho_I(1)
    VarsGhostFace_V(iP_I(1))   = BodyP_I(1)
    VarsGhostFace_V(iRho_I(3)) = BodyRho_I(3)
    VarsGhostFace_V(iP_I(3))   = BodyP_I(3)

    ! The planet absorbs stellar wind particles by default floating BC applied
    ! above, but no stellar wind particles back into the grid from the planet
    do iFluid = 1,nFluid

      if (iFluid == 2 .or. iFluid == 4) CYCLE

      call select_fluid

      if (VarsGhostFace_V(iUx) < 0 .and. iSide == 1) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif

      if (VarsGhostFace_V(iUx) > 0 .and. iSide == 2) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif

      if (VarsGhostFace_V(iUy) < 0 .and. iSide == 3) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif

      if (VarsGhostFace_V(iUy) > 0 .and. iSide == 4) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif

      if (VarsGhostFace_V(iUz) < 0 .and. iSide == 5) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif

      if (VarsGhostFace_V(iUz) > 0 .and. iSide == 6) then
        VarsGhostFace_V(iUx:iUz) = 0.0
      endif
    enddo

    !VarsGhostFace_V(iUx_I(1)) = 0.0
    !VarsGhostFace_V(iUy_I(1)) = 0.0
    !VarsGhostFace_V(iUz_I(1)) = 0.0
    !VarsGhostFace_V(iUx_I(3)) = 0.0
    !VarsGhostFace_V(iUy_I(3)) = 0.0
    !VarsGhostFace_V(iUz_I(3)) = 0.0

    call test_stop(NameSub, DoTest)

  end subroutine user_set_face_boundary

!==============================================================================
! Boundary conditions for the edges of the simulation domain. Active when in
! PARAM.in the #OUTERBOUNDARY is  'userfloat' or 'userInjectSW' for edges.
!==============================================================================
  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, found)

    use ModVarIndexes,   ONLY: nFluid
    use ModMultiFluid,   ONLY: iFluid, select_fluid, iRho, iRhoUx, iRhoUy, &
         iRhoUz, iUx, iUy, iUz, iP
    use ModAdvance,      ONLY: State_VGB
    use ModCellBoundary, ONLY: iMin, iMax, jMin, jMax, kMin, kMax
    use BATL_lib,        ONLY: nI, nJ, nK

    implicit none

    integer, intent(in)          :: iBlock, iSide
    logical, intent(out)         :: found
    character(len=*), intent(in) :: TypeBc

    ! Local variables
    integer :: i, j, k

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

          do iFluid = 1,nFluid
            ! NOTE: iFluid from ModMultifluid will be overwritten in order
            !       to correctly access select_fluid
            call select_fluid
            if (State_VGB(iRhoUx,1,j,k,iBlock) > 0.0) then
              State_VGB(iRhoUx,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo

        ! Inflow of stellar wind fluid
        call fix_stellar_wind

      case(2) ! The maximum X boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)

          do iFluid = 1,nFluid
            call select_fluid
            if (State_VGB(iRhoUx,nI,j,k,iBlock) < 0.0) then
              State_VGB(iRhoUx,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo

      case(3) ! The minimum Y boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,1,k,iBlock)

          do iFluid = 1,nFluid
            call select_fluid
            if (State_VGB(iRhoUy,i,1,k,iBlock) > 0.0) then
              State_VGB(iRhoUy,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo

      case(4) ! The maximum Y boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,nJ,k,iBlock)

          do iFluid = 1,nFluid
            call select_fluid
            if (State_VGB(iRhoUy,i,nJ,k,iBlock) < 0.0) then
              State_VGB(iRhoUy,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo

      case(5) ! The minimum Z boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,1,iBlock)

          do iFluid = 1,nFluid
            call select_fluid
            if (State_VGB(iRhoUz,i,j,1,iBlock) > 0.0) then
              State_VGB(iRhoUz,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo

      case(6) ! The maximum Z boundary
        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax
          State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,nK,iBlock)

          do iFluid = 1,nFluid
            call select_fluid
            if (State_VGB(iRhoUz,i,j,nK,iBlock) < 0.0) then
              State_VGB(iRhoUz,i,j,k,iBlock) = 0.0
            endif
          enddo
        enddo; enddo; enddo
      end select
    end select

    call test_stop(NameSub, DoTest, iBlock)

  contains

    subroutine fix_stellar_wind

      use ModPhysics,  ONLY: No2Si_V, Si2No_V, UnitX_, UnitTemperature_, rSun
      use ModGeometry, ONLY: Xyz_DGB
      use ModConst,    ONLY: cPi

      ! Local variables
      real :: x, y, z, SwVrad, rFrameStar, rFrameStarRstar

      FLUID: do iFluid = 1,nFluid

        call select_fluid

        if (iFluid == 2 .or. iFluid == 4) CYCLE

        do k = kMin,kMax; do j = jMin,jMax; do i = iMin,iMax

          x = Xyz_DGB(x_,i,j,k,iBlock)
          y = Xyz_DGB(y_,i,j,k,iBlock)
          z = Xyz_DGB(z_,i,j,k,iBlock)
          rFrameStar = sqrt( (Rorbit + x)**2 + y**2 + z**2 )

          ! Radial coordinate in stellar frame in units of Rstar
          rFrameStarRstar = rFrameStar * No2Si_V(UnitX_) / (Rstar*rSun)

          ! beta=1 velocity law
          SwVrad = SwVinf * (1.0 - 1.0/rFrameStarRstar)

          ! Mass density for each fluid
          State_VGB(iRho,i,j,k,iBlock) = &
               SwMdot/(4.0*cPi * rFrameStar * rFrameStar * SwVrad)

          if (iFluid == 3) &
               State_VGB(iRho,i,j,k,iBlock) = 1e-8*State_VGB(iRho,i,j,k,iBlock)

          ! Stellar wind velocities are purely radial
          State_VGB(iUx:iUz,i,j,k,iBlock) = &
               SwVrad * [ (Rorbit+x), y, z ]/rFrameStar

          ! Add phi component to account for the stellar wind being radial in
          ! the inertial (stellar) frame, but it has an azimuthal component
          ! in the rotating (planet) frame because of the orbiting planet
          State_VGB(iUx,i,j,k,iBlock) = State_VGB(iUx,i,j,k,iBlock) &
               + y * OmegaOrbit
          State_VGB(iUy,i,j,k,iBlock) = State_VGB(iUy,i,j,k,iBlock) &
               - (Rorbit+x) * OmegaOrbit

          ! Thermal pressure fixed at isothermal value
          State_VGB(iP,i,j,k,iBlock) = &
               State_VGB(iRho,i,j,k,iBlock) * 1e6*Si2No_V(UnitTemperature_)

          if (iFluid == 1) &
               State_VGB(iP,i,j,k,iBlock) = 2.0 * State_VGB(iP,i,j,k,iBlock)

          ! Convert velocities to momenta
          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
               State_VGB(iUx:iUz,i,j,k,iBlock) * State_VGB(iRho,i,j,k,iBlock)
        enddo; enddo; enddo
      enddo FLUID

    end subroutine fix_stellar_wind

  end subroutine user_set_cell_boundary

!==============================================================================
! Calculates the initial conditions of the grid for the planetary wind outflow.
! After a restart, this subroutine will not be accessed.
!==============================================================================
  subroutine user_set_ics(iBlock)

    use ModVarIndexes, ONLY: nFluid
    use ModMultiFluid, ONLY: iFluid, select_fluid, iRho, iRhoUx, iRhoUy, &
         iRhoUz, iP
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitX_, UnitTemperature_, &
         rBody, BodyRho_I, BodyP_I, Gamma_I, rSun
    use ModConst,      ONLY: cPi
    use ModGeometry,   ONLY: Xyz_DGB, true_cell, r_BLK
    use ModAdvance,    ONLY: State_VGB
    use BATL_lib,      ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    implicit none

    integer, intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k
    real    :: x0, y0, z0, rp, Vrad, SwVrad, rFrameStar, rFrameStarRstar

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    FLUID: do iFluid = 1,nFluid

      call select_fluid

      do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI

        x0 = Xyz_DGB(x_,i,j,k,iBlock)
        y0 = Xyz_DGB(y_,i,j,k,iBlock)
        z0 = Xyz_DGB(z_,i,j,k,iBlock)
        rp = r_BLK(i,j,k,iBlock)

        if (.not.true_cell(i,j,k,iBlock) .or. rp < rBody) CYCLE

        !****************
        ! Planetary wind
        !****************
        if (iFluid == 2 .or. iFluid == 4) then
          Vrad = VinfPw &
               * (1.0 - (1.0-(Vrad0Pw/VinfPw)**BetaIndex)*rBody/rp)**BetaIndex

          ! Mass density for each fluid from mass conservation
          State_VGB(iRho,i,j,k,iBlock) = BodyRho_I(iFluid) * Vrad0Pw/Vrad &
               * (rBody/rp)**2

          ! Thermal gas pressure for each fluid is a polytrope
          ! NOTE: ion fluid is corrected for electrons
          State_VGB(iP,i,j,k,iBlock) = BodyP_I(iFluid) &
               * ( State_VGB(iRho,i,j,k,iBlock) &
               /   BodyRho_I(iFluid) )**Gamma_I(iFluid)

          if (iFluid == 2) &
               State_VGB(iP,i,j,k,iBlock) = 2.0 * State_VGB(iP,i,j,k,iBlock)

          ! Momenta of each fluid (radial component, spherically symmetric)
          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
               State_VGB(iRho,i,j,k,iBlock) * Vrad * [ x0, y0, z0 ]/rp
        endif

        !**************
        ! Stellar wind
        !**************
        if (iFluid == 1 .or. iFluid == 3) then
          rFrameStar = sqrt( (Rorbit + x0)**2 + y0**2 + z0**2 )

          ! Radial coordinate in stellar frame in units of Rstar
          rFrameStarRstar = rFrameStar * No2Si_V(UnitX_) / (Rstar*rSun)

          ! beta=1 velocity law
          SwVrad = SwVinf * (1.0 - 1.0/rFrameStarRstar)

          ! Mass density for each fluid
          State_VGB(iRho,i,j,k,iBlock) = &
               SwMdot/(4.0*cPi * rFrameStar * rFrameStar * SwVrad)

          if (iFluid == 3) &
               State_VGB(iRho,i,j,k,iBlock) = 1e-8*State_VGB(iRho,i,j,k,iBlock)

          ! Stellar wind momenta are from purely radial velocity
          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = State_VGB(iRho,i,j,k,iBlock)&
               * SwVrad * [ (Rorbit+x0), y0, z0 ]/rFrameStar

          ! Thermal gas pressure for each fluid is initialised at isothermal
          ! NOTE: ion fluid is corrected for electrons
          State_VGB(iP,i,j,k,iBlock) = &
               State_VGB(iRho,i,j,k,iBlock) * 1e6*Si2No_V(UnitTemperature_)

          if (iFluid == 1) &
               State_VGB(iP,i,j,k,iBlock) = 2.0 * State_VGB(iP,i,j,k,iBlock)
        endif
      enddo; enddo; enddo
    enddo FLUID

    ! Initialising the uniform grid of densities for optical depth calculation
    if (iBlock == 1) call init_uniform_3d_grid_rho

    call test_stop(NameSub, DoTest, iBlock)

  contains
    !==========================================================================
    ! Initialisation of radiation grid with density and optical depth that will
    ! be used in the first instance of user_calc_sources.
    !==========================================================================
    subroutine init_uniform_3d_grid_rho

      use ModPhysics, ONLY: rBody, BodyRho_I

      ! Local variable
      real :: rRadGrid

      do k = 1,nRadGridZ; do j = 1,nRadGridY; do i = 1,nRadGridX
        rRadGrid = sqrt(RadGridX_C(i)**2 + RadGridY_C(j)**2 + RadGridZ_C(k)**2)

        if (rRadGrid <= rBody) then
          RhoRadGrid_C(i,j,k) = BodyRho_I(4)
          CYCLE
        end if

        Vrad = VinfPw * ( 1.0 - (1.0 - (Vrad0Pw/VinfPw)**BetaIndex) &
             *            rBody/rRadGrid )**BetaIndex

        RhoRadGrid_C(i,j,k) = BodyRho_I(4) * Vrad0Pw/Vrad * (rBody/rRadGrid)**2
      enddo; enddo; enddo

      ! This will create the initial array of tau
      call create_uniform_3d_grid_tau

      ! Reset for next timestep
      RhoRadGrid_C(1:nRadGridX,1:nRadGridY,1:nRadGridZ) = 0.0

    end subroutine init_uniform_3d_grid_rho

  end subroutine user_set_ics

!==============================================================================
! Set variables for point-implicit treatment and method specific settings.
!==============================================================================
  subroutine user_init_point_implicit

    use ModVarIndexes,    ONLY: HpRhoUx_, HpRhoUy_, HpRhoUz_, H1RhoUx_,      &
         H1RhoUy_, H1RhoUz_, HpP_, H1P_, SwHpRhoUx_, SwHpRhoUy_, SwHpRhoUz_, &
         SwH1RhoUx_, SwH1RhoUy_, SwH1RhoUz_, nFluid
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, &
         DoBalancePointImplicit, IsDynamicPointImplicit

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_init_point_implicit'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Allocate and set iVarPointImpl_I
    ! All fluid momenta and thermal pressure equations need to be implicit
    allocate(iVarPointImpl_I(4*nFluid-2))

    iVarPointImpl_I = [ &
         SwHpRhoUx_, SwHpRhoUy_, SwHpRhoUz_, &
         SwH1RhoUx_, SwH1RhoUy_, SwH1RhoUz_, &
         HpRhoUx_, HpRhoUy_, HpRhoUz_,       &
         H1RhoUx_, H1RhoUy_, H1RhoUz_,       &
         HpP_, H1P_ ]

    ! Whether Jacobian matric dS/dU will be set analytically. If true the
    ! DsDu_VVC matrix has to be programmed below
    IsPointImplMatrixSet = .false.

    ! Whether point implicit blocks should be load balanced. If set to true,
    ! UsePointImplicit_B should be set in user_calc_sources
    DoBalancePointImplicit = .false.

    ! Tell the load balancing scheme if the assignment of UsePointImplicit_B
    ! keeps changing so load balancing is needed every time step.
    IsDynamicPointImplicit = .false.

    call test_stop(NameSub, DoTest)

  end subroutine user_init_point_implicit

!==============================================================================
! Calculate extra source terms with (part-)explicit or point-implicit scheme.
!==============================================================================
  subroutine user_calc_sources(iBlock)

    use ModVarIndexes,    ONLY: HpRho_, H1Rho_, HpP_, H1P_, HpEnergy_,   &
         H1Energy_, MassFluid_I, HpRhoUx_, HpRhoUz_, H1RhoUx_, H1RhoUz_, &
         nVar, nFluid
    use ModMultiFluid,    ONLY: iFluid, select_fluid, iRho, iRhoUx, iRhoUy, &
         iRhoUz, iEnergy
    use ModPhysics,       ONLY: No2Si_V, No2Io_V, Io2No_V, Si2No_V, UnitRho_, &
         UnitU_, UnitX_, UnitN_, UnitTemperature_, UnitT_, UnitEnergyDens_,   &
         mSun, rBody, GammaMinus1_I
    use ModGeometry,      ONLY: Xyz_DGB, r_BLK, true_cell
    use ModAdvance,       ONLY: State_VGB, Source_VC
    use ModPointImplicit, ONLY: UsePointImplicit, IsPointImplSource, &
         UsePointImplicit_B
    use ModInterpolate,   ONLY: interpolate_scalar
    use ModConst,         ONLY: cGravitation, cBoltzmann, cPi
    use BATL_lib,         ONLY: nI, nJ, nK

    implicit none

    integer, intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k
    real :: x, y, z, rFrameStar, uHp_D(3), uH1_D(3)
    real :: SrhoSpecies, Smom_D(3), SEHp, SEH1, State_V(nVar)
    real :: Tion, Tneutral, Telectron, NumDensHp, NumDensH1, NumDensElectron
    real :: PhotoIonisationRate, CollisionRate, RecombinationRate
    real :: HeatPhotoIon, CoolLyalpha, CoolCollision
    real :: AlphaCoef, CollisionCoef, ExpMinTau
    real :: GravityStar_D(3), ForceCentrifugal_D(3)
    real :: ElasticCollisionCrossSection, SfacSq, CxCrossSection
    real :: CxRate, uThIon, uThNeutral, uCx, uCxStar

    logical :: DoTest, DoExtrapolate=.true.
    character(len=*), parameter :: NameSub = 'user_calc_source'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Evaluate explicit/part-implicit or point-implicit user source terms
    if (.not.UsePointImplicit) then
      ! Explicit scheme or part-implicit scheme
      call user_explicit_sources

      ! Part-implicit also uses UsePointimplicit_B. By default it will be the
      ! entire grid, but can be specified here to only do certain region
      if (UsePointImplicit_B(iBlock)) call user_implicit_sources
    elseif (IsPointImplSource) then
      ! Point-implicit sources are added; true when called in ModPointImplicit
      call user_implicit_sources
    else
      ! Explicit scheme
      call user_explicit_sources
    endif

    call test_stop(NameSub, DoTest, iBlock)

  contains
    !==========================================================================
    ! Sources of ionisation-recombination processes, heating/cooling, stellar
    ! gravity, centrifugal, and Coriolis forces.
    !==========================================================================
    subroutine user_explicit_sources

      do k = 1,nK; do j = 1,nJ; do i = 1,nI
        if (.not.true_cell(i,j,k,iBlock)) CYCLE

        ! Dump state vector into local state variables
        State_V = State_VGB(:,i,j,k,iBlock)

        ! Cell centred coordinates
        x = Xyz_DGB(x_,i,j,k,iBlock)
        y = Xyz_DGB(y_,i,j,k,iBlock)
        z = Xyz_DGB(z_,i,j,k,iBlock)

        NumDensHp       = State_V(HpRho_)/MassFluid_I(2)
        NumDensH1       = State_V(H1Rho_)/MassFluid_I(4)
        NumDensElectron = NumDensHp

        uHp_D = State_V(HpRhoUx_:HpRhoUz_) / State_V(HpRho_)
        uH1_D = State_V(H1RhoUx_:H1RhoUz_) / State_V(H1Rho_)

        ! Ion and neutral (code units), and electron temperature (Kelvin)
        ! NOTE: Tion = Te such that Tplasma = Tion + Te = 2*Tion
        Tion      = State_V(HpP_) / NumDensHp / 2.0
        Tneutral  = State_V(H1P_) / NumDensH1
        Telectron = Tion * No2Si_V(UnitTemperature_)

        ! To find attenuation in a particular cell, interpolate exp(-tau) from
        ! uniform grid to the BATSRUS AMR grid
        if (r_BLK(i,j,k,iBlock) <= rBody) then
          ExpMinTau = 0.0
        else
          ExpMinTau = &
               interpolate_scalar(ExpTauRadGrid_C(:,:,:), 3, [ 1, 1, 1 ],     &
                                  [ nRadGridX, nRadGridY, nRadGridZ ],        &
                                  [ x, y, z ],                                &
                                  RadGridX_C(:), RadGridY_C(:), RadGridZ_C(:),&
                                  DoExtrapolate)
        endif

        !******************
        ! Atomic processes
        !******************

        ! 1. Photo-ionisation rate (1/s) by radiation:
        !    H1 + hv -> H+ + e- + Eheat
        PhotoIonisationRate = SigmaNu0Cgs * FxuvCgs * ExpMinTau/EphotCgs
        PhotoIonisationRate = PhotoIonisationRate / Io2No_V(UnitT_)

        ! 1.1 Photo-ionising heat (erg/cm3/s) carried away by the electron
        HeatPhotoIon = EpsilonNu * SigmaNu0Cgs * NumDensH1*No2Io_V(UnitN_) &
             * FxuvCgs * ExpMinTau
        HeatPhotoIon = HeatPhotoIon * Io2No_V(UnitEnergyDens_)/Io2No_V(UnitT_)

        ! 2. Collisional ionisation by free electrons: H1 + e- -> H+ + 2e-
        ! Collisional coefficient (cm^3/s), Black (1981) eq. 3
        CollisionCoef = 5.83d-11 * sqrt(Telectron)* exp(-157822d0/Telectron)

        ! Collison rate (1/s) at which H+ is created
        CollisionRate = CollisionCoef * NumDensElectron * No2Io_V(UnitN_)
        CollisionRate = CollisionRate / Io2No_V(UnitT_)

        ! 2.1 Cooling from ejected bound electron carrying away heat
        ! NOTE: 2.18e-11 is 13.6eV in erg = ionization potential Hydrogen
        CoolCollision = CollisionCoef * 2.18d-11 &
             * NumDensElectron * NumDensH1 * No2Io_V(UnitN_)**2
        CoolCollision = CoolCollision &
             * Io2No_V(UnitEnergyDens_)/Io2No_V(UnitT_)

        ! 3. Radiative (ion-electron) recombination: H+ + e- -> H1 + hv
        ! Case-B coefficient (cm^3/s), Storey & Hummer (1995) table 1
        AlphaCoef = 2.70d-13 * (1.0d4/Telectron)**0.9

        ! Recombination rate (1/s) at which H1 is created
        RecombinationRate = AlphaCoef * NumDensElectron * No2Io_V(UnitN_)
        RecombinationRate = RecombinationRate / Io2No_V(UnitT_)

        ! 4. Lyman alpha cooling (erg/cm3/s), Black (1981) table 3
        CoolLyalpha = 7.5d-19 * exp(-118348d0/Telectron) &
             * NumDensElectron * NumDensH1 * No2Io_V(UnitN_)**2
        CoolLyalpha = CoolLyalpha * Io2No_V(UnitEnergyDens_)/Io2No_V(UnitT_)

        ! 5. Update hydro equations (see also Draine 1986, sec. 4)

        ! Rate per volume at which H+ mass is produced from H1 mass
        SrhoSpecies = (PhotoIonisationRate + CollisionRate) * State_V(H1Rho_) &
             - RecombinationRate * State_V(HpRho_)

        Source_VC(HpRho_,i,j,k) = Source_VC(HpRho_,i,j,k) + SrhoSpecies
        Source_VC(H1Rho_,i,j,k) = Source_VC(H1Rho_,i,j,k) - SrhoSpecies

        ! Rate per volume of momentum transfer from H+ to H1
        Smom_D = (PhotoIonisationRate + CollisionRate)*State_V(H1Rho_)*uH1_D &
             - RecombinationRate * State_V(HpRho_)*uHp_D

        Source_VC(HpRhoUx_:HpRhoUz_,i,j,k) = &
             Source_VC(HpRhoUx_:HpRhoUz_,i,j,k) + Smom_D
        Source_VC(H1RhoUx_:H1RhoUz_,i,j,k) = &
             Source_VC(H1RhoUx_:H1RhoUz_,i,j,k) - Smom_D

        ! Rate per volume at which TOTAL energy is added to the H+ fluid
        SEHp = (PhotoIonisationRate + CollisionRate)               &
             * (0.5*State_V(H1Rho_) * sum((uH1_D - uHp_D)**2)      &
             +  3.0/2.0 * State_V(H1P_))                           &
             - 3.0/2.0 * RecombinationRate * State_V(HpP_)/2.0     &
             + sum(Smom_D * uHp_D) - 0.5*SrhoSpecies*sum(uHp_D**2) &
             + HeatPhotoIon - CoolLyalpha - CoolCollision

        ! Rate per volume at which TOTAL energy is added to the H1 fluid
        SEH1 = RecombinationRate                                              &
             * (0.5*State_V(HpRho_) * sum((uH1_D - uHp_D)**2)                 &
             +  3.0/2.0 * State_V(HpP_)/2.0)                                  &
             - 3.0/2.0 * (PhotoIonisationRate + CollisionRate) * State_V(H1P_)&
             - sum(Smom_D * uH1_D) + 0.5*SrhoSpecies*sum(uH1_D**2)            &
             + HeatPhotoIon - CoolLyalpha - CoolCollision

        Source_VC(HpEnergy_,i,j,k) = Source_VC(HpEnergy_,i,j,k) + SEHp
        Source_VC(H1Energy_,i,j,k) = Source_VC(H1Energy_,i,j,k) + SEH1

        ! When doing non-conservative option
        Source_VC(HpP_,i,j,k) = Source_VC(HpP_,i,j,k) + GammaMinus1_I(2) &
             * (SEHp - sum(Smom_D * uHp_D) + 0.5*SrhoSpecies*sum(uHp_D**2))
        Source_VC(H1P_,i,j,k) = Source_VC(H1P_,i,j,k) + GammaMinus1_I(4) &
             * (SEH1 + sum(Smom_D * uH1_D) - 0.5*SrhoSpecies*sum(uH1_D**2))

        ! Store for printing variables
        PlotRateRec_CB(i,j,k,iBlock)      = RecombinationRate
        PlotRateCollIon_CB(i,j,k,iBlock)  = CollisionRate
        PlotRatePhotoIon_CB(i,j,k,iBlock) = PhotoIonisationRate
        PlotHeat_CB(i,j,k,iBlock)         = HeatPhotoIon
        PlotCoolLya_CB(i,j,k,iBlock)      = CoolLyalpha
        PlotCoolColl_CB(i,j,k,iBlock)     = CoolCollision
        PlotExp_mTau_CB(i,j,k,iBlock)     = ExpMinTau
      enddo; enddo; enddo

      ! Stellar gravity for stellar wind
      ! do iFluid = 1,nFluid
      !   call select_fluid
      !   if (iFluid==2 .or. iFluid==4) CYCLE
      !   do k = 1,nK; do j = 1,nJ; do i = 1,nI
      !     if (.not.true_cell(i,j,k,iBlock)) CYCLE

      !     x = Xyz_DGB(x_,i,j,k,iBlock)
      !     y = Xyz_DGB(y_,i,j,k,iBlock)
      !     z = Xyz_DGB(z_,i,j,k,iBlock)
      !     rFrameStar = sqrt( (Rorbit + x)**2 + y**2 + z**2 )

      !     ! Stellar gravity: -GMstar/R^2, acts in -r direction
      !     GravityStar_D = -gStar * [ (Rorbit + x), y, z ] / rFrameStar**3

      !     ! Update momenta and energy
      !     Source_VC(iRhoUx:iRhoUz,i,j,k) = Source_VC(iRhoUx:iRhoUz,i,j,k) &
      !          + State_VGB(iRho,i,j,k,iBlock) * GravityStar_D

      !     Source_VC(iEnergy,i,j,k) = Source_VC(iEnergy,i,j,k) &
      !          + sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) * GravityStar_D)
      !   enddo; enddo; enddo
      ! enddo

      if (.not.UseTidal .and. .not.UseCoriolis) RETURN

      ! Extra sources for all fluids if activated
      FLUID: do iFluid = 1,nFluid

        call select_fluid

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

            ! Update momenta and energy
            Source_VC(iRhoUx:iRhoUz,i,j,k) = Source_VC(iRhoUx:iRhoUz,i,j,k) &
                 + State_VGB(iRho,i,j,k,iBlock) * GravityStar_D

            Source_VC(iEnergy,i,j,k) = Source_VC(iEnergy,i,j,k) &
                 + sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) * GravityStar_D)

            ! Centrifugal force: -rho * Omega x (Omega x r), in +r direction
            ! NOTE: has only x,y component if rotation axis is along z
            ForceCentrifugal_D = OmegaOrbit**2 * [ (Rorbit + x), y, 0.0 ]

            ! Update momenta and energy
            Source_VC(iRhoUx:iRhoUz,i,j,k) = Source_VC(iRhoUx:iRhoUz,i,j,k) &
                 + State_VGB(iRho,i,j,k,iBlock) * ForceCentrifugal_D

            Source_VC(iEnergy,i,j,k) = Source_VC(iEnergy,i,j,k) &
                 + sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)    &
                 *     ForceCentrifugal_D)
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

            Source_VC(iRhoUx,i,j,k) = Source_VC(iRhoUx,i,j,k) &
                 + 2.0*OmegaOrbit * State_VGB(iRhoUy,i,j,k,iBlock)

            Source_VC(iRhoUy,i,j,k) = Source_VC(iRhoUy,i,j,k) &
                 - 2.0*OmegaOrbit * State_VGB(iRhoUx,i,j,k,iBlock)
          enddo; enddo; enddo
        endif
      enddo FLUID

    end subroutine user_explicit_sources
    !==========================================================================
    ! Mass-loading processes and elastic collisions between fluids.
    !==========================================================================
    subroutine user_implicit_sources

	! TO DO: implement SW-PW coupling
	return

    end subroutine user_implicit_sources

  end subroutine user_calc_sources

!==============================================================================
! Perform hydrodynamic update and set hydrogen mass density on radiation grid.
!==============================================================================
  subroutine user_update_states(iBlock)

    use ModUpdateState, ONLY: update_state_normal
    use ModMultiFluid,  ONLY: iRho_I, iP_I
    use ModPhysics,     ONLY: rBody
    use ModGeometry,    ONLY: r_BLK
    use ModAdvance,     ONLY: StateOld_VGB, State_VGB
    use ModFaceFlux,    ONLY: iLeft, jLeft, kLeft, iRight, jRight, kRight, &
         iBlockFace, FluxLeft_V, FluxRight_V

    implicit none

    integer, intent(in) :: iBlock

    ! Local variables
    integer :: i, j, k

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Kill the hydro fluxes inside the planet
    ! if(r_BLK(iLeft,jLeft,kLeft,iBlockFace) < rBody .and. &
    !      r_BLK(iRight,jRight,kRight,iBlockFace) < rBody) then
    !   FluxLeft_V(iRho_I(1):iP_I(1))  = 0.0
    !   FluxRight_V(iRho_I(3):iP_I(3)) = 0.0
    ! endif

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

    use ModPhysics, ONLY: No2Io_V, UnitRho_, UnitT_, UnitEnergyDens_, UnitU_
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
    	NameTecVar = 'exp(-`t)'

    case('heat') ! Heating rate
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotHeat_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
    	NameTecUnit = '(erg/cm^3/s)'
    	NameIdlUnit = 'erg/cm^3/s'
    	NameTecVar  = 'H'

    case('cool') ! Total cooling rate
    	PlotVar_G(1:nI,1:nJ,1:nK) = &
          PlotCoolLya_CB(:,:,:,iBlock) + PlotCoolColl_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
    	NameIdlUnit = 'erg/cm^3/s'
    	NameTecVar  = 'C'

    case('clya') ! Cooling rate Ly alpha
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotCoolLya_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
    	NameIdlUnit = 'erg/cm^3/s'
    	NameTecVar  = 'C_L_y_`a'

    case('ccoll') ! Cooling rate collisions
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotCoolColl_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G * No2Io_V(UnitEnergyDens_) / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(erg/cm^3/s)'
    	NameIdlUnit = 'erg/cm^3/s'
    	NameTecVar  = 'C_c_o_l_l'

    case('rion')
    	PlotVar_G(1:nI,1:nJ,1:nK) = &
          PlotRatePhotoIon_CB(:,:,:,iBlock) + PlotRateCollIon_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
    	NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
    	NameTecVar  = 'R_i_o_n'

    case('rphoto')
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotRatePhotoIon_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
    	NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
    	NameTecVar  = 'R_p_h_,_i_o_n'

    case('rcoll')
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotRateCollIon_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
    	NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
    	NameTecVar  = 'R_c_o_l_l_,_i_o_n'

    case('rrec')
    	PlotVar_G(1:nI,1:nJ,1:nK) = PlotRateRec_CB(:,:,:,iBlock)
    	if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
    	NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
    	NameTecVar  = 'R_r_e_c'

    case('rcx')
      PlotVar_G(1:nI,1:nJ,1:nK) = PlotChExcHpH1_CB(:,:,:,iBlock)
      if (IsDimensional) then
        PlotVar_G = PlotVar_G / No2Io_V(UnitT_)
      endif
      NameTecUnit = '(1/s)'
      NameIdlUnit = '1/s'
      NameTecVar  = 'R_Cx_Hp_H1'

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
        allocate(PlotChExcHpH1_CB(1:nI,1:nJ,1:nK,MaxBlock))
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
      if (allocated(PlotChExcHpH1_CB))    deallocate(PlotChExcHpH1_CB)

    case('initial condition done')
      ! Initialising grid of densities after startup OR restart
      do iBlock = 1,nBlock
        call create_uniform_3d_grid_rho(iBlock)
      enddo

      if (nProc > 1) then
        call mpi_reduce_real_array(RhoRadGrid_C, size(RhoRadGrid_C), &
                                   MPI_SUM, 0, iComm, iError)

        if (iProc == 0) call create_uniform_3d_grid_tau

        call MPI_BCAST(ExpTauRadGrid_C, size(ExpTauRadGrid_C), &
                       MPI_REAL, 0, iComm, iError)

        ! Reseting for next timestep
        RhoRadGrid_C = 0.0
      endif

    case('write progress')
      ! Output some extra model setup information
      call print_model_info
    end select

    call test_stop(NameSub, DoTest)

  end subroutine user_action

!==============================================================================
! Setup a uniform base grid for performing the ray tracing to compute the
! optical depth integral. To be used within the source term calculation.
!==============================================================================
  subroutine make_uniform_radiation_grid

    use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2

    implicit none

    ! Local variables
    integer :: idum, nPointLeft, nPointRight
    real    :: ResolutionMid, ResolutionRight, ResolutionLeft
    real    :: LengthMid, LengthLeft, LengthRight

    ! Grid is in normalised units
    ! === x-axis ===
  	LengthMid   = 2.0 - (-2.0)
  	LengthLeft  = -2.0 - x1
  	LengthRight = x2 - 2.0

    nPointRight = int(LengthRight * (nRadGridX+2-41)/(LengthLeft+LengthRight))
  	nPointLeft  = nRadGridX+2-41 - nPointRight

    ! 40 deltas within +/- 2, the same in all three axis
  	ResolutionMid   = LengthMid/(41-1)
  	ResolutionLeft  = LengthLeft/(nPointLeft-1)
  	ResolutionRight = LengthRight/(nPointRight-1)

    ! Generate x-grid
    RadGridX_C(1:nPointLeft) = &
        [( x1 + ResolutionLeft*(idum-1.0), idum=1,nPointLeft )]

    RadGridX_C(nPointLeft:nPointLeft+41) = &
        [( -2.0 + ResolutionMid*(idum-nPointLeft), &
            idum=nPointLeft,nPointLeft+41 )]

    RadGridX_C(nPointLeft+41:nRadGridX) = &
        [( 2.0 + ResolutionRight*(idum -(nPointLeft+40)), &
            idum=nPointLeft+41,nRadGridX )]

    ! === y-axis ====
    LengthLeft  = -2.0 - y1
    LengthRight = y2 - 2.0

    nPointRight = int(LengthRight * (nRadGridY+2-41)/(LengthLeft+LengthRight))
    nPointLeft  = nRadGridY+2-41 - nPointRight

    ResolutionLeft  = LengthLeft/(nPointLeft-1)
    ResolutionRight = LengthRight/(nPointRight-1)

    RadGridY_C(1:nPointLeft) = &
        [( y1 + ResolutionLeft*(idum-1.0), idum=1,nPointLeft )]

    RadGridY_C(nPointLeft:nPointLeft+41) = &
        [( -2.0 + ResolutionMid*(idum-nPointLeft), &
            idum=nPointLeft,nPointLeft+41 )]

    RadGridY_C(nPointLeft+41:nRadGridY) = &
        [( 2.0 + ResolutionRight*(idum -(nPointLeft+40)), &
            idum=nPointLeft+41,nRadGridY )]

    ! === z-axis ===
    LengthLeft  = -2.0 - z1
    LengthRight = z2 - 2.0

    nPointRight = int(LengthRight * (nRadGridZ+2-41)/(LengthLeft+LengthRight))
    nPointLeft  = nRadGridZ+2-41 - nPointRight

    ResolutionLeft  = LengthLeft/(nPointLeft-1)
    ResolutionRight = LengthRight/(nPointRight-1)

    RadGridZ_C(1:nPointLeft) = &
        [( z1 + ResolutionLeft*(idum-1.0), idum=1,nPointLeft )]

    RadGridZ_C(nPointLeft:nPointLeft+41) = &
        [( -2.0 + ResolutionMid*(idum-nPointLeft), &
            idum=nPointLeft,nPointLeft+41 )]

    RadGridZ_C(nPointLeft+41:nRadGridZ) = &
        [( 2.0 + ResolutionRight*(idum -(nPointLeft+40)), &
            idum=nPointLeft+41,nRadGridZ )]

  end subroutine make_uniform_radiation_grid

!==============================================================================
! Interpolate values of density from the AMR grid to the radiation grid. The
! computation is performed on the MHD grid spanned by each individual block.
! The densities for the entire grid are computed in ModBatsrusMethods.f90.
!==============================================================================
  subroutine create_uniform_3d_grid_rho(iBlock)

    use ModVarIndexes,  ONLY: H1Rho_
    use ModPhysics,     ONLY: rBody, BodyRho_I
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

    ! Minimum index where RadGridX_C > x_min_iblock
    iMinX = minloc(RadGridX_C, 1, mask = RadGridX_C > CoordMin_DB(x_, iBlock))

    ! Maximum index where RadGridX_C <= x_max_iblock
    iMaxX = maxloc(RadGridX_C, 1, mask = RadGridX_C <= CoordMax_DB(x_, iBlock))

    ! Same procedure for y,z direction
    iMinY = minloc(RadGridY_C, 1, mask = RadGridY_C > CoordMin_DB(y_, iBlock))
    iMaxY = maxloc(RadGridY_C, 1, mask = RadGridY_C <= CoordMax_DB(y_, iBlock))
    iMinZ = minloc(RadGridZ_C, 1, mask = RadGridZ_C > CoordMin_DB(z_, iBlock))
    iMaxZ = maxloc(RadGridZ_C, 1, mask = RadGridZ_C <= CoordMax_DB(z_, iBlock))

    do k = iMinZ,iMaxZ; do j = iMinY,iMaxY; do i = iMinX,iMaxX

      ! No interpolation inside planet
      if (sqrt(RadGridX_C(i)**2 + RadGridY_C(j)**2 + RadGridZ_C(k)**2) &
           <= rBody) then
        RhoRadGrid_C(i,j,k) = BodyRho_I(2)
        CYCLE
      endif

      ! There are two ways of interpolating: considering the ghost cells or not
      RhoRadGrid_C(i,j,k) = &
        interpolate_scalar(State_VGB(H1Rho_,:,:,:,iBlock),                   &
                           3, [ MinI, MinJ, MinK ], [ MaxI, MaxJ, MaxK ],    &
                           [ RadGridX_C(i), RadGridY_C(j), RadGridZ_C(k) ],  &
                           Xyz_DGB(x_,MinI:MaxI,1        ,1        ,iBlock), &
                           Xyz_DGB(y_,1        ,MinJ:MaxJ,1        ,iBlock), &
                           Xyz_DGB(z_,1        ,1        ,MinK:MaxK,iBlock), &
                           DoExtrapolate)
    enddo; enddo; enddo

  end subroutine create_uniform_3d_grid_rho

!==============================================================================
! Calculate the sum over n(Hydrogen) to get optical depth. This subroutine is
! called from ModBatsrusMethods.f90 after message passing to iProc=0. Thus,
! the value of ExpTauRadGrid_C is only complete for iProc=0.
!==============================================================================
  subroutine create_uniform_3d_grid_tau

    use ModVarIndexes, ONLY: MassFluid_I
    use ModPhysics,    ONLY: No2Si_V, No2Io_V, UnitX_, UnitN_

    implicit none

    ! Local variables
    integer :: idum
    real    :: CellSizeX, TauRadGrid_C(nRadGridX, nRadGridY, nRadGridZ), cTau

    ! Constant optical depth (factor 100 is to go from m -> cm)
    cTau = SigmaNu0Cgs * No2Io_V(UnitN_) * 100.*No2Si_V(UnitX_)

  	! Summing density in all cells along x to get optical depth
    ! This is zero in the first point, at the left edge of the grid
    TauRadGrid_C(1,:,:) = RhoRadGrid_C(1,:,:) * (RadGridX_C(2) - RadGridX_C(1))

  	do idum = 2,nRadGridX
      CellSizeX = RadGridX_C(idum) - RadGridX_C(idum-1)
      TauRadGrid_C(idum,:,:) = RhoRadGrid_C(idum,:,:) * CellSizeX &
           + TauRadGrid_C(idum-1,:,:)
  	enddo

    ! Normalisation of optical depth
    TauRadGrid_C    = TauRadGrid_C * cTau/MassFluid_I(2)
    ExpTauRadGrid_C = exp(-TauRadGrid_C)

  end subroutine create_uniform_3d_grid_tau

!==============================================================================
! Print some extra information of interest to the screen when called for.
!==============================================================================
  subroutine print_model_info

    use ModMain,        ONLY: IsStandAlone, UseRotatingFrame, UseRotatingBC, &
         UseRayTrace
    use ModPhysics,     ONLY: No2Si_V, No2Io_V, UnitX_, UnitRho_, UnitU_,  &
         UnitN_, UnitP_, UnitB_, UnitMass_, UnitT_, UnitRhoU_, UnitJ_,     &
         UnitDivB_, UnitAngle_, UnitTemperature_, UnitEnergyDens_,         &
         UnitElectric_, UnitPoynting_, BodyRhoSpecies_I, BodyRho_I,        &
         BodyNDim_I, BodyP_I, rPlanetSi, RotPeriodSi, OmegaBody
    use ModVarIndexes,  ONLY: IonFirst_, IonLast_, ScalarFirst_, ScalarLast_, &
         MassFluid_I, nFluid
    use ModGeometry,    ONLY: TypeGeometry
    use ModAdvance,     ONLY: UseNonConservative
    use ModPlanetConst, ONLY: mPlanet_I, rPlanet_I, Jupiter_
    use CON_planet,     ONLY: MassPlanet
    use ModConst,       ONLY: cSecondPerDay
    use ModIO,          ONLY: restart

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
    write(*,'(5X,A23,f10.3,A10)')  'Rorbit = ', RorbitAu, ' au'
    write(*,'(5X,A23,es10.3,A10)') 'OmegaOrbit = ', OmegaOrbitSi,' 1/s'
    write(*,'(5X,A23,es10.3,A10)') 'PeriodOrbit = ', PeriodOrbitDays,' days'
    write(*,'(5X,A23,es10.3,A10)') 'OmegaOrbit = ', OmegaOrbit,'normalised'
    write(*,'(5X,A23,es10.3,A10)') 'Rorbit = ', Rorbit,' normalised'
    write(*,'(2X, A)') ' '

    if (Mstar == 0.0) &
        call stop_mpi('Missing #USR_SYSTEMPROPERTIES in PARAM.in. &
        Maybe due to a restart of the code?')

    if (FxuvCgs == 0.0) &
        call stop_mpi('Missing #USR_FXUV in PARAM.in. &
        Maybe due to a restart of the code?')

    write(*,'(2X, A)') ' **** Stellar properties:'
    write(*,'(5X,A23,f10.3,A10)') 'Mstar = ', Mstar, ' Msun'
    write(*,'(5X,A23,f10.3,A10)') 'Rstar = ', Rstar, ' Rsun'
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Radiation properties:'
    write(*,'(5X,A23,es10.3,A10)') 'Fxuv orbit = ', FxuvCgs,' erg/cm2/s'
    write(*,'(5X,A23,es10.3,A10)') 'sigma_nu0 = ', SigmaNu0Cgs,' cm^2'
    write(*,'(5X,A23,es10.3,A10)') 'epsilon_nu = ', EpsilonNu,'  '
    write(*,'(5X,A23,es10.3,A10)') 'hnu_photon_eV = ', EphotEv,' eV'
    write(*,'(5X,A23,es10.3,A10)') 'hnu_photon = ', EphotCgs,' erg'
    write(*,'(2X, A)') ' '

    write(*,'(2X, A)') ' **** Other fluid properties:'
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyRho fluid 1  = ', BodyRho_I(1),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyRho fluid 2  = ', BodyRho_I(2),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyRho fluid 3  = ', BodyRho_I(3),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyRho fluid 4  = ', BodyRho_I(4),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyNDim fluid 1 = ', BodyNDim_I(1),' '
    write(*,'(5X,A23,es10.3,A10)') &
       'BodyNDim fluid 2  = ', BodyNDim_I(2),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyNDim fluid 3 = ', BodyNDim_I(3),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'BodyNDim fluid 4  = ', BodyNDim_I(4),' '
    write(*,'(5X,A23,es10.3,A10)') &
      'FluidMass 1 = ', MassFluid_I(1),'  '
    write(*,'(5X,A23,es10.3,A10)') &
      'FluidMass 2 = ', MassFluid_I(2),'  '

    write(*,'(5X,A23,i10)') 'nFluid = ', nFluid
    write(*,'(5X,A23,i10)') 'IonFirst_ = ', IonFirst_
    write(*,'(5X,A23,i10)') 'IonLast_ = ' , IonLast_
    write(*,'(5X,A23,i10)') 'ScalarFirst_ = ', ScalarFirst_
    write(*,'(5X,A23,i10)') 'ScalarLast_ = ', ScalarLast_

    write(*,'(5X,A23,i10)') 'nRadGridX = ', nRadGridX
    write(*,'(5X,A23,i10)') 'nRadGridY = ', nRadGridY
    write(*,'(5X,A23,i10)') 'nRadGridZ = ', nRadGridZ
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

  end subroutine print_model_info

end module ModUser
