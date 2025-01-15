module held_suarez_cam

  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  ! ----------------------------------------------------------------------------------
  ! Modified to inlude an option for the Polvani and Kushner (2002) ,
  ! GRL, 29, 7, 10.1029/2001GL014284 (PK02) equilibrium temperature profile.
  !
  ! If pkstrat = .True. then Polvani and Kushner relaxation temperature profile
  ! is used.
  !
  ! Namelist parameter: vgamma, sets the vortex strength (gamma parameter in PK02)
  !
  ! Modifications are denoted by
  !
  ! !PKSTRAT
  ! blah blah blah
  ! !END-PKSTRAT
  !
  ! Isla Simpson 8th June 2017
  !
  ! Updated for CESM2 release, Isla Simpson, 30th May 2018
  ! Modified for CESM2.2.3 alpha 17b release, Yi Wang, Aug 2024
  ! ----------------------------------------------------------------------------------
  ! Modified to allow user adding more namelist parameters for setting up PK02 forcings
  !
  ! Added namelist parameter:
  !      - noPV: whether to activate NO POLAR VOTEX case.
  !              set to .False. to generate polar vortex
  !                     .True.  for NO POLAR VORTEX case .i.e., w=0, vgamma and lat0 inactive)
  !      - pret: lower limit (Pa) for upper level damping (sponge layer)
  !      - lat0: (phi_0 parameter in PK02)
  !      - dellat: (delta_phi parameter in PK02)
  !      - dely: (delta_y parameter in PK02)
  !      - eps: (epsilon parameter in PK02)
  !      - delz: (delta_z parameter in PK02)
  !
  ! Modifications are denoted by:
  !
  ! ! PKSTRAT (Yi Wang, Sep 2024)
  ! blah blah blah
  ! ! END-PKSTRAT (Yi Wang, Sep 2024)
  !
  ! Modified for CESM2.2.3 alpha 17b release, Yi Wang, Sep 2024
  ! ----------------------------------------------------------------------------------



  use shr_kind_mod, only: r8 => shr_kind_r8

  ! PKSTAT
  use cam_abortutils  ,only:endrun
  use spmd_utils  ,only:masterproc
  use spmd_utils, only: mpicom,mstrid=>masterprocid, mpi_integer, mpi_real8, &
                        mpi_logical, mpi_character
  ! END-PKSTRAT

  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: held_suarez_init, held_suarez_tend

  ! PKSTRAT
  public :: pkstrat_readnl
  ! END-PKSTRAT

  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

  ! PKSTRAT
  logical :: pkstrat !.True. to use the PK02 TEQ
  real(r8) :: vgamma ! vortex gamma parameter (controlling vortex strength)
  ! END PKSTRAT

  ! PKSTRAT (Yi Wang, Sep 2024)
  logical :: noPV   ! whether to activate NO POLAR VOTEX case
  integer :: pret   ! lower limit (Pa) for upper level damping (sponge layer), p_sp parameter in PK02
  integer :: lat0   ! phi_0 parameter in PK02
  integer :: dellat ! delta_phi parameter in PK02
  integer :: dely   ! delta_y parameter in PK02
  integer :: eps    ! epsilon parameter in PK02
  integer :: delz   ! delta_z parameter in PK02
  ! END-PKSTRAT (Yi Wang, Sep 2024)

!=======================================================================
contains
!=======================================================================

  subroutine held_suarez_init()
    use physics_buffer,     only: physics_buffer_desc
    use cam_history,        only: addfld, add_default
    use ref_pres,           only: psurf_ref
    use held_suarez_1994,   only: held_suarez_1994_init

    ! Local variables
    character(len=512) :: errmsg
    integer            :: errflg

    ! Set model constant values
    call held_suarez_1994_init(psurf_ref, errmsg, errflg)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
    call add_default('QRS', 1, ' ')
 end subroutine held_suarez_init

  subroutine held_suarez_tend(state, ptend, ztodt)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !  algorithm 1: Held/Suarez IDEALIZED physics
    !  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
    !  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
    !
    ! Author: J. Olson
    !
    !-----------------------------------------------------------------------
    use air_composition,    only: cappav, cpairv
    use ref_pres,           only: pref_mid_norm
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use held_suarez_1994,   only: held_suarez_1994_run

    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    integer                            :: i, k             ! Longitude, level indices

    character(len=64)                  :: scheme_name      ! CCPP-required variables (not used in CAM)
    character(len=512)                 :: errmsg
    integer                            :: errflg

    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)

    ! PKSTRAT (Yi Wang, Sep 2024)
    ! call held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappav(1:ncol,:,lchnk), &
    !                           cpairv(1:ncol,:,lchnk), state%pmid(1:ncol,:),            &
    !                           state%u(1:ncol,:), state%v(1:ncol,:), state%t(1:ncol,:), &
    !                           ptend%u(1:ncol,:), ptend%v(1:ncol,:), ptend%s(1:ncol,:), &
    !                           scheme_name, errmsg, errflg)
    call held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappav(1:ncol,:,lchnk), &
                              cpairv(1:ncol,:,lchnk), state%pmid(1:ncol,:),            &
                              state%u(1:ncol,:), state%v(1:ncol,:), state%t(1:ncol,:), &
                              ptend%u(1:ncol,:), ptend%v(1:ncol,:), ptend%s(1:ncol,:), &
                              scheme_name, errmsg, errflg, pkstrat, vgamma, noPV, pret, lat0, dellat, dely, eps, delz)
    ! END-PKSTRAT (Yi Wang, Sep 2024)

    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)

  end subroutine held_suarez_tend

! PKSTRAT
  subroutine pkstrat_readnl(nlfile)
! ---------------------------------------------
! Read in namelist parameters to control the PK02 stratosphere
! logical: pkstrat.  If (pkstrat=.True.) then use the PK02
! stratosphere
!
! Isla Simpson (6th March 2017)
! Yi Wang (Sep 2024)

  use shr_kind_mod,only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use namelist_utils  ,only:find_group_name
  use units           ,only:getunit,freeunit
  character(len=*),intent(in)::nlfile
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'pkstrat_readnl'

  ! PKSTRAT (Yi Wang, Sep 2024)
  namelist /pkstrat_nl/ pkstrat, vgamma, noPV, pret, lat0, dellat, dely, eps, delz
  ! END-PKSTRAT (Yi Wang, Sep 2024)

  !Set default namelist parameters
  pkstrat=.False.
  vgamma=2._r8
  ! PKSTRAT (Yi Wang, Sep 2024)
  noPV=.False.
  pret=50
  lat0=-50
  dellat=10
  dely=60
  eps=-10
  delz=10
  ! END-PKSTRAT (Yi Wang, Sep 2024)

  !Read in namelist parameters
  if (masterproc) then
   unitn = getunit()
   open(unitn,file=trim(nlfile),status='old')
   call find_group_name(unitn,'pkstrat_nl',status=ierr)
   if (ierr.eq.0) then
     read(unitn,pkstrat_nl,iostat=ierr)
     if (ierr.ne.0) then
       call endrun('pkstrat_readnl:: ERROR reading namelist')
     endif
   endif
   close(unitn)
   call freeunit(unitn)
  endif

  call mpi_bcast(pkstrat   , 1, mpi_logical, mstrid, mpicom, ierr)
  call mpi_bcast(vgamma    , 1, mpi_real8, mstrid, mpicom, ierr)

  ! PKSTRAT (Yi Wang, Sep 2024)
  call mpi_bcast(noPV    , 1, mpi_logical, mstrid, mpicom, ierr)
  call mpi_bcast(pret    , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(lat0    , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(dellat    , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(dely    , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(eps    , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(delz    , 1, mpi_integer, mstrid, mpicom, ierr)
  ! PKSTRAT (Yi Wang, Sep 2024)

  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_1994_run")

  return
  end subroutine

! ENDPKSTRAT

end module held_suarez_cam
