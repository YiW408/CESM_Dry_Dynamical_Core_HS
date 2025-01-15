!< \section arg_table_held_suarez_1994
!! \htmlinclude held_suarez_1994.html
module held_suarez_1994

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
! Modified for CESM2 release, Isla Simpson, 30th May 2018
! Modified for CESM2.2.3 alpha 17b release, Yi Wang, Aug 2024
! ----------------------------------------------------------------------------------
! Modified to allow user adding more namelist parameters for setting up PK02 forcings
!
! 1. Added namelist parameter:
!      - noPV: whether to activate NO POLAR VOTEX case.
!              set to .False. to generate polar vortex
!                     .True.  for NO POLAR VORTEX case .i.e., w=0, vgamma and lat0 inactive)
!      - pret: lower limit (Pa) for upper level damping (sponge layer)
!      - lat0: (phi_0 parameter in PK02)
!      - dellat: (delta_phi parameter in PK02)
!      - dely: (delta_y parameter in PK02)
!      - eps: (epsilon parameter in PK02)
!      - delz: (delta_z parameter in PK02)
! 2. Add different weighting function W(phi) formula for positive lat0
!
! Modifications are denoted by:
!
! ! PKSTRAT (Yi Wang, Sep 2024)
! blah blah blah
! ! END-PKSTRAT (Yi Wang, Sep 2024)
!
! Modified for CESM2.2.3 alpha 17b release, Yi Wang, Sep 2024
! ----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  ! PKSTRAT
  use cam_logfile, only: iulog
  ! END-PKSRAT

  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994_run

  !!
  !! Forcing parameters
  !!
  real(kind_phys), parameter :: efoldf  =  1._kind_phys  ! efolding time for wind dissipation
  real(kind_phys), parameter :: efolda  = 40._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: efolds  =  4._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: sigmab  =  0.7_kind_phys ! threshold sigma level
  real(kind_phys), parameter :: t00     = 200._kind_phys ! minimum reference temperature
  real(kind_phys), parameter :: kf      = 1._kind_phys/(86400._kind_phys*efoldf) ! 1./efolding_time for wind dissipation

  real(kind_phys), parameter :: onemsig = 1._kind_phys - sigmab ! 1. - sigma_reference

  real(kind_phys), parameter :: ka      = 1._kind_phys/(86400._kind_phys * efolda) ! 1./efolding_time for temperature diss.
  real(kind_phys), parameter :: ks      = 1._kind_phys/(86400._kind_phys * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(kind_phys)            :: pref    = 0.0_kind_phys  ! Surface pressure

  ! PKSTRAT
  ! Forcing parameters for PK02 option
  real(kind_phys), parameter :: efoldstrat = 0.5_kind_phys
  real(kind_phys), parameter :: kfstrat     = 1._kind_phys/(86400._kind_phys*efoldstrat)
  ! END-PKSTRAT


!=======================================================================
contains
!=======================================================================

!> \section arg_table_held_suarez_1994_init Argument Table
!! \htmlinclude held_suarez_1994_init.html
  subroutine held_suarez_1994_init(pref_in, errmsg, errflg)

    !! Dummy arguments
    real(kind_phys),    intent(in)  :: pref_in

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    pref   = pref_in

  end subroutine held_suarez_1994_init

!> \section arg_table_held_suarez_1994_run Argument Table
!! \htmlinclude held_suarez_1994_run.html

  ! subroutine held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappa, &
  !       cpair, pmid, uwnd, vwnd, temp, du, dv, ds, scheme_name, errmsg, errflg)
  ! PKSTRAT (Yi Wang, Sep 2024)
  subroutine held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappa, &
       cpair, pmid, uwnd, vwnd, temp, du, dv, ds, scheme_name, errmsg, errflg, &
       pkstrat, vgamma, noPV, pret, lat0, dellat, dely, eps, delz)
  ! END-PKSTRAT (Yi Wang, Sep 2024)

    !-----------------------------------------------------------------------
    !
    ! 7th March 2017 (Isla Simpson) - modified to include the option of the
    ! Polvani and Kushner (2002), GRL, 29, 7, 10.1029/2001GL014284 stratospheric
    ! relaxation temperature profile
    !
    ! Modifications denoted by
    !
    ! !PKSTRAT
    ! blah blah blah
    ! !END-PKSTRAT
    !
    ! Updated for CESM2 release, Isla Simpson, May 30th
    ! Modified for CESM2.2.3 alpha 17b release, Yi Wang, Aug 2024
    !-----------------------------------------------------------------------

    ! PKSTRAT
    use physconst,      only: rair, gravit ! Gas constant and gravity for pkstrat
    ! END-PKSTRAT

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pver                    ! Num vertical levels
    integer,  intent(in)  :: ncol                    ! Num active columns
    real(kind_phys), intent(in)  :: pref_mid_norm(:) ! reference pressure normalized by surface pressure
    real(kind_phys), intent(in)  :: clat(:)          ! latitudes(radians) for columns
    real(kind_phys), intent(in)  :: cappa(:,:)       ! ratio of dry air gas constant to specific heat at constant pressure
    real(kind_phys), intent(in)  :: cpair(:,:)       ! specific heat of dry air at constant pressure
    real(kind_phys), intent(in)  :: pmid(:,:)        ! mid-point pressure
    real(kind_phys), intent(in)  :: uwnd(:,:)        ! Zonal wind (m/s)
    real(kind_phys), intent(in)  :: vwnd(:,:)        ! Meridional wind (m/s)
    real(kind_phys), intent(in)  :: temp(:,:)        ! Temperature (K)
    !
    ! Output arguments
    !
    real(kind_phys),   intent(out) :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(out) :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(out) :: ds(:,:)   ! Heating rate
    character(len=64), intent(out) :: scheme_name
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg

    ! PKSRAT
    logical, intent(in) :: pkstrat  !pkstrat=.True. to use the PK02 TEQ
    real(kind_phys), intent(in) :: vgamma  !gamma parameter in PK02 (controling vortex strength)
    ! END-PKSTRAT

    ! PKSTRAT (Yi Wang, Sep 2024)
    logical, intent(in) :: noPV       ! whether to activate NO POLAR VOTEX case
    integer, intent(in) :: pret       ! lower limit (Pa) for upper level damping (sponge layer), p_sp parameter in PK02
    integer, intent(in) :: lat0       ! phi_0 parameter in PK02
    integer, intent(in) :: dellat     ! delta_phi parameter in PK02
    integer, intent(in) :: dely       ! delta_y parameter in PK02
    integer, intent(in) :: eps        ! epsilon parameter in PK02
    integer, intent(in) :: delz       ! delta_z parameter in PK02
    ! END-PKSTRAT (Yi Wang, Sep 2024)

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(kind_phys) :: kv            ! 1./efolding_time (normalized) for wind
    real(kind_phys) :: kt            ! 1./efolding_time for temperature diss.
    real(kind_phys) :: trefa         ! "radiative equilibrium" T
    real(kind_phys) :: trefc         ! used in calc of "radiative equilibrium" T
    real(kind_phys) :: cossq(ncol)   ! coslat**2
    real(kind_phys) :: cossqsq(ncol) ! coslat**4
    real(kind_phys) :: sinsq(ncol)   ! sinlat**2
    real(kind_phys) :: coslat(ncol)  ! cosine(latitude)

    ! PKSTRAT
    real(kind_phys) :: sinlat(ncol) !sin(latitude)
    real(kind_phys) :: deglat(ncol) ! Latitude in degrees for columns
    real(kind_phys) :: pi !pi for computing latitude in degrees from clat
    ! END-PKSTRAT

    !
    !-----------------------------------------------------------------------
    !

    ! PKSRAT
    real(kind_phys) :: vgammaperm !gamma in K/m (as opposed to K/km)
    real(kind_phys) :: w, fac1, fac2, delt !weight function (Eq A2 in PK02)

      ! US standard atmosphere params
      integer, parameter :: nstd=7 ! # of levels specified in US standard atmosphere
      real(kind_phys) :: pstd_norm(nstd) ! pressures of US standard atmosphere,normalized by PS =101325
      real(kind_phys) :: tstd(nstd) ! Temperutre of standard levels
      real(kind_phys) :: lapse(nstd) ! Lapse rate of standard levels
      real(kind_phys) :: tvalstd ! Standard atmosphere T at current level
      real(kind_phys) :: tpv ! Polar vortex temperature (Eqn 3 of Kushner and Polvani 2004)
      real(kind_phys) :: pbase,lapsebase,tbase ! US standard parameters at base of layer
      integer :: k2
      real(kind_phys) :: tstd100 ! standard atmosphere temperature at 100hPa

      pstd_norm=(/ 1., 0.223361, 0.0540330, 0.00856668, 0.00109456, 0.000660636, 3.90468e-05/)
      tstd=(/ 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65 /)
      lapse=(/ -0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002 /)

      pi=4.*DATAN(1.D0)
    ! END-PKSTRAT

    errmsg = ' '
    errflg = 0
    scheme_name = "HELD_SUAREZ"

    ! print*,pkstrat,vgamma,noPV,pret,lat0,dellat,dely,eps,delz

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)

      ! PKSTRAT
      deglat(i) = (clat(i)/pi)*180.
      sinlat(i) = sin(clat(i))
      ! END-PKSTRAT

    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !

    ! PK-STRAT
    if (pkstrat) then
       tstd100=216.65_kind_phys*(100._kind_phys/110.906_kind_phys)**(-1._kind_phys*rair*0._kind_phys/9.81_kind_phys) !US standard atmosphere at 100hPa
       vgammaperm=vgamma/1000._kind_phys ! convert gamma param into K/m

       !Set up US standard atmosphere
       do k=1,pver
         ! Determine US standard atmosphere params at base of layer
         if (pref_mid_norm(k).lt.pstd_norm(nstd)) then
           pbase=pstd_norm(nstd)
           lapsebase=lapse(nstd)
           tbase=tstd(nstd)
         else
           do k2=1,nstd-1
            if (pref_mid_norm(k).gt.pstd_norm(k2+1)) then
             pbase=pstd_norm(k2)
             lapsebase=lapse(k2)
             tbase=tstd(k2)
             exit
            endif
           end do
         end if

         !US standard atmosphere at level
         tvalstd=tbase*(pref_mid_norm(k)/pbase)**(-1.*(rair*lapsebase/gravit))


         do i = 1,ncol
           if (pmid(i,k).lt.1e4_kind_phys) then !pre < 100hPa
             tpv=tstd100*(pmid(i,k)/1e4_kind_phys)**(rair*vgammaperm/gravit)
             ! PKSTRAT (Yi Wang, Sep 2024)
             if (noPV) then
               w=0.
             else
               if (lat0.lt.0) then
                 w=0.5*(1-tanh( (deglat(i)-lat0)/dellat))
               else
                 w=0.5*(1+tanh( (deglat(i)-lat0)/dellat))
               end if
             end if
             ! END-PKSTRAT (Yi Wang, Sep 2024)
             trefa=(1-w)*tvalstd + w*tpv
           else
             fac1=tstd100
             delt=dely*sinsq(i) + eps*sinlat(i) + delz*log(pmid(i,k)/1e5_kind_phys)*cossq(i)
             fac2=(315._kind_phys-delt)*(pmid(i,k)/1e5_kind_phys)**cappa(i,k)
             if (fac1.gt.fac2) then
               trefa=fac1
             else
               trefa=fac2
             end if
           end if !end pre < 100hPa

           if (pref_mid_norm(k) > sigmab) then
             kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
             ds(i,k) = (trefa - temp(i,k))*kt*cpair(i,k)
           else
             ds(i,k) = (trefa - temp(i,k))*ka*cpair(i,k)
           end if

         end do
       end do

    else  !using Held-Suarez TEQ
    ! END-PKSTRAT

    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt      = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._kind_phys - (60._kind_phys * sinsq(i))
          trefa   = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/pref)))*(pmid(i,k)/pref)**cappa(i,k)
          trefa   = max(t00,trefa)
          ds(i,k) = (trefa - temp(i,k))*kt*cpair(i,k)
        end do
      else
        do i = 1, ncol
          trefc   = 315._kind_phys - 60._kind_phys*sinsq(i)
          trefa   = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/pref)))*(pmid(i,k)/pref)**cappa(i,k)
          trefa   = max(t00,trefa)
          ds(i,k) = (trefa - temp(i,k))*ka*cpair(i,k)
        end do
      end if
    end do

    ! PK-STRAT
    end if !end if pkstrat
    ! END-PKSTRAT

    !
    ! Add diffusion near the surface for the wind fields
    !
    du(:,:) = 0._kind_phys
    dv(:,:) = 0._kind_phys

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*uwnd(i,k)
          dv(i,k) = -kv*vwnd(i,k)
        end do
      end if

      ! PKSTRAT
      !Add sponge layer in upper levels (see Appendix of PK02)
      if (pkstrat) then
        do i = 1,ncol
          if (pmid(i,k).lt.pret) then
            kv=kfstrat*((pret-pmid(i,k))/pret)**(2.)
            du(i,k) = -kv*uwnd(i,k)
            dv(i,k) = -kv*vwnd(i,k)
           endif !pmid < pret
        end do
      endif !pkstrat
      ! END-PKSTRAT
    end do

  end subroutine held_suarez_1994_run

end module held_suarez_1994
