program test
  use iso_fortran_env, only : error_unit

  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none
  ! constants
  real   , parameter :: grav = 9.81 ! m/s2
  real   , parameter :: vonkarm = 0.4
  integer, parameter :: max_iter = 20
  real   , parameter :: error=1.e-04, zeta_min=1.e-06, small=1.e-04
  real   , parameter :: RDGAS  = 287.04           !< Gas constant for dry air [J/kg/deg]
!   real   , parameter :: RVGAS  = 461.50           !< Gas constant for water vapor [J/kg/deg]
  real   , parameter :: KAPPA  = 2.0/7.0  !< RDGAS / CP_AIR [dimensionless]
  real   , parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]

  real   , parameter :: gust_zi = 1000.0 ! m, boundary layer depth for gustiness

  ! + namelists
  character(32) :: stable_option = '1'
  real    :: rich_crit      = 2.0
  real    :: drag_min_heat  = 1.e-05
  real    :: drag_min_moist = 1.e-05
  real    :: drag_min_mom   = 1.e-05
  real    :: zeta_trans     = 0.5

  character(32) :: rsl_option = 'none'
  real    :: rsl_mu_1 = 0.67 ! parameter of RSL correction
  real    :: rsl_mu_m = 2.59 ! parameter of RSL momentum correction
  real    :: rsl_mu_t = 0.95 ! parameter of RSL heat and tracer correction
  ! parameters of RSL integrals Im and It lookup tables
  logical :: use_RSL_lookup = .TRUE. !> use loookup tables to compute RSL integrals; otherwise
                                !! calculate integrals directly: this can be used for, say,
                                !! testing of quality of lookup in a single point runs, but would
                                !! very likely be prohibitively slow in global simulations
  real    :: a_min    = 1e-5    !> lower lookup table limit for parameter a of I_m and I_h RSL integrals: a_min > 0.
  real    :: a_max    = 100     !> upper lookup table limit for parameter a of I_m and I_h RSL integrals: a_max > a_min > 0.
  integer :: a_nsteps = 100     !> number of lookup table steps along the axis a.
  real    :: b_min    = -10.0   !> lower lookup table limit for parameter b of I_m and I_h RSL integrals.
  real    :: b_max    =  1000.0 !> upper lookup table limit for parameter b of I_m and I_h RSL integrals
  integer :: b_nsteps = 100     !> number of lookup table steps along the axis b.

  namelist /monin_obukhov_nml/ rich_crit, drag_min_heat, drag_min_moist, drag_min_mom, &
                               stable_option, zeta_trans, & !miz
                               rsl_option, rsl_mu_1, rsl_mu_m, rsl_mu_t, &
                               use_RSL_lookup, a_min, a_max, a_nsteps, b_min, b_max, b_nsteps

  real :: p_atm = 1e5
  real :: t_atm = 300.0
  real :: t_sfc = 295.0
  real :: u_atm = 5.0
  real :: z_atm = 17.5
  real :: z0m  = 0.1
  real :: k_over_B = 2.0
  real :: zR = 1.0 ! roughness sublayer thickness, m
  real :: gust_factor = 1.0
  ! sampling parameters
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /most_nml/ p_atm, t_atm, t_sfc, u_atm, z_atm, z0m, k_over_B, zR, gust_factor, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: n, i
  real    :: z0s
  class(most_functions_T), pointer :: most
  class(rsl_functions_T),  pointer :: rsl

  ! inputs
  real,    dimension(1) :: z, z0, zt, zq, zR1, x
  logical :: avail(1)
  ! outputs
  real,    dimension(1) :: &
      t_sfc1, t_atm1, u_atm1, rho, flux_t, flux_m, &
      cd_m, cd_t, cd_q, u_star, b_star, gust, rich, zeta, ga, ra
  ! for calculations of sensible flux derivatives:
  integer :: ier

  character(512) :: message

  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios, iomsg=message)
  if (ios/=0) then
     write (error_unit,'(a)')'Error reading monin_obukhov_nml : '//trim(message)
     stop 1
  endif
  read (701, most_nml, iostat=ios, iomsg=message)
  if (ios/=0) then
     write(error_unit,'(a)')'Error reading most_nml : '//trim(message)
     stop 1
  endif

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)
  write(*,most_nml)

  select case(trim(stable_option))
  case('1')
     most=>make_most1_functions(rich_crit)
  case('2')
     most=>make_most2_functions(rich_crit, zeta_trans)
  case('brutsaert')
     most=>make_brutsaert_functions(rich_crit)
  case default
     write (error_unit,*)'stable_option = "'//trim(stable_option)//'" is incorrect'
     stop 1
  end select

  ! set up roughness sublayer (RSL) corrections
  select case(trim(rsl_option))
  case('none')
     rsl=>NULL()
  case('ridder2010')
     rsl=>make_rsl_ridder2010_functions(rsl_mu_m,rsl_mu_t)
  case('ghannam2022')
     rsl=>make_rsl_ghannam2022_functions(rsl_mu_1,rsl_mu_m,rsl_mu_t)
  case default
     write (error_unit,*)'rsl_option = "'//trim(rsl_option)//'" is incorrect'
     stop 1
  end select

  write(error_unit, *)'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', 'Will set up RSL functions'
  call most%set_rsl_functions(rsl,use_RSL_lookup,a_min,a_max,a_nsteps,b_min,b_max,b_nsteps)
  write(error_unit, *)'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', 'Did set up RSL functions'

  write(*,*)'RESULTS:'
  write(*,'(a)')'zR,t_sfc,u_atm,flux_t,flux_m,cd_m,cd_t,cd_q,ga,ra,u_star,b_star,rich,zeta,gust,gust+u_atm,ier'
  z0s = z0m*exp(-k_over_B)
  n = 1; avail = .TRUE.
  u_atm1 = u_atm
  t_atm1 = t_atm
  t_sfc1 = t_sfc
  do i = 0, nsamples
     z = z_atm
     z0 = z0m
     zt = z0s
     zq = z0s
     zR1 = zR
     x = x0+i*(x1-x0)/nsamples
     select case(var)
     case ('u_atm')
        u_atm1 = x
     case ('t_sfc')
        t_sfc1 = x
     case ('zR')
        zR1 = x
     case default
        write (*,*) 'Incorrect sampling variable'
        stop 1
     end select

     rho = p_atm / (rdgas * t_atm) ! density

     call monin_obukhov_drag_1d(most, grav, vonkarm,                      &
          & error, zeta_min, max_iter, small,                             &
          & drag_min_heat, drag_min_moist, drag_min_mom,                  &
          & n, t_atm1, t_sfc1, z, z0, zt, zq, zR1, u_atm1, cd_m, cd_t,    &
          & cd_q, u_star, b_star, rich, zeta, ier, avail)
     ga     = cd_q * rho * abs(u_atm) ! conductance of constant flux layer for tracers
     ra     = 1.0/(max(ga,1e-6))      ! resistance of constant flux layer for tracers
     flux_t = cd_t * rho * abs(u_atm1) * cp_air * (t_sfc1 - t_atm1)  ! flux of sensible heat (W/m**2)
     flux_m = cd_m * rho * abs(u_atm1) * u_atm1                      ! flux of momentum (N/m2)

     ! calculate gustiness
     where (b_star > 0.)
         gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
     else where
         gust = 0.
     end where

    write(*,'(99(g14.5,:,","))') zR1, t_sfc1, u_atm1, flux_t, flux_m, cd_m, cd_t, cd_q, &
         ga, ra, u_star, b_star, rich, zeta, gust, sqrt(gust**2+u_atm1**2), ier
  enddo
end program test
