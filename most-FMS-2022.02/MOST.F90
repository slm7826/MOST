program test
  use monin_obukhov_inter

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

  real   , parameter :: delta_T = 0.01   ! increment for calculation of sensible flux derivative
  real   , parameter :: gust_zi = 1000.0 ! m, boundary layer depth for gustiness

  ! + namelists
  real    :: rich_crit      = 2.0
  real    :: drag_min_heat  = 1.e-05
  real    :: drag_min_moist = 1.e-05
  real    :: drag_min_mom   = 1.e-05
  logical :: neutral        = .false.
  integer :: stable_option  = 1
  real    :: zeta_trans     = 0.5
  logical :: new_mo_option  = .false.

  namelist /monin_obukhov_nml/ rich_crit, neutral, drag_min_heat, &
                               drag_min_moist, drag_min_mom,      &
                               stable_option, zeta_trans, new_mo_option !miz

  real :: p_atm = 1e5
  real :: t_atm = 300.0
  real :: t_sfc = 295.0
  real :: u_atm = 5.0
  real :: z_atm = 17.5
  real :: z0m  = 0.1
  real :: k_over_B = 2.0
  real :: gust_factor = 1.0
  ! sampling parameters
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /most_nml/ p_atm, t_atm, t_sfc, u_atm, z_atm, z0m, k_over_B, gust_factor, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: n, i
  real    :: z0s

  ! inputs
  real,    dimension(1) :: z, z0, zt, zq, x
  logical :: avail(1), lavail
  ! outputs
  real,    dimension(1) :: &
      t_sfc1, t_atm1, u_atm1, rho, flux_t, flux_m, &
      cd_m, cd_t, cd_q, u_star, b_star, gust, rich, zeta, ga, ra, deriv1
  ! for calculations of sensible flux derivatives:
  real, dimension(1) :: t_sfc0, flux_t0
  ! "deriv0" in the output is rho*Cp*Cd*|U|; "deriv1" is the difference in fluxes divided by
  ! differences in temperature. Therefore "deriv1" takes into account dependence of Cd on
  ! stability, while "deriv0" does not. Flux exchange actually uses "deriv0".
  integer :: ier

  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios)
  if (ios/=0) stop 'Error reading monin_obukhov_nml'
  read (701, most_nml,          iostat=ios)
  if (ios/=0) stop 'Error reading most_nml'

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)
  write(*,most_nml)
  write(*,*)'RESULTS:'
  write(*,'(a)') 't_sfc,u_atm,flux_t,flux_m,cd_m,cd_t,cd_q,ga,ra,u_star,b_star,rich,zeta,gust,gust+u_atm,deriv0,deriv1'

  z0s = z0m*exp(-k_over_B)
  n = 1; lavail=.TRUE.; avail = .TRUE.
  u_atm1 = u_atm
  t_atm1 = t_atm
  t_sfc1 = t_sfc
  do i = 0, nsamples
     z = z_atm
     z0 = z0m
     zt = z0s
     zq = z0s
     x = x0+i*(x1-x0)/nsamples
     select case(var)
     case ('u_atm')
        u_atm1 = x
     case ('t_sfc')
        t_sfc1 = x
     case default
        write (*,*) 'Incorrect sampling variable'
        stop 1
     end select

     rho = p_atm / (rdgas * t_atm) ! density

     ! calculate parameters for T_sfc - delta_T
     t_sfc0 = t_sfc1 - delta_T
     call monin_obukhov_drag_1d(grav, vonkarm,                            &
          & error, zeta_min, max_iter, small,                             &
          & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &!miz
          & drag_min_heat, drag_min_moist, drag_min_mom,                  &
          & n, t_atm1, t_sfc0, z, z0, zt, zq, u_atm1, cd_m, cd_t,         &
          & cd_q, u_star, b_star, rich, zeta, lavail, avail, ier)
     flux_t0 = cd_t * rho * abs(u_atm1) * cp_air * (t_sfc0 - t_atm1)  ! flux of sensible heat (W/m**2)

     call monin_obukhov_drag_1d(grav, vonkarm,                            &
          & error, zeta_min, max_iter, small,                             &
          & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &!miz
          & drag_min_heat, drag_min_moist, drag_min_mom,                  &
          & n, t_atm1, t_sfc1, z, z0, zt, zq, u_atm1, cd_m, cd_t,         &
          & cd_q, u_star, b_star, rich, zeta, lavail, avail, ier)
     ga     = cd_q * rho * abs(u_atm1) ! conductance of constant flux layer for tracers
     ra     = 1.0/(max(ga,1e-6))       ! resistance of constant flux layer for tracers
     flux_t = cd_t * rho * abs(u_atm1) * cp_air * (t_sfc1 - t_atm1)  ! flux of sensible heat (W/m**2)
     flux_m = cd_m * rho * abs(u_atm1) * u_atm1                      ! flux of momentum (N/m2)
     deriv1 = (flux_t - flux_t0)/delta_T ! derivative of sensible heat flux wrt surface temperature

     ! calculate gustiness
     where (b_star > 0.)
         gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
     else where
         gust = 0.
     end where

    write(*,'(99(g14.5,:,","))') t_sfc1, u_atm1, flux_t, flux_m, cd_m, cd_t, cd_q, &
         ga, ra, u_star, b_star, rich, zeta, gust, sqrt(gust**2+u_atm1**2), &
         cd_t * rho * abs(u_atm1) * cp_air, deriv1
  enddo
end program test
