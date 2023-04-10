program test
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
  real   , parameter :: STEFAN  = 5.6734e-8 !< Stefan-Boltzmann constant [W/m^2/deg^4]

  real   , parameter :: day = 86400.0 ! seconds in day
  real   , parameter :: pi  = 3.14159265358979
  real   , parameter :: gust_zi = 1000.0 ! m, boundary layer depth for gustiness


  ! + namelists
  character(32) :: stable_option = '1'
  real    :: rich_crit      = 2.0
  real    :: drag_min_heat  = 1.e-05
  real    :: drag_min_moist = 1.e-05
  real    :: drag_min_mom   = 1.e-05
  logical :: neutral        = .false.
  real    :: zeta_trans     = 0.5

  namelist /monin_obukhov_nml/ stable_option, rich_crit, neutral, drag_min_heat, &
                               drag_min_moist, drag_min_mom, zeta_trans

  real :: p_atm     = 1.0e5   ! surface pressure, Pa
  real :: t_atm_ave   = 300.0 ! mean atmos T, K
  real :: t_atm_range = 10.0  ! range of atmos , K
  real :: t_atm_shift = 3.0   ! atmos T phase shift, hrs
  real :: wind_atm  = 5.0     ! wind in the atmosphere, m/s
  real :: z_atm     = 17.5    ! height of the atmosphere, m
  real :: z0m       = 0.1     ! roughness length, m
  real :: k_over_B  = 2.0     ! ln(z0m/z0s)
  real :: dt        = 1800.0  ! time step, s
  real :: swnet_max = 200.0   ! max of downward short-wave, W/m2
  real :: lwdn      = 200.0   ! downward long-wave
  real :: cap       = 0.0     ! surface heat capacity, J/(m2 K)
  real :: gust_factor = 0.0   ! gustiness factor

  namelist /idealized_nml/ p_atm, t_atm_ave, t_atm_range, t_atm_shift, &
      wind_atm, z_atm, z0m, k_over_B, &
      dt, swnet_max, lwdn, cap, gust_factor
  ! - end of namelists

  class(most_functions_T), pointer :: most

  integer :: io
  integer :: n
  real    :: z0s
  real, dimension(1) :: rho, rho_drag, gust

  ! inputs
  real, dimension(1) :: Ta, Ts0, z, z0, zt, zq, wind
  logical :: avail(1), lavail
  ! outputs
  real,    dimension(1) :: drag_m, drag_t, drag_q, u_star, b_star, rich, zeta
  integer :: ier

  real :: time ! seconds
  real :: Ts, delta_Ts ! surface temperature and its time step tendency
  real :: lwup0, rnet0, shflx0 ! fluxes before implicit time step
  real :: lwup,  rnet,  shflx  ! final vales of the fluxes
  real :: DRDT ! derivalive of lwup wrt surface temperature
  real :: swnet ! net shortwave
  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=io)
  read (701, idealized_nml,     iostat=io)

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)
  write(*,idealized_nml)

  z0s = z0m*exp(-k_over_B)
  z = z_atm
  z0 = z0m
  zt = z0s
  zq = z0s
  n = size(Ta); lavail=.TRUE.; avail = .TRUE.
  gust = 0.0

  write(*,*) 'RESULTS:'
  write(*,'(99(a14,:,","))') 'Time','Ts','Ta','rnet','swnet','lwdn','lwup','lwnet','shflx',&
            'rnet0','lwup0','lwnet0','shflx0','rho','CD_t','CD_m','rho_CD_U','gust','wind',&
            'ustar','bstar','rich','zeta'

  select case(trim(stable_option))
  case('1')
     most=>make_most1_functions(rich_crit)
  case('2')
     most=>make_most2_functions(rich_crit, zeta_trans)
  case('brutsaert')
     most=>make_brutsaert_functions(rich_crit)
  case default
     write (*,*)'stable_option = "'//trim(stable_option)//'" is incorrect'
     stop 1
  end select

  time = 0.0
  Ts = t_atm_ave
  do
     lwup0 = STEFAN*Ts**4
     DRDT  = 4*STEFAN*Ts**3
     swnet = swnet_max*max(-sin(2*pi*time/day),0.0)
     rnet0 = swnet + lwdn - lwup0

     if (time > 5*day) exit

     Ts0 = Ts
     Ta  = t_atm_ave - 0.5*t_atm_range*sin(2*pi*(time/day-t_atm_shift/24.0))

     wind = sqrt(wind_atm**2 + gust**2)

     call monin_obukhov_drag_1d(most, grav, vonkarm,                  &
          & error, zeta_min, max_iter, small,                         &
          & drag_min_heat, drag_min_moist, drag_min_mom,              &
          & n, Ta, Ts0, z, z0, zt, zq, wind, drag_m, drag_t,          &
          & drag_q, u_star, b_star, rich, zeta, ier, avail)
     rho = p_atm / (rdgas * Ta(1)) ! density
     rho_drag = cp_air * drag_t * rho * wind_atm
     ! solve the linearized energy balance implicitly
     shflx0 = rho_drag(1) * (Ts - Ta(1))  ! flux of sensible heat (W/m**2)
     delta_Ts = (rnet0 - shflx0)/(cap/dt+rho_drag(1)+DRDT)
     ! updated values of the fluxes
     shflx = shflx0 + rho_drag(1) * delta_Ts
     lwup  = lwup0  + DRDT        * delta_Ts
     rnet  = swnet  + lwdn - lwup
     Ts    = Ts + delta_Ts
     time  = time+dt
     write(*,'(99(g14.5,:,","))') time/day, Ts, Ta, rnet, swnet, lwdn, lwup, lwdn-lwup, shflx, &
            rnet0,lwup0,lwdn-lwup0,shflx0, &
            rho, drag_t, drag_m, rho_drag, gust, wind, u_star, b_star, rich, zeta

     ! calculate gustiness
     where (b_star > 0.)
         gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
     else where
         gust = 0.
     end where

  enddo
end program test
