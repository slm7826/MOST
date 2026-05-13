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
  real   , parameter :: STEFAN  = 5.6734e-8 !< Stefan-Boltzmann constant [W/m^2/deg^4]

  real   , parameter :: day = 86400.0 ! seconds in day
  real   , parameter :: pi  = 3.14159265358979
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
  real :: cap_g     = 0.0     ! surface heat capacity, J/(m2 K)
  real :: cap_c     = 0.0     ! canopy air heat capacity, J/(m2 K)
  real :: gust_factor = 0.0   ! gustiness factor
  real :: con_g_h_factor = 0.1 ! ratio of ground-to-cana conductance to rho_drag
  integer :: n_gust_iter = 1  ! number of gustiness iterations

  namelist /idealized_nml/ p_atm, t_atm_ave, t_atm_range, t_atm_shift, &
      wind_atm, z_atm, z0m, k_over_B, &
      dt, swnet_max, lwdn, cap_g, cap_c, gust_factor, n_gust_iter, con_g_h_factor
  ! - end of namelists

  integer :: io
  integer :: n
  real    :: z0s
  real, dimension(1) :: rho, rho_drag, gust

  ! inputs
  real, dimension(1) :: Ta, Tc0, Tg0, z, z0, zt, zq, wind
  logical :: avail(1), lavail
  ! outputs
  real,    dimension(1) :: drag_m, drag_t, drag_q, u_star, b_star, rich, zeta
  integer :: ier, i

  real :: time ! seconds
  real :: Tc, delta_Tc ! canopy air temperature and its time step tendency
  real :: Tg, delta_Tg ! surface temperature and its time step tendency
  real :: lwup0, rnet0, shflx0 ! fluxes before implicit time step
  real :: lwup,  rnet,  shflx  ! final vales of the fluxes
  real :: DRDT ! derivalive of lwup wrt surface temperature
  real :: Ha0, DHaDTc
  real :: Hg0, DHgDTc, DHgDTg
  real :: swnet ! net shortwave

  integer,parameter :: iTc = 1, iTg = 2
  real :: A(2,2),B(2),det, con_g_h
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
  write(*,'(a)') 'Time,Tg,Tc,Ta,rnet,swnet,lwdn,lwup,lwnet,shflx,rnet0,lwup0,lwnet0,shflx0,rho,CD_t,CD_m,rho_CD_U,gust,wind,ustar,bstar,rich,zeta'

  time = 0.0
  Tg = t_atm_ave
  Tc = t_atm_ave
  do
     lwup0 = STEFAN*Tg**4
     DRDT  = 4*STEFAN*Tg**3
     swnet = swnet_max*max(-sin(2*pi*time/day),0.0)
     rnet0 = swnet + lwdn - lwup0

     if (time > 5*day) exit

     Tc0 = Tc; Tg0 = Tg
     Ta  = t_atm_ave - 0.5*t_atm_range*sin(2*pi*(time/day-t_atm_shift/24.0))

     b_star = 0.0; u_star = 0.0
     do i = 1,n_gust_iter
        ! calculate gustiness for the next iteration, using u_star and b_star that was
        ! just calculated
        where (b_star > 0.)
            gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
        else where
            gust = 0.
        end where
        wind = sqrt(wind_atm**2 + gust**2)

        call monin_obukhov_drag_1d(grav, vonkarm,               &
             & error, zeta_min, max_iter, small,                         &
             & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &!miz
             & drag_min_heat, drag_min_moist, drag_min_mom,              &
             & n, Ta, Tc0, z, z0, zt, zq, wind, drag_m, drag_t,         &
             & drag_q, u_star, b_star, rich, zeta, lavail, avail, ier)
     end do

     rho = p_atm / (rdgas * Ta(1)) ! density
     rho_drag = cp_air * drag_t * rho * wind_atm
     con_g_h  = rho_drag(1) * con_g_h_factor

     ! solve the linearized energy balance implicitly
     Ha0    = rho_drag(1) * (Tc0(1) - Ta(1))  ! sensible heat from canopy air to atmos, W/m2
     DHaDTc = rho_drag(1)
     Hg0    = con_g_h * (Tg0(1)-Tc0(1)) ! sensible heat flux from ground to canopy air, W/m2
     DHgDTg = con_g_h ; DHgDTc = -con_g_h

     ! form the matrix
     A(iTc,iTc) = cap_c/dt + DHaDTc - DHgDTc
     A(iTc,iTg) = -DHgDTg
     B(iTc)     = Hg0 - Ha0
     A(iTg,iTc) = DHgDTc
     A(iTg,iTg) = cap_g/dt + DHgDTg + DRDT
     B(iTc)     = rnet0 - Hg0
     ! solve the system
     det      = A(iTc,iTc)*A(iTg,iTg) - A(iTc,iTg)*A(iTg,iTc)
     delta_Tc = (B(iTc)*A(iTg,iTg) - B(iTg)*A(iTg,iTc))/det
     delta_Tg = (A(iTc,iTc)*B(iTg) - A(iTc,iTg)*B(iTc))/det

     ! updated values of the fluxes
     shflx = shflx0 + rho_drag(1) * delta_Tc
     lwup  = lwup0  + DRDT        * delta_Tg
     rnet  = swnet  + lwdn - lwup
     Tc    = Tc + delta_Tc
     Tg    = Tg + delta_Tg
     time  = time+dt
     write(*,'(99(g14.5,:,","))') time/day, Tg, Tc, Ta, rnet, swnet, lwdn, lwup, lwdn-lwup, shflx, &
            rnet0,lwup0,lwdn-lwup0,shflx0, &
            rho, drag_t, drag_m, rho_drag, gust, wind, u_star, b_star, rich, zeta

  enddo
end program test
