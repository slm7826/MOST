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


  ! Ca dTa/dt = Sa + Ha
  ! Cg dTg/dT = Sg - Lg - Ha
  ! where:
  !    Ta, Tg - (prognostic) temperatures of the atm and ground
  !    Ca, Cg are the heat capacities of lowest atm layer and the ground,
  !    Ha - calculated sensible heat flux between ground and the atm, positive upward
  !    Sa(Ta) - forcing applied to the atm layer
  !    Sg - net SW radiation to the ground
  !    Lg(Tg) - net LW radiation to the ground

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
  real :: tau_atm   = 1800.0  ! relaxation time for the atm layer, s
  real :: wind_atm  = 5.0     ! wind in the atmosphere, m/s
  real :: z_atm     = 17.5    ! height of the atmosphere, m
  real :: z0m       = 0.1     ! roughness length, m
  real :: k_over_B  = 2.0     ! ln(z0m/z0s)
  real :: dt        = 1800.0  ! time step, s
  real :: swnet_max = 200.0   ! max of downward short-wave, W/m2
  real :: lwdn      = 200.0   ! downward long-wave
  real :: Ca        = 4200.0  ! heat capacity of the atm layer, J/(m2 K)
  real :: Cg        = 4200.0  ! heat capacity of the ground, J/(m2 K)
  real :: gust_factor = 0.0   ! gustiness scaling factor, unitless
  integer :: n_gust_iter = 1  ! number of gustiness iterations

  namelist /idealized_nml/ p_atm, t_atm_ave, t_atm_range, t_atm_shift, tau_atm, &
      wind_atm, z_atm, z0m, k_over_B, &
      dt, swnet_max, lwdn, Ca, Cg, gust_factor, n_gust_iter
  ! - end of namelists

  integer :: io
  integer :: n
  real    :: z0s
  real, dimension(1) :: rho, rho_drag, gust

  ! inputs
  real, dimension(1) :: Ta, Ts0, z, z0, zt, zq, wind, gust0
  logical :: avail(1), lavail
  ! outputs
  real,    dimension(1) :: drag_m, drag_t, drag_q, u_star, b_star, rich, zeta
  integer :: ier, i

  real :: time ! seconds
  real :: Tg   ! ground surface temperature, degK
  real :: delta_Ta, delta_Tg ! tendencies of atm and ground surface temperature
  real :: Ta0 ! prescribed temperature of the atmosphere, degK
  real :: Sa0, DSaDTa ! forcing to the atmosphere (W/m2) and its derivative
  real :: Ha, Ha0, DHaDTa, DHaDTg ! sensible heat flux (W/m2) and its derivaties
  real :: Sg0 ! net short-wave to the ground, W/m2
  real :: Lg0, DLgDTg ! net LW to the ground (W/m2) and its derivative
  real :: Gamma, e, f ! elimination coefficients

  real :: lwup0, rnet0 ! fluxes before implicit time step
  real :: lwup  ! final vales of the fluxes
  ! read the namelists
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
  write(*,'(a)') 'Time,Tg,Ta,lwdn,lwup,lwnet,rnet0,lwup0,lwnet0,shflx0,rho,CD_t,CD_m,rho_CD_U,gust,wind,ustar,bstar,rich,zeta'

  time = 0.0
  Tg = t_atm_ave
  Ta = t_atm_ave
  do
     if (time > 5*day) exit

     Ts0 = Tg
     Ta0 = t_atm_ave - 0.5*t_atm_range*sin(2*pi*(time/day-t_atm_shift/24.0))

     do i = 1,n_gust_iter
        gust0 = gust
        wind = sqrt(wind_atm**2 + gust**2)

        call monin_obukhov_drag_1d(grav, vonkarm,               &
             & error, zeta_min, max_iter, small,                         &
             & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &!miz
             & drag_min_heat, drag_min_moist, drag_min_mom,              &
             & n, Ta, Ts0, z, z0, zt, zq, wind, drag_m, drag_t,         &
             & drag_q, u_star, b_star, rich, zeta, lavail, avail, ier)
        ! re-calculate gustiness for the next iteration, using u_star and b_star that was
        ! just calculated
        where (b_star > 0.)
            gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
        else where
            gust = 0.
        end where
     end do

     rho = p_atm / (rdgas * Ta(1)) ! density
     rho_drag = cp_air * drag_t * rho * wind_atm

     ! forcing to the atm
     Sa0    = Ca/tau_atm*(Ta0 - Ta(1))
     DSaDTa = - Ca/tau_atm
!      write(*,*) 'Sa0=', Sa0, 'Ta0=',Ta0, 'Ta=', Ta

     ! net ground SW, W/m2, positive to ground
     Sg0     = swnet_max*max(-sin(2*pi*time/day),0.0)

     ! net ground LW, W/m2, positive to ground
     lwup0   = STEFAN*Tg**4
     Lg0     = lwdn - lwup0 ! ner long-wave
     DLgDTg  = -4*STEFAN*Tg**3

!      write(*,*) 'Sa0=', Sa0, 'DSaDTa=', DSaDTa, 'Sg0=',Sg0, 'Lg0=',Lg0,'DLgDTg=',DLgDTg

     rnet0 = Sg0 + Lg0

    ! flux of sensible heat, W/m2, positive upwards
     Ha0     = rho_drag(1) * (Tg - Ta(1))
     DHaDTg  = + rho_drag(1)
     DHaDTa  = - rho_drag(1)
!      write(*,*) 'Ha0=',Ha0,'DHaDTg=',DHaDTg,'DHaDTa=',DHaDTa

     ! solve the linearized energy balance equations implicitly
     Gamma = Ca/dt - DSaDTa - DHaDTa
     f = (Sa0 + Ha0)/Gamma
     e = DHaDTg/Gamma

     delta_Tg = (Sg0 + Lg0 - Ha0 - DHaDTa*f)/ \
                (Cg/dt - DLgDTg + DHaDTg + DHaDTa*e)
     delta_Ta = f + e*delta_Tg

     ! updated values of the fluxes
     Ha    = Ha0 + DHaDTg * delta_Tg + DHaDTa * delta_Ta
     lwup  = lwup0  - DLgDTg * delta_Tg

     ! updated values of prognostic variables
     Ta = Ta + delta_Ta
     Tg = Tg + delta_Tg

     time  = time+dt
     write(*,'(99(g14.5,:,","))') time/day, Tg, Ta, lwdn, lwup, lwdn-lwup, &
            rnet0,lwup0,lwdn-lwup0,Ha, &
            rho, drag_t, drag_m, rho_drag, gust0, wind, u_star, b_star, rich, zeta

  enddo
end program test
