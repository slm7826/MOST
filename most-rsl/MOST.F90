program test
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod
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

  call monin_obukhov_init()
  open (701, file='input.nml')
  read (701, most_nml, iostat=ios, iomsg=message)
  if (ios/=0) then
     write(error_unit,'(a)')'Error reading most_nml : '//trim(message)
     stop 1
  endif
  close(701)
  write(*,most_nml)

  write(*,*)'RESULTS:'
  write(*,'(99(a14,:,","))')'zR','t_sfc','u_atm','flux_t','flux_m','cd_m','cd_t','cd_q','ga','ra',&
                'u_star','b_star','rich','zeta','gust','gust+u_atm','ier'
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
