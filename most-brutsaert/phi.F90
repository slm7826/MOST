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
  real   , parameter :: RVGAS  = 461.50           !< Gas constant for water vapor [J/kg/deg]
  real   , parameter :: KAPPA  = 2.0/7.0  !< RDGAS / CP_AIR [dimensionless]
  real   , parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]

  real   , parameter :: delta_T = 0.01
  real   , parameter :: gust_zi = 1000.0 ! m, boundary layer depth for gustiness

  ! + namelists
  character(32) :: stable_option = '1'
  real    :: rich_crit      = 2.0
  real    :: drag_min_heat  = 1.e-05
  real    :: drag_min_moist = 1.e-05
  real    :: drag_min_mom   = 1.e-05
  logical :: neutral        = .false.
  real    :: zeta_trans     = 0.5
  logical :: new_mo_option  = .false.

  namelist /monin_obukhov_nml/ stable_option, rich_crit, neutral, drag_min_heat, &
                               drag_min_moist, drag_min_mom, zeta_trans
  ! sampling parameters
  real    :: x0 = 0.1       ! start of x axis
  real    :: x1 = 10.0      ! end of x axis
  integer :: nsamples = 200 ! number of intervals
  namelist /most_nml/ x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: n, i
  class(most_functions_T), pointer :: most

  ! inputs
  logical :: mask(1)
  ! outputs
  real, dimension(1) :: zeta, phi_m, phi_h
  real :: log0, log1, logz, sign
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

  select case(trim(stable_option))
  case('1')
     most=>make_most1_functions(rich_crit)
     write(*,*)'label = stable_option_1'
  case('2')
     most=>make_most2_functions(rich_crit, zeta_trans)
     write(*,*)'label = stable_option_2'
  case('brutsaert')
     most=>make_brutsaert_functions(rich_crit)
     write(*,*)'label = Brutsaert(2005)'
  case default
     write (*,*)'stable_option = "'//trim(stable_option)//'" is incorrect'
     stop 1
  end select
  if (x0*x1<0) then
     write(*,*) 'x0 and x1 must be of the same sign'
     stop 1
  endif

  if (x1>0) then
     sign = 1.0
  else
     sign = -1.0
  endif
  x0 = abs(x0); x1 = abs(x1)
  write(*,*)'RESULTS:'
  write(*,'(a)') 'zeta,-zeta,phi_m,phi_h,phi_m-1,phi_h-1'
  n = 1; mask = .TRUE.
  log0 = log(x0); log1=log(x1)
  do i = 0, nsamples
     logz = log0+i*(log1-log0)/nsamples
     zeta = sign*exp(logz)
     call most%derivative_m(n,mask,zeta,phi_m,ier)
     call most%derivative_t(n,mask,zeta,phi_h,ier)

     write(*,'(99(g14.5,:,","))') zeta, -zeta, phi_m, phi_h, phi_m-1,phi_h-1
  enddo
end program test
