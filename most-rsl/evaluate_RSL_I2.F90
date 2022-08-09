program test
  use iso_fortran_env, only : error_unit

  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

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

  ! sampling parameters
  character(8) :: var = '' ! variable to sample
  real    :: a = 1.0, b = 1.0
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /evaluate2_nml/ a, b, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: i, ierr_m, ierr_t
  character(512) :: msg
  class(most_functions_T), pointer :: most
  class(rsl_functions_T),  pointer :: rsl
  real  :: x,s_m,s_t,a_inv,b_inv

  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading monin_obukhov_nml ::'//trim(msg)
     stop 1
  endif
  read (701, evaluate2_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading evaluate2_nml ::'//trim(msg)
     stop 1
  endif

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)
  write(*,evaluate2_nml)

  ! set up stability corrections
  select case(trim(stable_option))
  case('1')
     most=>make_most1_functions(rich_crit)
  case('2')
     most=>make_most2_functions(rich_crit, zeta_trans)
  case('brutsaert')
     most=>make_brutsaert_functions(rich_crit)
  case('neutral')
     most => make_neutral_functions(rich_crit)
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

!   write(*,'(a20," = ",g14.5)') "a",a
  write(*,'(a20," = ",g14.5)') "1/a",1/a
!   write(*,'(a20," = ",g14.5)') "b",b
  write(*,'(a20," = ",g14.5)') "1/b",1/b

  write(*,*)'RESULTS:'
  write(*,'(99(a14,:,","))')'a','b','1/a','1/b','I_m','I_h','ierr_m','ierr_t'
  do i = 0,nsamples
     x = x0+i*(x1-x0)/nsamples
     select case(var)
     case ('a')
        a = x
     case ('a2')
        a = x**2
     case ('1/a')
        a = 1/x
     case ('b')
        b = x
     case ('b2')
        if (x<0) then
           b = -x**2
        else
           b =  x**2
        endif
     case default
        write (*,*) 'Incorrect sampling variable'
        stop 1
     end select
     call RSL_integral_I_m(most,a,b,s_m,ierr_m)
     call RSL_integral_I_t(most,a,b,s_t,ierr_t)
     a_inv = 0.0;  if (a/=0.0) a_inv = 1/a
     b_inv = 0.0;  if (b/=0.0) b_inv = 1/b
     write(*,'(99(g14.5,:,","))') a, b, a_inv, b_inv, s_m, s_t, ierr_m, ierr_t
  enddo
end program test
