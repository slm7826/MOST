program test
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
  logical :: neutral        = .false.
  real    :: zeta_trans     = 0.5
!   logical :: new_mo_option  = .false.

  character(32) :: rsl_option = 'none'
  real    :: rsl_mu_1 = 0.67 ! parameter of RSL correction
  real    :: rsl_mu_m = 2.59 ! parameter of RSL momentum correction
  real    :: rsl_mu_t = 0.95 ! parameter of RSL heat and tracer correction


  namelist /monin_obukhov_nml/ stable_option, rich_crit, neutral, drag_min_heat, &
                               drag_min_moist, drag_min_mom, zeta_trans, &
                               rsl_option, rsl_mu_1, rsl_mu_m, rsl_mu_t

  ! sampling parameters
  real    :: z_a           = 1.0 ! top of the constant flux layer above displacement height, m
  real    :: z_R           = 1.0 ! roughness sublayer thickness, m
  real    :: L_inv = 0.0 ! 1/L, reciprocal of Monin-Obukhov length
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /input_nml/ z_a, z_R, L_inv, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: i, ierr_m, ierr_t
  character(512) :: msg
  class(most_functions_T), pointer :: most
  class(rsl_functions_T),  pointer :: rsl
  real  :: x,a,b,s_m,s_t

  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading monin_obukhov_nml ::'//trim(msg)
     stop 1
  endif
  read (701, input_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading most_nml ::'//trim(msg)
     stop 1
  endif

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)
  write(*,input_nml)

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
     write (*,*)'stable_option = "'//trim(stable_option)//'" is incorrect'
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
     write (*,*)'rsl_option = "'//trim(rsl_option)//'" is incorrect'
     stop 1
  end select
  call most%set_rsl_functions(rsl)

  a = z_a/z_R
  b = z_R*L_inv

  write(*,'(a20," = ",g14.5)') "z_a/z_R",a
  write(*,'(a20," = ",g14.5)') "z_R/L",b
  write(*,'(a20," = ",g14.5)') "1/L", L_inv

  write(*,*)'RESULTS:'
  write(*,'(a)')'z_a,z_R,1/L,z_a/z_R,z_R/L,z_R/z_a,rsl_integral_m,rsl_integral_t,ierr_m,ierr_t'
  do i = 0,nsamples
     x = x0+i*(x1-x0)/nsamples
     select case(var)
     case ('z_a/z_R')
        a = x
        z_R = z_a/a
     case ('z_R/L')
        b = x
        L_inv = b/z_R
     case ('z_R')
        z_R = x
        a = z_a/z_R
        b = z_R*L_inv
     case ('1/L')
        a = z_a/z_R
        L_inv = x
        b = z_R*L_inv
     case default
        write (*,*) 'Incorrect sampling variable'
        stop 1
     end select
     call RSL_integral_I_m(most,a,b,s_m,ierr_m)
     call RSL_integral_I_t(most,a,b,s_t,ierr_t)
     write(*,'(99(g14.5,:,","))') z_a, z_R, L_inv, a, b, z_R/z_a, s_m, s_t, ierr_m, ierr_t
  enddo
end program test
