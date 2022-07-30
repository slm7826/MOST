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
  integer, parameter :: n = 1
  real    :: z_a(n)      = 1.0 ! top of the constant flux layer above displacement height, m
  real    :: z0m(n)      = 1.0 ! roughness length for momentum, m
  real    :: k_over_B(n) = 2.0 !
  real    :: z_R(n)      = 1.0 ! roughness sublayer thickness, m
  real    :: L_inv(n)    = 0.0 ! 1/L, reciprocal of Monin-Obukhov length
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /input_nml/ z_a, z0m, k_over_B, z_R, L_inv, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: i
  real    :: z0s(n) ! roughness length for heat, m
  real    :: z1m(n), z1s(n), z2(n) ! limits of RSL integral
  real    :: zRL(n) ! z_R/L
  real    :: zR_inv(n) ! 1/z_R
  real    :: ln_z_z0(n), ln_z_zs(n), zeta0m(n), zeta0s(n), zeta_a(n) ! for MOST integral
  logical, dimension(1) :: mask ! for MOST integral
  class(most_functions_T), pointer :: most
  class(rsl_functions_T),  pointer :: rsl

  character(512) :: msg
  real :: F_m(n), FF_m(n)
  real :: F_t(n), FF_t(n)
  integer :: ierr
  ! inputs
  real  :: x

  ! read namelists
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading monin_obukhov_nml ::'//trim(msg)
     stop 1
  endif
  read (701, input_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading input_nml ::'//trim(msg)
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

  z0s = z0m*exp(-k_over_B) ! roughness length for heat, m
  ln_z_z0 = log(z_a/z0m)
  ln_z_zs = log(z_a/z0s)

  write(*,'(a20," = ",g14.5)') "z0s",z0s
  write(*,'(a20," = ",g14.5)') "1/L", L_inv
  write(*,'(a20," = ",g14.5)') "zeta", z_a*L_inv
  zR_inv = 0.0
  where (z_R>0) &
     zR_inv = 1/z_R
  write(*,'(a20," = ",g14.5)') "z_a/z_R", z_a*zR_inv
  write(*,'(a20," = ",g14.5)') "z_R/z_a", z_R/z_a

  mask(:) = .TRUE.
  write(*,*)'RESULTS:'
  write(*,'(a)')'z_a,z0m,z0s,z_R,1/L,zeta,z_a/z_R,z_R/z_a,MO_integral_m,full_integral_m,rsl_m_ratio,MO_integral_t,full_integral_t,rsl_t_ratio'
  do i = 0,nsamples
     x = x0+i*(x1-x0)/nsamples
     select case(var)
     case ('1/L')
        L_inv = x
     case ('z_R')
        z_R = x
     case default
        write (*,*) 'Incorrect sampling variable'
        stop 1
     end select
     ! assign derived parameter values
     zR_inv = 0.0
     where (z_R>0) &
        zR_inv = 1/z_R
     z1m = z0m*zR_inv; z2 = z_a*zR_inv ! boundaries of RSL integral I_1
     z1s = z0s*zR_inv
     zRL = z_R*L_inv      ! z_R/L, parameter of RSL integral I_1
     zeta0m = z0m*L_inv   ! z0m/L, lower limit of MOST integral F
     zeta0s = z0s*L_inv   ! z0s/L, lower limit of MOST integral F
     zeta_a = z_a*L_inv   ! z_a/L, upper limit of MOST integral F

     call most%integral_m(n, mask, zeta_a, zeta0m, ln_z_z0, F_m, ierr)
     call most%integral_t(n, mask, zeta_a, zeta0s, ln_z_zs, F_t, ierr)

     FF_m = F_m
     FF_t = F_t
     call most%add_rsl_integral_m(n, mask, L_inv, z0m, z_a, z_R, FF_m, ierr=ierr)
     call most%add_rsl_integral_t(n, mask, L_inv, z0s, z_a, z_R, FF_t, ierr=ierr)

     write(*,'(99(g14.5,:,","))') &
         z_a, z0m, z0s, z_R, L_inv, z_a*L_inv, z_a*zR_inv, z_R/z_a,&
         F_m, FF_m, FF_m/F_m, &
         F_t, FF_t, FF_t/F_t
  enddo

end program test
