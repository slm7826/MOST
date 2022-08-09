program test
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod
  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

  ! sampling parameters
  real    :: z_a           = 1.0 ! top of the constant flux layer above displacement height, m
  real    :: z_R           = 1.0 ! roughness sublayer thickness, m
  real    :: L_inv = 0.0 ! 1/L, reciprocal of Monin-Obukhov length
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /evaluateI_nml/ z_a, z_R, L_inv, &
       var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: i, ierr_m, ierr_t
  character(512) :: msg
  real  :: x,a,b,s_m,s_t

  call monin_obukhov_init()

  open (701, file='input.nml')
  read (701, evaluateI_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(error_unit,*)'Error reading evaluateI_nml ::'//trim(msg)
     stop 1
  endif
  close(701)
  write(*,evaluateI_nml)

  a = z_a/z_R
  b = z_R*L_inv

  write(*,'(a20," = ",g14.5)') "z_a/z_R",a
  write(*,'(a20," = ",g14.5)') "z_R/L",b
  write(*,'(a20," = ",g14.5)') "1/L", L_inv

  write(*,*)'RESULTS:'
  write(*,'(99(a14,:,","))')'z_a','z_R','1/L','z_a/z_R','z_R/L','z_R/z_a','rsl_integral_m','rsl_integral_t','ierr_m','ierr_t'
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
