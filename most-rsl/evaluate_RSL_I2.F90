program test
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod
  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

  ! + namelists
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
  real  :: x,s_m,s_t,a_inv,b_inv


  call monin_obukhov_init()

  ! read namelists
  open (701, file='input.nml')
  read (701, evaluate2_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(*,*)'Error reading evaluate2_nml ::'//trim(msg)
     stop 1
  endif
  close(701)
  write(*,evaluate2_nml)


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
