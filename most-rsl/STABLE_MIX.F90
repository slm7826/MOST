program test
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod
  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

  ! - namelist
  real :: rich  = 1.0   ! Richardson number
  real :: zR    = 0.0   ! roughness sublayer scale (m)
  real :: z_atm = 17.5  ! atmospheric height (m)
  ! sampling parameters
  character(8) :: var = '' ! variable to sample
  real    :: x0 ! start of x axis
  real    :: x1 ! end of x axis
  integer :: nsamples ! number of intervals
  namelist /stable_mix_nml/  z_atm, rich, zR, &
    var,x0,x1,nsamples
  ! - end of namelists
  integer :: ios
  integer :: i
  real    :: x
  character(512) :: message

  ! inputs
  real, dimension(1,1,1) :: Ri3, z3
  real, dimension(1,1)   :: zR2
  ! outputs
  real, dimension(1,1,1) :: mix3


  call monin_obukhov_init()
  open (701, file='input.nml')
  read (701, stable_mix_nml, iostat=ios, iomsg=message)
  if (ios/=0) then
     write(error_unit,'(a)')'Error reading stable_mix_nml : '//trim(message)
     stop 1
  endif
  close(701)
  write(*,stable_mix_nml)

  write(*,*)'RESULTS:'
  write(*,'(99(a14,:,","))')'zR','rich','z','mix'
  do i = 0, nsamples
    z3   = z_atm
    zR2  = zR
    Ri3  = rich
    x = x0+i*(x1-x0)/nsamples
    select case(var)
    case ('z_atm')
       z3 = x
    case ('rich')
       Ri3 = x
    case ('zR')
       zR2 = x
    case default
       write (*,*) 'Incorrect sampling variable'
       stop 1
    end select

    call stable_mix(Ri3, z3, zR2, mix3)

    write(*,'(99(g14.5,:,","))') zR2, Ri3, z3, mix3
  enddo
end program test
