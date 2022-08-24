program test
  use netcdf
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod, only : monin_obukhov_init, most
  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

  integer :: ncid
  integer :: a_dimid, b_dimid
  integer :: a_varid, b_varid, im0_varid, it0_varid, im1_varid, it1_varid
  integer :: a0_varid, b0_varid, it_varid, im_varid, pm_varid, pt_varid

  integer :: ios
  character(512) :: msg

  real, allocatable :: a(:), b(:), Im0(:,:), Im1(:,:), It0(:,:),  It1(:,:)
  real :: x, y, x0, x1, y0, y1
  integer :: i,j
  integer :: ierr_m, ierr_t

  real    :: a_min    = 1e-5    !> lower lookup table limit for parameter a of I_m and I_h RSL integrals: a_min > 0.
  real    :: a_max    = 100     !> upper lookup table limit for parameter a of I_m and I_h RSL integrals: a_max > a_min > 0.
  integer :: a_nsteps = 500     !> number of lookup table steps along the axis a.
  real    :: b_min    = -1000.0 !> lower lookup table limit for parameter b of I_m and I_h RSL integrals.
  real    :: b_max    =  1000.0 !> upper lookup table limit for parameter b of I_m and I_h RSL integrals
  integer :: b_nsteps = 500     !> number of lookup table steps along the axis b.
  namelist /lookup_test_2D_nml/ a_min, a_max, a_nsteps, b_min, b_max, b_nsteps

  ! read namelists
  open (701, file='input.nml')
  read (701, lookup_test_2D_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(error_unit,*)'Error reading lookup_test_2D_nml :: '//trim(msg)
     stop 1
  endif
  close(701)
  write(*,lookup_test_2D_nml)

  call monin_obukhov_init()

  allocate(a(a_nsteps+1))
  x0 = log(a_min); x1 = log(a_max)
  do i = 1,a_nsteps+1
     x = x0+(x1-x0)/a_nsteps*(i-1)
     a(i) = exp(x)
  enddo

  allocate(b(b_nsteps+1))
  ! NOTE: sign (a, b) returns the absolute value of a times the sign of b
  y0 = sign(abs(b_min)**(1./3.),b_min); y1 = sign(abs(b_max)**(1./3.),b_max)
  do j = 1,b_nsteps+1
     y = y0+(y1-y0)/b_nsteps*(j-1)
     b(j) = y**3
  enddo

  allocate(Im0(size(a),size(b)),Im1(size(a),size(b)))
  allocate(It0(size(a),size(b)),It1(size(a),size(b)))
  do j = 1, size(b)
  write(*,'(i5.4)',advance='NO') j
  do i = 1, size(a)
     call RSL_integral_R_m (most, a(i), HUGE(1.0), b(j), Im0(i,j), ierr_m)
     call RSL_integral_R_t (most, a(i), HUGE(1.0), b(j), It0(i,j), ierr_t)

     call RSL_lookup_I (most, a(i), b(j), most%Im, most%pm, Im1(i,j), ierr_m)
     call RSL_lookup_I (most, a(i), b(j), most%It, most%pt, It1(i,j), ierr_t)
  enddo
  enddo

  call check(nf90_create('lookup_test.nc',nf90_clobber,ncid))

  ! define the variables for the lookup table
  call check(nf90_def_dim(ncid,'lookup_a',size(most%a),a_dimid))
  call check(nf90_def_dim(ncid,'lookup_b',size(most%b),b_dimid))

  call check(nf90_def_var(ncid,'lookup_a',nf90_double,[a_dimid],a0_varid))
  call check(nf90_put_att(ncid,a0_varid,'long_name','parameter a = z1/zR'))

  call check(nf90_def_var(ncid,'lookup_b',nf90_double,[b_dimid],b0_varid))
  call check(nf90_put_att(ncid,b0_varid,'long_name','parameter b = zR/L'))

  call check(nf90_def_var(ncid,'lookup_Im',nf90_double,[a_dimid,b_dimid],im_varid))
  call check(nf90_put_att(ncid,im_varid,'long_name','lookup table for momentum'))
  call check(nf90_def_var(ncid,'lookup_It',nf90_double,[a_dimid,b_dimid],it_varid))
  call check(nf90_put_att(ncid,it_varid,'long_name','lookup table for heat'))

  call check(nf90_def_var(ncid,'lookup_pm',nf90_double,[a_dimid],pm_varid))
  call check(nf90_put_att(ncid,pm_varid,'long_name','power for extrapolation of Im beyond lower limit of b'))
  call check(nf90_def_var(ncid,'lookup_pt',nf90_double,[a_dimid],pt_varid))
  call check(nf90_put_att(ncid,pt_varid,'long_name','power for extrapolation of It beyond lower limit of b'))

  ! define variables for lookup quality evaluation
  call check(nf90_def_dim(ncid,'a',size(a),a_dimid))
  call check(nf90_def_dim(ncid,'b',size(b),b_dimid))

  call check(nf90_def_var(ncid,'a',nf90_double,[a_dimid],a_varid))
  call check(nf90_put_att(ncid,a_varid,'long_name','parameter a = z1/zR'))

  call check(nf90_def_var(ncid,'b',nf90_double,[b_dimid],b_varid))
  call check(nf90_put_att(ncid,b_varid,'long_name','parameter b = zR/L'))

  call check(nf90_def_var(ncid,'Im0',nf90_double,[a_dimid,b_dimid],im0_varid))
  call check(nf90_put_att(ncid,im0_varid,'long_name','integrated I for momentum'))
  call check(nf90_def_var(ncid,'Im1',nf90_double,[a_dimid,b_dimid],im1_varid))
  call check(nf90_put_att(ncid,im1_varid,'long_name','looked-up I for momentum'))

  call check(nf90_def_var(ncid,'It0',nf90_double,[a_dimid,b_dimid],it0_varid))
  call check(nf90_put_att(ncid,it0_varid,'long_name','integrated I for heat'))
  call check(nf90_def_var(ncid,'It1',nf90_double,[a_dimid,b_dimid],it1_varid))
  call check(nf90_put_att(ncid,it1_varid,'long_name','looked-up I for heat'))

  call check(nf90_enddef(ncid))

  ! write lookup table
  call check(nf90_put_var(ncid,a0_varid,most%a))
  call check(nf90_put_var(ncid,b0_varid,most%b))

  call check(nf90_put_var(ncid,im_varid,most%Im))
  call check(nf90_put_var(ncid,it_varid,most%It))

  call check(nf90_put_var(ncid,pm_varid,most%pm))
  call check(nf90_put_var(ncid,pt_varid,most%pt))

  ! write interpolation test
  call check(nf90_put_var(ncid,a_varid,a))
  call check(nf90_put_var(ncid,b_varid,b))

  call check(nf90_put_var(ncid,im0_varid,Im0))
  call check(nf90_put_var(ncid,it0_varid,It0))

  call check(nf90_put_var(ncid,im1_varid,Im1))
  call check(nf90_put_var(ncid,it1_varid,It1))

  call check(nf90_close(ncid))
  write(*,*)

contains
  subroutine check(ncerr)
    integer, intent(in) :: ncerr

    if (ncerr.eq.nf90_noerr) return
    write(error_unit,'(a)') nf90_strerror(ncerr)
    stop 1
  end subroutine check
end program test
