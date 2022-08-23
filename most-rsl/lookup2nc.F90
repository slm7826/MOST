program test
  use netcdf
  use iso_fortran_env, only : error_unit

  use monin_obukhov_mod
  use rsl_functions_mod
  use monin_obukhov_functions_mod
  use monin_obukhov_kernel

  implicit none

  integer :: ncid
  integer :: a_dimid, b_dimid
  integer :: a_varid, b_varid, im_varid, it_varid, pm_varid, pt_varid

  call monin_obukhov_init()

  call check(nf90_create('lookup.nc',nf90_clobber,ncid))
  call check(nf90_def_dim(ncid,'a',size(most%a),a_dimid))
  call check(nf90_def_dim(ncid,'b',size(most%b),b_dimid))

  call check(nf90_def_var(ncid,'a',nf90_double,[a_dimid],a_varid))
  call check(nf90_put_att(ncid,a_varid,'long_name','parameter a = z1/zR'))

  call check(nf90_def_var(ncid,'b',nf90_double,[b_dimid],b_varid))
  call check(nf90_put_att(ncid,b_varid,'long_name','parameter b = zR/L'))

  call check(nf90_def_var(ncid,'Im',nf90_double,[a_dimid,b_dimid],im_varid))
  call check(nf90_put_att(ncid,im_varid,'long_name','integral I for momentum'))
  call check(nf90_def_var(ncid,'It',nf90_double,[a_dimid,b_dimid],it_varid))
  call check(nf90_put_att(ncid,it_varid,'long_name','integral I for heat'))

  call check(nf90_def_var(ncid,'pm',nf90_double,[a_dimid],pm_varid))
  call check(nf90_put_att(ncid,pm_varid,'long_name','power for extrapolation of Im beyond lower limit of b'))
  call check(nf90_def_var(ncid,'pt',nf90_double,[a_dimid],pt_varid))
  call check(nf90_put_att(ncid,pt_varid,'long_name','power for extrapolation of It beyond lower limit of b'))

  call check(nf90_enddef(ncid))

  call check(nf90_put_var(ncid,a_varid,most%a))
  call check(nf90_put_var(ncid,b_varid,most%b))

  call check(nf90_put_var(ncid,im_varid,most%Im))
  call check(nf90_put_var(ncid,it_varid,most%It))

  call check(nf90_put_var(ncid,pm_varid,most%pm))
  call check(nf90_put_var(ncid,pt_varid,most%pt))

  call check(nf90_close(ncid))
contains
  subroutine check(ncerr)
    integer, intent(in) :: ncerr

    if (ncerr.eq.nf90_noerr) return
    write(error_unit,'(a)') nf90_strerror(ncerr)
    stop 1
  end subroutine check
end program test
