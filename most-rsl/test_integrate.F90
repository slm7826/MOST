program test
  use integrate_mod, only: integrand, &
     integrate_trapezoid, integrate_midpoint, integrate_simpson, &
     integrate_romberg_trapezoid, integrate_romberg_midpoint, &
     integrate_romberg_midpoint_inv, integrate_romberg_midpoint_exp
  implicit none

  real :: a, b, rtol, exact

  rtol=1e-8

  a = 1.0; b=2.0; exact = 1.0
  call test_integrate(f0,'f = 1',a,b,rtol,exact)

  a = 1.0; b=2.0; exact = (b**3-a**3)/3.0
  call test_integrate(f2,'f = x**2',a,b,rtol,exact)

  a = 1.0; b=2.0; exact = (b**4-a**4)/4.0
  call test_integrate(f3,'f = x**3',a,b,rtol,exact)

  a = 1.0; b=10.0; exact = exp(-a)-exp(-b)
  call test_integrate(fexp,'f = exp(-x)',a,b,rtol,exact)

  a = 1.0; b=HUGE(1.0)*1e-4; exact = exp(-a)-exp(-b)
  call test_integrate(fexp,'f = exp(-x)',a,b,rtol,exact)

  a = 1.0; b=999.0; exact = exp(-a)-exp(-b)
  call test_integrate(fexp,'f = exp(-x)',a,b,rtol,exact)

  a = 1.0; b=999; exact = 0.2193839343955202736771637754601216490310472934069082075779
  call test_integrate(fexp1,'f = exp(-x)/x',a,b,rtol,exact)

  a = 1.0; b=HUGE(1.0)*1e-4; exact = 0.2193839343955202736771637754601216490310472934069082075779
  call test_integrate(fexp1,'f = exp(-x)/x',a,b,rtol,exact)

contains

subroutine test_integrate(f,description,a,b,rtol,exact)
  procedure(integrand) :: f
  character(*), intent(in) :: description
  real, intent(in) :: a,b   ! limits of integration
  real, intent(in) :: rtol  ! relative tolerance
  real, intent(in) :: exact ! exact value of integral

  integer :: ierr, nDeg, nSteps
  real :: s

  write(*,*)
  write(*,'("Integrating ",a," from ",g9.3," to ",g9.3,"; exact value=",g23.16)')trim(description),a,b, exact
  call integrate_trapezoid(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Trapezoid = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_midpoint(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Midpoint = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_simpson(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Simpson = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_romberg_trapezoid(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Romberg Trapezoid = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_romberg_midpoint(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Romberg Midpoint = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_romberg_midpoint_inv(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Romberg Midpoint Inv = ',s,abs(s-exact),ierr,nSteps,nDeg
  call integrate_romberg_midpoint_exp(f,a,b,rtol,s,ierr,nDeg=nDeg,nSteps=nSteps)
  write(*,100)'Romberg Midpoint Exp = ',s,abs(s-exact),ierr,nSteps,nDeg

100 format(a24,g23.16,'  error=',g11.4,'  ierr=',i1,'  N steps=',i9,'  Refinement degrees=',i2)
end subroutine test_integrate

real function f0(x) result(f)
    real, intent(in) :: x
    f = 1
end function f0

real function f1(x) result(f)
    real, intent(in) :: x
    f = x
end function f1

real function f2(x) result(f)
    real, intent(in) :: x
    f = x**2
end function f2

real function f3(x) result(f)
    real, intent(in) :: x
    f = x**3
end function f3

real function fexp(x) result(f)
    real, intent(in) :: x
    f = exp(-x)
end function fexp

real function fexp1(x) result(f)
    real, intent(in) :: x
    f = exp(-x)/x
end function fexp1

end program test

