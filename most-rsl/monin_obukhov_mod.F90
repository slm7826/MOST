module monin_obukhov_mod

use iso_fortran_env, only : error_unit

use rsl_functions_mod
use monin_obukhov_functions_mod
use monin_obukhov_kernel

implicit none

!---- namelist
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
!---- end of namelist

!---- module variables
class(most_functions_T), pointer :: most
class(rsl_functions_T),  pointer :: rsl

contains ! ==============================================================================

subroutine monin_obukhov_init()

  integer        :: ios ! i/o status
  character(512) :: msg ! error message

  ! read namelist
  open (701, file='input.nml')
  read (701, monin_obukhov_nml, iostat=ios, iomsg=msg)
  if (ios/=0) then
     write(error_unit,*)'Error reading monin_obukhov_nml ::'//trim(msg)
     stop 1
  endif
  close(701)

  write(*,*)'SETTINGS:'
  write(*,monin_obukhov_nml)

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
  write(error_unit, *)'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD :: Setting up RSL functions'
  call most%set_rsl_functions(rsl,use_RSL_lookup,a_min,a_max,a_nsteps,b_min,b_max,b_nsteps)
  write(error_unit, *)'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD :: Did set up RSL functions'
end subroutine monin_obukhov_init

!=======================================================================
subroutine stable_mix(rich, z_ag, rsl_scale, mix)

real, intent(in) , dimension(:,:,:)  :: rich      ! Richardson number
real, intent(in) , dimension(:,:,:)  :: z_ag      ! height above ground, m
real, intent(in) , dimension(:,:)    :: rsl_scale ! roughness sublayer scale, m
real, intent(out), dimension(:,:,:)  :: mix       !

integer :: i,j,k,n
integer :: ier ! error code returned by most%stable_mix
real, dimension(size(rich,1),size(rich,2),size(rich,3)) :: &
    pm2, & ! RSL correction for momentum squared
    Ri     ! Richardson number scaled with RSL corrections
real :: pt ! RSL correction for heat

! if(size(rich,3).ne.size(z_ag,3)) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
!      'vertical sizes of "rich" ('//string(size(rich,3))//') and "z_ag" (' &
!      //string(size(z_ag,3))//') are inconsistent', FATAL)

n = size(rich,1)*size(rich,2)*size(rich,3)

if(associated(most%rsl)) then
  ! scale Richardson number with roughness sublayer corrections, where necessary
  do j = 1,size(rich,2)
  do i = 1,size(rich,1)
    if (rsl_scale(i,j)>0.0) then
      do k = 1,size(rich,3)
        pm2(i,j,k) = most%rsl%rsl_m(z_ag(i,j,k)/rsl_scale(i,j))**2
        pt         = most%rsl%rsl_t(z_ag(i,j,k)/rsl_scale(i,j))
        Ri (i,j,k) = rich(i,j,k) * pm2(i,j,k)/pt
      enddo
    else
      pm2(i,j,:) = 1.0
      Ri (i,j,:) = rich(i,j,:)
    endif
  enddo
  enddo

  call most%stable_mix(n, Ri, mix, ier)

  ! scale mixing factor with roughness sublayer correction
  mix(:,:,:) = mix(:,:,:)/pm2(:,:,:)
else
  call most%stable_mix(n, rich, mix, ier)
endif

! if (ier.ne.0) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
!      'stable_mix calculations for stable_option "'//trim(stable_option)//'" returned an error', FATAL)

end subroutine stable_mix
end module
