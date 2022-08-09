#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"

pushd "$codeDir"
make PROFILE.x
popd

outroot=output/profile/unstable
# t_atm=299.7
# t_sfc=294.6
#   z0m=3.8891
#    zR=63.799
#  wind=1.3591

# t_atm=269.30
# t_sfc=255.36
#   z0m=1.8156
#    zR=62.365
#  wind=3.0584

t_atm=300.35
t_sfc=303.12
  z0m=2.4474
   zR=85.361
 wind=2.6518

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 2
       rich_crit = 1.0
       zeta_trans =  0.5
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       ! interpolation parameters
!       a_min=1e-6, a_max=100,   a_nsteps=100
       b_min=-100,  b_max=1000,  b_nsteps=200
/
 &profile_nml
       z0m = $z0m
       k_over_B = 2.0
       u_atm = $wind
       t_sfc = $t_sfc
       t_atm = $t_atm
       zR    = $zR
       nsamples = 100
/
EOF
$codeDir/PROFILE.x
