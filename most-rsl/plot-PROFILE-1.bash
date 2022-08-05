#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"

pushd "$codeDir"
make
popd

outroot=output/profile/unstable
t_atm=299.7
t_sfc=294.6
  z0m=3.8891
   zR=63.799
 wind=1.3591

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5,
       rsl_option = "ghannam2022"
/
 &most_nml
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
