#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
pushd $codeDir
make
popd

# variables available in the output:
# Time,Ts,Ta,rnet,swnet,lwdn,lwup,lwnet,shflx,rnet0,lwup0,lwnet0,shflx0,rho,CD_t,CD_m,rho_CD_U,gust,wind,ustar,bstar,rich,zeta

outdir=output/SIM
tmpdir=$outdir

mkdir -p $outdir $tmpdir
rm -f $tmpdir/*.csv

cat <<EOF > input.nml
&monin_obukhov_nml
    stable_option = 2
    rich_crit     = 1.0
    zeta_trans    = 0.5
/
&idealized_nml
    p_atm       = 100000.0
    t_atm_ave   = 300.0000
    t_atm_range = 10.00000
    t_atm_shift = 0.000000
    wind_atm    = 2.000000
    z_atm       = 17.50000
    z0m         = 2.000000
    k_over_b    = 2.0
    dt          = 1800.000
    swnet_max   = 1000.000
    lwdn        = 400.0000
    cap         = 4200.0
    gust_factor = 0.0
/
EOF
$codeDir/SIM.x > $tmpdir/SIM.csv
commonFlags="--x=Time --label=var-name --xlab=None --xlim=0.5:2.5 --width=12 --aspect=0.2"
$tooldir/plot.py $commonFlags --var=Ts,Ta             --ylab=degK --save=$outdir/temp.pdf   $tmpdir/SIM.csv
$tooldir/plot.py $commonFlags --var=swnet,lwnet,shflx --ylab=W/m2 --save=$outdir/fluxes.pdf $tmpdir/SIM.csv
$tooldir/plot.py $commonFlags --var=CD_m,CD_t                     --save=$outdir/CD.pdf     $tmpdir/SIM.csv
# $tooldir/plot.py $commonFlags --var=ustar             --ylab=m/s  --save=$outdir/ustar.pdf  $tmpdir/SIM.csv
# $tooldir/plot.py $commonFlags --var=bstar                         --save=$outdir/bstar.pdf  $tmpdir/SIM.csv
$tooldir/plot.py $commonFlags --var=ustar,bstar       --ylab=ustar,bstar  --save=$outdir/ustar.pdf  $tmpdir/SIM.csv
$tooldir/plot.py $commonFlags --var=rich,zeta         --ylab=rich,zeta    --save=$outdir/rich.pdf   $tmpdir/SIM.csv
$tooldir/plot.py $commonFlags --var=wind,gust         --ylab=m/s    --save=$outdir/wind.pdf   $tmpdir/SIM.csv
