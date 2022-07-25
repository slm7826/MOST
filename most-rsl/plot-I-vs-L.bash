#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
 outdir=output/integral_I

pushd "$codeDir"
make
popd

mkdir -p $outdir

echo 'Calculating integralI vs 1/L ...'
tmpdir=$outdir/vs_L
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
# for x in 1.0 2.0 4.0 8.75 17.5
for x in 17.5 8.75 4.0 2.0 1.0
do
   cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5
!       rsl_option = 'ridder2010'
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
/
 &input_nml
       z_a=17.5, z_R=$x, L_inv = 1.0
       var='1/L', x0 = -10, x1=10, nsamples = 500
/
EOF
    echo z_R = $x
    time $codeDir/evaluate_RSL_I.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done
echo 'Plotting integralI vs 1/L ...'

commonFlags="--x=1/L --label=z_R/z_a,z_a/z_R"
$tooldir/plot.py $commonFlags --y=rsl_integral_m        --save=$outdir/I_m_vs_L.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_m --ylog --save=$outdir/I_m_vs_L_log.pdf  $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t        --save=$outdir/I_t_vs_L.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t --ylog --save=$outdir/I_t_vs_L_log.pdf  $tmpdir/*.csv
