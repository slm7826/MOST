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

echo 'Calculating integral I vs z_R ...'
tmpdir=$outdir/vs_z_R
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
# for x in -10.0 -5.0 -1.0 0.0 1.0 5.0 10.0
for x in 10.0 5.0 1.0 0.0 -1.0 -5.0 -10.0
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
       z_a=17.5, L_inv = $x
       var='z_R', x0 = 0.01, x1=35, nsamples = 200
/
EOF
    time $codeDir/evaluate_RSL_I.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done
echo 'Plotting integral I vs z_R/z_a ...'
commonFlags="--x=z_R/z_a --xlim=0.2:2 --label=1/L"
$tooldir/plot.py $commonFlags --y=rsl_integral_m                        --save=$outdir/I_m_vs_zR.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_m --ylog --ylim=1e-7:1e3 --save=$outdir/I_m_vs_zR_log.pdf  $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t                        --save=$outdir/I_t_vs_zR.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t --ylog --ylim=1e-7:1e3 --save=$outdir/I_t_vs_zR_log.pdf  $tmpdir/*.csv

echo 'Plotting integral I vs z_a/z_R ...'
commonFlags="--x=z_a/z_R --xlim=0:3 --label=1/L --xlim=0:5"
$tooldir/plot.py $commonFlags --y=rsl_integral_m                         --save=$outdir/I_m_vs_zzR.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_m  --ylog --ylim=1e-7:1e3 --save=$outdir/I_m_vs_zzR_log.pdf  $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t                         --save=$outdir/I_t_vs_zzR.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_integral_t  --ylog --ylim=1e-7:1e3 --save=$outdir/I_t_vs_zzR_log.pdf  $tmpdir/*.csv
