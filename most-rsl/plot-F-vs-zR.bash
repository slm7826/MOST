#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars


codeDir="."
tooldir="../tools"
 outdir=output/integral_F

pushd "$codeDir"
make
popd

mkdir -p $outdir

echo 'Calculating integralR vs z_R ...'
tmpdir=$outdir/vs_zR
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
for x in -10.0 -5.0 -1.0 0.0 1.0 5.0 10.0
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
       z_a=17.5, z_R=$x, L_inv = $x, k_over_B = 0.0
!       var='z_R', x0 = 0.01, x1=30, nsamples = 200
       var='z_R', x0 = 0.01, x1=200, nsamples = 200
/
EOF
    time $codeDir/evaluate_RSL_F.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done

echo 'Plotting integralR vs z_R/z_a ...'
flags="--x=z_R/z_a --label=1/L"
$tooldir/plot.py $flags --y=full_integral_m        --save=$outdir/F_m_vs_zR.pdf       $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_m --ylog --save=$outdir/F_m_vs_zR_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $flags --y=rsl_m_ratio            --save=$outdir/F_m_ratio_vs_zR.pdf $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_t        --save=$outdir/F_t_vs_zR.pdf       $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_t --ylog --save=$outdir/F_t_vs_zR_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $flags --y=rsl_t_ratio            --save=$outdir/F_t_ratio_vs_zR.pdf $tmpdir/*.csv

echo 'Plotting integralR vs z_a/z_R ...'
flags="--x=z_R/z_a --label=1/L --xlim=0:12"
$tooldir/plot.py $flags --y=full_integral_m        --save=$outdir/F_m_vs_zzR.pdf       $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_m --ylog --save=$outdir/F_m_vs_zzR_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $flags --y=rsl_m_ratio            --save=$outdir/F_m_ratio_vs_zzR.pdf $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_t        --save=$outdir/F_t_vs_zzR.pdf       $tmpdir/*.csv
$tooldir/plot.py $flags --y=full_integral_t --ylog --save=$outdir/F_t_vs_zzR_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $flags --y=rsl_t_ratio            --save=$outdir/F_t_ratio_vs_zzR.pdf $tmpdir/*.csv
