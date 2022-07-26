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

echo 'Calculating integralR vs L ...'
tmpdir=$outdir/vs_L
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
# for x in 0.1 0.5 1.0 2.0 5.0 10.0 20.0
for x in 17.5 8.75 4.0 2.0 1.0
do
   echo z_R = $x
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
       var='1/L', x0 = -10, x1=10, nsamples = 200
/
EOF
    time $codeDir/evaluate_RSL_F.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done

echo 'Plotting integralR vs L ...'
commonFlags="--x=1/L --label=z_R,z_R/z_a,z_a/z_R"
$tooldir/plot.py $commonFlags --y=full_integral_m         --save=$outdir/F_m_vs_L.pdf       $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=full_integral_m  --ylog --save=$outdir/F_m_vs_L_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_m_ratio             --save=$outdir/F_m_ratio_vs_L.pdf $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=full_integral_t         --save=$outdir/F_t_vs_L.pdf       $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=full_integral_t  --ylog --save=$outdir/F_t_vs_L_log.pdf   $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=rsl_t_ratio             --save=$outdir/F_t_ratio_vs_L.pdf $tmpdir/*.csv
