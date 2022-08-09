#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
 outdir=output/integral_I2

pushd "$codeDir"
make
popd

mkdir -p $outdir

echo 'Calculating integralI vs a ...'
tmpdir=$outdir/vs_a
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
for x in 10.0 5.0 1.0 0.0 -1.0 -5.0 -10.0
do
   cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5
!       rsl_option = 'ridder2010'
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       use_RSL_lookup = F
/
 &evaluate2_nml
!       var='1/a', b = $x, x0 = 0.1, x1 = 100, nsamples = 200
!       var='a', b = $x, x0 = 0.1, x1 = 100, nsamples = 500
       var='a2', b = $x, x0 = 0.1, x1 = 3.162, nsamples = 100
/
EOF
    echo z_R = $x
    time $codeDir/evaluate_RSL_I2.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done
echo 'Plotting integralI vs a ...'

commonFlags="--x=a --label=b"
$tooldir/plot.py $commonFlags --y=I_m        --save=$outdir/Im_vs_a.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=I_h        --save=$outdir/Ih_vs_a.pdf      $tmpdir/*.csv

commonFlags="--x=1/a --label=b"
$tooldir/plot.py $commonFlags --y=I_m        --save=$outdir/Im_vs_aInv.pdf   $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=I_h        --save=$outdir/Ih_vs_aInv.pdf   $tmpdir/*.csv
