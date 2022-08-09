#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
 outdir=output/integral_IB

pushd "$codeDir"
make
popd

mkdir -p $outdir

echo 'Calculating integralI vs a ...'
tmpdir=$outdir/vs_b
mkdir -p $tmpdir
rm -f $tmpdir/*.csv
i=1
for x in 0.1 0.2 0.5 1.0 2.0 5.0
do
   cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  'brutsaert',
       rich_crit = 1.0,
       zeta_trans =  0.5
!       rsl_option = 'ridder2010'
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       use_RSL_lookup = F
/
 &evaluate2_nml
       var='b2', a = $x, x0 = -3.162, x1 = +3.162, nsamples = 100
/
EOF
    echo z_R = $x
    time $codeDir/evaluate_RSL_I2.x > $tmpdir/`printf "%2.2d" $i`.csv
    (( i++ ))
done
echo 'Plotting integralI vs b ...'

commonFlags="--x=b --label=a"
$tooldir/plot.py $commonFlags --y=I_m        --save=$outdir/Im_vs_b.pdf      $tmpdir/*.csv
$tooldir/plot.py $commonFlags --y=I_h        --save=$outdir/Ih_vs_b.pdf      $tmpdir/*.csv

# commonFlags="--x=b --label=b"
# $tooldir/plot.py $commonFlags --y=I_m        --save=$outdir/Im_vs_bInv.pdf   $tmpdir/*.csv
# $tooldir/plot.py $commonFlags --y=I_h        --save=$outdir/Ih_vs_bInv.pdf   $tmpdir/*.csv
