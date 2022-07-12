#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

outdir="./output/phi"
codeDir="."
tooldir="../tools"
pushd $codeDir
make
popd

mkdir -p $outdir/tmp

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 'brutsaert'
       rich_crit     = 1.0
       zeta_trans    = 0.5
/
 &most_nml
       x0 = 0.01, x1=100.0, nsamples = 200
/
EOF
$codeDir/phi.x > $outdir/tmp/stable_brutsaert.csv

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 2
       rich_crit     = 1.0
       zeta_trans    = 0.5
/
 &most_nml
       x0 = 0.01, x1=100.0, nsamples = 200
/
EOF
$codeDir/phi.x > $outdir/tmp/stable_2.csv

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 'brutsaert'
       rich_crit     = 1.0
       zeta_trans    = 0.5
/
 &most_nml
       x0 = -0.01, x1=-100.0, nsamples = 200
/
EOF
$codeDir/phi.x > $outdir/tmp/unstable_brutsaert.csv

cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 2
       rich_crit     = 1.0
       zeta_trans    = 0.5
/
 &most_nml
       x0 = -0.01, x1=-100.0, nsamples = 200
/
EOF
$codeDir/phi.x > $outdir/tmp/unstable_2.csv

    commonFlags="--label=label --xlog --xlim=0.01:100 --width=5 --aspect=0.5"

$tooldir/plot.py $commonFlags --x=zeta  --y='phi_m-1' --ylog $outdir/tmp/stable*.csv --save $outdir/stable_m.pdf
$tooldir/plot.py $commonFlags --x=zeta  --y='phi_h-1' --ylog $outdir/tmp/stable*.csv --save $outdir/stable_h.pdf
$tooldir/plot.py $commonFlags --x=-zeta --y='phi_m'   --ylim=-0.05:1.05  $outdir/tmp/unstable*.csv --save $outdir/unstable_m.pdf
$tooldir/plot.py $commonFlags --x=-zeta --y='phi_h'   --ylim=-0.05:1.05  $outdir/tmp/unstable*.csv --save $outdir/unstable_h.pdf
