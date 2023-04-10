#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"

pushd "$codeDir"
make
popd

outroot=output/stable_mix
rich=1.0
# for z in 17.5, 100, 500, 1000; do
for z in 17.5 50 100; do
    outdir=$outroot/vs-rich/zAtm${z}; tmpdir=$outdir/tmp

    echo "Calculating ${outdir}..."
    mkdir -p $outdir $tmpdir; rm -f $tmpdir/*.csv
    i=1
#     for zR in 0 1.0 2.0 5.0 10.0 20.0
    for zR in 0 5.0 20.0 50.0
    do
      cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5,
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       use_RSL_lookup = F
/
 &stable_mix_nml
       rich  = $rich
       zR    = $zR
       z_atm = $z
       var='rich', x0=0.1, x1=1.0, nsamples = 100
/
EOF
        $codeDir/STABLE_MIX.x > $tmpdir/`printf "%2.2d" $i`.csv
        (( i++ ))
    done

    echo "Plotting ${outdir}..."
    commonFlags="--width=4.5 --aspect=1.0 --xlim=0.0:1.0 --ylim=0.0:0.3 --x=rich --y=mix --label=zR,z_atm"
    $tooldir/plot.py $commonFlags --save "$outdir/mix-vs-rich.pdf" $tmpdir/*.csv
done

for zR in 5.0 20.0 50.0
do
    outdir=$outroot/vs-z/zR${zR}; tmpdir=$outdir/tmp

    echo "Calculating ${outdir}..."
    mkdir -p $outdir $tmpdir; rm -f $tmpdir/*.csv
    i=1
    for rich in 0.2 0.5 0.75 1.0
    do
      cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5,
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       use_RSL_lookup = F
/
 &stable_mix_nml
       rich  = $rich
       zR    = $zR
       z_atm = $z
       var='z_atm', x0=10, x1=100, nsamples = 100
/
EOF
        $codeDir/STABLE_MIX.x > $tmpdir/`printf "%2.2d" $i`.csv
        (( i++ ))
    done

    echo "Plotting ${outdir}..."
    commonFlags="--width=4.5 --aspect=1.0 --ylim=0.0:100.0 --xlim=0.0:0.3 --x=mix --y=z --label=rich,zR"
    $tooldir/plot.py $commonFlags --save "$outdir/mix-vs-z.pdf" $tmpdir/*.csv
done
