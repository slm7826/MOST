#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"

pushd "$codeDir"
make
popd

outroot=output/profile/unstable
t_sfc=301.0
for z0m in 1.0 0.1; do
    for wind in 1.0 5.0; do
        outdir=$outroot/z0m${z0m}w${wind}; tmpdir=$outdir/tmp

        echo "Calculating ${outdir}..."
        mkdir -p $outdir $tmpdir; rm -f $tmpdir/*.csv
        i=1
        for zR in 0 1.0 2.0 5.0 10.0 20.0
        do
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
       t_atm = 300.0
       zR    = $zR
       nsamples = 100
/
EOF
            $codeDir/PROFILE.x > $tmpdir/`printf "%2.2d" $i`.csv
            (( i++ ))
        done

        echo "Plotting ${outdir}..."
        for var in t u
        do
            commonFlags="--width=5 --aspect=1.2 --x=$var --y=z --ylim=0:17.5 --label=zR,zeta,flux_t"
            $tooldir/plot.py $commonFlags --save "$outdir/$var.pdf" $tmpdir/*.csv
        done
    done
done

outroot=output/profile/stable
t_sfc=299.0
for z0m in 1.0 0.1; do
    for wind in 1.0 5.0; do
        outdir=$outroot/z0m${z0m}w${wind}; tmpdir=$outdir/tmp

        echo "Calculating ${outdir}..."
        mkdir -p $outdir $tmpdir; rm -f $tmpdir/*.csv
        i=1
        for zR in 0 1.0 2.0 5.0 10.0 20.0
        do
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
       t_atm = 300.0
       zR    = $zR
       nsamples = 100
/
EOF
            $codeDir/PROFILE.x > $tmpdir/`printf "%2.2d" $i`.csv
            (( i++ ))
        done

        echo "Plotting ${outdir}..."
        for var in t u
        do
            commonFlags="--width=5 --aspect=1.2 --x=$var --y=z --ylim=0:17.5 --label=zR,zeta,flux_t"
            $tooldir/plot.py $commonFlags --save "$outdir/$var.pdf" $tmpdir/*.csv
        done
    done
done
