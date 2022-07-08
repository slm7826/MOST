#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
pushd $codeDir
make
popd

# variables available in the output:
# t_sfc,u_atm,flux_t,flux_m,cd_m,cd_t,cd_q,ga,ra,u_star,b_star,rich,zeta,gust,gust+u_atm

# for z0m in 0.1 2; do
for z0m in 0.1; do
    for k_over_B in 2 0; do
        outdir=output/z0m${z0m}kb${k_over_B}
        tmpdir=$outdir/tmp

        mkdir -p $outdir $tmpdir
        rm -f $tmpdir/*.csv

        i=1
        for wind in 0.5 1.0 2.0 5.0 10.0
        do
            cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  'brutsaert'
       rich_crit = 1.0
       zeta_trans =  0.5
/
 &most_nml
       z0m = $z0m
       k_over_B = $k_over_B
       u_atm = $wind
       t_atm = 300.0
       var='t_sfc', x0 = 295.0, x1=305.0, nsamples = 200
/
EOF
            $codeDir/MOST.x > $outdir/tmp/`printf "%2.2d" $i`.csv
            (( i++ ))
        done

        commonFlags="--x=t_sfc --label=u_atm --xlim=295:305 --aspect=0.5"
        var=flux_t
        $tooldir/plot.py $commonFlags --y=$var --ylim=-300:450 --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=b_star
        $tooldir/plot.py $commonFlags --y=$var --ylim=-0.02:0.16 --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=flux_m
        $tooldir/plot.py $commonFlags --y=$var --ylim=-0.05:0.85  --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=u_star
        $tooldir/plot.py $commonFlags --y=$var --ylim=-0.05:0.9 --save "$outdir/$var.pdf" $outdir/tmp/*.csv
    done
done
