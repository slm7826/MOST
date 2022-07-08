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

for z0m in 2; do
    for k_over_B in 2 0; do
        outdir=output/z0m${z0m}kb${k_over_B}-vsU
        tmpdir=$outdir/tmp

        mkdir -p $outdir $tmpdir
        rm -f $tmpdir/*.csv

        i=1
        for t_sfc in 300.5 301.0 302.5
        do
            cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option = 2
       rich_crit     = 1.0
       zeta_trans    = 0.5
/
 &most_nml
       z0m = $z0m
       k_over_B = $k_over_B

       t_atm = 300.0
       t_sfc = $t_sfc
       var='u_atm', x0 = 0.01, x1=5.0, nsamples = 200
/
EOF
            $codeDir/MOST.x > $outdir/tmp/`printf "%2.2d" $i`.csv
            (( i++ ))
        done

        commonFlags="--x=u_atm --label=t_sfc --xlim=0:5 --aspect=0.5"
        var=flux_t
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:1000 --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=flux_m
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:1.3  --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=b_star
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:0.5  --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=u_star
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:1.2  --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=gust
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:5    --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        var=gust+u_atm
        $tooldir/plot.py $commonFlags --y=$var --ylim=0:7    --save "$outdir/$var.pdf" $outdir/tmp/*.csv

        for var in rich zeta
        do
            $tooldir/plot.py $commonFlags --y=$var --save "$outdir/$var.pdf" $outdir/tmp/*.csv
        done
    done
done