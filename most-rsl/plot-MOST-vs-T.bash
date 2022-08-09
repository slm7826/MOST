#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"

pushd "$codeDir"
make
popd

for z0m in 1.0 0.1; do
    for wind in 1.0 5.0; do
        outdir=output/fluxes/z0m${z0m}w${wind}; tmpdir=$outdir/tmp

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
       use_RSL_lookup = F
/
 &most_nml
       z0m = $z0m
       k_over_B = 2.0
       u_atm = $wind
       t_atm = 300.0
       zR    = $zR
       var='t_sfc', x0 = 295.0, x1=305.0, nsamples = 200
/
EOF
            $codeDir/MOST.x > $tmpdir/`printf "%2.2d" $i`.csv
            (( i++ ))
        done

# var=flux_t
# ./plot.py --x=t_sfc --y=$var --label=u_atm --xlim=295:305 --ylim=-500:2500 --save "$outdir/$var.pdf" $tmpdir/*.csv
#
# var=b_star
# ./plot.py --x=t_sfc --y=$var --label=u_atm --xlim=295:305 --ylim=-0.05:0.5 --save "$outdir/$var.pdf" $tmpdir/*.csv
#
# var=flux_m
# ./plot.py --x=t_sfc --y=$var --label=u_atm --xlim=295:305 --ylim=-0.1:4.5 --save "$outdir/$var.pdf" $tmpdir/*.csv
#
# var=u_star
# ./plot.py --x=t_sfc --y=$var --label=u_atm --xlim=295:305 --ylim=-0.01:2.0 --save "$outdir/$var.pdf" $tmpdir/*.csv
#
# var=zeta
# ./plot.py --x=t_sfc --y=$var --label=u_atm --xlim=295:305 --save "$outdir/$var.pdf" $tmpdir/*.csv

# var=zeta
# ./plot.py --x=t_sfc --y=$var --label=zR --xlim=295:305 --ylim=-0.01:0.01 $tmpdir/*.csv  --save "$outdir/$var.pdf"

#t_sfc,u_atm,flux_t,flux_m,cd_m,cd_t,cd_q,ga,ra,u_star,b_star,rich,zeta,gust,gust+u_atm

        echo "Plotting ${outdir}..."
        for var in flux_t flux_m u_star zeta b_star cd_m cd_t gust gust+u_atm
        do
        #    $tooldir/plot.py --x=t_sfc --y=$var --title="$var: z0m=$z0m, kB=$k_over_B" --label=u_atm \
            $tooldir/plot.py --aspect=0.5 --x=t_sfc --y=$var --label=zR \
            --xlim=295:305 --save "$outdir/$var.pdf" $tmpdir/*.csv
        done
    done
done
