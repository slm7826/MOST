set -e # exit on errors
set -u # exit on undefined vars

codeDir="."
tooldir="../tools"
 outdir=output/test_lookup

pushd "$codeDir"
make
popd

mkdir -p $outdir


   cat <<EOF > input.nml
 &monin_obukhov_nml
       stable_option =  2,
       rich_crit = 1.0,
       zeta_trans =  0.5
       rsl_option = 'ghannam2022', rsl_mu_1 = 0.67, rsl_mu_m = 2.0, rsl_mu_t = 1.0
       a_min=  1.0e-5, a_max=  100.0,  a_nsteps=100,
       b_min= -1000.0, b_max=  1000.0, b_nsteps=300,
/
 &lookup_test_2D_nml
       a_min = 1e-5,  a_max = 100, a_nsteps = 400
       b_min = -1000, b_max = 1000, b_nsteps = 900
/
EOF

./test_lookup_2D.x
