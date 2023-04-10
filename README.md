# MOST: exploration of Monin-Obukhov parameterization

## Setting up environment

To set up environment on GFDL workstation, do:
```bash
source env-workstation
```

On a laptop with gfortran installed, I do:
```bash
source env-mac
```
## Files
- `most-FMS-2022.02/`: code from 2022.02 release of FMS
- `most-brutsaert/`: code from rev.f5ebc26b (2022-07-01 on branch `user/slm/sfcx`), with implementation of Brutsaert stability functions
- `most-rsl/` : code with roughness sublayer parameterization
- `tools/` : plotting tools
- `MOST.tex` : Notes and description of the parameterizations

## How to build figures

```bash
cd most-brutsaert
./plot-phi.bash
./plot-vs-T1.bash
./plot-vs-T2.bash
./plot-vs-U.bash

cd ../most-FMS-2022.02
./plot-vs-T1.bash
./plot-vs-T2.bash
./plot-vs-U.bash

cd ../most-rsl
./plot-STABLE_MIX.bash
./plot-F-vs-L.bash
./plot-F-vs-zR.bash
./plot-I-vs-A.bash
./plot-I-vs-B.bash
./plot-I-vs-A-Bru.bash
./plot-I-vs-B-Bru.bash
./plot-MOST-vs-T.bash
./plot-PROFILE.bash

# this takes a long time
./prepare-lookup2D.bash
../tools/interp2D.py lookup_test.nc -s interp2D.png
../tools/interp1D.py lookup_test.nc -s interp1d.pdf
```
