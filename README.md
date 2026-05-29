# MOST: exploration of Monin-Obukhov parameterization

## Setting up environment

To set up environment on GFDL workstation, do:
```bash
source env-workstation
```

On a laptop with gfortran installed, I do:
```bash
source env-workstation
```
## Files
`most-FMS-2022.02/`: this directory contains the code that uses the implementation of MOST
from the 2022.02 release of FMS.

- `SIM.F90`: toy model simulating fluxes and temperatures given
idealized input of short-wave and long-wave radiation, prescribed wind, and
prescribed temperature in the lowest layer of the atmosphere. The ground is a
simple layer with constant heat capacity (which can be zero), directly connected
to the surface fluxes.

- `plot-SIM.bash`: script that compiles the code, runs the
SIM.F90 simulation with the parameter it sets, and plots output sving the
figures into output directory. Note that it overrides previous results in the
output dir, if any.

- `SIMA.F90`: same as `SIM.F90`, except temperature of the
atmosphere is also prognostic.

- `plot-SIMA.bash`: as `plot-SIM.bash`, except using the
executable with prognostic temperature of the atmospheric later `SIMA.F90`.

`most-brutsaert/`: this directory contains the tests that use MOST
implementation from rev.f5ebc26b (2022-07-01 on branch `user/slm/sfcx`), with
modified (0d) interfaces and ability to use Brutsaert stability functions