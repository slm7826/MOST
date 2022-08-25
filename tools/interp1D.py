#!/usr/bin/env python
from __future__ import print_function

import numpy      as np
import netCDF4    as nc
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='plot something')
parser.add_argument('-v','--verbose', dest='verb', default=0,
    help='increase verbosity', action='count')
parser.add_argument('-s','--save', nargs=1, metavar='FILENAME',
    help='save figure instead of plotting it on screen.')
parser.add_argument('input', metavar='FILENAME',
    help='input file')
args=parser.parse_args()

if args.verb > 0:
    print('using pakages:')
    for package in np,nc,mpl:
        print('    {:>20} : {}'.format(package.__name__,package.__version__))

src=nc.Dataset(args.input,'r')
a = src.variables['a'][...]
Im0 = src.variables['Im0'][...]
Im1 = src.variables['Im1'][...]
It0 = src.variables['It0'][...]
It1 = src.variables['It1'][...]

fig, axs = plt.subplots(1,2,facecolor='w',figsize=(12,5))

axs[0].plot(a,Im0[-1,:],label='$I_m$, integrated')
axs[0].plot(a,Im1[-1,:],label='$I_m$, interpolated',ls='--')
axs[0].plot(a,It0[-1,:],label='$I_t$, integrated')
axs[0].plot(a,It1[-1,:],label='$I_m$, interpolated',ls='--')

axs[1].plot(a,Im1[-1,:]-Im0[-1,:],label=r'Interpolated - Integrated, $I_m$')
axs[1].plot(a,It1[-1,:]-It0[-1,:],label=r'Interpolated - Integrated, $I_t$')
for ax in [axs[0], axs[1]]:
    ax.set_xlim(0,5.0)
    ax.set_xlabel(r'$a=z_1/z_R$')
    ax.grid()
    ax.legend(loc='best')

#     fig.colorbar(m1,ax=ax)
# axs[1].set_xlabel(r'$b=z_R/L$')
# axs[1].sharex(axs[0])
# axs[0].set_title(r'Interpolated - integrated momentum integral $I_m$',y=1.0,pad=-14)
# axs[1].set_title(r'Interpolated - integrated heat integral $I_t$',y=1.0,pad=-14)

if args.save:
   fig.savefig(args.save[0], transparent=True, bbox_inches='tight')
else:
   plt.show()

