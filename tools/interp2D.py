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

dIm = src.variables['Im1'][...]-src.variables['Im0'][...]
dIt = src.variables['It1'][...]-src.variables['It0'][...]

fig, axs = plt.subplots(2,1,facecolor='w',figsize=(10,8))
Y, X = np.meshgrid(src.variables['a'][...], src.variables['b'][...])

# plt.rcParams['axes.titley'] = 1.0

m1 = axs[0].pcolormesh(X, Y, dIm, cmap='seismic',vmin=-0.5,vmax=0.5)
m2 = axs[1].pcolormesh(X, Y, dIt, cmap='seismic',vmin=-0.5,vmax=0.5)
for ax in [axs[0], axs[1]]:
    ax.set_ylim(0.0,10.0)
    ax.set_ylabel(r'$a=z_1/z_R$')
    fig.colorbar(m1,ax=ax)
axs[1].set_xlabel(r'$b=z_R/L$')
axs[1].sharex(axs[0])
axs[0].set_title(r'Interpolated - integrated momentum integral $I_m$',y=1.0,pad=-14)
axs[1].set_title(r'Interpolated - integrated heat integral $I_t$',y=1.0,pad=-14)
if args.save:
   fig.savefig(args.save[0], transparent=True, bbox_inches='tight')
else:
   plt.show()

