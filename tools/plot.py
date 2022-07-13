#!/usr/bin/env python
import csv
import re
import numpy      as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict

# filter out text before the separator
def decomment(csvfile):
    pastHeader = False
    for row in csvfile:
        # print(pastHeader,row[0:2])
        if pastHeader:
            yield row
        if re.match(r'\s*RESULTS:',row):
            pastHeader = True

# parse limit string
def parse_limits(str):
    #regular expresion for float: r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
    try :
       ys,ye = str.split(':')
       return (float(ys),float(ye))
    except :
       die(f'limit spec "{str}" is incorrect, must be "float:float"')

# parse arguments
parser = argparse.ArgumentParser(description='plot something')
parser.add_argument('-v','--verbose', dest='verb',
    help='increase verbosity', action='count', default=0)
parser.add_argument('-x','--x', action='store', default='t_sfc',
    help='plot x axis')
parser.add_argument('-y','--y', action='store', default='flux_t',
    help='variable to plot')
parser.add_argument('--label', action='store', default=None,
    help='list of settings to include in legend')
parser.add_argument('--xlimits', action='store', metavar='XS:XE',
    help='X limits for the plot')
parser.add_argument('--ylimits', action='store', metavar='YS:YE',
    help='Y limits for the plot')
parser.add_argument('--xlog', action='store_true', default=False,
    help='logarithmic X')
parser.add_argument('--ylog', action='store_true', default=False,
    help='logarithmic Y')
parser.add_argument('--width',type=float, default=8.0,
    help='width of the plot')
parser.add_argument('--aspect',type=float, default=1.0/1.618,
    help='aspect (height/width) of the plot')
parser.add_argument('--title', default=None,
    help='plot title')
parser.add_argument('-s','--save', nargs=1, metavar='FILENAME',
    help='save figure instead of plotting it on screen.')
parser.add_argument('input', nargs='+', metavar='INPUT',
    help='input files')
args=parser.parse_args()

if args.verb > 0:
    print('using pakages:')
    for package in np,mpl:
        print('    {:>20} : {}'.format(package.__name__,package.__version__))

figW=args.width; figH=figW*args.aspect
# fig,ax = plt.subplots(1,1,facecolor='w',figsize=(figW,figH))
fig = plt.figure(facecolor='w',figsize=(figW,figH))
ax  = fig.add_axes([0.1,0.13,0.85,0.85])
if (args.title): ax.set_title(args.title)
ax.set_xlabel(args.x)
ax.set_ylabel(args.y)
if args.xlimits:
    xs,xe=parse_limits(args.xlimits)
    ax.set_xlim(xs,xe)
if args.ylimits:
    ys,ye=parse_limits(args.ylimits)
    ax.set_ylim(ys,ye)
if args.xlog:
    ax.set_xscale('log')
if args.ylog:
    ax.set_yscale('log')
ax.grid(True)

for infile in args.input:
    # read CSV data
    columns = defaultdict(list) # each value in each column is appended to a list
    with open(infile) as f:
        reader = csv.DictReader(decomment(f)) # read rows into a dictionary format
        for row in reader: # read a row as {column1: value1, column2: value2,...}
            for (k,v) in row.items(): # go over each column name and value
                columns[k].append(float(v)) # append the value into the appropriate list
                                      # based on column name k
    # read namelist data for the legend
    if args.label:
        ll=args.label.split(',')
        label=[]
        for l in ll:
            with open(infile) as f:
                for row in f:
                    m = re.match(r'\s*'+l+r'\s*=\s*([^\s,]+)',row,flags=re.I)
                    if m:
                        try:
                            value = float(m.group(1))
                            label.append('{}={:g}'.format(l,value))
                        except:
                            label.append(m.group(1))
        label=', '.join(label)
    else:
        label=''
    ax.plot(columns[args.x],columns[args.y],label=label)

ax.legend(loc='best')

# print bounding box in some coordinates
# bb = ax.get_tightbbox(fig.canvas.get_renderer())
# print(bb.xmin,bb.xmax,bb.ymin,bb.ymax)

if args.save:
   fig.savefig(args.save[0], transparent=True, bbox_inches='tight')
#    fig.savefig(args.save[0], transparent=True)
else:
   plt.show()
