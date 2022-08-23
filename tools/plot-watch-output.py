#!/usr/bin/env python
from __future__ import print_function

import numpy      as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
#import matplotlib.ticker as ticker
import argparse
import datetime as dt
import re

parser = argparse.ArgumentParser(description='''
Extract relevant lines from watch point output and plot the extracted values.
If there are several matching values per time step (as may happen if the water
consistency steps or fog steps are triggered, then the last value is plotted)
''')
parser.add_argument('-v', '--verbose', dest='verb', action='count', default=0,
    help='''
    increase verbosity. Can be used multiple times. On first verbosity level, the names
    and patterns are printed as they are extracted; on the second level second
    the extracted extracted lines are printed, with corresponding dates.
    ''')
parser.add_argument('-e','--extract', nargs='+', action='append',
    help='''
    text to extract, with options for plotting. The plotting option can include
    any valid keyword arguments to matplotlib plot util, in "key=value" format, for example:
    "-e drag_q lw=0.5 ls=dashed marker=. color=k".
    In addition, ke/value pair var=pattern can be used to specify a value to plot from
    extracted line, if different from the patter used to search for matching lines, for
    example: ""
    ''')
parser.add_argument('--ylimits', action='store', metavar='YS:YE',
    help='Y limits for the plot')
parser.add_argument(
    '-s','--save', nargs=1, metavar='FILENAME',
    help='save figure instead of plotting it on screen.')
#
parser.add_argument('-i','--input', required=True,
    help='input text file with watch point output')
# example:
# plot-watch-output.py -v -v -e drag_q lw=0.5 ls=dashed marker=. color=k  -e atmos_wind  \
#     -i /lustre/f2/dev/Sergey.Malyshev/lm4p2/xanadu/lm4p2-GSWP3-potveg/ncrc3.intel18-prod/stdout/run/lm4p2-GSWP3-potveg.o205274178

args=parser.parse_args()

def parse_limits(str):
    #regular expresion for float: r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
    try :
       ys,ye = str.split(':')
       return (float(ys),float(ye))
    except :
       die(f'limit spec "{str}" is incorrect, must be "float:float"')

# parse y limits, if present
if args.ylimits: ys,ye=parse_limits(args.ylimits)

if args.verb > 0:
    print('using pakages:')
    for package in np,mpl:
        print('    {:>20} : {}'.format(package.__name__,package.__version__))

# pattern for model date string extraction
dateP = re.compile(r'.*update_land_model_fast_0d begins:(.*)')

# initialize matplotlib figure
figW=10.0; figH=figW/2
fig,ax = plt.subplots(1,1,facecolor='w',figsize=(figW,figH))

# loop over requested variables
for inp in args.extract:
    # parse list of individual plot parameters and files
    opts={} # empty dictionary for plot options

    # parse input arguments
    linePattern = valuePattern = inp[0] # patterns for line matching and value matching, initially they are identical
    # see if user wants to plot a different patterns for line match and plotted values
    for f in inp[1:]:
        m = re.match(r'(\w+)=(.*)',f)
        if m:
            if m.group(1) in ['plot','var']:
                # replace pattern for value extraction
                valuePattern = m.group(2)
            else:
                # add key/value pair to plot options
                opts[m.group(1)] = m.group(2)
    if args.verb>0:
        print('Extracting "{}" from lines containing "{}"'.format(valuePattern,linePattern))
    # compile string and value extraction pattern
    lineRe  = re.compile(r'\b'+re.escape(linePattern) +r'\b') # pattern for extracted string matching
    valueRe = re.compile(r'\b'+re.escape(valuePattern)+r'\b\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)') # pattern for value extraction

    # extract dates and relevant values
    dates=[];  lines=[]; values=[]; i = -1
    with open(args.input) as textFile:
        for line in textFile :
            m = dateP.match(line) # search for the beginning of the time step
            if m:
                # this if the starting line of the land model time time step
                dates.append(dt.datetime.strptime(m[1],r'%Y-%m-%d %H:%M:%S'))
                lines.append('')
                values.append(np.nan)
                i += 1
                continue # nothing to do with this line any more
            if i<0 :
                continue # we haven't fond the first date yet, skip the lines till then
            # now we got the line after the date, search it for the line match pattern
            if lineRe.search(line) :
                lines[i] = line.strip()
                # try to find
                m = valueRe.search(line)
                try:
                    values[i] = float(m[1])
                except:
                    print('CANNOT FIND VALUE PATTERN "{}" IN: "{} {}"'.
                           format(valuePattern,dates[i],lines[i]))

    if args.verb > 1:
        for l1,l2 in zip(dates,lines):
            print(l1,l2)

    # print(opts)
    ax.plot(dates,values,label=valuePattern,**opts)
#        ax.plot(values,label=valueS)

if args.ylimits :
    ax.set_ylim(ys,ye)
# ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(12))
locator   = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)

# ax.set_title('my title')
# ax.set_xlabel('x label')
ax.grid()
ax.legend(loc='best')

if args.save:
   fig.savefig(args.save[0], transparent=True, bbox_inches='tight')
else:
   plt.show()
