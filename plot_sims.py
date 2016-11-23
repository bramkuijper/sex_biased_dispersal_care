#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib
from matplotlib.ticker import AutoMinorLocator

import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^patch.*",line) != None:
        parline = idx - 1;
        break;

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    histdat = pd.read_csv(filename, sep=";")

# generate the figure
minorLocator = AutoMinorLocator(5)

# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 2

# add first subplot
ax1 = plt.subplot(num_rows,1,1)
ax1.yaxis.set_minor_locator(minorLocator)
ax1.plot(
        histdat["generation"],histdat["uf"],'b',
        histdat["generation"],histdat["um"],'r',
        )
ax1.set_ylim([0,.5])
ax1.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
ax1.set_ylabel(r'Parental effort, $u_{x}$')
ax1.legend((r'$u_{\mathrm{f}}$',r'$u_{\mathrm{m}}$'))

# add subplot for variances
ax2 = plt.subplot(num_rows,1,2)
ax2.plot(histdat["generation"],histdat["varuf"],'b',
            histdat["generation"],histdat["varum"],'r')
ax2.set_ylabel(r'Variances')
ax2.legend((r'$\sigma_{u_{\mathrm{f}}}^{2}$',r'$\sigma_{u_{\mathrm{m}}}^{2}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
