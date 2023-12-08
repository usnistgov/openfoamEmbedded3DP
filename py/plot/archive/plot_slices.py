#!/usr/bin/env python
'''Functions for plotting cross-sections'''

# external packages
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='10.0')


#-------------------------------------------

def plotSlices(data:pd.DataFrame) -> plt.Figure:
    '''use this to plot all slices within a pandas dataframe containing the points on the interface
    data is a pandas dataframe, e.g. one imported from importpoints'''
    d = data
    ticks = [-1, -0.5, 0, 0.5, 1]
    size = 8
    if len(xpts(data))>1:
        cmap = 'RdBu'
        color = 'x'
        sidelegend = True
    else:
        cmap = None
        color = 'Black'
        sidelegend = False
    p = d.plot.scatter(x='y', y='z', c=color, xlim=[-1, 1], ylim=[-1, 1], legend=False,\
               colormap = cmap, figsize=[size, size], \
               xticks=ticks, yticks=ticks, sharex=False, s=0.1)
    p.set_aspect('equal', adjustable='box')
    p.set_xlabel('y (mm)')
    p.set_ylabel('z (mm)')
    if sidelegend:
        f = plt.gcf()
        cax = f.get_axes()[1]
        cax.set_ylabel('x (mm)')
        cax.figsize=[0.5, size]
    return p

# def plotXSOnAx(pts:pd.DataFrame, ax, color) -> None:
#     '''sort the points by radial position and plot as line'''
#     ptlists = [[]]
    
#     # split the points into discrete objects
#     for i,row in pts.iterrows():
#         if i==0:
#             ptlists[0].append(dict(row))
#         else:
#             j = 0
#             while j<len(ptlists) and np.sqrt((row['y']-ptlists[j][-1]['y'])**2+(row['z']-ptlists[j][-1]['z'])**2)>0.1:
#                 # point too far away
#                 j+=1
#             if j<len(ptlists):
#                 ptlists[j].append(row)
#             else:
#                 ptlists.append([row])
                
#     # plot each object as a line
#     for ptsi in ptlists:
#         ptsidf = pd.DataFrame(ptsi)
#         ax.plot(ptsidf['y'], ptsidf['z'], color=color, linewidth=1, marker=None)
