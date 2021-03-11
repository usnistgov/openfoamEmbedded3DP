import os
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
import statistics as st
from PIL import Image
import interfacemetrics as intm
import scipy as sp
from shapely.geometry import Polygon
import re
import folderparser as fp
import random
import math
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10
from typing import List, Dict, Tuple, Union, Any, TextIO
from interfacemetricsplots import *

#### plot slices

# use this to plot all slices within a pandas dataframe containing the points on the interface
    # data is a pandas dataframe, e.g. one imported from importpoints
def plotSlices(data:pd.DataFrame) -> plt.Figure:
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



#### cross-section plots

# plot cross-sections for whole folder
    # xs is a pandas DataFrame holding points
    # folder is a full pathname
    # cp is a comboPlot object
def XSPlot(xs:pd.DataFrame, folder:str, cp:comboPlot) -> None:
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except:
        return
    xlist = list(xs['y']+x0)
    ylist = list(xs['z']+y0)
    try:
        x,y = intm.pdCentroid(xs)
    except:
        return
    else:
        xmid = x + x0
        ymid = y + y0
    if cp.split:
        ax = cp.axs[sigmapos]
    else:
        ax = cp.axs[0]
    size = 0.01*6/(len(cp.xmlist)*(cp.ncol))
    ax.scatter(xlist, ylist, c=color, s=0.01, rasterized=True)
    ax.arrow(x0, y0, xmid-x0, ymid-y0, ls='--', head_width=0.05, head_length=0.1, fc=color, ec=color, length_includes_head=True)
    
# this plots an ideal cross-section on the cross-section plot
# fs is a folderStats object
def XSPlotIdeal(cp:comboPlot, fs) -> None:
    x0 = cp.xrtot[0]+cp.dx/2
    y0 = cp.yrtot[-1]-cp.dy/2
    color='Black'
    plotCircle(cp.axs[0], x0, y0, fs.niw/2, 'Ideal', color)

# plots cross-section from one file
def XSPlotf(folder:str, time:float, x:float, cp:comboPlot) -> None:
    ptsx = intm.importPtsSlice(folder, time, x)
    if len(ptsx)>0:
        XSPlot(ptsx, folder, cp)

# plot all cross-sections together
    # topFolder is the folder that holds all the files
    # time is the time at which to take the cross-section
    # x is the distance behind the center of the nozzle to take the cross-section
    # sigmalist0 is the list of surface tensions to include in this plot
def XSPlots0(topFolder:str, exportFolder:str, time:float, x:float, sigmalist0:List[float], overwrite:bool=False, **kwargs) -> None:
    label = 'xs_'+str(x)+'_t_'+str(time)
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    dx = 0.7
    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 6.5, sigmalist=sigmalist0, **kwargs)
    cp.sigmalist = sigmalist0
    fs = intm.folderStats(cp.flist[0])

    xvalue = fs.ncx + x
    for folder in cp.flist:
        XSPlotf(folder, time, xvalue, cp)
    cp.figtitle = 'Cross sections, '+str(x)+' mm behind nozzle, t = '+str(time)+' s'
    cp.clean()
    for ax in cp.axs:
        ax.set_ylim([cp.yrtot[0], cp.yrtot[1]+1])
    XSPlotIdeal(cp,fs)
    for ax in cp.axs:
        ax.grid(linestyle='-', linewidth='0.25', color='#949494')
        
    intm.exportIm(fn, cp.fig)