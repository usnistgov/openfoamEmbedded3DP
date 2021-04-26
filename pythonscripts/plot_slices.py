#!/usr/bin/env python
'''Functions for plotting cross-sections'''

# external packages
import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages
import interfacemetrics as intm
from plot_general import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# plotting
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------

#### plot slices


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



#### cross-section plots



def XSPlot(xs:pd.DataFrame, folder:str, cp:comboPlot) -> None:
    '''plot cross-sections for whole folder
    xs is a pandas DataFrame holding points
    folder is a full pathname
    cp is a comboPlot object'''
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
    

def XSPlotIdeal(cp:comboPlot, fs:intm.folderStats) -> None:
    '''this plots an ideal cross-section on the cross-section plot
    fs is a folderStats object'''
    x0 = cp.xrtot[0]+cp.dx/2
    y0 = cp.yrtot[-1]-cp.dy/2
    color='Black'
    plotCircle(cp.axs[0], x0, y0, fs.niw/2, 'Ideal', color)


def XSPlotf(folder:str, time:float, x:float, cp:comboPlot) -> None:
    '''plots cross-section from one file'''
    ptsx = intm.importPtsSlice(folder, time, x)
    if len(ptsx)>0:
        XSPlot(ptsx, folder, cp)


def XSPlots0(topFolder:str, exportFolder:str, time:float, x:float, sigmalist0:List[float], overwrite:bool=False, **kwargs) -> None:
    '''plot all cross-sections together
    topFolder is the folder that holds all the files
    time is the time at which to take the cross-section
    x is the distance behind the center of the nozzle to take the cross-section
    sigmalist0 is the list of surface tensions to include in this plot'''
    label = 'xs_'+str(x)+'_t_'+str(time)
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    dx = 0.7
    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 6.5, sigmalist=sigmalist0, **kwargs)
    cp.sigmalist = sigmalist0
    fs = intm.folderStats(cp.flist[0])
    (cp.flist).sort(key=lambda folder:extractTP(folder)['sigma']) # sort folders by sigma value so they are stacked in the right order
    

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
        
#     intm.exportIm(fn, cp.fig)