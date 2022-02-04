#!/usr/bin/env python
'''Functions for plotting cross-sections'''

# external packages
import sys
import os
# currentdir = os.path.dirname(os.path.realpath(__file__))
# parentdir = os.path.dirname(currentdir)
# sys.path.append(parentdir)
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
import interfacemetrics as intm
from plot_general import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='10.0')

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
        traceback.print_exc()
        return
#     xlist = list(xs['y']+x0)
#     ylist = list(xs['z']+y0)
    xs['y'] = xs['y']+x0
    xs['z'] = xs['z']+y0
    try:
        xmid,ymid = intm.pdCentroid(xs)
    except:
        return
#     else:
#         xmid = x + x0
#         ymid = y + y0
    if cp.split:
        ax = cp.axs[sigmapos]
    else:
        ax = cp.axs[0]
#     size = 0.01*6/(len(cp.xmlist)*(cp.ncol))
    
    bottom = xs.z.min()
    height= xs.z.max() - bottom
    left = xs.y.min()
    width = xs.y.max() - left
    
    leftpts = xs.copy()
    leftpts = leftpts[(leftpts.y<left+0.25*width)]
    leftpts.sort_values(by='z', inplace=True)
    toppts = xs.copy()
    toppts = toppts[(toppts.y>left+0.25*width)&(toppts.y<left+0.75*width)&(toppts.z>leftpts.z.median())]
    toppts.sort_values(by='y', inplace=True)
    botpts = xs.copy()
    botpts = botpts[(botpts.y>left+0.25*width)&(botpts.y<left+0.75*width)&(botpts.z<leftpts.z.median())]
    botpts.sort_values(by='y', inplace=True)
    rightpts = xs.copy()
    rightpts = rightpts[(rightpts.y>left+0.55*width)]
    rightpts.sort_values(by='z', inplace=True)
    
    for pts in [leftpts, toppts, botpts, rightpts]:
        ax.plot(pts['y'], pts['z'], color=color, linewidth=1, marker=None)
    
    
#     ax.scatter(xlist, ylist, color=color, s=0.01, marker='.', facecolor=color, edgecolor=None, rasterized=True)
    ax.arrow(x0, y0, xmid-x0, ymid-y0,  head_width=0.05, head_length=0.1, fc=color, ec=color, length_includes_head=True)
    

def XSPlotIdeal(cp:comboPlot, fs:intm.folderStats) -> None:
    '''this plots an ideal cross-section on the cross-section plot
    fs is a folderStats object'''
    x0 = cp.xrtot[0]+cp.dx/2
    y0 = cp.yrtot[-1]+cp.dy/2
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
    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 12, sigmalist=sigmalist0, **kwargs)
    cp.sigmalist = sigmalist0
    
    fs = intm.folderStats(cp.flist[0])
    (cp.flist).sort(key=lambda folder:extractTP(folder)['sigma']) # sort folders by sigma value so they are stacked in the right order
    xvalue = fs.ncx + x
    for folder in cp.flist:
        XSPlotf(folder, time, xvalue, cp)
    cp.figtitle = f'Cross sections, {x} mm behind nozzle, t = {time} s'
    cp.clean()
    for ax in cp.axs:
        if len(cp.yrtot)>1: # RG
            ax.set_ylim([cp.yrtot[0], cp.yrtot[1]+1])
    XSPlotIdeal(cp,fs)
    for ax in cp.axs:
        ax.grid(linestyle='-', linewidth='0.25', color='#949494')
        
    intm.exportIm(fn, cp.fig, **kwargs)
