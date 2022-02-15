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
    

def XSPlotIdeal(cp:comboPlot, le:dict) -> None:
    '''this plots an ideal cross-section on the cross-section plot
    fs is a dictionary holding metadata, from legendUnique'''
    xind = int(cp.indicesreal.x.min())
    x0 = cp.xmlist[xind]
    yind = int(cp.indicesreal[cp.indicesreal.x==xind].y.max())
    y0 = cp.ymlist[yind]+cp.dy
    color='Black'
    plotCircle(cp.axs[0], x0, y0, float(le['nozzle_inner_width'])/2, 'Ideal', color)
    cp.indicesreal = cp.indicesreal.append({'x':xind, 'y':yind+1}, ignore_index=True)


def XSPlotf(folder:str, time:float, xbehind:float, cp:comboPlot, xunits:str='mm', ref:dict={'nozzle_inner_width':0, 'ink_velocity':0, 'bath_velocity':0}) -> None:
    '''plots cross-section from one file'''
    ptsx = intm.importPtsSlice(folder, time, xbehind, xunits=xunits)
    if len(ptsx)==0:
        return
    if float(ref['nozzle_inner_width'])>0:
        # rescale the points so all cross-sections reference the same 
        le = fp.legendUnique(folder)
        dEst = float(le['nozzle_inner_width'])*np.sqrt(float(le['ink_velocity'])/float(le['bath_velocity'])) # ideal final diameter
        dEst0 = float(ref['nozzle_inner_width'])*np.sqrt(float(ref['ink_velocity'])/float(ref['bath_velocity'])) # ideal final diameter
        for s in ['y', 'z']:
            ptsx[s] = [dEst0/dEst*i for i in ptsx[s]]
    if len(ptsx)>0:
        XSPlot(ptsx, folder, cp)


def XSPlots0(topFolder:str, exportFolder:str, time:float, xbehind:float, xunits:str='mm', overwrite:bool=False, dx:float=0.5, **kwargs) -> None:
    '''plot all cross-sections together
    topFolder is the folder that holds all the files
    time is the time at which to take the cross-section
    x is the distance behind the center of the nozzle to take the cross-section
    sigmalist0 is the list of surface tensions to include in this plot'''
    label = f'xs_{xbehind}{xunits}_t_{time}'
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return

    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 6.5, **kwargs)
    
    (cp.flist).sort(key=lambda folder:extractTP(folder)['sigma']) # sort folders by sigma value so they are stacked in the right order
    le0 = fp.legendUnique(cp.flist[0]) # use first file as reference dimensions
    for folder in cp.flist:
        XSPlotf(folder, time, xbehind, cp, xunits=xunits, ref=le0)
    xunname = xunits.replace('nozzle_inner_width', '$d_i$')
    cp.figtitle = f'Cross sections, {xbehind} {xunname} behind nozzle, t = {time} s'
    XSPlotIdeal(cp,le0)
    cp.clean()

    for ax in cp.axs:
        ax.grid(linestyle='-', linewidth='0.25', color='#949494')
        
    intm.exportIm(fn, cp.fig, **kwargs)
