#!/usr/bin/env python
'''Functions for plotting cross-sections'''

# external packages
import sys
import os
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

def angle360(row:pd.Series, center:List) -> float:
    '''get the angle between the vector from center to row[y,z], relative to the x axis'''
    v = [row['y']-center[0],row['z']-center[1]]
    angle = np.arccos(np.dot(v, [1,0])/(np.linalg.norm(v)))
    if v[1]<0:
        angle = 2*np.pi-angle
    return angle

def plotXSOnAx(pts:pd.DataFrame, ax, color) -> None:
    '''sort the points by radial position and plot as line'''
    left = pts.y.min()
    right = pts.y.max()
    bottom = pts.z.min()
    top = pts.z.max()
    center = [(left+right)/2, (top+bottom)/2]
    pts['theta'] = [angle360(row, center) for i,row in pts.iterrows()]
    pts.sort_values(by='theta', inplace=True, ignore_index=True)
    ax.plot(pts['y'], pts['z'], color=color, linewidth=1, marker=None)

def XSPlot(xs:pd.DataFrame, folder:str, cp:comboPlot) -> None:
    '''plot cross-sections for whole folder
    xs is a pandas DataFrame holding points
    folder is a full pathname
    cp is a comboPlot object'''
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except:
        logging.info(folder)
        traceback.print_exc()
        return
    if xs.z.max()>cp.dy/2:
        # xs is going to fall outside of the plot bounds
        xind = cp.xmlist.index(x0)
        yind = cp.ymlist.index(y0)
        zout = (xs.z.max()-(cp.dy/2))/cp.dy + cp.dy/10
        cp.indicesreal = cp.indicesreal.append({'x':xind, 'y':yind+zout}, ignore_index=True)

    xs['y'] = xs['y']+x0  # shift the points relative to the point on the plot where this sim's origin should be
    xs['z'] = xs['z']+y0
    try:
        xmid,ymid = intm.pdCentroid(xs)
    except:
        return
    if cp.split:
        ax = cp.axs[sigmapos]
    else:
        ax = cp.axs[0]
    plotXSOnAx(xs, ax, color)
    ax.arrow(x0, y0, xmid-x0, ymid-y0,  head_width=0.05, head_length=0.1, fc=color, ec=color, length_includes_head=True)
    
    

def XSPlotIdeal(cp:comboPlot, le:dict) -> None:
    '''this plots an ideal cross-section on the cross-section plot
    le is a dictionary holding metadata, from legendUnique'''
    if len(cp.ylistreal)==1:
        # put to the left
        xind = int(cp.indicesreal.x.min())
        x0 = cp.xmlist[xind]-cp.dx
        xind = xind-1
        yind = int(cp.indicesreal.y.max())
        y0 = cp.ymlist[yind]
    else:
        # put above
        xind = int(cp.indicesreal.x.min())
        x0 = cp.xmlist[xind]
        yind = int(cp.indicesreal[cp.indicesreal.x==xind].y.max())
        y0 = cp.ymlist[yind]+cp.dy
        yind = yind+1
    color='Black'
    plotCircle(cp.axs[0], x0, y0, float(le['nozzle_inner_width'])/2, 'Ideal', color)
    cp.indicesreal = cp.indicesreal.append({'x':xind, 'y':yind}, ignore_index=True)


def XSPlotf(folder:str, time:float, xbehind:float, cp:comboPlot, xunits:str='mm', ref:dict={'nozzle_inner_width':0, 'ink_velocity':0, 'bath_velocity':0}) -> None:
    '''plots cross-section from one file
    time is the time of the cross-section in seconds
    xbehind is the distance behind the nozzle to take the cross-section, in xunits
    cp is the comboPlot object to plot on
    xunits is 'mm' or 'nozzle_inner_width'
    ref holds metadata about a reference simulation. leave values at 0 to not do any scaling. Otherwise, scale the size of the cross-section relative to the reference values
    '''
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
    xbehind is the distance behind the center of the nozzle to take the cross-section
    xunits is 'mm' or 'nozzle_inner_width'
    overwrite True to overwrite values
    dx is the spacing between cross-sections, in mm
    
    '''
    label = f'xs_{xbehind}{xunits}_t_{time}'
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return

    if 'dy' in kwargs:
        dy = kwargs['dy']
    else:
        dy = dx
    if 'imsize' in kwargs:
        imsize = kwargs['imsize']
        kwargs.pop('imsize')
    else:
        imsize = 6.5
    cp = comboPlot(topFolder, [-dx, dx], [-dy, dy], imsize, **kwargs)
    
    if len(cp.flist)==0:
        logging.warning("No files in list")
        return
    
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
