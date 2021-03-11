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


#################### STEADY STATE PLOTS  


# stsp gets the steady times and steady positions for a folder
    # folder is a full path name
    # outputs st and sp are pandas dataframes
def stsp(folder:str):
    stfn = os.path.join(folder, 'steadytimes.csv')
    spfn = os.path.join(folder, 'steadypositions.csv')
    if not os.path.exists(stfn):
        raise stfn+' does not exist'
    if not os.path.exists(spfn):
        raise spfn+' does not exist'
    st = intm.plainIm(stfn, 0)
    sp = intm.plainIm(spfn, 0)
    return st,sp


# plot steady state conditions in time and space. This produces two regions, one where the filament is steady in space, and another where the filament is steady in time
    # st is a dataframe where the filament is steady in time, as produced by steadytime
    # sp is a dataframe where the filament is steady in position, as produced by steadypos
def plotSteadyMetrics(st:pd.DataFrame, sp:pd.DataFrame):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    if len(sp)>0:
        ax.plot(sp['t'], sp['xf'], color='royalblue', label='steady in position')
        ax.plot(sp['t'], sp['x0'], color='royalblue')
        ax.fill_between(sp['t'], sp['x0'], sp['xf'], interpolate=True, alpha=0.2, facecolor='royalblue')
    if len(st)>0:
        ax.plot(st['t0'], st['x'], color='maroon', label='steady in time')
        ax.plot(st['tf'], st['x'], color='maroon')
        ax.fill_betweenx(st['x'], st['t0'], st['tf'], interpolate=True, alpha=0.2, facecolor='maroon')
    plt.xlim(0,2.5)
    plt.ylim(0,8)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Distance behind nozzle (mm)')
    if len(sp)>0 or len(st)>0:
        ax.legend()
    return ax

        
# given a folder name, plot the steady regions
    # folder is a full path name
def plotSteadyFromFolder(folder:str):
    st,sp=stsp(folder)
    return plotSteadyMetrics(st,sp)

# plot steady state conditions in time and space. This produces two regions, one where the filament is steady in space, and another where the filament is steady in time
    # ax is the axis to plot the results on
    # st is a dataframe where the filament is steady in time, as produced by steadytime
    # sp is a dataframe where the filament is steady in position, as produced by steadypos
    # c is a color
def plotSteadyColor(ax, st:pd.DataFrame, sp:pd.DataFrame, c, labels:bool=False):
    # plot steady in position
    c = '#940f0f'
    if len(sp)>0:
        ax.plot(sp['t'], sp['xf'], color=c, linewidth=0.5)
        ax.plot(sp['t'], sp['x0'], color=c, linewidth=0.5)
        
        if labels:
            ax.fill_between(sp['t'], sp['x0'], sp['xf'], interpolate=True, alpha=0.5, facecolor=c, edgecolor=c, hatch="\\\\\\", label='steady in position', linewidth=0.5)
        else:
            ax.fill_between(sp['t'], sp['x0'], sp['xf'], interpolate=True, alpha=0.5, facecolor=c, edgecolor=c, hatch="\\\\\\", linewidth=0.5)
    
    # plot steady in time
    c = '#356577'
    if len(st)>0:
        ax.plot(st['t0'], st['x'], color=c, linewidth=0.5)
        ax.plot(st['tf'], st['x'], color=c, linewidth=0.5)
        if labels:
            ax.fill_betweenx(st['x'], st['t0'], st['tf'], interpolate=True, alpha=0.2, facecolor=c, label='steady in time', linewidth=0.5)
        else:
            ax.fill_betweenx(st['x'], st['t0'], st['tf'], interpolate=True, alpha=0.2, facecolor=c, linewidth=0.5)
    return ax

    
def plotSteadyAx(ax, folder:str, color, xlim:float=5.1, ylim:float=7.2, labels:bool=False) -> None:
    try:
        st,sp = stsp(folder)
        
    except:
        return
    if len(st)<2 or len(sp)<2:
        return
    l = intm.importLegend(folder)
    finaltime = float(l.loc[6, 'val'])
    st.replace(1000, finaltime)
    sp.replace(1000, ylim)
    color='gray'
    ax.plot([finaltime, finaltime],[-100,100], color=color, linestyle='--', linewidth='1')
    st['tf'] = [finaltime for i in range(len(st))]
    plotSteadyColor(ax, st, sp, color, labels=labels)  
    ax.set_xlim(0,xlim)
    ax.set_ylim(0,ylim)
    

# plot steady regions for one folder in a grid of plots
    # gp is a gridOfPlots object
    # folder is a path name
def plotSteadyInGrid(gp:gridOfPlots, folder:str) -> None:
    try:
        color, x0, y0, sigmapos = vvplot(folder, gp)
        # in the grid of plots, the row# is y, and the col# is x
    except:
        print('vvplot error')
        return
    ax = gp.axs[len(gp.axs)-y0-1, x0]
    plotSteadyAx(ax, folder, color)

# steadyPlots plots all of the folders on one grid fo plots
    # topFolder is the full path name to the folder that holds all of the folders
    # imsize is the size of each plot
def steadyPlots(topFolder:str, imsize:int, exportFolder:str, sigmalist:List[float], overwrite:bool=False, **kwargs) -> None:
    labeli = 'steady'
    for s in sigmalist:
        labeli = labeli+'_'+str(s)
    fn = intm.imFn(exportFolder, labeli, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    gp = gridOfPlots(topFolder, imsize, sigmalist=sigmalist, **kwargs)
    for folder in gp.flist:
        tp = extractTP(folder)
        sig = tp['sigma']
        if sig in sigmalist:
            plotSteadyInGrid(gp, folder)
    
    gp.clean()
    
    # in the grid of plots, the row# is y, and the col# is x
    for i in range(len(gp.axs)):
        # going down the rows
        gp.axs[i][0].set_ylabel('Distance (mm)')
    for j in range(len(gp.axs[0])):
        # going across the columns
        gp.axs[-1][j].set_xlabel('Time (s)')
 
    intm.exportIm(fn, gp.fig)
    
    

###### single file stability plots

# cmap is a colormap (matplotlib.colors.ListedColormap)
def stab1(folder:str, t:float, x:float, i:int, num:int, ax:plt.Axes, dx:float, cmap) -> None:
    pts = intm.importPtsSlice(folder, t, x)
    if len(pts)>0:
        xlist = list(pts['y']+dx*i)
        ylist = list(pts['z'])
        color = cmap(i/(num-1))
        ax.scatter(xlist, ylist, color=color, s=0.01)

# go back to the stability plot and fix it after the fact
def adjustStabilityPlot(ax:plt.Axes, numx:float, xlabel:str, xlist:List[float], dx:float) -> None:
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('z (mm)')
    ax.set_xticks(dx*np.linspace(0, numx-1, num=numx), minor=False)
    ax.set_xticklabels(xlist, fontname="Arial", fontsize=10) 
    ax.axhline(y=0, color='k')

# a stability plot plots cross-sections across time and space
# it will have two rows: one for a constant time but shifting position
# another for a constant position but shifting time
# this is all from the same simulation
def stabilityPlot(folder:str, exportFolder:str, tconst:float, xconst:float, export:bool=True) -> None:
    fs = intm.folderStats(folder)
    fig, axs = plt.subplots(1,1, figsize=(6.5, 2.4))
    axs0 = plt.subplot2grid((2, 3), (0, 0), colspan=2, frameon=False)
    axs1 = plt.subplot2grid((2, 3), (1, 0), colspan=2, frameon=False)
    axs2 = plt.subplot2grid((2, 3), (0, 2), rowspan=2)
        
    # color map
    cmap = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)

    # plot xs as a function of x
    dx = fs.niw
    nxs = 10 # number of cross-sections per plot\
    xmin = 0.5
    xmax = 6.5
    numx = round((xmax-xmin)/0.5)+1
    xlist = np.linspace(xmin, xmax, num=numx)
    for i,x in enumerate(xlist): # nozzle right edge to bath right x
        stab1(folder, tconst, x+fs.ncx, i, numx, axs0, dx, cmap)
        
    # plot xs as a function of t
    tmin = 0
    tmax = 5
    numt = round((tmax-tmin)/0.25)+1
    tlist = np.linspace(tmin, tmax, num=numt)
    for i,t in enumerate(tlist): # nozzle right edge to bath right x
        stab1(folder, t, xconst+fs.ncx, i, numt, axs1, dx, cmap)
        
    # format the plots
    adjustStabilityPlot(axs0, numx, 'slice position (mm)', xlist, dx)
    adjustStabilityPlot(axs1, numt, 'slice time (s)', tlist, dx)
    
    # add labels
    tx = -0.1
    ty = -0.5
    axs0.text(tx, ty, 'time = '+str(tconst)+' s', transform=axs0.transAxes, horizontalalignment='left', verticalalignment='top')
    axs1.text(tx, ty, 'position = '+str(xconst)+' mm', transform=axs1.transAxes, horizontalalignment='left', verticalalignment='top')
    
    # plot steady metrics
    plotSteadyAx(axs2, folder, '#000000', labels=True)
    axs2.set_xlabel('Time (s)')
    axs2.set_ylabel('Position (mm)')
    axs2.legend()
    
    fig.tight_layout()
    fbase = os.path.basename(folder)
    fname = 'stability_'+fbase+'_t_'+str(tconst)+'_x_'+str(xconst)
    if export:
        intm.exportIm(os.path.join(exportFolder, fname), fig)
        