#!/usr/bin/env python
'''Functions for plotting steady state metrics'''

# external packages
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
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


#################### STEADY STATE PLOTS  



def stsp(folder:str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''stsp gets the steady times and steady positions for a folder
    folder is a full path name
    outputs st and sp are pandas dataframes'''
    stfn = os.path.join(folder, 'steadytimes.csv')
    spfn = os.path.join(folder, 'steadypositions.csv')
    if not os.path.exists(stfn):
        raise stfn+' does not exist'
    if not os.path.exists(spfn):
        raise spfn+' does not exist'
    st, stunits = intm.plainIm(stfn, 0)
    sp, spunits = intm.plainIm(spfn, 0)
    return st,sp

#----------------------
# time and position different colors

def plotSteadyMetrics(st:pd.DataFrame, sp:pd.DataFrame) -> plt.Axes:
    '''plot steady state conditions in time and space. This produces two regions, one where the filament is steady in space, and another where the filament is steady in time. Time and position are different colors.
    st is a dataframe where the filament is steady in time, as produced by steadytime
    sp is a dataframe where the filament is steady in position, as produced by steadypos'''
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

        

def plotSteadyFromFolder(folder:str) -> plt.Axes:
    '''given a folder name, plot the steady regions
    folder is a full path name'''
    st,sp=stsp(folder)
    return plotSteadyMetrics(st,sp)


#-----------------------
# time and position same color

def plotSteadyColor(ax, st:pd.DataFrame, sp:pd.DataFrame, c, labels:bool=False) -> plt.Axes:
    '''plot steady state conditions in time and space. This produces two regions, one where the filament is steady in space, and another where the filament is steady in time, where both time and position are the same color but have different hatching
    ax is the axis to plot the results on
    st is a dataframe where the filament is steady in time, as produced by steadytime
    sp is a dataframe where the filament is steady in position, as produced by steadypos
    c is a color'''
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

    
def plotSteadyAx(ax:plt.Axes, folder:str, color, xlim:float=5.1, ylim:float=7.2, labels:bool=False) -> None:
    '''Plot steady regions on a plot for a single folder, where both time and position are the same color but have different hatching
    ax is the axis to plot on
    folder is the simulation folder path
    color is the color to plot both time and position
    xlim is the x limit of the axis
    ylim is the y limit of the axis
    labels False to not include axis labels'''
    try:
        st,sp = stsp(folder)
    except:
        return
    if len(st)<2 or len(sp)<2:
        return

    finaltime = fp.currentTime(folder)[0]
    
    for i,row in st.iterrows():
        if float(row['tf'])>finaltime:
            st.loc[i,'tf'] = finaltime
    ylims = str(ylim) # RG
    for i,row in sp.iterrows():
        if float(row['xf'])>ylim:
            sp.loc[i,'xf'] = ylim 
    color='gray'
    ax.plot([finaltime, finaltime],[-100,100], color=color, linestyle='--', linewidth='1')
    plotSteadyColor(ax, st, sp, color, labels=labels)
    ax.set_ylim(0,ylim)
    ax.set_xticks([0,2.5])
    

def plotSteadyInGrid(gp:gridOfPlots, folder:str) -> None:
    '''plot steady regions for one folder in a grid of plots, where color depends on surface tension
    gp is a gridOfPlots object
    folder is a path name'''
    try:
        color, x0, y0, sigmapos = vvplot(folder, gp)
        # in the grid of plots, the row# is y, and the col# is x
    except:
        print('vvplot error')
        return
    dim = np.ndim(gp.axs)
    ax = gp.axs[y0, x0]
    plotSteadyAx(ax, folder, color)


def steadyPlots(topFolder:str, imsize:int, exportFolder:str, sigmalist:List[float], overwrite:bool=False, **kwargs) -> None:
    '''steadyPlots plots all of the folders on one grid of plots, where color depends on surface tension
    topFolder is the full path name to the folder that holds all of the folders
    imsize is the size of each plot
    exportFolder is the folder to export the plot to
    sigmalist is a list of sigma values to include. This is most legible with only one sigma value.'''
    
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
    if np.ndim(gp.axs)!=1: # RG
        for i in range(len(gp.axs)):
            # going down the rows
            gp.axs[i][0].set_ylabel('Distance (mm)')

        for j in range(len(gp.axs[0])):
            # going across the columns
            gp.axs[-1][j].set_xlabel('Time (s)')
    else:
        gp.axs[0].set_ylabel('Distance (mm)')
        mid = len(gp.axs)//2 # put y label at the middle plot
        gp.axs[mid].set_xlabel('Time (s)')

#     intm.exportIm(fn, gp.fig)
    
    
#-------------------------------------------------------------------
###### single file stability plots


def stab1(folder:str, t:float, x:float, i:int, num:int, ax:plt.Axes, dx:float, cmap) -> None:
    '''Plot one slice
    folder is the path name for the simulation
    t is the time
    x is the position in the bath (absolute, not relative to nozzle)
    i is the index of the slice, used to determine color
    num is the total number of slices we will plot
    ax is the axis to plot on
    dx is the spacing between slices on the plot
    cmap is a colormap (matplotlib.colors.ListedColormap)'''
    pts = intm.importPtsSlice(folder, t, x)
    if len(pts)>0:
        xlist = list(pts['y']+dx*i)
        ylist = list(pts['z'])
        color = cmap(i/(num-1))
        ax.scatter(xlist, ylist, color=color, s=0.01, rasterized=True)


def adjustStabilityPlot(ax:plt.Axes, numx:float, xlabel:str, xlist:List[float], dx:float) -> None:
    '''go back to the stability plot and clean it up.
    ax is the axes
    numx is the number of ticks on the x axis
    xlabel is the x axis label
    xlist is the list of x ticks
    dx is the spacing between ticks'''
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('z (mm)')
    ax.set_xticks(dx*np.linspace(0, numx-1, num=numx), minor=False)
    ax.set_xticklabels(xlist, fontname="Arial", fontsize=10) 
    ax.axhline(y=0, color='k')


def stabilityPlot(folder:str, exportFolder:str, tconst:float, xconst:float, export:bool=True, **kwargs) -> None:
    '''a stability plot plots cross-sections across time and space and a region plot of the steady regions
    the cross-section plot section will have two rows: one for a constant time but shifting position
    another for a constant position but shifting time
    this is all from the same simulation
    folder is the simulation folder
    exportFolder is the location to save the iamge to
    tconst is the time that stays constant in the position sweep
    xconst is the position that stays constant in the time sweep
    export true to export plot. Otherwise, just show the figure'''
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
    numx = round((xmax-xmin)/0.5)+1  # number of ticks
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
    
    if export:
        fbase = os.path.basename(folder)
        fname = 'stability_'+fbase+'_t_'+str(tconst)+'_x_'+str(xconst)
        intm.exportIm(os.path.join(exportFolder, fname), fig, **kwargs)
        