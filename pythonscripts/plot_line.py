#!/usr/bin/env python
'''Plotting line traces from Paraview'''

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

from plot_general import *
import interfacemetrics as intm

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------

############ line plots

  
    

def linePlot(folder:str, time:float, ax:plt.Axes, colorf, rang:List[float], mode:int=0) -> None:
    '''plot the result of a line trace collected with paraview
        colorf is the function to use to determne the color of the plotted line, e.g. sigfunc
        rang is the total list of values that you get when you evaluate colorf over the whole folder
        mode is 0 to plot velocities, 1 to plot viscosities'''
    t1,units = intm.importLine(folder, time)
    if len(t1)==0:
        print('File is missing in ', folder)
        return 

    if mode==0:
        # velocity mode
        ystrink = 'vx'
        ystrsup = 'vx'
    else:
        tp = extractTP(folder)
        ystrsup = 'nu_sup'
        ystrink = 'nu_ink'
        # viscosity mode
        if 'nu_sup' in t1:
            t1['nu_sup']=t1['nu_sup']*10**3 # this is specifically when the density is 10**3 kg/m^3
        else:
            nu_suplist = [tp['nusup'] for i in range(len(t1))]
            t1['nu_sup'] = nu_suplist
        if 'nu_ink' in t1:
            t1['nu_ink']=t1['nu_ink']*10**3
        else:
            nu_inklist = [tp['nuink'] for i in range(len(t1))]
            t1['nu_ink'] = nu_inklist

    inkPts = t1[t1['alpha']>0.5]
    zname = 'z'
    minx = min(inkPts['z'])
    maxx = max(inkPts['z'])
    supPtsLeft = t1[t1['z']<minx]
    supPtsRight = t1[t1['z']>maxx]
    
    val = folderToFunc(folder, colorf)
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    fracval = decideRatio(val, rang)
    if fracval==0.5:
        color = 'gray'
    else:
        color = cmap(fracval)

    for suppts in [supPtsLeft, supPtsRight]:
        ax.plot(suppts['z'], suppts[ystrsup], color=color)
    ax.plot(inkPts['z'], inkPts[ystrink], color=color, linestyle='--')
    pts = t1[(t1['z']==minx) | (t1['z']==maxx)]
    ax.scatter(pts['z'], pts[ystrink], color=color, label=decideFormat(val))
    ax.scatter(pts['z'], pts[ystrsup], color=color)
    
    
    
    

def linePlots(folders:List[str], func, time:float, imsize:int, mode:int=0) -> plt.Figure:
    '''files is a list of folders (e.g. nb16, nb17) to include in the plot
    func is the function to use for deciding plot colors. func should be as a function of transport properties dictionary
    e.g. func could be multfunc'''
    funcvals = unqListFolders(folders, func)
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False, figsize=(imsize, imsize))
    for f in folders:
        linePlot(f, time, ax, func, funcvals, mode)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    ax.set_xlabel('z (mm)')
    if mode==0:
        ax.set_ylabel('x velocity (mm/s)')
    else:
        ax.set_ylabel('Viscosity (Pa.s)')
        ax.set_yscale('log')
    return fig