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


############ line plots

# files is a list of folders (e.g. nb16, nb17) to include in the plot
# func is the function to use for deciding plot colors. func should be as a function of inkvisc, supvisc, sigma
# e.g. func could be multfunc
def linePlots(folders:List[str], func, time:float, imsize:int, mode:int=0) -> plt.Figure:
    funcvals = unqListFolders(folders, func)
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False, figsize=(imsize, imsize))
    for f in folders:
        linePlot(f, time, ax, func, funcvals, mode)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    ax.set_xlabel('z (mm)')
    if mode==0:
        ax.set_ylabel('x velocity (mm/s)')
    else:
        ax.set_ylabel('Viscosity (Pa s)')
        ax.set_yscale('log')
    return fig
    
    
# plot the result of a line trace collected with paraview
# colorf is the function to use to determne the color of the plotted line, e.g. sigfunc
# rang is the total list of values that you get when you evaluate colorf over the whole folder
def linePlot(folder:str, time:float, ax:plt.Axes, colorf, rang:List[float], mode:int=0) -> None:
    t1 = intm.importLine(folder, time)
    if len(t1)==0:
        print('File is missing in ', folder)
        return 

    if mode==0:
        # velocity mode
        ystrink = 'u:0'
        ystrsup = 'u:0'
    else:
        tp = extractTP(folder)
        ystrsup = 'nu2'
        ystrink = 'nu1'
        # viscosity mode
        if 'nu2' in t1:
            t1['nu2']=t1['nu2']*10**3
        else:
            nu2list = [tp['nusup'] for i in range(len(t1))]
            t1['nu2'] = nu2list
        if 'nu1' in t1:
            t1['nu1']=t1['nu1']*10**3
        else:
            nu1list = [tp['nuink'] for i in range(len(t1))]
            t1['nu1'] = nu1list

    inkPts = t1[t1['alpha.ink']>0.5]
    minx = min(inkPts['points:2'])
    maxx = max(inkPts['points:2'])
    supPtsLeft = t1[t1['points:2']<minx]
    supPtsRight = t1[t1['points:2']>maxx]
    
    val = folderToFunc(folder, colorf)
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    fracval = decideRatio(val, rang)
    if fracval==0.5:
        color = 'gray'
    else:
        color = cmap(fracval)

    for suppts in [supPtsLeft, supPtsRight]:
        ax.plot(suppts['points:2'], suppts[ystrsup], color=color)
    ax.plot(inkPts['points:2'], inkPts[ystrink], color=color, linestyle='--')
    pts = t1[(t1['points:2']==minx) | (t1['points:2']==maxx)]
    ax.scatter(pts['points:2'], pts[ystrink], color=color, label=decideFormat(val))
    ax.scatter(pts['points:2'], pts[ystrsup], color=color)