#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

# external packages
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
import itertools
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import tools.strings as st
from folder_stats import folderStats
from plot.sizes import sizes
from plot.var_plots import varPlots

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='10.0')

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

#-------------------------------------------

class folderPlots(varPlots):
    '''A generic class used for plotting many folders at once. Subclasses are comboPlot, which puts everything on one plot, and gridOfPlots, which puts everything in separate plots based on viscosity.'''
    
    def __init__(self, topFolder:str
                 , adjustBounds:bool=True
                 , **kwargs):
        super().__init__(topFolder, **kwargs)
        self.ab = adjustBounds
        self.getIndices()
            

    
    #-------------------------------------------------------------------
    
    def getIndices(self):
        '''get the position of each folder in each list'''
        
        # fill in indices
        for s in self.varlist:
            l = getattr(self, f'{s}list')
            self.filedf[f'{s}index'] = [l.index(i) for i in self.filedf[f'{s}var']]
    
    def getPos(self, fs:folderStats) -> dict:
        '''get the position in the plot for a folder'''
        folder = fs.folder
        row = self.filedf[self.filedf.folder==folder].iloc[0]
        return dict([[f'{s}index', row[f'{s}index']] for s in ['x', 'y', 'c', 'splitx', 'splity']])
    
    def getXYRow(self, row:pd.Series) -> dict:
        '''get the x,y,color, and split of the folder from the row in self.filedf'''
        x0 = self.xmlist[row['xindex']]
        y0 = self.ymlist[row['yindex']]
        color = self.colors.getColor(row['cvar'])
        axcol = row['splitxindex']
        axrow = row['splityindex'] 
        ax = self.axs[axrow][axcol]
        self.indicesreal = self.indicesreal.append({'x':row['xindex'], 'y':row['yindex'], 'color':row['cindex'], 'axx':axcol, 'axy':axrow}, ignore_index=True)
        self.xlistreal.append(x0)
        self.ylistreal.append(y0)
        self.clistreal.append(color)
        return {'x0':x0, 'y0':y0, 'color':color, 'ax':ax}
    
    #-------------------------------------------------------------------    
    
    def putAbove(self, axcol:int, axrow:int=0) -> Tuple[int,int,float,float]:
        '''get the indices and positions for the ideal plot to put the plot above the top left corner'''
        ir = self.indicesreal[(self.indicesreal.axx==axcol)&(self.indicesreal.axy==axrow)]  # select the rows in this axis
        xind = int(ir.x.min())
        x0 = self.xmlist[xind]
        yind = int(ir[ir.x==xind].y.max())
        y0 = self.ymlist[yind]+self.dy
        yind = yind+1
        return xind, yind, x0, y0

    def putLeft(self, axcol:int, axrow:int=0) -> Tuple[int,int,float,float]:
        '''get the indices to put the plot left of the bottom left corner'''
        ir = self.indicesreal[(self.indicesreal.axx==axcol)&(self.indicesreal.axy==axrow)]  # select the rows in this axis
        xind = int(ir.x.min())
        x0 = self.xmlist[xind]-self.dx
        xind = xind-1
        yind = int(ir.y.max())
        y0 = self.ymlist[yind]
        return xind, yind, x0, y0