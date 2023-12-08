#!/usr/bin/env python
'''Plotting tools for plotting everything on one axis'''

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
from plot.multi_plot import multiPlot
import plot.colors as co
from points.folder_points import folderPoints

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

class tracePlot(multiPlot):
    '''for plotting yvar vs. x position at different times or different spacings, for various folders'''
    
    def __init__(self, topFolder:Union[List[str], str], exportFolder:str, yvar:str, cvar:str, xunits:str='niw', **kwargs):
        if cvar=='time':
            cc = 'const'
        else:
            cc = cvar
        self.cvarreal = cvar
        super().__init__(topFolder, exportFolder, 'xbehind', yvar, cvar=cc, **kwargs)  # feed initialization dummy variables because xvar and yvar might not be folderStats attributes
        self.units = {}
        self.cvar = cvar
        if cvar=='time':
            self.colors = co.plotColors([round(0.1*x,2) for x in range(26)], 'time', 'time (s)', **kwargs)
            self.legend.colors = self.colors
        self.xunits = xunits
        self.plot()
        self.clean()
        if self.export:
            self.exportIm(**kwargs) 
        
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'trace_{self.cvarreal}_{self.yvarreal}'
        
    def xbStr(self) -> str:
        '''convert xbehind to a readable string'''
        if type(self.xbehind) is dict:
            xbstr =''
            for key,val in self.xbehind.items():
                xbstr = f'{xbstr}{os.path.basename(key)}{val}_'
        else:
            xbstr = str(self.xbehind)
        return xbstr

        
    def plotFolderAx(self, ax, row:pd.Series) -> None:
        '''plot all lines for one folder on the axis'''
        fs = self.fstats[row['folder']]
        fp = folderPoints(fs)
        df, u = fp.importSummary()
        if len(df)==0:
            return
        df = df[(df.xbehind>-2)&(df.xbehind<6)] # throw out the areas too close to the boundaries
        if self.xunits=='niw' and self.xvarreal in ['x', 'xbehind']:
            df, u = fp.convertXunits(df, u, self.xunits, xvar=self.xvarreal)
        self.units = {**self.units, **u}
        
        if not self.cvar=='time':
            df = df[df.time==2.5]  # just select the slices at 2.5 seconds
            
        if self.cvar in df:
            for val in df[self.cvar].unique():
                color = self.colors.getColor(val)
                df1 = df[df[self.cvar]==val]
                ax.plot(df1[self.xvarreal], df1[self.yvarreal], color=color, marker=None, lw=0.75)
        else:
            color = self.colors.getColor(row['cvar'])
            ax.plot(df[self.xvarreal], df[self.yvarreal], color=color, marker=None, lw=0.75)
                
    def plotAx(self, row:int, col:int):
        '''plot values on the axis'''
        ax = self.axs[row][col]
        folders = self.selectFiles(row,col)
        for i,row in folders.iterrows():
            self.plotFolderAx(ax, row)
        self.addIdeal(ax, self.yvarreal)
        ax.axvline(0, ls='--', c='gray', lw=0.75)  # add a line at the nozzle
        
        
    def plot(self):
        '''plot the xvar and yvar on each axis'''
        for row in range(self.nrow):
            for col in range(self.ncol):
                self.plotAx(row, col)
        if self.cvar=='time':
            self.legend.colorBar(self.fig)