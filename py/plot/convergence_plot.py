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
from summarize.log_reader import logReader

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

class convergencePlot(multiPlot):
    '''for plotting residuals vs. time, for various folders'''
    
    def __init__(self, topFolder:Union[List[str], str], exportFolder:str, yvar:str, **kwargs):
        super().__init__(topFolder, exportFolder, 'simTime', yvar, **kwargs)  # feed initialization dummy variables because xvar and yvar might not be folderStats attributes
        self.units = {}
        self.plot()
        self.clean()
        
        if self.export:
            self.exportIm(**kwargs) 
        
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'convergence_{self.yvarreal}'
        
    def plotFolderAx(self, ax, row:pd.Series) -> None:
        '''plot all lines for one folder on the axis'''
        lr = logReader(row['folder'])
        df = lr.df
        self.units = {**self.units, **lr.u}

        color = self.colors.getColor(row['cvar'])
        df=df[df[self.yvarreal]>0]
        ylist = df[self.yvarreal]
        if not self.scaling==1:
            ylist = ylist/self.scaling
        ax.plot(df['simtime'], ylist, color=color, marker=None, lw=0.5)
                
    def plotAx(self, row:int, col:int):
        '''plot values on the axis'''
        ax = self.axs[row][col]
        folders = self.selectFiles(row,col)
        for i,row in folders.iterrows():
            self.plotFolderAx(ax, row)
        
    def plot(self):
        '''plot the xvar and yvar on each axis'''
        for row in range(self.nrow):
            for col in range(self.ncol):
                self.plotAx(row, col)
