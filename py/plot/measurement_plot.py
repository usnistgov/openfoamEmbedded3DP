#!/usr/bin/env python
'''Functions for plotting values from simulations'''

# external packages
import sys
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot.colors import plotColors

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot.value_plot import valuePlot
from folder_stats import folderStats
from points.folder_points import folderPoints
from plainIm import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------


class measurementPlot(valuePlot):
    '''plot a circle or square representing a value from the folderStats
    topFolder is the folder that holds all the files
    overwrite True to overwrite plots
    logValues to plot the values on a log scale
    shape can be square or circle
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , var:str
                 , time:float
                 , xbehind:float
                 , xunits:str='mm'
                 , **kwargs):
        self.time = time
        self.xbehind = xbehind
        self.xunits = xunits
        self.var = var
        self.fPoints = {}   
            # dictionary of folderPoints objects
        super().__init__(topFolder, exportFolder, timeplot=False, shape='square', gridlines=False, **kwargs)
        
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'meas_{self.var}_{self.xbehind}{self.xunits}_t_{self.time}'
        
        
    def metaItem(self, fs:folderStats) -> float:
        '''get the value to plot from the folderStats object. holder text, replace this for subclasses'''
        if fs.folder in self.fPoints:
            fp = self.fPoints[fs.folder]
        else:
            fp = folderPoints(fs)
            self.fPoints[fs.folder] = fp
        
        row, u = fp.importSummarySlice(self.time, self.xbehind, self.xunits)
        if len(row)==0:
            return np.nan
        if not hasattr(self, 'units'):
            self.units = u
        if not 'behind' in self.fnc.figTitle:
            xunname = self.xunits.replace('nozzle_inner_width', '$d_i$')
            l, style = self.getLabel(self.var)
            self.fnc.figTitle = f'{l}, {self.xbehind} {xunname} behind nozzle, t = {self.time} s'
        row = row.iloc[0]
        return row[self.var]
        