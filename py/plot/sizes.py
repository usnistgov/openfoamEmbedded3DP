#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

# external packages
import sys
import os
import numpy as np
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class sizes:
    '''for setting sizes of figures, fonts, and markers'''
    
    def __init__(self, rows:int, cols:int, plotType:str='notebook'):
        self.rows = rows
        self.cols = cols
        self.plotType = plotType
        if self.plotType=='ppt':
            self.fs = 18
            self.getFigSize(14, 7)
            self.markersize=100
            self.linewidth = 2
        elif self.plotType=='paper':
            self.fs = 8
            self.getFigSize(6.5, 8.5)
            self.markersize=20
            self.linewidth = 1
        elif self.plotType=='notebook':
            self.fs = 10
            self.getFigSize(10, 10)
            self.markersize = 40
            self.linewidth = 2
        else:
            raise ValueError(f'Unknown plot type {self.plotType}')
            
    def values(self):
        return self.fs, self.figsize, self.markersize, self.linewidth
            
            
    def getFigSize(self, wmax:float, hmax:float) -> None:
        self.ar = self.rows/self.cols
        wider = [wmax, wmax*self.ar]
        if wider[1]>hmax:
            wider = [w*hmax/wider[1] for w in wider]
        self.figsize = tuple(wider)