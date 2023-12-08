#!/usr/bin/env python
'''Functions for plotting values from simulations'''

# external packages
import sys
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot.value_plot import valuePlot
from folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class ratePlot(valuePlot):
    '''plot the simulation rates as circles'''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , **kwargs):
        super().__init__(topFolder, exportFolder, timeplot=True, shape='circle', logValues=False, minval=0, maxval=300, **kwargs)
        self.figtitle = f'Simulation rate (real hr/sim s)'
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'rate'
        
        
    def metaItem(self, fs:folderStats) -> float:
        '''get the value to plot from the folderStats object. holder text, replace this for subclasses'''
        sr = fs.simulation_rate
        if sr=='':
            return np.nan
        else:
            return float(sr)