#!/usr/bin/env python
'''Functions for plotting number of cells from simulations'''

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
from summarize.log_reader import logReader
from folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class cellsPlot(valuePlot):
    '''plot the number of cells as circles'''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , **kwargs):
        super().__init__(topFolder, exportFolder, timeplot=True, shape='circle', logValues=False, minVal=0, **kwargs)
        self.figtitle = f'Number of cells'
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'cells'
        
        
    def metaItem(self, fs:folderStats) -> float:
        '''get the value to plot from the folderStats object. import the number of cells on the last step from the logReader table'''
        lr = logReader(fs.folder)
        if hasattr(lr, 'df'):
            row = lr.df.iloc[-1]
            return row['cells']