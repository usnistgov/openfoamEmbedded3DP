#!/usr/bin/env python
'''Functions for plotting names of simulations for easy referencing'''

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
from plot.meta_plot import metaPlot
from folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class txtPlot(metaPlot):
    '''plot the simulation folder basenames
    topFolder is the folder that holds all the files
    overwrite True to overwrite plots
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , xr:List[float] = [-0.5, 0.5]
                 , yr:List[float] = [-0.9, 0.9]
                 , **kwargs):
        super().__init__(topFolder, exportFolder=exportFolder, xr=xr, yr=yr, **kwargs)

    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'name'
        
    def metaItem(self, fs:folderStats) -> str:
        '''get the text to plot from the folderStats object. holder text, replace this for subclasses'''
        return fs.bn
