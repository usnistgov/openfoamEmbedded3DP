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
from plot.combo_plot import comboPlot
from folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class metaPlot(comboPlot):
    '''plot as text a single piece of metadata from the folderStats
    topFolder is the folder that holds all the files
    overwrite True to overwrite plots
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , **kwargs):
        super().__init__(topFolder, exportFolder=exportFolder, **kwargs)
        self.getFN(**kwargs)
        if not self.checkOverwrite(export=self.export, overwrite=self.overwrite):
            return
        
        # iterate through folders and plot data
        for i,row in self.filedf.iterrows():
            self.plotFolder(row)
            
        self.clean()
            
        if self.export:
            self.exportIm(**kwargs) 
            
        #-------------------
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'meta'
        
        
    def plotFolder(self, row) -> None:
        '''given a row in the pandas dataframe, plot the slices'''
        # identify if we need to add a new ideal plot
        fs = self.fstats[row['folder']]
        # plot basename
        pos = self.getXYRow(row)
        xm = pos['x0']
        ym = pos['y0'] + ((row['cindex']+1)/(len(self.clist)+1) - 0.5)*(self.yr[1]-self.yr[0])   # shift vertically if we have multiple colors
        pos['ax'].text(xm, ym, self.metaItem(fs), horizontalalignment='center', verticalalignment='center', c=pos['color'])
        
    def metaItem(self, fs:folderStats) -> str:
        '''get the text to plot from the folderStats object. holder text, replace this for subclasses'''
        return ''
