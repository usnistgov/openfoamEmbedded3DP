#!/usr/bin/env python
'''Summarize slices from all folders in a single table'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from points.folder_points import folderPoints
from folder_stats import folderStats
import file.file_handling as fh
from plainIm import *


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#------------------------------------------------------

class superSummary:
    '''given a time and position, collect slices from all of the folders and compile them into a single summary table'''
    
    def __init__(self, topFolder:Union[List[str], str]
                 , exportFolder:str
                 , time:float
                 , xbehind:float
                 , xunits:str='mm'
                 , **kwargs):
        self.topFolder = topFolder
        self.exportFolder = exportFolder
        self.time = time
        self.xbehind = xbehind
        self.xunits = xunits
        self.fileName()
        
    def fileName(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'summary_{self.xbehind}{self.xunits}_t_{self.time}.csv'
        self.fn = os.path.join(self.exportFolder, self.label)

    def importFile(self):
        '''import the summary from file'''
        if os.path.exists(self.fn):
            self.df, self.units = plainIm(self.fn, 0)
        else:
            self.getTable()
            
    def addRow(self, f:str) -> None:
        '''add the folder to the dataframe'''
        fs = folderStats(f)
        d,u = fs.metaRow()
        fp = folderPoints(fs)
        row, u2 = fp.importSummarySlice(self.time, self.xbehind, self.xunits)
        if len(row)>0:
            d = {**d, **dict(row.iloc[0])}
            u = {**u, **u2}
            self.rowlist.append(d)
            self.units = {**self.units, **u}
    
    def getTable(self):
        '''create the table by scraping data from folders'''
        flist = fh.simFolders(self.topFolder)
        self.rowlist = []
        self.units = {}
        for f in flist:
            self.addRow(f)
        self.df = pd.DataFrame(self.rowlist)
        plainExp(self.fn, self.df, self.units)
        