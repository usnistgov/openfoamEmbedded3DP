#!/usr/bin/env python
'''Analyzing simulated single filaments'''

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
from file.file_handling import folderHandler
from plainIm import *


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#------------------------------------------------------

class summarizer:
    '''given a simulation, get critical statistics on the filament shape
    folder is the full path name to the folder holding all the files for the simulation
        there must be an interfacePoints folder holding csvs of the interface points
    overwrite true to overwrite existing files, false to only write new files
    a return value of 0 indicates success
    a return value of 1 indicates failure'''
    
    def __init__(self, folder:str, overwrite:bool=False):
        self.folder = folder
        self.success = False
        if not os.path.exists(folder):
            return
        self.fn = self.ssFile()
        if os.path.exists(self.fn) and not overwrite:
            self.success = True
            return
        self.fh = folderHandler(self.folder)
        self.cf = self.fh.caseFolder()
        if not os.path.exists(self.cf):
            return
        
        self.fs = folderStats(self.folder)
        self.fp = folderPoints(self.fs)
        self.summarize()
        
    def importSS(self) -> Tuple[pd.DataFrame, dict]:
        '''import slice summaries. folder is full path name'''
        self.df, self.units = plainIm(self.fn, 0)
        if len(d)==0:
            return [], []
        try:
            for s in self.df:
                self.df[s] = pd.to_numeric(self.df[s], errors='coerce') # remove non-numeric times
        except Exception as e:
            pass
        self.df = self.df.dropna() # drop NA values
        
        return self.df, self.units    
        
    def ssFile(self) -> str:
        '''slice Summaries file name'''
        return os.path.join(self.folder, 'sliceSummaries.csv')

    def addFile(self, f:str) -> None:
        '''summarize the file and add it to the stack'''
        print(f)
        data, self.units = self.fp.importPointsFile(os.path.join(self.ipfolder, f))
        if len(data)==0:
            return
        xlist = self.xlist(data)
        for x in xlist:
            sli = data[data['x']==x]
            if len(sli)>9:
                ss1 = self.sliceSummary(sli)
                self.s.append(ss1)
                    
    def getTable(self) -> pd.DataFrame:
        '''go through all the interface points files and summarize each x and time slice
        fs should be a folderStats object
        outputs a pandas DataFrame'''
        self.ipfolder = os.path.join(self.folder, 'interfacePoints')
        if not os.path.exists(self.ipfolder):
            raise ValueError('No interface points')
        ipfiles = os.listdir(self.ipfolder)
        if len(ipfiles)==0:
            logging.info(f'No slices recorded in {self.folder}')
            header = self.defaultHeader()
            self.df = pd.DataFrame([], columns=header)  
            self.units = {}
            return
        self.s = []
        
        for f in ipfiles:
            self.addFile(f)
            
        self.df = pd.DataFrame(self.s, dtype=np.float64)
        self.df.dropna(inplace=True)
        
    def saveSummary(self):
        '''save the summary points'''
        if not hasattr(self, 'df'):
            return
        self.df.sort_values(by=['time', 'x'], inplace=True)
        plainExp(self.fn, self.df, self.sliceUnits())
        
    def xlist(self, data) -> list:
        '''get a list of x positions to probe'''
        xlist = np.sort(data.x.unique())
        return xlist
        
            
    def summarize(self) -> None:
        '''summarize the slices'''
        
        # summarizeSlices will raise an exception if there are no interface
        # points files
        self.getTable()

        # if there are points in the summary, save them
        if len(self.df)>0:
            self.saveSummary()
            self.success = True
