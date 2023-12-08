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
from plainIm import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#-------------------------------------------

class steadyMetrics:
    '''for determining the point at which metrics remain steady.
    steadyList determines a list of times and positions which have reached steady state, for either steady in time or steady in position.
        folder is a full path name
        dt is the size of the chunk of time in s over which we want to evaluate if the chunk is "steady", e.g. 1
        vdcrit is the maximum range of the variable of interest to be considered "steady", e.g. 0.01
        col is the column name of the variable we're watching to see if it's steady, e.g. 'vertdispn'
        mode is 'xbehind' or 'time': 
            'xbehind' means that for each x, we're finding a range of steady times. 
            'time' means that for each time, we're finding a range of steady xbehinds'''
    
    def __init__(self, folder:str, overwrite:bool=False, mode:str='time', dother:float=1, vdcrit:float=0.01, col:str='vertdispn'):
        self.folder = folder
        self.success = False
        if not os.path.exists(folder):
            return
        self.mode = mode
        self.dother = dother
        self.vdcrit = vdcrit
        self.col = col
        self.fn = self.fileName()
        if os.path.exists(self.fn) and not overwrite:
            self.success = True
            return
        self.findSteady()
        self.export()
    
    def fileName(self):
        if self.mode=='time':
            return os.path.join(self.folder, f'steadyTimes_{self.dother}_{self.col}.csv')
        elif self.mode=='xbehind':
            return os.path.join(self.folder, f'steadyPositions_{self.dother}_{self.col}.csv')
        
    def findSliceVal(self, l1:pd.DataFrame, sliceval:float) -> None:
        '''for a single slice in time or position, find the point where it has flattened out'''
        l1 = l1.sort_values(by=self.other) # sort the slice by time if mode is x, x if mode is time
        list2 = l1[self.other].unique()
        vdrange = 100
        i = -1
        while vdrange>self.vdcrit and i+1<len(list2):
            i+=1
            xtother = list2[i]
            l2 = l1[(l1[self.other]>=xtother-self.dother/2) & (l1[self.other]<=xtother+self.dother/2)] 
                # get a chunk of size dt centered around this time
            vdrange = l2[self.col].max()-l2[self.col].min()
        if i+1<len(list2):
            other0 = xtother
            while vdrange<=self.vdcrit and i+1<len(list2):
                i+=1
                xtother = list2[i]
                l2 = l1[(l1[self.other]>=xtother-self.dother/2) & (l1[self.other]<=xtother+self.dother/2)] 
                    # get a chunk of size dt centered around this time
                vdrange = l2[self.col].max()-l2[self.col].min()
            if i+1<len(list2):
                otherf = xtother
            else:
                otherf = 1000
            if self.mode=='xbehind':
                sliceval=round(sliceval,3)
            else:
                other0 = round(other0,3)
                otherf = round(otherf,3)
            self.flatlist.append({self.var0():sliceval, self.var1():other0, self.var2():otherf}) # stable within this margin
                                  
    def var0(self):
        return self.mode[0]
    
    def var1(self):
        return f'{self.o0}0'
    
    def var2(self):
        return f'{self.o0}f'
        
    def findSteady(self) -> pd.DataFrame:
        '''go through all the slices and find the range where the values are steady'''
        
        sfn = os.path.join(self.folder, 'sliceSummaries.csv')
        if not os.path.exists(sfn):
            return
        ss, ssunits = plainIm(sfn)

        self.flatlist = []
        if self.mode=='xbehind': # mode is the variable that we use to split into groups
            self.other='time' # other is the variable that we scan across
        else:
            self.other='xbehind'
                
        self.o0 = self.other[0]
        # copy units from the summary table
        self.units = {self.var0():ssunits[self.mode], self.var1():ssunits[self.other], self.var2():ssunits[self.other]} 

        if len(ss)<2:
            return pd.DataFrame([], columns=self.units.keys())

        slicevals = ss[self.mode].unique() # this gets the list of unique values for the mode variable
        for sliceval in slicevals:
            self.findSliceVal(ss[ss[self.mode]==sliceval], sliceval)
        self.df = pd.DataFrame(self.flatlist)
    
    def export(self) -> None:
        '''export the file'''
        if len(self.df)>0:
            plainExp(self.fn, self.df, self.units)
        