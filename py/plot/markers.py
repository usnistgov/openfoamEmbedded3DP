#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

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
import tools.strings as st

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

class plotMarkers:
    '''for deciding on marker values for plots'''
    
    def __init__(self, vallist:list, mvar:str, mlabel:str
                 , markerSize:int, line:bool
                 , markerList:list=['o','v','s','^','X','D', 'P', '<', '>', '8', 'p', 'h', 'H']
                 , lineList:list = ['solid', 'dotted', 'dashed', 'dashdot']
                 , **kwargs):
        self.vallist = vallist  # list of values used to determine color
        self.mvar = mvar
        self.mlabel = mlabel
        self.markerSize = markerSize
        self.markerList = markerList
        self.lineList = lineList
        self.line = line
        
        if mvar=='const':
            if not 'marker' in kwargs:
                kwargs['marker'] = self.markerList[0]
            if not 'line' in kwargs:
                kwargs['line'] = self.lineList[0]
        if 'marker' in kwargs:
            self.mfunc = self.constMarker
            self.mvalFunc = self.constFunc
            self.marker = kwargs['marker']
        else:
            self.mfunc = self.listMarker
            self.mvalFunc = self.indexFunc
            
        if 'markerDict' in kwargs:
            self.mfunc = self.dictMarker
            self.mvalFunc = self.exactFunc
            self.mDict = kwargs['markerDict']
            
        if not line:
            self.lfunc = self.constLine
            self.lvalFunc = self.constFunc
            self.line0 = 'None'
        else:
            if 'lineStyle' in kwargs:
                self.lfunc = self.constLine
                self.lvalFunc = self.constFunc
                self.line = kwargs['lineStyle']
            elif 'lineDict' in kwargs:
                self.lfunc = self.dictLine
                self.lvalFunc = self.exactFunc
                self.lDict = kwargs['lineDict']
            else:
                self.lfunc = self.listLine
                self.lvalFunc = self.indexFunc
            
    #---------------------------
            
    def indexFunc(self, val:Any) -> float:
        '''get the index of this value in the list'''
        if not val in self.vallist:
            return 0
        else:
            return self.vallist.index(val)
        
    def constFunc(self, val:Any) -> float:
        return 0
    
    def exactFunc(self, val:Any) -> Any:
        return val
    
    #----------------------------
            
    def listMarker(self, val):
        '''get the marker from a list'''
        return self.markerList[val]

    def constMarker(self, val):
        '''always return the same marker'''
        return self.marker
    
    def dictMarker(self, val):
        '''get the marker from a dictionary'''
        return self.mDict[val]
    
    #-----------------------------
                                    
    def listLine(self, val):
        '''get the line from a list'''
        return self.lineList[val]
    
    def constLine(self, val):
        '''always return the same line'''
        return self.line  
    
    def dictLine(self, val):
        '''get the line style from a dictionary'''
        return self.lDict[val]
    
    #-----------------------------
    
    def getMarker(self, val):
        '''get the marker from a value'''
        return self.mfunc(self.mvalFunc(val))
    
    def getLine(self, val):
        '''get the line style from a value'''
        if not self.line:
            return 'None'
        return self.lfunc(self.lvalFunc(val))
