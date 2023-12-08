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
def sigmaVelocityFunc(sigma:int, velocity:int) -> str:
    '''get a string representing this value'''
    if sigma==0:
        if velocity==0:
            return '#ad5555'
        else:
            return '#852113'
    else:
        if velocity==0:
            return '#52abcc'
        else:
            return '#1a5f78'

class plotColors:
    '''for deciding on color values for plots'''
    
    def __init__(self, vallist:list, cvar:str, clabel:str, byIndices:bool=True, logScale:bool=False, defaultColor:str='#000000', **kwargs):
        self.vallist = vallist  # list of values used to determine color
        self.cvar = cvar
        self.clabel = clabel
        self.defaultColor = defaultColor
        
        # explicitly set the bounds of the range used for value scaling
        if 'minval' in kwargs:
            self.minval = kwargs['minval']
        else:
            self.minval = min(self.vallist)
        if 'maxval' in kwargs:
            self.maxval = kwargs['maxval']
        else:
            self.maxval = max(self.vallist)
            
        # values are strings. no fractional scaling allowed
        if type(list(self.vallist)[0]) is str:
            byIndices=True
        
        if byIndices:
            # select based on the index in the list
            self.valFunc = self.indexFunc
        else:
            if logScale:
                # select based on log-scaled fractional value within a range
                self.valFunc = self.logFracFunc
            else:
                # select based on fractional value within a range
                self.valFunc = self.fracFunc
           
        if 'color' in kwargs:
            # always one color
            self.cfunc = self.oneColor
            self.valfunc = self.exactFunc
            self.color = kwargs['color']
        elif 'colorList' in kwargs:
            # select value from list of colors
            self.colorList = kwargs['colorList']
            self.cfunc = self.listFunc
            self.valFunc = self.indexFunc
        elif 'colorDict' in kwargs:
            # select value from dictionary of colors
            self.colorDict = kwargs['colorDict']
            self.cfunc = self.dictFunc
            self.valFunc = self.exactFunc
        elif 'cname' in kwargs:
            # select value from color palette
            self.cname = kwargs['cname']
            if self.cname=='cubeHelix':
                self.makeCubeHelix()
            elif self.cname=='diverging':
                self.makeDiverging()
            else:
                self.makePalette()
        elif self.cvar=='sigma':
            # use specific colors for interfacial tension
            self.colorDict = {0:'#940f0f', 40:'#61bab0', 70:'#1e5a85'}
            self.cfunc = self.dictFunc
            self.valFunc = self.exactFunc
        else:
            # use one color
            if len(self.vallist)==1:
                self.color = 'black'
                self.cfunc = self.oneColor
                self.valFunc = self.exactFunc
            else:
                self.cname = 'viridis'
                self.makePalette()
                
        

                
    def fracFunc(self, val:Any) -> float:
        '''get the position of the value scaled by value as a fraction from 0 to 1'''
        return (val-self.minval)/(self.maxval-self.minval)
        
    def logFracFunc(self, val:Any) -> float:
        return (np.log10(val)-np.log10(self.minval))/(np.log10(self.maxval)-np.log10(self.minval))
        
    def indexFunc(self, val:Any) -> float:
        '''get the position of the value scaled by order in the list as a fraction from 0 to 1'''
        if not val in self.vallist:
            raise ValueError(f'Color value {val} is not in color value list {self.vallist}')
        
        return self.vallist.index(val)/len(self.vallist)
    
    def exactFunc(self, val:Any) -> Any:
        return val
    
    #-------------------------
    
    def oneColor(self, val:Any) -> str:
        '''always return the same color'''
        return self.color
                
    def listFunc(self, val:Any) -> str:
        ''' get a color from a list of values '''
        i = int(val*len(self.vallist))
        if i<0 or i>len(self.colorList):
            return self.defaultColor
        return self.colorList[i]
    
    def dictFunc(self, val:Any) -> str:
        '''get a color from a dictionary'''
        if not val in self.colorDict:
            return self.defaultColor
        return self.colorDict[val]
    
    def makeCubeHelix(self):
        self.cmap = sns.cubehelix_palette(as_cmap=True, rot=-0.4)
        self.cfunc = self.cmapFunc
    
    def makeDiverging(self):
        self.cmap = sns.diverging_palette(220, 20, as_cmap=True)
        self.cfunc = self.cmapFunc
    
    def makePalette(self):
        self.cmap = sns.color_palette(self.cname, as_cmap=True)
        self.cfunc = self.cmapFunc
    
    def cmapFunc(self, val:Any) -> str:
        if val<0 or val>1:
            return self.defaultColor
        cmap = self.cmap
        return cmap(val)
    
    #-------------------------
    
    def getColor(self, val):
        return self.cfunc(self.valFunc(val))
             