#!/usr/bin/env python
'''Handling file names'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from tools.strings import varNicknames

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#################################################################


class fnCreator:
    '''Construct an image file name with no extension. 
    Exportfolder is the folder to export to. 
    Label is any given label. 
    Topfolder is the folder this image refers to, e.g. HBHBsweep. Insert any extra values in kwargs as keywords'''
    
    def __init__(self, exportFolder:str
                 , labels:Union[list, str]
                 , topFolder:str
                 , titleVars:dict={}
                 , **kwargs) -> str:    
        self.vn = varNicknames()  
        self.figTitle = ''
        if type(labels) is list:
            for i,label in enumerate(labels):
                if i==0:
                    self.s = f'{label}'
                else:
                    self.s = f'{self.s}_{label}'
        else:
            self.s = f'{labels}'
            
        if type(topFolder) is list:
            self.bn = ','.join([os.path.basename(tf) for tf in topFolder])
        else:
            self.bn = os.path.basename(topFolder)
        self.s = f'{self.s}_{self.bn}'
        for key,val in kwargs.items():
            self.addToStack(key,val)
            self.addToFigTitle(key,val)
        for key,val in titleVars.items():
            self.addToFigTitle(key,val)
        
        self.s = self.s.replace('*', 'x')
        self.s = self.s.replace('/', 'div')
        self.s = self.vn.shorten(self.s)
        self.fn = os.path.join(exportFolder, self.bn, 'plots', self.s)
        if len(self.fn)>252:
            self.fn = self.fn[:252]
            
    def addToFigTitle(self, key:str, val:Any) -> str:
        '''get a figure title to describe the restrictions placed on this plot'''
        if '_list' in key:
            ki = self.vn.shortSymbol(key[:-5])
            vi = ', '.join([self.vn.shortSymbol(str(i)) for i in val])
            if len(self.figTitle)>0:
                self.figTitle = f'{self.figTitle}; {ki}: {vi}'
            else:
                self.figTitle = f'{ki}: {vi}'
        if key=='restrictions':
            for keyi, vali in val.items():
                self.addToFigTitle(f'{keyi}_list', vali)
        

    def addToStack(self, key:str, val:Any) -> str:
        '''add a variable definition to the name'''
        if key in ['adjustBounds'
                   , 'cname', 'colorDict', 'colorList', 'crops'
                   , 'display','dispUnits'
                   , 'ef', 'eps', 'export'
                   , 'gridlines'
                   , 'horizLabels'
                   , 'insideLabels'
                   , 'legendCols','line', 'lineDict'
                   , 'makeLegend', 'markerDict', 'markerList'
                   , 'overwrite'
                   ,'plotType', 'png'
                   , 'subLabels', 'split', 'svg'
                   ,  'xr', 'yr', 'xticks', 'yticks']:
            return ''
        if key=='restrictions':
            self.s = f'{self.s}_'
        else:
            self.s = f'{self.s}_{key}_'
        if type(val) is list:
            for ki in val:
                self.s = f'{self.s}{ki}-'
            self.s = self.s[0:-1]  # remove last dash
        elif type(val) is dict:
            for key0,val0 in val.items():
                if os.path.isdir(key0):
                    self.s = f'{self.s}{os.path.basename(key0)}_'
                else:
                    self.s = f'{self.s}{key0}_'
                if type(val0) is list:
                    for ki in val0:
                        self.s = f'{self.s}{ki}-'
                    self.s = self.s[0:-1]  # remove last dash
                else:
                    if os.path.isdir(val0):
                        self.s = f'{self.s}{os.path.basename(val0)}'
                    else:
                        self.s = f'{self.s}{val0}'
                self.s = f'{self.s}_'
            self.s = self.s[:-1]
        else:
            self.s = f'{self.s}{val}'
        self.s = self.s.replace(' ', '')
    
    def makeNew(self, export:bool=False, overwrite:bool=False) -> bool:
        '''determine if we should make a new file'''
        if overwrite:
            return True
        if not export:
            return True
        if os.path.exists(self.png):
            return
    
    def png(self):
        return f'{self.fn}.png'
    
    def svg(self):
        return f'{self.fn}.svg'
    
    def eps(self):
        return f'{self.fn}.eps'
    
    def saveFig(self, fig, fn:str) -> None:
        fig.savefig(fn, bbox_inches='tight', dpi=300, transparent=True)
    
    def export(self, fig, svg:bool=True, png:bool=True, eps:bool=False, **kwargs):
        '''export the png and/or svg'''
        if not png and not svg:
            return
        create = self.fn
        clist = []
        while not os.path.exists(os.path.dirname(create)):
            clist.append(os.path.basename(create))
            create = os.path.dirname(create)
        while len(clist)>0:
            os.mkdir(create)
            logging.info(f'Created directory {create}')
            create = os.path.join(create, clist.pop(-1))
        if svg:
            self.saveFig(fig, self.svg())
        if png:
            self.saveFig(fig, self.png())
        if eps:
            self.saveFig(fig, self.eps())
        logging.info(f'Exported {self.fn}')
