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
import matplotlib
import colorsys
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot.multi_plot import multiPlot
from summarize.super_summary import superSummary


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#------------------------------------------------------

class superSummaryPlot(multiPlot):
    '''plot data from the supersummary'''
    
    def __init__(self, topFolder:Union[List[str], str]
                 , exportFolder:str
                 , ssFolder:str
                 , xvar:Union[str, List[str]]
                 , yvar:Union[str, List[str]]
                 , time:float
                 , xbehind:float
                 , xunits:str='niw'
                 , reference:bool=True
                 , referenceStyle:str='scatter'
                 , refLoc:str='left'
                 , **kwargs):
        self.ssFolder = ssFolder
        self.time = time
        self.xbehind = xbehind
        self.xunits = xunits
        self.reference = reference
        self.referenceStyle=referenceStyle
        self.refLoc = refLoc
        super().__init__(topFolder, exportFolder, xvar, yvar, **kwargs)  
        self.ss = superSummary(topFolder, ssFolder, time, xbehind, xunits)
        self.ssRef = superSummary(topFolder, ssFolder, time, -3, 'niw')   # get a reference slice for the initial values
        
        self.ss.importFile()   # import or create values
        self.ssRef.importFile()
        self.units = self.ss.units
        
        self.plot()
        self.clean()
        if self.export:
            self.exportIm(**kwargs) 
        
    def combine(self, var:Union[list, str]) -> None:
        if type(var) is list:
            return ','.join(var)
        else:
            return var
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'summary_{self.combine(self.xvarreal)}_{self.combine(self.yvarreal)}_{self.xbStr()}{self.xunits}_t_{self.time}'
        

    def plotAx(self, row:int, col:int) -> None:
        '''plot values on a single axis'''
        v = self.varTableRow(row, col)  # get the names of the x and y variables
        xvar = v['x']
        yvar = v['y']
        
        # filter the list to just the folders in this splitval
        folders = self.selectFiles(row,col)

        # split by color and plot scatter
        for cval in self.clist:
            for mval in self.mlist:
                foldersi = folders[(folders.cvar==cval)&(folders.mvar==mval)]    # narrow to the folders for this color
                if len(foldersi)>0:
                    ax = self.axs[row][col]                   # determine which axis to plot on
                    color = self.colors.getColor(cval)        # the color for this color value
                    marker = self.markers.getMarker(mval)
                    self.plotPoints(ax, foldersi, xvar, yvar, color, marker, mval)
                    # add a reference value ahead of the nozzle  
                    if self.reference:
                        self.addReference(ax, foldersi, xvar, yvar, color, marker, cval, mval)
        self.addIdeal(ax, yvar)
        self.separateRef(ax, xvar)
        
    def separateRef(self, ax, xvar:str):
        '''add a vertical line to separate the reference from the actual points'''
        if self.reference and self.referenceStyle=='scatter':
            ax.axvline(self.xout(xvar, 0.5), ls='--', c='gray', lw=0.75)
    
    def plotPoints(self, ax, foldersi:pd.DataFrame, xvar:str, yvar:str, color, marker, mval):
        '''plot the points and line on the plot'''
        df = self.ss.df
        ss = df[df.folder.isin(foldersi.folder)]  # get the data for these folders
        ax.scatter(ss[xvar], ss[yvar], edgecolor=color, facecolor='None', marker=marker, s=self.markerSize, linewidths=1)
        if self.line:
            ss = ss.sort_values(by=xvar)
            ls = self.markers.getLine(mval)
            ax.plot(ss[xvar], ss[yvar], color=color, marker=None, linestyle=ls)
                
    def xout(self, xvar:str, xfrac:float):
        '''get a value on the right edge of the plot'''
        df = self.ss.df
        xl = sorted(list(df[xvar].unique()))
        if len(xl)==1:
            dx = 0.1
        else:
            dx = min(np.diff(xl))
        if self.refLoc=='right':
            return xl[-1]+dx*xfrac 
        else:
            return xl[0]-dx*xfrac
        
    
    def addReference(self, ax, foldersi:pd.DataFrame, xvar:str, yvar:str, color, marker, cval, mval):
        '''add a reference point at the left to compare these values to the values ahead of the nozzle'''
        f0 = foldersi.iloc[0]['folder']
        if not self.fstats[f0].ink.dynamic['v']==0:
            # only draw a reference line if there is no flow
            return
        dfref = self.ssRef.df
        ssref = dfref[dfref.folder==f0]
        if len(ssref)==0:
            return
        if self.referenceStyle=='line':
            # plot a horizontal line across the whole plot
            ax.axhline(ssref.iloc[0][yvar], color=color, linestyle='dashed', lw=1)
        else:
            # plot one point on the left or right edge
            ax.scatter([self.xout(xvar, 1.4)], [ssref.iloc[0][yvar]], color=color, edgecolor='None', marker='X', s=self.markerSize)
            self.legend.addReferenceHandle(color, cval, mval)
            
    def makeGrayer(self, color):
        '''lower the saturation of the given color'''
        h,l,s = colorsys.rgb_to_hls(*matplotlib.colors.ColorConverter.to_rgb(color))
        return colorsys.hls_to_rgb(h,l,s*0.5)
        
    def plot(self):
        '''plot the xvar and yvar on each axis'''
        for row in range(self.nrow):
            for col in range(self.ncol):
                self.plotAx(row, col)
 