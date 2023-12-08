#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

# external packages
import sys
import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback
import string

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import tools.strings as st
from folder_stats import folderStats
from file_names import fnCreator
import plot.colors as co
import plot.markers as ma
import plot.legend as le
import file_handling as fh
from plot.sizes import sizes

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

class varPlots:
    '''a general plotting tool for getting lists of files, setting up variables
    topFolder is the folder we're plotting
    exportFolder is the folder to export iamges to
    overwrite is true to overwrite existing images
    export is true to export images
    display is true to display the image in the jupyter notebook
    plotType can be notebook, paper, ppt
       '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str = ''
                 , overwrite:bool = False
                 , export:bool = False
                 , display:bool = True
                 , plotType:str='notebook'
                 , line:bool=False
                 , horizLabels:bool = True
                 , insideLabels:bool = True
                 , subLabels:bool = True
                 , legendCols:int=4
                , makeLegend:bool=True
                , varlist:List[str] = ['c', 'splitx', 'splity', 'x', 'y', 'm']
                 , **kwargs):
        self.topFolder = topFolder
        self.exportFolder = exportFolder
        self.overwrite=overwrite
        self.export=export
        self.display = display
        self.plotType = plotType
        self.varlist = varlist
        self.line = line
        self.horizLabels = horizLabels
        self.insideLabels = insideLabels
        self.legendCols = legendCols
        self.subLabels = subLabels
        self.makeLegend = makeLegend
        self.vn = st.varNicknames()
        self.kwargs = kwargs
        self.getFN(**kwargs)
        self.getSize(**kwargs)
        self.getVars()
        self.getRestrictions(**kwargs)
        self.getFiles()
        self.getVarLists()
        self.legendList(**kwargs)
        
    def getFileLabel(self):
        '''placeholder. overwrite this for subclasses'''
        return ''
        
    def getFN(self, **kwargs) -> None:
        '''get a file name creator'''
        self.getFileLabel()
        self.fnc = fnCreator(self.exportFolder, self.label, self.topFolder, **kwargs)
        
    def checkOverwrite(self, overwrite:bool=False, export:bool=True, **kwargs) -> bool:
        '''check if we should make a new plot'''
        if not export:
            return True
        if overwrite:
            return True
        return os.path.exists(self.fnc.png())
        
    def exportIm(self, **kwargs):
        '''export images'''
        self.fnc.export(self.fig, **kwargs)
        
    def getSize(self, plotType:str='notebook', **kwargs):
        '''determine the plot size, marker style, font size, etc.'''
        
        # load automatic variables
        if self.plotType=='ppt':
            self.fs = 18
            self.imsize=14
            self.markerSize=100
            self.linewidth = 2
            self.fontName = 'Avenir'
        elif self.plotType=='paper':
            self.fs = 8
            self.imsize = 6.5
            self.markerSize=30
            self.linewidth = 1
            self.fontName = 'Arial'
        elif self.plotType=='paperhalf':
            self.fs = 8
            self.imsize = 3.25
            self.markerSize = 30
            self.linewidth = 1
            self.fontName = 'Arial'
        elif self.plotType=='notebook':
            self.fs = 10
            self.imsize = 10
            self.markerSize = 40
            self.linewidth = 2
            self.fontName = 'Arial'
        else:
            raise ValueError(f'Unknown plot type {self.plotType}')
            
        # override values if explicitly given
        if 'fontsize' in kwargs:
            self.fs = kwargs['fontsize']
        if 'fontName' in kwargs:
            self.fontName = kwargs['fontName']
        if 'imsize' in kwargs:
            self.imsize = kwargs['imsize']
        if 'markerSize' in kwargs:
            self.markerSize = kwargs['markerSize']
        if 'linewidth' in kwargs:
            self.linewidth = kwargs['linewidth']
            
        plt.rc('font', size=self.fs) 
        
    def getVars(self):
        '''get the variables to plot on'''
        # get the variables or expressions we want to operate on
        for var in ['x', 'y', 'splitx', 'splity', 'c', 'm']:
            if f'{var}var' in self.kwargs:
                setattr(self, f'{var}var', self.kwargs.pop(f'{var}var'))
            else:
                setattr(self, f'{var}var', 'const')
                
    def getRestrictions(self, **kwargs):
        '''find the variables which have a specific list of values'''
        self.restrictions = {}
        for key,val in kwargs.items():
            if '_list' in key:
                k = key.replace('_list', '')
                self.restrictions[k] = val
        if 'restrictions' in kwargs:
            for key,val in kwargs['restrictions'].items():
                self.restrictions[key] = val
                
    def getFileList(self) -> list:
        '''get a list of all the simulation files in the topfolders'''
        return fh.simFolders(self.topFolder) # list of all folders in the top folder
    
    def getFileDF(self) -> None:
        '''get a dataframe with file metadata useful for plotting'''
        flist = self.getFileList()
        self.fstats = {}
        stats = []
        # get folderStats objects for each folder and compile critical information into a dataframe
        for folder in flist:
            fs = folderStats(folder)
            self.fstats[folder] = fs
            stats.append(fs.meta(self.xvar, self.yvar, self.cvar, self.splitxvar, self.splityvar, self.mvar, self.restrictions))  # get a list of variables relevant to the plot
        self.filedf = pd.DataFrame(stats)
        
    def filterFiles(self) -> None:
        '''filter down the filedf to only files that meet restrictions'''
        for k,val in self.restrictions.items():
            # iterate through keys and lists for restricted values and winnow the list of files
            self.filedf[k], prec = st.expFormatList(list(self.filedf[k]), returnPrecision=True)  # round floats
            l0 = st.expFormatList(val, forceFormat=((prec==1000)), returnPrecision=False, prec=prec)  # make sure all values are using the same format
            self.filedf = self.filedf[self.filedf[k].isin(l0)] # only take rows that are in the list
        if len(self.filedf)==0:
            raise FileNotFoundError('No files found that fit restrictions')
            
    def roundValues(self) -> None:
        '''round the values to eliminate rounding errors where values that should be the same are interpreted as different'''
        for var in self.varlist:
            vname = f'{var}var'
            if not getattr(self, vname) in ['xvar', 'yvar', 'var']:
                formatted, prec = st.expFormatList(list(self.filedf[vname]), returnPrecision=True)
                if prec==1000:
                    self.filedf[vname]=formatted
                elif prec>-20:
                    self.filedf[vname] = [round(v, prec) for v in self.filedf[vname]]
                
    def getFiles(self): 
        '''get a list of files that meet restrictions. varlist is a list of variables, e.g. ['x', 'y', 'c', 'splitx']'''
        self.getFileDF()  # create the dataframe
        self.filterFiles() # filter out files that don't meet restrictions
        # round the values to eliminate rounding errors
        self.roundValues()
              
        # sort by color so they're stacked in the same order at every x,y point
        if 'c' in self.varlist:
            if 'm' in self.varlist:
                l = ['cvar', 'mvar']
            else:
                l = ['cvar']
        elif 'm' in self.varlist:
            l = ['mvar']
        else:
            return
        self.filedf.sort_values(by=l, inplace=True)
            
    def mixedSort(self, l1:list) -> list:
        '''sort a list that may have mixed data types'''
        try:
            return sorted(l1)
        except:
            l2 = [0 if i=='' else i for i in l1]
            try:
                return sorted(l1)
            except:
                l1 = [str(i) for i in l1]
                return sorted(l1)
            
    def getVarLists(self):
        '''get a list of unique values that will be plotted'''
        d = {}
        if 'x' in self.varlist:
            self.xlist = self.mixedSort(self.filedf.xvar.unique())
            self.xlistreal = []
            d['x'] = []
        if 'y' in self.varlist:
            self.ylist = self.mixedSort(self.filedf.yvar.unique())
            self.ylistreal = []
            d['y'] = []
        if 'c' in self.varlist:
            self.clist = self.mixedSort(self.filedf.cvar.unique())
            self.clistreal = []
            d['color'] = []
        if 'm' in self.varlist:
            self.mlist = self.mixedSort(self.filedf.mvar.unique())
            self.mlistreal = []
            d['marker'] = []
        if 'splitx' in self.varlist:
            self.splitxlist = self.mixedSort(self.filedf.splitxvar.unique())
            self.ncol = len(self.splitxlist)
            d['axx'] = []
        else:
            self.ncol = 1
        if 'splity' in self.varlist:
            self.splitylist = self.mixedSort(self.filedf.splityvar.unique())
            self.nrow = len(self.splitylist)
            d['axy'] = []
        else:
            self.nrow = 1
        self.indicesreal = pd.DataFrame(d)

            
    def legendList(self, **kwargs):
        '''Make a legend from the list of color values and store it for later'''
        if hasattr(self, 'clist'):
            txt, txtstyle=self.getLabel(self.cvar, short=True)
            self.colors = co.plotColors(self.clist, self.cvar, txt, **kwargs)
        else:
            self.colors = co.plotColors([0], 'const', '', **kwargs)
        if hasattr(self, 'mlist'):
            txt, txtstyle=self.getLabel(self.mvar, short=True)
            self.markers = ma.plotMarkers(self.mlist, self.mvar, txt, self.markerSize, self.line, **kwargs)
        else:
            self.markers = ma.plotMarkers([0], 'const', '', self.markerSize, self.line, **kwargs)
        self.legend = le.plotLegend(self.colors, self.markers, self.filedf, **kwargs)
        return 
    
    #--------------------------------------------------------
    
    def firstFS(self) -> folderStats:
        '''get the first folderStats object'''
        if not hasattr(self, 'fs1'):
            self.fs1 = self.fstats[self.filedf.iloc[0]['folder']]
        return self.fs1
    
    def axIterate(self, func) -> None:
        '''generic function for applying action over all axes'''
        for axrow in self.axs:
            for ax in axrow:
                func(ax)
                
    def ax0(self):
        '''get the first axis'''
        return self.axs[0][0]
    
    def addGridlines(self):
        '''add gridlines to each plot'''
        self.axIterate(lambda ax:ax.grid(linestyle='-', linewidth='0.25', color='#949494'))
        
    def subFigureLabel(self, ax, label:str) -> None:
        '''add a subfigure label to the top left corner. 
        inside=True to put inside the frame, inside=False to put outside the frame'''
        if self.insideLabels:
            x=0.05
            y = 0.95
            ha = 'left'
            va = 'top'
        else:
            x = -0.05
            y = 1.05
            ha = 'right'
            va = 'bottom'

        ax.text(x, y, label
                , transform=ax.transAxes
                , horizontalalignment=ha, verticalalignment=va
                , fontsize=self.fs+2 , fontname=self.fontName)

    def subFigureLabels(self) -> None:
        '''add subfigure labels to all axes
        horiz=True to order within rows, then columns. False to order by columns, then rows
        inside=True to put inside the frame, inside=False to put outside the frame'''
        if self.nrow==1 and self.ncol==1:
            return
        if not self.subLabels:
            return
        alphabet_string = string.ascii_uppercase
        alphabet_list = list(alphabet_string)
        if self.horizLabels:
            # 2d array
            for axrow in self.axs:
                for ax in axrow:
                    self.subFigureLabel(ax, alphabet_list.pop(0))
        else:
            w = len(self.axs[0])
            h = len(self.axs)
            for i in range(w):
                for j in range(h):
                    self.subFigureLabel(self.axs[j][i], alphabet_list.pop(0))
                    
    #--------------------------------
    
    def getLabelList(self, var:str) -> Tuple[list,dict]:
        '''get a list of tick labels given the var x or y'''
        varname = getattr(self, f'{var}var')
        if varname=='sigma,ink_velocity':
            labels = []
            styles = []
            for y in getattr(self, f'{var}list'):
                l,style = self.sigmaInkLabel(y)
                labels.append(l)
                styles.append(style)
        elif 'transportModel' in varname:
            labels = [self.vn.shorten(y) for y in getattr(self, f'{var}list')]
            styles = [{} for y in getattr(self, f'{var}list')]
        else:
            labels = st.expFormatList(getattr(self, f'{var}list'))
            styles = [{} for y in getattr(self, f'{var}list')]
        return labels, styles
    
    def sigmaInkLabel(self, valy:str) -> Tuple[str, dict]:
        '''convert a sigma,ink_velocity value to a string'''
        kwargs = {'color':'Black', 'rotation':0, 'fontname':self.fontName, 'fontsize':self.fs}
        if valy=='Ideal':
            return valy, kwargs
        spl = re.split(', ', valy)
        kwargs['color'] = co.sigmaVelocityFunc(int(spl[0]), int(spl[1]))
        txt = self.vn.sigmaVelocity(valy)
        return txt, kwargs
            
    def getLabel(self, varname:str, short:bool=True, name:bool=True, units:bool=True, **kwargs) -> Tuple[str, dict]:
        '''get an axis label and style'''
        styling = {'color':'Black', 'fontname':self.fontName, 'fontsize':self.fs}
        if varname=='sigma,ink_velocity' and 'val' in kwargs:
            return self.sigmaInkLabel(kwargs['val'])
        if hasattr(self, 'units') and varname in self.units:
            # this is in the units stored in this object
            u = self.units[varname]
            if len(u)>0 and units:
                if name:
                    if 'val' in kwargs:
                        v = kwargs['val']
                        s = f'{varname}={v}{u}'
                    else:
                        s = f'{varname} ({u})'
                else:
                    if 'val' in kwargs:
                        v = kwargs['val']
                        s = f'{v}{u}'
                    else:
                        s = u
            else:
                if name:
                    if 'val' in kwargs:
                        v = kwargs['val']
                        s = f'{varname}={v}'
                    else:
                        s = varname
                else:
                    s = ''
            if short:
                txt = self.vn.symbolic(self.vn.shorten(self.vn.toEnglish(s)))
            else:
                txt = s
        else:
            # this is in a fileStats object
            txt = self.firstFS().getLabel(varname, short, name=name, units=units, **kwargs)
        return txt,styling