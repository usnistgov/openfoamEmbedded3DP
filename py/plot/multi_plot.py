#!/usr/bin/env python
'''Plotting tools for plotting everything on one axis'''

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

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import tools.strings as st
from folder_stats import folderStats
from plot.var_plots import varPlots
from summarize.ideals import ideals
import plot.colors as co

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

class multiPlot(varPlots):
    '''plot values across multiple axes. plots directly using regular values, no conversions like comboPlot. splitxvar and splityvar can be for different variables, or '''
    
    def __init__(self, topFolder:Union[List[str], str]
                 , exportFolder:str
                 , xvar:Union[str, List[str]]
                 , yvar:Union[str, List[str]]
                 , gridlines:bool=False
                 , titley:float=1.1
                 , short:bool=True
                 , dispUnits:bool=True
                 , scaling:float=1
                 , **kwargs):
        self.gridlines = gridlines  
        self.xvarreal = xvar
        self.yvarreal = yvar
        self.titley = titley
        self.short = short
        self.dispUnits = dispUnits
        self.referenceHandles = []
        self.ideals = ideals()
        self.scaling = scaling
        super().__init__(topFolder, exportFolder, xvar='bn', yvar='bn', **kwargs)  # feed initialization dummy variables because xvar and yvar might not be folderStats attributes
        self.setupVars()
        self.setupAxes()
        self.figtitle = ''
        
            
    def xbStr(self) -> str:
        '''convert xbehind to a readable string'''
        if type(self.xbehind) is dict:
            xbstr =''
            for key,val in self.xbehind.items():
                xbstr = f'{xbstr}{os.path.basename(key)}{val}_'
        else:
            xbstr = str(self.xbehind)
        return xbstr
    
    def selectFiles(self, row:int, col:int) -> pd.DataFrame:
        '''get the list of files that should be plotted in this axis'''
        folders = self.filedf.copy()
        if not self.splitxvar in ['xvar', 'yvar']:
            splitxval = self.splitxlist[col]
            folders = folders[folders.splitxvar==splitxval]
        if not self.splityvar in ['xvar', 'yvar']:
            splityval = self.splitylist[row]
            folders = folders[folders.splityvar==splityval]
        return folders

        
    def setupVars(self):
        '''set up the variables to plot. in multiPlots, each axis can have a different xvar and yvar. i is column and j are the position in the axes array, and x and y are the variables to plot'''
        if type(self.xvarreal) is list:
            if self.splitxvar=='xvar':
                # plot a different xvar on each column
                if type(self.yvarreal) is list:
                    if self.splityvar=='yvar':
                        # plot a different yvar on each row
                        self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                       , 'x':self.xvarreal[col]
                                                       , 'y':self.yvarreal[row]} 
                                                      for col in range(len(self.xvarreal)) 
                                                      for row in range(len(self.yvarreal))])
                    else:
                        # we have a list but nothing to plot it over
                        raise ValueError('splityvar or splitxvar must be yvar to plot list of yvars')
                else:
                    # plot same yvar on each plot in the column
                    self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                   , 'x':self.xvarreal[col]
                                                   , 'y':self.yvarreal} 
                                                  for col in range(len(self.xvarreal)) 
                                                  for row in range(len(self.splitylist))])
            elif self.splityvar=='xvar':
                # plot a different xvar on each row
                if type(self.yvarreal) is list:
                    if self.splitxvar=='yvar':
                        # plot a different yvar on each column
                        self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                       , 'x':self.xvarreal[row]
                                                       , 'y':self.yvarreal[col]} 
                                                      for col in range(len(self.yvarreal)) 
                                                      for row in range(len(self.xvarreal))])
                    else:
                        # we have a list but nothing to plot it over
                        raise ValueError('splityvar or splitxvar must be yvar to plot list of yvars')
                else:
                    # plot same yvar on each plot
                    self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                   , 'x':self.xvarreal[row]
                                                   , 'y':self.yvarreal} 
                                                  for row in range(len(self.xvarreal)) 
                                                  for col in range(len(self.splitylist))])
        else:
            # plot same xvar on each plot
            if type(self.yvarreal) is list:
                if self.splitxvar=='yvar':
                    # plot a different yvar on each column
                    self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                   , 'x':self.xvarreal
                                                   , 'y':self.yvarreal[col]} 
                                                  for row in range(len(self.splitylist)) 
                                                  for col in range(len(self.yvarreal))])
                elif self.splityvar=='yvar':
                    # plot a different yvar on each row
                    self.vartable = pd.DataFrame([{'row':row, 'col':col
                                                   , 'x':self.xvarreal
                                                   , 'y':self.yvarreal[row]} 
                                                  for row in range(len(self.yvarreal)) 
                                                  for col in range(len(self.splitxlist))])
                else:
                    # we have a list but nothing to plot it over
                    raise ValueError('splityvar or splitxvar must be yvar to plot list of yvars')
            else:
                # plot same yvar on each plot
                self.vartable = pd.DataFrame([{'row':row, 'col':col
                                               , 'x':self.xvarreal
                                               , 'y':self.yvarreal} 
                                              for col in range(len(self.splitxlist)) 
                                              for row in range(len(self.splitylist))])
                
        self.nrow = len(self.vartable.row.unique())
        self.ncol = len(self.vartable.col.unique())

    def sharey(self):
        '''determine if all plots share the same y axis'''
        if self.splityvar=='yvar':
            return 'row'
        if self.splitxvar=='yvar':
            return 'col'
        return True
    
    def sharex(self):
        '''determine if all plots share the same x axis'''
        if self.splityvar=='xvar':
            return 'row'
        if self.splitxvar=='xvar':
            return 'col'
        return True
        
    def setupAxes(self):
        '''set up the figure and axes'''
        width = self.imsize
        ar = 1
        height = width*ar
        self.fig, self.axs = plt.subplots(nrows=self.nrow
                                          , ncols=self.ncol
                                          , figsize=(width,height)
                                          , sharey=self.sharey()
                                          , sharex=self.sharex())
        
        if self.sharey():
            wspace = 0
        else:
            wspace = 0.1
        if self.sharex():
            hspace = 0
        else:
            hspace = 0.1
        if self.nrow>1:
            hspace = hspace+0.1
        self.fig.subplots_adjust(wspace=wspace, hspace=hspace)

        if self.ncol==1:
            self.axs = [self.axs]
        if self.nrow==1:
            self.axs = [self.axs]
            
        # vert/horizontal grid
        if self.gridlines:
            self.addGridlines()
        
    def varTableRow(self, row:int, col:int) -> pd.Series:
        '''get the row in the variable name table for this row and column in the figure'''
        row = self.vartable[(self.vartable.row==row)&(self.vartable.col==col)]
        return row.iloc[0]
            
    def labelAx(self, row:int, col:int):
        '''add axis labels and tick labels to the plot'''
        v = self.varTableRow(row, col)
        xvar = v['x']
        yvar = v['y']
        if row==self.nrow-1 or not self.sharex():
            # only set x label if there is a unique value
            txt, txtkwargs = self.getLabel(xvar, short=True, name=True, units=True)
            self.axs[row][col].set_xlabel(txt, **txtkwargs)
        if col==0 or not self.sharey():   
            # only set y label if there is a unique value
            yname, txtkwargs = self.getLabel(yvar, short=True, name=True, units=True)
            if not self.scaling==1:
                yname = f'{yname}$\\times${self.scaling}'
            self.axs[row][col].set_ylabel(yname, **txtkwargs)

        

    def subplotTitles(self) -> None:
        '''add titles to the top and right edges of the plot'''
        if not self.splitxvar in ['xvar', 'yvar'] and self.ncol>1:
            # label the columns
            for i,valx in enumerate(self.splitxlist):
                txt, txtkwargs = self.getLabel(self.splitxvar, short=True, name=True, units=(i==len(self.splitxlist)-1), val=valx)
                self.axs[0][i].set_title(txt, **txtkwargs)
        if not self.splityvar in ['xvar', 'yvar'] and self.nrow>1:
            # label the rows
            for j,valy in enumerate(self.splitylist):
                txt, txtkwargs = self.getLabel(self.splityvar, short=True, name=True, units=(j==len(self.splitylist)-1), val=valy)
                if not 'rotation' in txtkwargs:
                    txtkwargs['rotation'] = 270
                ax = self.axs[j][-1]
                ax.text(1.03, 0.5, txt, transform=ax.transAxes, ha='left', va='center', **txtkwargs)
                
            
    def setSquare(self, ax):
        '''make the axis square'''
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
        
    def fixTicks(self, ax):
        '''manually set ticks if they were fed into the function call in the first place'''
        if 'xticks' in self.kwargs:
            ax.set_xticks(self.kwargs['xticks'], minor=False)
        if 'yticks' in self.kwargs:
            ax.set_yticks(self.kwargs['yticks'], minor=False)
            
    def clean(self):
        '''post-processes the plot to add components after all plots are added '''
        # put x labels on all plots
        for row in range(self.nrow):
            for col in range(self.ncol):
                self.labelAx(row, col)
        self.addLegend()
        self.subplotTitles()            
        self.fig.suptitle(self.figtitle, y=self.titley, fontname=self.fontName, fontsize=self.fs)
        self.axIterate(self.setSquare)  # make all axes square
        self.axIterate(self.fixTicks)
        self.subFigureLabels()    # add subfigure labels if requested
        width = self.imsize
        height = width*(self.nrow)/(self.ncol)
        self.fig.set_size_inches(width,h=height, forward=True)
        if self.display:
            display(self.fig)
        plt.close()

    def addLegend(self):
        '''Add a color legend'''
        if not self.makeLegend:
            return
        if len(self.clist)>1:
            self.ax0().legend(handles=self.legend.patches(), loc='center', ncol=self.legendCols, bbox_to_anchor=(1.3, self.titley+0.2), frameon=False)
            
    def addIdeal(self, ax, yvar):
        '''add an ideal value as a horizontal line on the axis'''
        ideal = self.ideals.getValue(yvar)
        if not ideal=='':
            ax.axhline(ideal, ls='--', c='gray', lw=0.75)
      #      ax.text(min(df[xvar]), ideal, 'ideal', ha='left', va='bottom', color='gray') # key for ideal horizontal lines
        
