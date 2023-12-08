#!/usr/bin/env python
'''Functions for plotting values from simulations'''

# external packages
import sys
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot.colors import plotColors

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

class valuePlot(comboPlot):
    '''plot a circle or square representing a value from the folderStats
    topFolder is the folder that holds all the files
    overwrite True to overwrite plots
    logValues to plot the values on a log scale
    shape can be square or circle
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , cname:str='diverging'
                 , logValues:bool=False
                 , shape:str='square'
                 , dx:float=0.6
                 , timeplot:bool=False
                 , labelFreq:int=1
                 , labelLocs:str='auto'
                 , **kwargs):
        super().__init__(topFolder, exportFolder=exportFolder, xr=[-dx,dx], yr=[-dx,dx], **kwargs)
        self.getFN(**kwargs)
        if not self.checkOverwrite(export=self.export, overwrite=self.overwrite):
            return
        self.logValues = logValues
        self.shape = shape
        self.dx = dx
        self.cname = cname
        self.timeplot = timeplot
        self.labelLocs = labelLocs
        self.findValues(**kwargs) # find all of the values we are going to plot
        
        
        # plot circles in order of value
        self.filedf.sort_values(by=['value'])   
        self.filedf.reset_index(drop=True)
        
        # determine how often to place a text label
        if self.shape=='square':
            self.setColors()
        spacing = labelFreq 
        # iterate through folders and plot data
        for i,row in self.filedf.iterrows():
            self.plotFolder(row, caption=i%spacing==0)
            
        if self.shape=='square':
            self.valueLegend(**kwargs)
            
        self.clean()
            
        if self.export:
            self.exportIm(**kwargs) 
            
        #-------------------
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names. placeholder'''
        self.label = f'value'
            
    def metaItem(self, fs:folderStats) -> float:
        '''get the value to plot from the folderStats object. holder text, replace this for subclasses'''
        return np.nan
    
    def findValues(self, **kwargs) -> str:
        '''find all of the values we are going to plot'''
        for i,row in self.filedf.iterrows():
            fs = self.fstats[row['folder']]
            val = self.metaItem(fs)
            self.filedf.loc[i, 'value'] = val
        if 'minval' in kwargs:
            self.minval = kwargs['minval']
        else:
            self.minval = self.filedf.value.dropna().min()
        if 'maxval' in kwargs:
            self.maxval = kwargs['maxval']
        else:
            self.maxval = self.filedf.value.dropna().max()   
        if self.shape=='circle':
            self.rmax = self.dx # maximum radius of circle is the spacing between points/2
                
    def setColors(self):
        '''set a new color scale based on value'''
        if hasattr(self, 'var'):
            self.cvar = self.var
        else:
            self.cvar = 'value'
        l,s = self.getLabel(self.cvar)
        self.colors = plotColors(self.filedf.value.dropna()
                                 , self.cvar, l, cname=self.cname
                                 , minval=self.minval, maxval=self.maxval, byIndices=False)
        self.legend.colors = self.colors
        self.filedf.loc[:, 'cvar'] = self.filedf.value
        
        
    def scale(self, row:pd.Series) -> float:
        '''find the value to plot, as a fraction from 0 to 1'''
        val = row['value']
        if self.logValues:
            # log scale the values 
            val = np.log10(val)
            minval = np.log10(self.minval)
            maxval = np.log10(self.maxval)
        else:
            minval = self.minval
            maxval = self.maxval
        return (val-minval)/(maxval - minval)    
    
    def caption(self, row:pd.Series) -> str:
        '''get the text to display, representing this value'''
        val = row['value']
        if np.isnan(val):
            return ''
        if val>1000:
            caption='%.2e'%(val)
        elif val>10 or self.maxval>100:
            caption='%2.0f'%(val)
        else:
            caption = '%1.2f'%(val)
        return caption

        
    def plotFolder(self, row:pd.Series, caption:bool) -> None:
        '''given a row in the pandas dataframe, plot the slices'''
        # identify if we need to add a new ideal plot
        val = row['value']
        if np.isnan(val):
            return
        pos = self.getXYRow(row)  # get position in the plot
        scale = self.scale(row)   # get scaled value to plot
        
        if self.shape=='square':
            self.plotSquare(row, pos, caption, scale)
        elif self.shape=='circle':
            self.plotCircle(row, pos, caption, scale)
        else:
            raise ValueError(f'Unexpected shape value {self.shape}')
            
    def plotSquare(self, row:pd.Series, pos:dict, caption:bool, scale:float) -> None:
        '''plot a square colored by value, with a caption'''
        color = self.colors.getColor(row['value'])
        box = plt.Rectangle([pos['x0']-self.dx,pos['y0']-self.dx], self.dx*2, self.dx*2, color=color, ec=None)
        pos['ax'].add_artist(box)
        if caption:
            txt = self.caption(row)
            # calculate brightness of color
            l = 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2]
            if l>0.4:
                txtcolor = 'black'
            else:
                txtcolor = 'white'
            pos['ax'].text(pos['x0'], pos['y0'], txt, horizontalalignment='center', verticalalignment='center', c=txtcolor)
            
    def getRadius(self, scale:float) -> float:
        '''get the radius of the circle to plot'''
        return np.sqrt(scale)*self.rmax
    
    def labelInside(self, pos:dict, row:pd.Series, txt:str, color) -> None:
        '''put a text label inside the circle'''
        txtx = pos['x0']                         
        txty = pos['y0'] + 0.5*((row['cindex']+1)/(len(self.clist)+1) - 0.5)*(self.yr[1]-self.yr[0])   # shift vertically if we have multiple colors
        pos['ax'].text(txtx, txty, txt, horizontalalignment='center', verticalalignment='center', color=color)
        
    def labelOutside(self, pos:dict, row:pd.Series, radius:float, txt:str, color) -> None:
        '''put a text label outside the circle'''
        angle=(90-30*row['cindex'])*(2*np.pi/360)  # angle to put the label at, in rad                            
        dar = self.dx/3                             # length of arrow
        arrowx = pos['x0']+radius*np.cos(angle)      # where the arrow points   
        arrowy = pos['y0']+radius*np.sin(angle)
        txtx = arrowx+dar*np.cos(angle)       # where the label is
        txty = arrowy+dar*np.sin(angle)
        pos['ax'].annotate(txt, (arrowx, arrowy), color=color, xytext=(txtx, txty), ha='center', arrowprops={'arrowstyle':'->', 'color':color})
            
    def plotCircle(self, row:pd.Series, pos:dict, caption:bool, scale:float) -> None:
        '''plot a circle with a label pointing to the circle'''
        color = pos['color']
        radius = self.getRadius(scale)
        circle = plt.Circle([pos['x0'], pos['y0']], radius, color=color, fill=False) # create the circle
        pos['ax'].add_artist(circle)                                         # put the circle on the plot
        if caption:
            txt = self.caption(row)
            if self.labelLocs=='inside':
                self.labelInside(pos, row, txt, color)
            elif self.labelLocs=='outside':
                self.labelOutside(pos, row, radius, txt, color)
            elif self.labelLocs == 'auto':
                if radius>self.dx/2:
                    # if the circle is big, put the label inside
                    self.labelInside(pos, row, txt, color)
                else:
                    # if the circle is small, put the label outside
                    self.labelOutside(pos, row, radius, txt, color)
            else:
                raise ValueError(f'Unexpected value for labelLocs: {self.labelLocs}. Valid options are inside, outside, or auto.')

    def valueLegend(self, **kwargs) -> None:
        '''Put a color legend for the gradient plot on the bottom'''
        self.legend.colorBar(self.fig, **kwargs)
        
