#!/usr/bin/env python
'''Plotting tools for plotting everything on one axis'''

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
from folder_stats import folderStats
import plot.colors as co
from plot.folder_plots import folderPlots

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

class comboPlot(folderPlots):
    '''stores variables needed to create big combined plots across several folders 
    topFolder is a full path name. topFolder contains all of the folders we want to plot
        xr is the min and max x value for each section of the plot, e.g. [-0.7, 0.7]
        yr is the min and max y value for each section of the plot
        gridlines true to show gridlines'''
    
    def __init__(self, topFolder:str
                 , xr:List[float]=[-0.5, 0.5]
                 , yr:List[float]=[-0.5, 0.5]
                 , gridlines:bool=True
                 , titley:float=1.1
                 , makeXLabels:bool=True
                 , makeYLabels:bool=True
                 , makeTitle:bool=True
                 , **kwargs):
        self.makeXLabels = makeXLabels
        self.makeYLabels = makeYLabels
        self.makeTitle = makeTitle
        super().__init__(topFolder, **kwargs)
        self.type="comboPlot"
        self.gridlines = gridlines
        self.titley = titley   
        self.setupXY(xr, yr)
        if len(self.xmlist)==0:
            return
        self.setupAxes()
        self.addLegend()
        
    def setupXY(self, xr:List[float], yr:List[float]) -> None:
        '''set up the list of x and y positions to plot at'''
        self.xr = xr # x bounds for each plot chunk
        self.yr = yr
        self.dx = xr[1]-xr[0] # size of each plot chunk
        self.dy = yr[1]-yr[0]
        self.legdy = 0
        self.xrtot = [xr[0], xr[0]+(len(self.xlist)+1)*self.dx] # total bounds of the whole plot
        self.yrtot = [yr[0], yr[0]+(len(self.ylist)+1)*self.dy]
        self.xmlist = [xr[0]+(i+1/2)*self.dx for i in range(len(self.xlist))] 
            # x displacement list. this is the midpoint of each section of the plot
        self.ymlist = [yr[0]+(i+1/2)*self.dy for i in range(len(self.ylist))] # y displacement list
        
    def setupAxes(self):
        '''set up the figure and axes'''
        width = self.imsize
        ar = ((self.yrtot[1]-self.yrtot[0])*self.nrow)/((self.xrtot[1]-self.xrtot[0])*self.ncol)
        height = width*ar
        self.fig, self.axs = plt.subplots(nrows=self.nrow, ncols=self.ncol, figsize=(width,height), sharey=True, sharex=True)
        if self.nrow==1:
            self.fig.subplots_adjust(wspace=0, hspace=0)
        else:
            self.fig.subplots_adjust(wspace=0, hspace=0.05)

        if self.nrow==1:
            self.axs = [self.axs]
        if self.ncol==1:
            self.axs = [[ax] for ax in self.axs]
        
            
        # vert/horizontal grid
        if self.gridlines:
            self.addGridlines()
                
                
    def addLegend(self):
        '''Add a color legend'''
        if not self.makeLegend:
            return
        if len(self.clist)>1:
            self.ax0().legend(handles=self.legend.patches(), loc='center', ncol=self.legendCols, bbox_to_anchor=(0.7, self.titley+0.1), frameon=False)
            
    def adjustBounds(self, indices:List[int], xr:List[float], legdy:float):
        '''adjust the bounds of the plot.
        indices is a list of indices which were included in the plot
        xr is the [min, max] position of each segment, e.g. [-0.7, 0.7]
        legdy is the block height for the legend'''
        if len(indices)>0:
            pos1 = min(indices)
            pos2 = max(indices)
            dx = (xr[1]-xr[0])
            if legdy>0:
                x0 = pos1-legdy/2
            else:
                x0 = xr[0]+pos1*dx
            xrtot = [x0, xr[0]+pos2*dx+dx]
        else:
            xrtot = xr
        return xrtot
            
    def bounds(self):
        ''' adjust the bounds of the plot to only include plotted data
        each time we added a folder to the plot, we added the 
        x and y values of the centers to xlistreal and ylistreal
        this is particularly useful if a big section of the plot 
        is unplottable. e.g., very low support viscosity/ink viscosity
        values produce filaments which curl up on the nozzle, so they don't produce cross-sections.
        This step cuts out the space we set out for those folders that didn't end up in the final plot
        if we were given adjustBounds=False during initialization, don't adjust the bounds
        '''        
        if self.ab:
            self.xrtot = self.adjustBounds(self.indicesreal.x, self.xr, 0)
            self.yrtot = self.adjustBounds(self.indicesreal.y, self.yr, self.legdy)
        else:
            self.xrtot[1] = self.xrtot[1]-self.dx
            self.yrtot[1] = self.yrtot[1]-self.dy
            self.yrtot[0] = self.yrtot[0]/2+self.indicesreal.y.min()*self.dy
            
            
    def labelAx(self, ax):
        '''add axis labels and tick labels to the plot'''
        if len(self.xlistreal)>1 and self.makeXLabels:
            txt,style = self.getLabel(self.xvar, short=True)
            ax.set_xlabel(txt, **style)
            ax.set_xticks(self.xmlist, minor=False)
            labels, styles = self.getLabelList('x')
            ax.set_xticklabels(labels, fontname=self.fontName, fontsize=self.fs)
            if 'color' in styles[0]:
                for i, xtick in enumerate(axrow[0].get_xticklabels()):
                    xtick.set_color(styles[i]['color'])
        else:
            ax.set_xticks([])

        # the way comboPlots is set up, it has one big plot, 
        # and each folder is plotted in a different section of the plot
        # because the viscosity inputs are on a log scale, 
        # it is more convenient to make these plots in some real space
        # ,e.g. if we're plotting cross-sections, make the whole
        # plot go from 0-8 mm, and give the sections centers at 1, 2, 3... mm
        # then go back in later and relabel 1,2,3 to the actual viscosities, 10, 100, 1000... Pa s
        # this is the relabeling step

        if len(self.ylistreal)>1:
            ax.set_yticks(self.ymlist, minor=False)
        # make each section of the plot square
        ax.set_aspect('equal', adjustable='box')

        if len(self.xrtot)==2:
            ax.set_xlim(self.xrtot) # set the limits to the whole bounds
        if len(self.yrtot)==2:
            ax.set_ylim(self.yrtot)
            
    def removey(self, ax) -> None:
        '''remove y information from the axis'''
        ax.set_ylabel("")
        ax.set_yticks([])
        
    def subplotTitles(self) -> None:
        '''add titles to each subplot'''
        if self.ncol>1:
            labelx,lxstyle = self.getLabel(self.splitxvar, short=True, name=True, units=False)
            unitsx,uxstyle = self.getLabel(self.splitxvar, short=True, name=False, units=True)
        if self.nrow>1:
            labely,lystyle = self.getLabel(self.splityvar, short=True, name=True, units=False)
            unitsy,uystyle = self.getLabel(self.splityvar, short=True, name=False, units=True)
        if self.ncol>1:
            if self.nrow>1:
                for i,valx in enumerate(self.splitxlist):
                    for j,valy in enumerate(self.splitylist):
                        txt = self.vn.shortSymbol(f'{labelx}={valx}{unitsx}, {labely}={valy}{unitsy}')
                        self.axs[j][i].set_title(txt, **lxstyle)
            else:
                for i,valx in enumerate(self.splitxlist):
                    txt = self.vn.shortSymbol(f'{labelx}={valx}{unitsx}')
                    self.axs[0][i].set_title(txt, **lxstyle)
        else:
            if self.nrow>1:
                for j,valy in enumerate(self.splitylist):
                    txt = self.vn.shortSymbol(f'{labely}={valy}{unitsy}')
                    self.axs[j][0].set_title(txt, **lystyle)
            
    def clean(self):
        '''post-processes the plot to add components after all plots are added '''
        self.bounds()  # adjust the bounds of the plot to remove unused space
        self.axIterate(self.labelAx)  # put x labels on all plots
        self.subplotTitles()   # add subplot titles
                    

        if len(self.ylistreal)<2 or not self.makeYLabels: 
            # since there is only one variable, remove y information
            self.axIterate(self.removey)
        else:
            for axrow in self.axs:
                if not self.yvar=='sigma,ink_velocity':
                    txt,style= self.getLabel(self.yvar, True)
                    axrow[0].set_ylabel(txt, **style)
                labels, styles = self.getLabelList('y')
                axrow[0].set_yticklabels(labels, fontname=self.fontName, fontsize=self.fs)
                if 'color' in styles[0]:
                    for i, ytick in enumerate(axrow[0].get_yticklabels()):
                        ytick.set_color(styles[i]['color'])
           
        if self.makeTitle:
            self.fig.suptitle(self.fnc.figTitle, y=self.titley, fontname=self.fontName, fontsize=self.fs)
        
        if self.ab:
            # reset the figure size so the title is in the right place
            if len(self.xlistreal)>0 and len(self.ylistreal)>0:
                width = self.imsize
                height = width*((self.yrtot[1]-self.yrtot[0])*self.nrow)/((self.xrtot[1]-self.xrtot[0])*self.ncol)
                self.fig.set_size_inches(width,h=height, forward=True)
                if self.nrow==1:
                    hspace=0
                else:
                    hspace = 0.2
                self.fig.subplots_adjust(wspace=0, hspace=hspace, top=0.95, bottom=0.05, left=0.05, right=0.99)
                
        self.subFigureLabels()    # add subfigure labels if requested
        
        if not self.makeFrame:
            self.axIterate(self.removeFrame)
                
        if self.display:
            display(self.fig)
        plt.close()
