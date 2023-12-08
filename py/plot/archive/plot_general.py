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
from figureLabels import *

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


#-----------------------------------------
# ranges

def logRatio(val:float, rang:List[float]) -> float:
    '''give the value's position within range as a fraction, on a log scale. Useful for making legends.
    rang is total list of values, or just the min and max
    val is a value within the list'''
    val2 = np.log10(val)
    expmax = np.log10(max(rang))
    expmin = np.log10(min(rang))
    if expmax-expmin==0:
        return 0
    return (val2-expmin)/(expmax-expmin)


def linRatio(val:float, rang:List[float]) -> float:
    '''give the value's position within range as a fraction. Useful for making legends.
    rang is total list of values, or just the min and max
    val is a value within the list'''
    maxv = max(rang) 
    minv = min(rang)
    if maxv-minv==0:
        return 0
    return (val-minv)/(maxv-minv)


def decideRatio(val:float, rang:List[float]) -> float:
    '''give the value's position within range as a fraction
    decide whether to use a log scale or linear scale depending on the size of the range
    rang is total list of values, or just the min and max
    val is a value within the list'''
    mr = min(rang)
    if mr==0 or max(rang)-min(rang)<max(rang)/10:
        return linRatio(val, rang)
    else:
        return logRatio(val, rang)

    
def logRatioFunc(tp:Dict, func, rang:List[float]) -> float:
    '''give the value's position within range as a fraction
    tp is transport properties dictionary
    func is the function to apply to transport properties to get one value
    rang is total list of values, or just the min and max'''
    val = func(tp)
    return logRatio(val, rang)



    
#-------------------------------------------------


class gridOfPlots(folderPlots):
    '''a grid of several plots'''
    
    def __init__(self, topFolder, imsize, **kwargs):
        '''topFolder is a full path name. topFolder contains all of the folders we want to plot
        imsize is the size of EACH plot'''
        super().__init__(topFolder, imsize, **kwargs)
        self.type = 'gridOfPlots'
        self.ylist.reverse() # we reverse the rows so values go upwards up the side of the plot

        # create figure
        # in the grid of plots, the row# is y, and the col# is x
        nrows = len(self.ylist)
        ncols = len(self.xlist)
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols,\
                                sharex='col', sharey='row',\
                                figsize=(imsize*len(self.xlist), imsize*len(self.ylist)))
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        
        if nrows==1 and ncols==1:
            axs = np.array([[axs]])
        elif nrows==1 or ncols==1:
            axs = np.array([axs])
                 
        # store axes and figure in object
        self.axs = axs
        self.fig = fig
    
    
    def clean(self):
        '''post-processes the plot to add components after all plots are added'''
        if self.ab:
            self.xlistreal.sort()
            self.ylistreal.sort()
            self.xlistreal = expFormatList(self.xlistreal)
            self.ylistreal = expFormatList(self.ylistreal)
            self.xlist = expFormatList(self.xlist)
            self.ylist = expFormatList(self.ylist)
            self.ylistreal.reverse() # we reverse the rows so values go upwards up the side of the plot
            xindices = []
            yindices = []
            # in the grid of plots, the row# is y, and the col# is x

            # go through the list of original axes, and remove the axis if it wasn't used
            for i, xval in enumerate(self.xlist):
                for j, yval in enumerate(self.ylist):
                    if ((xval in self.xlistreal) and (yval in self.ylistreal)):
                        xindices.append(i)
                        yindices.append(j)
                    else:
                        self.fig.delaxes(self.axs[j,i])
            firstx = min(xindices)
            firsty = min(yindices)
            lastx = max(xindices)
            lasty = max(yindices)

            # eliminate extra axes
            ax2 = []
            dim = np.ndim(self.axs)
            if dim==2:
                for i in range(firsty, lasty+1):
                    ax2.append(self.axs[i][firstx:lastx+1])
            elif dim==1:
                ax2 = [self.axs[firstx:lastx+1]]
            else:
                ax2 = [[self.axs]]
            self.axs = ax2
        else:
            self.xlistreal = self.xlist
            self.ylistreal = self.ylist
            firstx = 0
            firsty = 0
            lastx = len(self.xlist)
            lasty = len(self.ylist)

        fs = self.fontsize
        
        # plot labels
        # in the grid of plots, the row# is y, and the col# is x
        strvals = expFormatList(self.xlistreal)
        for j, xval in enumerate(self.xlistreal):
            # going across columns
            if np.ndim(self.axs)!=1: # RG
                strval = str(strvals[j])
                self.axs[0][j].set_title(strval, fontsize=fs) # put the title on top
            else:
                self.axs[j].set_title(xval, fontsize=fs) # put the title on top for 1D plot
        strvals = expFormatList(self.ylistreal)
        for i, yval in enumerate(self.ylistreal):
            # going down rows
            strval = str(strvals[i])
            if np.ndim(self.axs)!=1: # RG
                self.axs[i][-1].text(5.1, 4, strval, verticalalignment='center', rotation=270, fontsize=fs) # put the title on the right side
                
        # reset figure size
     #   self.fig.set_size_inches(self.imsize*len(self.xlistreal), self.imsize*len(self.ylistreal))
        
        # top level axis labels
        self.fig.suptitle(self.getLabel('x', False), y=1-(firsty/len(self.ylist)), fontsize=fs) # was 0.92 instead of 1.25 RG
        if np.ndim(self.axs)!=1: # RG
            self.fig.text((lastx/len(self.xlist))+0.2, 0.5, self.getLabel('y', True),\
                      verticalalignment='center', rotation=270, fontsize=fs)
        
        #### legends
        if np.ndim(self.axs)!=1: # RG
            midleftax = self.axs[-1][0]
            midrightax = self.axs[-1][-1]
        else:
            midleftax = self.axs[0]
            midrightax = self.axs[-1]
        spcolor = '#940f0f'                 
        hatch1 = mpatches.Patch(facecolor=spcolor,alpha=0.5,hatch="\\\\\\",label='Steady in position')
        stcolor = '#356577'
        hatch2 = mpatches.Patch(facecolor=stcolor,alpha=0.2,label='Steady in time')
        midleftax.legend(handles=[hatch1, hatch2], loc='center left', bbox_to_anchor=(0, -0.6)) # was (0, -0.5) but overlapped times RG

        midrightax.legend(handles=self.plist, loc='center right', bbox_to_anchor=(1, -0.6), ncol=4) # was (0, -0.5) but overlapped times RG
        
#         self.fig.tight_layout()
        
#-------------------------------------------------

#---------------------------------------
# plotting tools

 
def addDots(ax:plt.Axes, xlist:List[float], ylist:List[float]):
    '''adds a grid of dots at intersections to the viscosity map
    ax is the axis to add dots to
    xlist is a list of x points
    ylist is a list of y pointss'''
    xl = []
    yl = []
    for x in xlist:
        for y in ylist:
            xl.append(x)
            yl.append(y)
    ax.scatter(xl, yl, color='#969696', s=10)
    return


# def adjustBounds(xlistreal:List[float], xr:List[float], xlist:List[float], dx:float):
# def adjustBounds(indices:List[int], xr:List[float], legdy:float):
#     '''adjust the bounds of the plot.
#     indices is a list of indices which were included in the plot
#     xr is the [min, max] position of each segment, e.g. [-0.7, 0.7]
#     legdy is the block height for the legend'''
#     if len(indices)>0:
#         pos1 = min(indices)
#         pos2 = max(indices)
#         dx = (xr[1]-xr[0])
#         if legdy>0:
#             x0 = pos1-legdy/2
#         else:
#             x0 = xr[0]+pos1*dx
#         xrtot = [x0, xr[0]+pos2*dx+dx]
#     else:
#         xrtot = xr
#     return xrtot

def emptyYLabels(ax:plt.Axes):
    '''Leave the y tick labels empty. Useful for side-by-side plots.'''
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_yticklabels(empty_string_labels)


    
# def multiPlots(nvars:int, imsize:float=6.5, sharey:bool=False, sharex:bool=False):
#     '''set up a plot, given nvars number of plots'''
#     if sharey:
#         minw = 6.5/4 # plots for a full width image
#     else:
#         minw = 6.5/3 # plots for a full width image
        
#     cols = int(np.floor(imsize/minw))
#     rows = int(np.ceil(nvars/cols))

#     fig, axs = plt.subplots(rows, cols, sharex=sharex, sharey=sharey, figsize=(imsize, imsize*rows/cols))
#     if rows==1:
#         if cols>1:
#             axs = np.array([axs])
#         else:
#             axs = np.array([[axs]])
#     return fig,axs
            
            
#--------------------------------------
