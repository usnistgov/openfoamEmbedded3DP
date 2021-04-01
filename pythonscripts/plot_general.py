#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
import itertools
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

import interfacemetrics as intm
import folderparser as fp

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10


__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------

###############################################################################################################
###############################################################################################################
################# PLOTS #######################################################################################
###############################################################################################################
###############################################################################################################

# converts kinematic viscosity to dynamic viscosity
# can be a single float or a list of floats
def kinToDyn(v:Union[List[float], float]) -> Union[List[float], float]:
    if type(v) is list:
        v2 = []
        for vi in v:
            v2.append(kinToDyn(vi))
        return v2
    else:
        v2 = v*1000
        if v2>=1:
            v2 = int(round(v2))
        return v2

# preset colors for surface tensions
def sigmaColor(sigma:int) -> str:
    if sigma==0:
        return '#940f0f'
    elif sigma==40:
        return '#61bab0'
    elif sigma==70:
        return '#1e5a85'

# put a number in exponential format
def expFormat(x:float) -> str:
    return r'$10^{{{}}}$'.format(int(round(np.log10(x))))

# determine what kind of format to put the number in
# might return float or string
def decideFormat(x:float) -> Union[float, str]:
    if x==0:
        return x
    l = np.log10(x)
    if not l==round(l):
        return x
    else:
        return expFormat(x)

# put the whole list in exponential format
def expFormatList(xlist:List[float]) -> List[Any]:
    xout = []
    for x in xlist:
        xout.append(decideFormat(x))
    return xout  


#{'nuink':nuink, 'tau0ink':tau0ink, 'kink':kink, 'nink':nink, 'nusup':nusup, 'tau0sup':tau0sup, 'ksup':ksup, 'nsup':nsup, 'sigma':sigma}

# use this function to just pluck a value out of the transport properties
def tpFunc(tp:Dict, stri:str) -> float:
    return tp[stri]

# x,y functions
    # inkv is an ink viscosity
    # supv is a support viscosity
    # sigma is a surface tension
def multfunc(tp:Dict) -> float:
    return round(tp['nuink']*tp['nusup'], 10)

def divfunc(tp:Dict) -> float:
    return round(tp['nusup']/tp['nuink'], 10)

# color functions
def sigfuncc(tp:Dict) -> str:
    return sigmaColor(tp['sigma'])

# val should be 0-1
# returns a color
def cubehelix1(val:float):
    cm = sns.cubehelix_palette(as_cmap=True, rot=-0.4)
    return cm(val)

# give the value's position within range as a fraction, on a log scale
# rang is total list of values, or just the min and max
# val is a value within the list
def logRatio(val:float, rang:List[float]) -> float:
    val2 = np.log10(val)
    expmax = np.log10(max(rang))
    expmin = np.log10(min(rang))
    if expmax-expmin==0:
        return 0
    return (val2-expmin)/(expmax-expmin)

# give the value's position within range as a fraction
# rang is total list of values, or just the min and max
# val is a value within the list
def linRatio(val:float, rang:List[float]) -> float:
    maxv = max(rang) 
    minv = min(rang)
    if maxv-minv==0:
        return 0
    return (val-minv)/(maxv-minv)

# give the value's position within range as a fraction
# decide whether to use a log scale or linear scale depending on the size of the range
# rang is total list of values, or just the min and max
# val is a value within the list
def decideRatio(val:float, rang:List[float]) -> float:
    mr = min(rang)
    if mr==0 or max(rang)-min(rang)<max(rang)/10:
        return linRatio(val, rang)
    else:
        return logRatio(val, rang)

# give the value's position within range as a fraction
# inkv is ink viscosity
# supv is support viscosity
# sigma is surface tension
# func is the function to apply to inkv, supv, sigma to get one value
# rang is total list of values, or just the min and max
def logRatioFunc(tp:Dict, func, rang:List[float]) -> float:
    val = func(tp)
    logRatio(val, rang)

# given a folder, get one value to represent that folder
# func is the function to apply to inkv, supv, sigma to get one value
def folderToFunc(folder:str, func) -> float:
    tp = extractTP(folder)
    return func(tp)

def tpCombos(tplists:Dict) -> List:
    vallists = []
    for l in tplists:
        vallist = [[l[0:-4],i] for i in tplists[l]] # e.g. [['nuink', 5], ['nuink', 6]]
        vallists.append(vallist)
    l0 = list(itertools.product(*vallists))
    l0 = [dict(l) for l in l0]
    return l0

# unqList finds the unique x or y positions
    # f is the function to apply to ink viscosity, support viscosity, sigma
    # nuinklist is a list of ink viscosities
    # nusuplist is a list of support viscosities
    # sigmalist is a list of surface tensions
def unqList(f, tplists:Dict) -> List[float]:
    lout = []
    combos = tpCombos(tplists)
    for tp in combos:
        l = f(tp)
        if l not in lout:
            lout.append(l)
    lout.sort()
    return lout

# given a folder of many simulations, find unique x or y positions in comboPlot or gridOfPlots
def unqListFolders(folders:List[str], func) -> List[float]:
    funcvals = [] # outputs for the func
    for f in folders:
        val = folderToFunc(f, func)
        if val not in funcvals:
            funcvals.append(val)
    return funcvals


    
    
#### WORKING WITH SINGLE FILES




    #############################################################
#### WORKING WITH MULTIPLE PLOTS
    
# find the position of v in list l
    # l is a list
    # v is a value
def findPos(l, v):
    try:
        p = l.index(v)
    except ValueError:
        return -1
    return p

# find values needed for plotting in grids
    # nuink is the viscosity of the ink
    # nusup is the viscosity of the support
    # sigma is the surface tension
    # xpv is a comboPlot or gridOfPlots object
def vv(tp:Dict, xpv):
    x = xpv.xfunc(tp)
    y = xpv.yfunc(tp)
    xpos = findPos(xpv.xlist, x)
    ypos = findPos(xpv.ylist, y)
    sigmapos = findPos(xpv.tplists['sigmalist'], tp['sigma'])
    if xpos<0 or ypos<0 or sigmapos<0:
        raise ValueError
    if xpv.type=="comboPlot":
        x0 = xpv.xmlist[xpos]
        y0 = xpv.ymlist[ypos]
    else:
        x0 = xpos
        y0 = len(xpv.ylist)-ypos-1
    if x not in xpv.xlistreal:
        xpv.xlistreal.append(x)
    if y not in xpv.ylistreal:
        xpv.ylistreal.append(y)
    color = xpv.cfunc(tp)
    return color, x0, y0, sigmapos

# le is the legend
# nui is the index of the viscosity
# getHB is true to get herschel-bulkley parameters
def extractTPfluid(le:pd.DataFrame, nui:int, getHB:bool) -> Tuple[float]:
    nu = kinToDyn(float(le.loc[nui, 'val'])) # get the actual viscosity in Pa s
    if getHB:
        tau0 = kinToDyn(float(le.loc[nui+1, 'val']))
        k = kinToDyn(float(le.loc[nui+2, 'val']))
        n = float(le.loc[nui+3, 'val'])
    else:
        tau0 = 0
        k = 0
        n = 0
    return nu, tau0, k, n

# extractTP gets the transport properties (viscosities and surface tension) for a folder
    # folder is a full path name
def extractTP(folder:str) -> Tuple[float, float, float]:
    le = intm.importLegend(folder)
    if len(le)==0:
        raise Exception('No values in legend')
    i0 = le[le['title']=='transportProperties'].index[0] # find the index where the transport properties start
    # express viscosities as a dynamic viscosity 
    # (multiply it by 1000 because all of our simulations were done with a density of 1000)
    inkmode = le.loc[i0+2, 'val']
    # find the row that contains the ink model.
    if inkmode=='Newtonian':
        # if it's newtonian, the viscosity is stored 3 rows below 'transportProperties'
        nuinki = i0+3
        getHB = False
    else:
        # otherwise the zero shear viscosity is 4 rows below
        nuinki = i0+4
        getHB = True
    nuink, tau0ink, kink, nink = extractTPfluid(le, nuinki, getHB)
    
    # find the support transport properties
    # in some legends, there are no Herschel Bulkley boxes, and in some there are
    # we need to find where the support section starts
    nusupi = nuinki 
    
    # iterate through rows until we find the support transportmodel
    while not le.loc[nusupi, 'title']=='transportModel':
        nusupi+=1
        
    if le.loc[nusupi, 'val']=='Newtonian':
        nusupi = nusupi+1
        getHB = False
    else:
        nusupi = nusupi+2
        getHB = True
    nusup, tau0sup, ksup, nsup = extractTPfluid(le, nusupi, getHB)    
    
    # likewise we need to find the row where the surface tension is stored 
    # because we don't know if the legend contains HB rows or not
    sigmai = nusupi
    while not le.loc[nusupi, 'title']=='sigma':
        nusupi+=1
        
    # convert surface tension to mJ/m^2
    sigma = int(round(1000*float(le.loc[nusupi, 'val'])))
#     return nuink, nusup, sigma
    return {'nuink':nuink, 'tau0ink':tau0ink, 'kink':kink, 'nink':nink, 'nusup':nusup, 'tau0sup':tau0sup, 'ksup':ksup, 'nsup':nsup, 'sigma':sigma}

# vvplot finds variables used for plotting in grids
    # folder is a folder name
    # xpv is a comboPlot or gridOfPlots object
def vvplot(folder:str, xpv):
#     nuink, nusup, sigma = extractTP(folder)
    tp = extractTP(folder)
    return vv(tp, xpv)

# get a list of all of the viscosity and surface tension values in this list of folders
    # flist is a list of folders (full path names)
def listTPvalues(flist, **kwargs):
    flist2 = flist
    ap0 = False # ap0 tells us if we should automatically include or exclude files
    tp = extractTP(flist[0]) # this gives us an initial list of strings that extractTP pulls out
    lists = dict([[s+'list', []] for s in tp])
    
    # if we've already defined sup viscosity, ink viscosity, or surface tension list, only use the files in that list
    for l in lists:
        if l in kwargs:
            flist2 = []
            ap0 = True 
            lists[l] = kwargs[l]

    # for each folder, extract transport properties and add to lists
    for folder in flist:
        tp = extractTP(folder)
        append = ap0
        
        for s in tp:
            l = s+'list'
            if l in kwargs:
                # if we've already designated a list, don't include this file if it's not in the list
                if not tp[s] in lists[l]:
                    append = False
            else:
                # if we haven't designated a list, include this file
                if not tp[s] in lists[l]:
                    lists[l].append(tp[s])
        
        if append:
            flist2.append(folder)
    
    # sort all of the lists
    for l in lists:
        lists[l].sort()
        
    return flist2, lists

################# 
### CLASSES

#### WORKING WITH GRIDS




#-------------------------------------------------


# topFolder is the folder we're plotting
# imsize is the size of the total image in inches
# split is true to split by surface tension
class folderPlots:
    def __init__(self, topFolder:str, imsize:float, split:bool=False, **kwargs):
        self.kwargs = kwargs
        self.ab = not 'adjustBounds' in self.kwargs or self.kwargs['adjustBounds']==True
        self.topFolder = topFolder
        self.imsize = imsize
        self.split = split
        self.plotsLists(**kwargs)
        
    def convertFunc(self, var):
        if len(var)>0:
            if var in self.strings:
                func = lambda tp: tpFunc(tp, var)
                return func
            else:
                for s in self.strings:
                    var = var.replace(s, "tp[\'"+s+"\']")
                func = lambda tp: eval(var)
                return func
                

                
        
    # plotsLists initializes gridOfPlots and comboPlots objects
    # obj is a gridOfPlots or comboPlots object
    def plotsLists(self, **kwargs):
        self.flist = fp.caseFolders(self.topFolder) # list of all folders in the top folder
        self.flist, self.tplists = listTPvalues(self.flist, **kwargs) # list of transport property lists
        self.cfunc = sigfuncc
        
        self.strings = [s[0:-4] for s in list(self.tplists.keys())]
        
        if 'xvar' in kwargs and 'yvar' in kwargs:
            self.xvar = kwargs['xvar']
            self.yvar = kwargs['yvar']
        if 'mode' in kwargs:
            self.xvar = ''
            self.yvar = ''
            self.mode = kwargs['mode']
            if self.mode==0:
                # plot ink viscosity vs. support viscosity
                self.xvar = 'nuink'
                self.yvar = 'nusup'
            elif self.mode==1:
                # plot ink viscosity/support viscosity vs. inkviscosity*supportviscosity
                self.xfunc = multfunc
                self.yfunc = divfunc
            elif self.mode==2:
                self.xfunc = divfunc
                self.yfunc = sigfunc
            elif self.mode==3:
                self.xfunc = supfunc
                self.yfunc = divfunc
            elif self.mode==4:
                self.xfunc = supfunc
                self.yfunc = divfunc
            elif self.mode==5:
                self.xvar = 'nuink'
                self.yvar = 'nusup'
            else:
                raise ValueError('Invalid mode: '+self.mode)    
        else:
            if not 'xvar' in kwargs and not 'yvar' in kwargs:
                raise ValueError('Invalid function: set mode or set xvar and yvar')
            if not 'xvar' in kwargs:
                raise ValueError('Invalid function: missing xvar') 
            if not 'yvar' in kwargs:
                raise ValueError('Invalid function: missing yvar')  
        
        self.xfunc = self.convertFunc(self.xvar)
        self.yfunc = self.convertFunc(self.yvar)
        
        try:
            self.xlist = unqList(self.xfunc, self.tplists)
            self.ylist = unqList(self.yfunc, self.tplists)
            self.xlistreal = []
            self.ylistreal = []
            self.legendList()
        except:
            raise ValueError('Failed to identify x and y variables')
        return self
    
    def legendList(self):
        if not self.split:
            sigmalist = self.tplists['sigmalist']
            plist = [mpatches.Patch(color=sigmaColor(sigmalist[i]), label=sigmalist[i]) for i in range(len(sigmalist))]
            ph = [plt.plot([],marker="", ls="", label='\u03C3 (mJ/m$^2$)')[0]]; # Canvas
            self.plist = ph + plist
            plt.close()
        return 
    
    def getLabel(self, var, short):
        
        # determine which axis we're trying to name
        if var=='x':
            func = self.xfunc
        else:
            func = self.yfunc
         
        namedefs = {'nuink':'Ink viscosity (Pa$\cdot$s)', 'tau0ink':'Ink yield stress (Pa)','kink':'Ink k (Pa*s^n)', 'nink':'Ink n',\
                    'nusup':'Support viscosity (Pa$\cdot$s)', 'tau0sup':'Support yield stress (Pa)', 'ksup':'Support k (Pa*s^n)', 'nsup':'Support n',\
                    'sigma':'Surface tension (mJ/m$^2$)', 'product':'Support viscosity \u00d7 ink viscosity (Pa$^2\cdot$s$^2$)', 'ratio':'Support viscosity / ink viscosity'}
        
        # convert the function to a label
        if func==multfunc:
            varstr = 'product'
        elif func==divfunc:
            label = 'ratio'
        else:
            # input was xvar yvar, not mode
            if var=='x':
                varstr = self.kwargs['xvar']
            else:
                varstr = self.kwargs['yvar']
             
        if varstr in namedefs:
            label = namedefs[varstr]
        else:
            for s in namedefs:
                varstr = varstr.replace(s, namedefs[s])
            label = varstr

        if short:
            label = label.replace('viscosity', '\u03B7')
            label = label.replace('Surface tension', '\u03C3')
        return label
    
#-------------------------------------------------

# this creates a grid of several plots
    # topFolder is a full path name. topFolder contains all of the folders we want to plot
    # imsize is the size of each image

class gridOfPlots(folderPlots):
    def __init__(self, topFolder, imsize, **kwargs):
        super().__init__(topFolder, imsize, **kwargs)
        self.type = 'gridOfPlots'
        self.ylist.reverse() # we reverse the rows so values go upwards up the side of the plot

        # create figure
        # in the grid of plots, the row# is y, and the col# is x
        fig, axs = plt.subplots(nrows=len(self.ylist), ncols=len(self.xlist),\
                                sharex='col', sharey='row',\
                                figsize=(imsize*len(self.xlist), imsize*len(self.ylist)))
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
                 
        # store axes and figure in object
        self.axs = axs
        self.fig = fig 
    
    # clean post-processes the plot to add components after all plots are added
    def clean(self):
        
        if self.ab:
            self.xlistreal.sort()
            self.ylistreal.sort()
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
            for i in range(firsty, lasty+1):
                ax2.append(self.axs[i][firstx:lastx+1])
            self.axs = ax2
        else:
            self.xlistreal = self.xlist
            self.ylistreal = self.ylist
            firstx = 0
            firsty = 0
            lastx = len(self.xlist)
            lasty = len(self.ylist)

        fs = 10
        
        # plot labels
        # in the grid of plots, the row# is y, and the col# is x
        for j, xval in enumerate(self.xlistreal):
            # going across columns
            strval = str(expFormat(xval)) # title
            self.axs[0][j].set_title(strval, fontsize=fs) # put the title on top
        for i, yval in enumerate(self.ylistreal):
            # going down rows
            strval = str(expFormat(yval)) # title
            self.axs[i][-1].text(5.1, 4, strval, verticalalignment='center', rotation=270, fontsize=fs)
                # put the title on the right side
                
        # reset figure size
     #   self.fig.set_size_inches(self.imsize*len(self.xlistreal), self.imsize*len(self.ylistreal))
        
        # top level axis labels
        self.fig.suptitle(self.getLabel('x', False), y=0.92-(firsty/len(self.ylist)), fontsize=fs)
        self.fig.text((lastx/len(self.xlist))-0.1, 0.5, self.getLabel('y', False),\
                      verticalalignment='center', rotation=270, fontsize=fs)
        
        #### legends
        midleftax = self.axs[-1][0]
        midrightax = self.axs[-1][-1]
        spcolor = '#940f0f'                 
        hatch1 = mpatches.Patch(facecolor=spcolor,alpha=0.5,hatch="\\\\\\",label='Steady in position')
        stcolor = '#356577'
        hatch2 = mpatches.Patch(facecolor=stcolor,alpha=0.2,label='Steady in time')
        midleftax.legend(handles=[hatch1, hatch2], loc='center left', bbox_to_anchor=(0, -0.5))

        midrightax.legend(handles=self.plist, loc='center right', bbox_to_anchor=(1, -0.5), ncol=4)
        
#         self.fig.tight_layout()
        
#-------------------------------------------------
        
### this class stores variables needed to create big combined plots across several folders 
    # topFolder is a full path name. topFolder contains all of the folders we want to plot
    # xr is the min and max x value for each section of the plot, e.g. [-0.7, 0.7]
    # yr is the min and max y value for each section of the plot
    # imsize is the size of the whole image in inches
    # gridlines true to show gridlines
class comboPlot(folderPlots):
    def __init__(self, topFolder:str, xr:List[float], yr:List[float], imsize:float, gridlines:bool=True, **kwargs):
        super().__init__(topFolder, imsize, **kwargs)
        self.type="comboPlot"
        
        self.figtitle = ""
        self.xr = xr # x bounds for each plot
        self.yr = yr
        self.dx = xr[1]-xr[0] # size of each plot chunk
        self.dy = yr[1]-yr[0]
        self.xrtot = [xr[0], xr[0]+(len(self.xlist)+1)*self.dx] # total bounds of the whole plot
        self.yrtot = [yr[0], yr[0]+(len(self.ylist)+1)*self.dy]
        self.xmlist = [xr[0]+(i+1/2)*self.dx for i in range(len(self.xlist))] 
            # x displacement list. this is the midpoint of each section of the plot
        self.ymlist = [yr[0]+(i+1/2)*self.dy for i in range(len(self.ylist))] # y displacement list

        # if split, make a row of 3 plots. If not split, make one plot
        if self.split:
            ncol = len(self.tplists['sigmalist'])
        else:
            ncol = 1
        self.ncol = ncol
        self.imsize = imsize
        fig, axs = plt.subplots(nrows=1, ncols=ncol, figsize=(imsize,imsize*len(self.ylist)/len(self.xlist)), sharey=True)
        fig.subplots_adjust(wspace=0)
        
        if not self.split:
            axs = [axs]
        
        # vert/horizontal grid
        if gridlines:
            for ax in axs:
                ax.grid(linestyle='-', linewidth='0.25', color='#949494')

        # set position of titley
        if not self.split:
            self.titley = 1
        else:
            self.titley = 1
        # store variables
        self.axs = axs
        self.fig = fig 
        
        self.addLegend()
        
    def addLegend(self):
        if not self.split:
            self.axs[0].legend(handles=self.plist, loc='upper center', ncol=4, bbox_to_anchor=(0.5, self.titley+0.1))
        
    # clean post-processes the plot to add components after all plots are added    
    def clean(self):
        # adjust the bounds of the plot to only include plotted data
        # each time we added a folder to the plot, we added the 
        # x and y values of the centers to xlistreal and ylistreal
        # this is particularly useful if a big section of the plot 
        # is unplottable. e.g., very low support viscosity/ink viscosity
        # values produce filaments which curl up on the nozzle, so they don't produce cross-sections.
        # This step cuts out the space we set out for those folders that didn't end up in the final plot
        
        
        
        if self.ab:
            self.xrtot = adjustBounds(self.xlistreal, self.xr, self.xlist)
            self.yrtot = adjustBounds(self.ylistreal, self.yr, self.ylist)
        else:
            self.xrtot[1] = self.xrtot[1]-self.dx
            self.yrtot[1] = self.yrtot[1]-self.dy
        
        
        # put x labels on all plots
        for ax in self.axs:
            ax.set_xlabel(self.getLabel('x', len(self.axs)>1), fontname="Arial", fontsize=10)

            # the way comboPlots is set up, it has one big plot, 
            # and each folder is plotted in a different section of the plot
            # because the viscosity inputs are on a log scale, 
            # it is more convenient to make these plots in some real space
            # ,e.g. if we're plotting cross-sections, make the whole
            # plot go from 0-8 mm, and give the sections centers at 1, 2, 3... mm
            # then go back in later and relabel 1,2,3 to the actual viscosities, 10, 100, 1000... Pa s
            # this is the relabeling step
            ax.set_xticks(self.xmlist, minor=False)
            ax.set_yticks(self.ymlist, minor=False)
            
            ax.set_xticklabels(expFormatList(self.xlist), fontname="Arial", fontsize=10)  
            #emptyYLabels(ax)
            if len(self.xrtot)==2:
                ax.set_xlim(self.xrtot) # set the limits to the whole bounds
            if len(self.yrtot)==2:
                ax.set_ylim(self.yrtot)

            # make each section of the plot square
            ax.set_aspect('equal', adjustable='box')
            
        if self.split:
            sigmalist = self.tplists['sigmalist']
            for i in range(len(sigmalist)):
                self.axs[i].set_title('\u03C3='+str(sigmalist[i])+' mJ/m$^2$', fontname="Arial", fontsize=10)

        self.axs[0].set_ylabel(self.getLabel('y', False), fontname="Arial", fontsize=10)
        
        yticklabels = expFormatList(self.ylist)
        self.axs[0].set_yticklabels(yticklabels, fontname="Arial", fontsize=10)
        
        # reset the figure size so the title is in the right place
        if self.ab and len(self.xlistreal)>0 and len(self.ylistreal)>0:
            width = self.imsize
            height = width*len(self.ylistreal)/(len(self.xlistreal)*len(self.axs))
            self.fig.set_size_inches(width, height)
       
        self.fig.suptitle(self.figtitle, y=self.titley, fontname="Arial", fontsize=10)
        
#         self.fig.tight_layout()
        
        return

# adds a grid of dots at intersections to the viscosity maps
    # ax is th axis to add dots to
    # xlist is a list of x points
    # ylist is a list of y points
def addDots(ax, xlist, ylist):
    xl = []
    yl = []
    for x in xlist:
        for y in ylist:
            xl.append(x)
            yl.append(y)
    ax.scatter(xl, yl, color='#969696', s=10)
    return
    
# adjust the bounds of the plot
    # xlistreal is a list of x points to be included in the plot
    # xr is the [min, max] position of each segment, e.g. [-0.7, 0.7]
    # xlist is the initial list of x points we included in the plot
def adjustBounds(xlistreal, xr, xlist):
    if len(xlistreal)>1:
        xmin = min(xlistreal)
        xmax = max(xlistreal)
        pos1 = xlist.index(min(xlistreal))
        pos2 = xlist.index(max(xlistreal))+1
        dx = xr[1]-xr[0]
        xrtot = [xr[0]+pos1*dx, xr[0]+pos2*dx]
    else:
        xrtot = [0]
    return xrtot

def emptyYLabels(ax):
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_yticklabels(empty_string_labels)

    
   
 # plotCircle plots a circle, with optional caption
    # ax is axis to plot on
    # x0 is the x position
    # y0 is the y position
    # radius is the radius of the circle, in plot coords
    # caption is the label. use '' to have no caption
    # color is the color of the circle
def plotCircle(ax:plt.Axes, x0:float, y0:float, radius:float, caption:str, color, sigmapos:int=0) -> None:
    circle = plt.Circle([x0,y0], radius, color=color, fill=False)
    ax.add_artist(circle)
    angle=(90-30*sigmapos)*(2*np.pi/360)
    dar = 0.2 # length of arrow
    if len(caption)>0:
        if radius>0.3:
            txtx = x0
            txty = y0+0.2*sigmapos
            ax.text(txtx, txty, caption, horizontalalignment='center', verticalalignment='center', color=color)
        else:
            arrowx = x0+radius*np.cos(angle)
            arrowy = y0+radius*np.sin(angle)
            txtx = arrowx+dar*np.cos(angle)
            txty = arrowy+dar*np.sin(angle)
            ax.annotate(caption, (arrowx, arrowy), color=color, xytext=(txtx, txty), ha='center', arrowprops={'arrowstyle':'->', 'color':color})
        
        

    
    