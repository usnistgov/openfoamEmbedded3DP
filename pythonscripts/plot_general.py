#!/usr/bin/env python
'''Plotting tools for analyzing OpenFOAM single filaments'''

# external packages
import sys
import os
import numpy as np
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
import interfacemetrics as intm
import folderparser as fp

# plotting
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------


# formatting units

def kinToDyn(v:Union[List[float], float], density:float=1000) -> Union[List[float], float]:
    '''converts kinematic viscosity to dynamic viscosity
    v can be a single float or a list of floats. The default assumes a density of 1000 kg/m^3'''
    if type(v) is list:
        v2 = []
        for vi in v:
            v2.append(kinToDyn(vi))
        return v2
    if type(v) is str:
        try:
            v2 = float(v)
        except:
            return ''
    if type(density) is str:
        try:
            density = float(density)
        except:
            return ''
    v2 = v*density
    if abs(round(v2)-v2)<10**-20:
        v2 = int(v2)
    return v2

def expFormat(x:float) -> str:
    '''display a number in exponential format'''
    try:
        p0 = np.log10(x)
        p = round(p0)
        if abs(p-p0)<10**-8:
            return r'$10^{{{}}}$'.format(int(p))
        else:
            return x
    except:
        return x


def decideFormat(x:float) -> Union[float, str]:
    '''determine what kind of format to put the number in
    might return float or string'''
    if x==0:
        return x
    try:
        l = np.log10(x)
    except:
        return x
    if not l==round(l):
        return x
    else:
        return expFormat(x)


def expFormatList(xlist:List[float]) -> List[Any]:
    '''put the whole list in exponential or appropriate format'''
    xout = []
    for x in xlist:
        xout.append(decideFormat(x))
    return xout  

#-----------------------------------------
# x, y functions

def tpFunc(tp:Dict, stri:str) -> float:
    '''use this function to pluck a value out of the transport properties. Candidates include: nuink (ink viscosity), tau0ink (ink yield stress), kink (ink consistency index), nink (ink power law index), nusup (support viscosity), tau0sup, ksup, nsup, sigma (surface tension)}'''
    return tp[stri]
 
def multfunc(tp:Dict) -> float:
    '''ink viscosity * support viscosity. tp is a transportProperties dict'''
    return round(tp['nuink']*tp['nusup'], 10)

def divfunc(tp:Dict) -> float:
    '''support viscosity / ink viscosity. tp is a transportProperties dict'''
    return round(tp['nusup']/tp['nuink'], 10)

#-----------------------------------------
# color functions

def sigmaColor(sigma:int) -> str:
    '''preset colors for surface tensions in mJ/m^2'''
    if sigma==0:
        return '#940f0f'
    elif sigma==40:
        return '#61bab0'
    elif sigma==70:
        return '#1e5a85'

def sigfuncc(tp:Dict) -> str:
    '''color for this surface tension. tp is a transportProperties dict'''
    return sigmaColor(tp['sigma'])


def cubehelix1(val:float):
    '''val should be 0-1. returns a color'''
    cm = sns.cubehelix_palette(as_cmap=True, rot=-0.4)
    return cm(val)

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




#--------- 
# get transport properties

def extractTPfluid(le:pd.DataFrame, nui:int, getHB:bool) -> Tuple[float]:
    '''Extract transport properties from a legend. 
        le is the legend
        nui is the index of the viscosity
        getHB is true to get herschel-bulkley parameters
        used by extractTP'''
    nuj = nui
    while not le.loc[nuj,'title']=='rho' and nuj<len(le):
        nuj+=1
    rho = float(le.loc[nuj, 'val']) # kg/m^3
    nu = kinToDyn(float(le.loc[nui, 'val']), density=rho) # get the actual viscosity in Pa s
    if getHB:
        tau0 = kinToDyn(float(le.loc[nui+1, 'val']), density=rho)
        k = kinToDyn(float(le.loc[nui+2, 'val']), density=rho)
        n = float(le.loc[nui+3, 'val'])
    else:
        tau0 = 0
        k = 0
        n = 0
    return nu, tau0, k, n, rho

def extractTP(folder:str, units:bool=False) -> Tuple[dict,dict]:
    '''extractTP gets the metadata for a folder
    folder is a full path name
    used by vvplot, listTPvalues, folderToFunc'''
    try:
        le, u = fp.legendUnique(folder, units=True) # get dictionary of legend values
    except:
        traceback.print_exc()
        print(folder)
        raise Exception('No values in legend') 
    if len(le)==0:
        raise Exception('No values in legend')
    if not 'ink_velocity' in le:
        fs.addUnitsToLegend(folder)
        le, u = fp.legendUnique(folder, units=True) # get dictionary of legend values
    for key in le:
        try:
            le[key] = float(le[key])
        except:
            pass
    le['folder'] = folder
        
    if u['ink_velocity'] == 'm/s':
        le['vink'] = 1000*le['ink_velocity']
        u['vink'] = 'mm/s'
        le['vsup'] = 1000*le['bath_velocity']
        u['vsup'] = 'mm/s'
    if u['ink_rho'] == 'kg/m^3':
        for fluid in ['ink', 'sup']:
            le[f'rho{fluid}'] = le[f'{fluid}_rho']/1000
            u[f'rho{fluid}'] = 'g/mL'
    for s in ['nu', 'nu0', 'tau0', 'k', 'n']:
        for fluid in ['ink', 'sup']:
            for dic in [le, u]:
                dic[f'{s}{fluid}'] = dic[f'{fluid}_{s}']
            # convert kinematic to dynamic
            if s in ['tau0', 'k', 'nu', 'nu0'] and type(le[f'{s}{fluid}']) is float:
                try:
                    le[f'{s}{fluid}'] = kinToDyn(le[f'{s}{fluid}'], le[f'{fluid}_rho'])
                except:
                    traceback.print_exc()
                    pass
                else:
                    u[f'{s}{fluid}'] = {'m^2/s':'Pa*s', 'm^2/s^2':'Pa', 'm^2*s^(n-2)':'Pa*s^n'}[u[f'{s}{fluid}']]
    
    if u['sigma'] =='J/m^2':
        le['sigma'] = 1000*float(le['sigma'])
        u['sigma'] = 'mJ/m^2'
    le['nozzleAngle'] = float(le['nozzle_angle'])
    le['dink'] = float(le['nozzle_inner_width'])
    le['dsup'] = float(le['nozzle_inner_width']) + 2*float(le['nozzle_thickness'])
    if units:
        return le,u
    else:
        return le


#---

def listTPvalues(flist, **kwargs) -> Tuple[List[str], Dict]:
    '''get a list of all of the transport property values in this list of folders
    flist is a list of folders (full path names)
    Returns a list of files and Dict of properties. 
        If we use kwargs to say we only want files with, e.g. nuink=10, then it will only return those files.
        The dictionary lists all of the values in the list of files for each transport property variable. e.g. {'nuinklist':[10,100], 'tau0inklist':[0], 'kinklist':[0], ...}
    '''
    
    tab = [extractTP(folder) for folder in flist]
    tab = pd.DataFrame(tab)
    
    for k in tab.keys():
        if k+'list' in kwargs:
            tab = tab[tab[k].isin(kwargs[k+'list'])]
    flist2 = list(tab['folder'])
    lists = dict([[k+'list', sorted(list(tab[k].unique()))] for k in tab.keys()[1:]])
    return flist2, lists


#---------------------------------------
# using transport properties to decide how to plot

def folderToFunc(folder:str, func) -> float:
    '''given a folder, get one value to represent that folder
    func is the function to apply to transport properties to get one value
    used by unqListFolders'''
    tp = extractTP(folder)
    return round(func(tp), 15)

def tpCombos(tplists:Dict) -> List:
    '''List of combinations of transportProperties values. tplists comes from listTPvalues
    used by unqList, and will look like {'nuinklist':[10, 100], 'tau0inklist':[10], ...}'''
    vallists = []
    for l in tplists:
        vallist = [[l[0:-4],i] for i in tplists[l]] # e.g. [['nuink', 10^5], ['nuink', 10^6]]
        vallists.append(vallist)
    l0 = list(itertools.product(*vallists))
    return l0



def unqList(f, tplists:Dict) -> List[float]:
    '''unqList finds the unique x or y positions in the plot
    f is the function to apply to transport properties
    tplists is a list of transport properties
    used by folderplot'''
    lout = []
    combos = tpCombos(tplists)
    for tp in combos:
        l = f(tp)
        if l not in lout:
            lout.append(l)
    lout.sort()
    return lout


def unqListFolders(folders:List[str], func) -> List[float]:
    '''given a folder of many simulations, find unique x or y positions in comboPlot or gridOfPlots'''
    funcvals = [] # outputs for the func
    for f in folders:
        val = folderToFunc(f, func)
        if val not in funcvals:
            funcvals.append(val)
    funcvals.sort()
    return funcvals


#--------------------------------

def findPos(l:List, v:Any) -> Any:
    '''find the position of v in list l. l is a list. v is a value in the list.
    used by vv'''
    try:
        p = l.index(v)
    except ValueError:
        return -1
    return p

#---

def vv(tp:Dict, xpv) -> Tuple[Any, float, float, int]:
    '''find values needed for plotting in grids
    tp holds transport properties for a single simulation
    xpv is a comboPlot or gridOfPlots object
    Returns the color (could be a string or other), x plot position, y plot position, and position of sigma in list of sigmas
    used by vvplot'''
    x = xpv.xfunc(tp)
    y = xpv.yfunc(tp)
    xpos = findPos(xpv.xlist, x)
    ypos = findPos(xpv.ylist, y)
    sigmapos = findPos(xpv.tplists['sigmalist'], tp['sigma'])
    
    # find the position in the plot x0,y0 for this simulation. Not the real value! Just a placeholder so we don't have to deal with logs.
    if xpos<0 or ypos<0 or sigmapos<0:
        raise ValueError
    if xpv.type=="comboPlot":
        x0 = xpv.xmlist[xpos]
        y0 = xpv.ymlist[ypos]
    else:
        x0 = xpos
        y0 = len(xpv.ylist)-ypos-1  # we need to flip y upside down for a grid of plots
        
    # expand the list of real x and y values to include this one
    if x not in xpv.xlistreal:
        xpv.xlistreal.append(x)
    if y not in xpv.ylistreal:
        xpv.ylistreal.append(y)
    color = xpv.cfunc(tp)
    return color, x0, y0, sigmapos

def vvplot(folder:str, xpv):
    '''vvplot finds variables used for plotting in grids
    folder is a folder name
    xpv is a comboPlot or gridOfPlots object'''
    tp = extractTP(folder)
    return vv(tp, xpv)


#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
# PLOT CLASSES
#-------------------------------------------------
#-------------------------------------------------

class folderPlots:
    '''A generic class used for plotting many folders at once. Subclasses are comboPlot, which puts everything on one plot, and gridOfPlots, which puts everything in separate plots based on viscosity.'''
    
    def __init__(self, topFolder:str, imsize:float, split:bool=False, **kwargs):
        '''topFolder is the folder we're plotting
            imsize is the size of the total image in inches
            split is true to split into separate plots by surface tension'''
        self.kwargs = kwargs
        self.ab = not 'adjustBounds' in self.kwargs or self.kwargs['adjustBounds']==True
        self.topFolder = topFolder
        self.imsize = imsize
        self.split = split
        self.plotsLists(**kwargs)
        
    def convertFunc(self, var):
        '''Convert a variable name or expression, e.g. 'nuink' or 'nusup/nuink' into a lambda function to be used on transport properties dict'''
        strs = self.strings.copy()
        strs.remove('sup')
        strs.remove('ink')
        if len(var)>0:
            if var in strs:
                # just one variable, can just pick the variable from transport properties
                func = lambda tp: tpFunc(tp, var)
                return func
            else:
                for s in strs:
                    var = var.replace(s, "tp[\'"+s+"\']")
                func = lambda tp: eval(var)
                return func
       
        
    def plotsLists(self, **kwargs):
        '''plotsLists initializes gridOfPlots and comboPlots objects, creating the initial figure'''
        self.flist = fp.caseFolders(self.topFolder) # list of all folders in the top folder
        self.flist, self.tplists = listTPvalues(self.flist, **kwargs) # list of transport property lists
        self.cfunc = sigfuncc
        self.strings = [s[0:-4] for s in list(self.tplists.keys())]
        
        if self.split:
            ncol = len(self.tplists['sigmalist'])
        else:
            ncol = 1
        self.ncol = ncol
        
        # get the variables or expressions we want to operate on
        if 'xvar' in kwargs and 'yvar' in kwargs:
            self.xvar = kwargs['xvar']
            self.yvar = kwargs['yvar']
        elif 'mode' in kwargs:
            # older versions of the code used preset plotting modes. Defining 'xvar' and 'yvar' directly allows for more flexibility.
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
        
        # convert those variables into functions that we can use on transport properties dicts
        self.xfunc = self.convertFunc(self.xvar)
        self.yfunc = self.convertFunc(self.yvar)
        try:
            # find lists of unique x values and y values
            self.xlist = unqListFolders(self.flist, self.xfunc)
            self.ylist = unqListFolders(self.flist, self.yfunc)
            self.xlistreal = []
            self.ylistreal = []
            self.legendList()
        except Exception as e:
            print(e)
            raise ValueError('Failed to identify x and y variables')
        return self
    
    
    def legendList(self):
        '''Make a legend from the list of sigma values and store it for later'''
        if self.ncol==1:
            sigmalist = self.tplists['sigmalist']
            plist = [mpatches.Patch(color=sigmaColor(sigmalist[i]), label=sigmalist[i]) for i in range(len(sigmalist))]
            ph = [plt.plot([],marker="", ls="", label='\u03C3 (mJ/m$^2$)')[0]]; # Canvas
            self.plist = ph + plist
            plt.close()
        return 
    
    def getLabel(self, var, short):
        '''Get label for the x or y axis'''
        # determine which axis we're trying to name
        if var=='x':
            func = self.xfunc
        else:
            func = self.yfunc
        #RG 
        namedefs = {'nuink':'Ink viscosity (Pa$\cdot$s)', 'tau0ink':'Ink yield stress (Pa)','kink':'Ink k (Pa*s^n)', 'nink':'Ink n',\
                    'nusup':'Support viscosity (Pa$\cdot$s)', 'tau0sup':'Support yield stress (Pa)', 'ksup':'Support k (Pa*s^n)', 'nsup':'Support n',\
                    'sigma':'Surface tension (mJ/m$^2$)', 'product':'Support viscosity \u00d7 ink viscosity (Pa$^2\cdot$s$^2$)', 'ratio':'Support viscosity / ink viscosity', 'nozzleAngle':'Nozzle angle (degrees)'} 
        
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
        fig, axs = plt.subplots(nrows=len(self.ylist), ncols=len(self.xlist),\
                                sharex='col', sharey='row',\
                                figsize=(imsize*len(self.xlist), imsize*len(self.ylist)))
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
                 
        # store axes and figure in object
        self.axs = axs
        self.fig = fig
    
    
    def clean(self):
        '''post-processes the plot to add components after all plots are added'''
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
            if np.ndim(self.axs)!=1: # RG
                strval = str(expFormat(xval)) # title
                self.axs[0][j].set_title(strval, fontsize=fs) # put the title on top
            else:
                self.axs[j].set_title(xval, fontsize=fs) # put the title on top for 1D plot
        for i, yval in enumerate(self.ylistreal):
            # going down rows
            strval = str(expFormat(yval)) # title
            if np.ndim(self.axs)!=1: # RG
                self.axs[i][-1].text(5.1, 4, strval, verticalalignment='center', rotation=270, fontsize=fs) # put the title on the right side
                
        # reset figure size
     #   self.fig.set_size_inches(self.imsize*len(self.xlistreal), self.imsize*len(self.ylistreal))
        
        # top level axis labels
        self.fig.suptitle(self.getLabel('x', False), y=1.25-(firsty/len(self.ylist)), fontsize=fs) # was 0.92 instead of 1.25 RG
        if np.ndim(self.axs)!=1: # RG
            self.fig.text((lastx/len(self.xlist))-0.2, 0.5, self.getLabel('y', False),\
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

class comboPlot(folderPlots):
    '''stores variables needed to create big combined plots across several folders '''
    
    def __init__(self, topFolder:str, xr:List[float], yr:List[float], imsize:float, gridlines:bool=True, **kwargs):
        '''topFolder is a full path name. topFolder contains all of the folders we want to plot
        xr is the min and max x value for each section of the plot, e.g. [-0.7, 0.7]
        yr is the min and max y value for each section of the plot
        imsize is the size of the whole image in inches
        gridlines true to show gridlines'''
        
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
        
        self.imsize = imsize
        fig, axs = plt.subplots(nrows=1, ncols=self.ncol, figsize=(imsize,imsize*len(self.ylist)/len(self.xlist)), sharey=True)
        fig.subplots_adjust(wspace=0)
        
        if self.ncol==1:
            axs = [axs]
        
        # vert/horizontal grid
        if gridlines:
            for ax in axs:
                ax.grid(linestyle='-', linewidth='0.25', color='#949494')


        # store variables
        self.axs = axs
        self.fig = fig
        
        self.addLegend()
        
    def addLegend(self):
        '''Add a surface tension legend'''
        if not self.ncol==1 and not self.xvar=='nozzleAngle': # surface tension legend was not needed for conical nozzles bc only sigma=0 was used RG
            self.axs[0].legend(handles=self.plist, loc='upper center', ncol=4, bbox_to_anchor=(0.5, self.titley+0.1))
        
    def clean(self):
        '''post-processes the plot to add components after all plots are added '''
        
        # adjust the bounds of the plot to only include plotted data
        # each time we added a folder to the plot, we added the 
        # x and y values of the centers to xlistreal and ylistreal
        # this is particularly useful if a big section of the plot 
        # is unplottable. e.g., very low support viscosity/ink viscosity
        # values produce filaments which curl up on the nozzle, so they don't produce cross-sections.
        # This step cuts out the space we set out for those folders that didn't end up in the final plot
        # if we were given adjustBounds=False during initialization, don't adjust the bounds
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
            
        if self.ncol>1:
            sigmalist = self.tplists['sigmalist']
            for i in range(len(sigmalist)):
                self.axs[i].set_title('\u03C3='+str(sigmalist[i])+' mJ/m$^2$', fontname="Arial", fontsize=10)

        self.axs[0].set_ylabel(self.getLabel('y', False), fontname="Arial", fontsize=10)
        
        yticklabels = expFormatList(self.ylist)
        self.axs[0].set_yticklabels(yticklabels, fontname="Arial", fontsize=10)
        
        # reset the figure size so the title is in the right place
        if len(self.xlistreal)>0 and len(self.ylistreal)>0:
            width = self.imsize
            height = width*len(self.ylistreal)/(len(self.xlistreal)*len(self.axs))
            self.fig.set_size_inches(width, height)
       
        # set position of titley
        ylim = self.axs[0].get_ylim()
        if ylim[1]-ylim[0]>2.5:
            self.titley = 1.1
        elif ylim[1]-ylim[0]>1.3:
            self.titley = 0.9
        else:
            self.titley = 0.8
        self.fig.suptitle(self.figtitle, y=self.titley, fontname="Arial", fontsize=10)

        if self.xvar=='nozzleAngle': # since there is only one variable, remove y information RG
            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_xticklabels(self.xlist, fontname="Arial", fontsize=10) # 10 degree angle should not be 10^1

#         self.fig.tight_layout()
        
        return
    
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


def adjustBounds(xlistreal:List[float], xr:List[float], xlist:List[float]):
    '''adjust the bounds of the plot.
    xlistreal is a list of x points to be included in the plot
    xr is the [min, max] position of each segment, e.g. [-0.7, 0.7]
    xlist is the initial list of x points we included in the plot'''
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

def emptyYLabels(ax:plt.Axes):
    '''Leave the y tick labels empty. Useful for side-by-side plots.'''
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_yticklabels(empty_string_labels)

def plotCircle(ax:plt.Axes, x0:float, y0:float, radius:float, caption:str, color, sigmapos:int=0) -> None:
    '''plotCircle plots a circle, with optional caption.
    ax is axis to plot on
    x0 is the x position
    y0 is the y position
    radius is the radius of the circle, in plot coords
    caption is the label. use '' to have no caption
    color is the color of the circle
    sigmapos is the position in the sigma list for this circle. Useful for timeplots, so labels don't stack on top of each other. If we're using this to plot an ideal cross-section or some other sigma-unaffiliated value, sigmapos=0 will put the label inside the circle or just above it.'''
    
    circle = plt.Circle([x0,y0], radius, color=color, fill=False) # create the circle
    ax.add_artist(circle)                                         # put the circle on the plot

    if len(caption)>0:
        if radius>0.3:
            # if the circle is big, put the label inside
            txtx = x0                         
            txty = y0+0.2*sigmapos
            ax.text(txtx, txty, caption, horizontalalignment='center', verticalalignment='center', color=color)
        else:
            # if the circle is small, put the label outside
            angle=(90-30*sigmapos)*(2*np.pi/360)  # angle to put the label at, in rad                            
            dar = 0.2                             # length of arrow
            arrowx = x0+radius*np.cos(angle)      # where the arrow points   
            arrowy = y0+radius*np.sin(angle)
            txtx = arrowx+dar*np.cos(angle)       # where the label is
            txty = arrowy+dar*np.sin(angle)
            ax.annotate(caption, (arrowx, arrowy), color=color, xytext=(txtx, txty), ha='center', arrowprops={'arrowstyle':'->', 'color':color})
            
            
#--------------------------------------
