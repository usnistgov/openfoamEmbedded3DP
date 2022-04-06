#!/usr/bin/env python
'''Functions for plotting survival metrics and stress metrics...'''

# external packages
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import math
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import folderparser as fp
import interfacemetrics as intm
from plot_general import *
from plainIm import *
from figureLabels import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='8.0')

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "NIST"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------
#--------------------------------

        
def unitName(zunits:str) -> str:
    '''convert the units name to something more readable and compact'''
    if zunits=='nozzle_inner_width':
        return '$d_i$'
    elif zunits=='vsup':
        return 'Translation speed'
    else:
        return zunits.replace('_', ' ')
    
    
def weights(a1:float=10**0, b1:float=0, c1:float=1, 
            a2:float=10**-2, b2:float=0.5, c2:float=0.5, 
            a3:float=10**-3, b3:float=1, c3:float=0, **kwargs):
    '''standard weights for the survival function'''
    w1 = {'a':a1, 'b':b1, 'c':c1}
    w2 = {'a':a2, 'b':b2, 'c':c2}
    w3 = {'a':a3, 'b':b3, 'c':c3}
    return [w1, w2, w3]

def survivalEqLabel(a,b,c):
    '''get a label for the survival equation'''
    a1 = '{-'+str(a)+'}'
    b1 = '{'+str(b)+'}'
    c1 = '{'+str(c)+'}'
    return '$S=exp('+a1+'\\tau^'+b1+'t^'+c1+')$'

class survival:
    '''holds info about survival over the length of the nozzle'''
    
    def __init__(self, rbar:float, a:float=10**-4, b:float=0.5, c:float=0.5):
        '''rbar is the normalized radius
        a,b,c are model parameters'''
        self.rbar = rbar
        self.zlist = [] # this must be in mm
        self.xlist = [] # this must be in mm
        self.totalI = 1 # current survival fraction
        self.Ilist = [1]
        self.tlist = [0] # this is in s
        self.dtlist = [0]
        self.taulist = [0]
        self.vlist = [0]
        self.a = a
        self.b = b
        self.c = c
        
    def addStep(self, row:pd.Series) -> None:
        '''add the point to the survival'''
        if len(self.zlist)==0:
            # start at 100% survival
            self.zlist.append(row['z'])
            self.xlist.append(row['x'])
            return
        else:
            v = row['magu'] # velocity magnitude
            traveled = np.sqrt((self.zlist[-1] - row['z'])**2 + (self.xlist[-1]-row['x'])**2)
            dt = traveled/v # difference in time since last step
            tau = row['shearstressmag'] 
            if v<0.1*self.vlist[-1]:
                # velocity discontinuity, abort
                return
            
            Ii = np.exp(-self.a*tau**self.b*dt**self.c) # survival during this step
            Ii = min(1, Ii)
            self.totalI = self.totalI*Ii # survival overall
            self.Ilist.append(self.totalI) # keep track of survival at each step
            self.tlist.append(self.tlist[-1]+dt)
            self.taulist.append(tau)
            self.dtlist.append(dt)
            self.xlist.append(row['x'])
            self.zlist.append(row['z'])
            self.vlist.append(v)
                
    def addUnits(self, normval:float) -> None:
        '''divide the units of the z variable by the given value'''
        self.zlist = [z/normval for z in self.zlist]

        
def survivalCalc(folder:str, time:float=2.5, a:float=10**-2, b:float=0.5, c:float=0.5, zunits:str='mm', dr:float=0.05, fcrit:float=0.9, volume:bool=True, **kwargs):
    '''calculate what cells will survive the process, if S=exp(-a*tau^b*t^c). 
    time is the time at which to collect stress data
    zunits is a string, usually mm or nozzle_inner_width
    dr is the spacing between relative radial positions to group by, as a fraction
    fcrit is the min fraction of points required to get returned. Otherwise too many points were skipped, and the measurement is not valid
    volume = True to use the whole volume of the nozzle, if there is a file. otherwise, use a slice collected from the center of the nozzle
    '''
    
    if volume:
        df,units = intm.importPtsNoz(folder, time) # get points in nozzle
    else:
        df,units = intm.importSliceNoz(folder, time) # get points in nozzle
    if len(df)==0:
        return [] 
    
    df = intm.takePlane(df, folder, dr=dr, xhalf=True) # take just the middle plane

    if 'rbarlist' in kwargs:
        rbarlist = kwargs['rbarlist']
    else:
        rbarlist = np.arange(0, 1, dr)                            # normalized radius evenly spaced from 0 to 1
        rbarlist = [round(rbar,5) for rbar in rbarlist]           # round to avoid floating point error
    rz = dict([[rbar,survival(rbar, a, b, c)] for rbar in rbarlist]) # table of survival as a function of r/r0 and z
    zlist = list(df.z.unique())
    zlist.sort()                  # put in order, where negative values are at top
    for z in zlist:
        df0 = df[df.z==z]         # select points in this plane
        for rbar in rbarlist:
            df1 = df0[(df0['rbar']==rbar)&(df0['magu']>0)]
            if len(df1)>0:
                row = df1.mean() # average all points
                rz[rbar].addStep(row)
                
    # remove series that are too short
    remlist = []
    for key in rz.keys():
#         if len(rz[key].zlist)<0.9*len(zlist):
        if len(rz[key].zlist)<0.75*len(zlist):
            remlist.append(key)
    for key in remlist:
        rz.pop(key)
            
    if not zunits=='mm':
        le = fp.legendUnique(folder)
        if zunits in le:
            normval = float(le[zunits])
            for key in rz:
                rz[key].addUnits(normval)
            
    return rz

def survivalRateRZ(rz:dict, dr) -> float:
    '''get the survival rate, given a dictionary rz that holds survival objects
    dr is the spacing between relative radial positions to group by, as a fraction
    '''
    weightedsum = 0
    weight = 0
    for rbar in rz:
        area = (rbar**2 - (rbar-dr)**2)
        weightedsum = weightedsum + rz[rbar].totalI*area # survival * area of ring
        weight = weight + area
    if weight>0:
        return weightedsum/weight
    else:
        return 0
    
    
def survivalRate(folder:str, time:float=2.5, a:float=10**-3, b:float=0.5, c:float=0.5, dr:float=0.05, **kwargs) -> float:
    '''get the percentage of surviving cells at the end of the nozzle
    time is the time at which to collect stress data
    a,b,c are survival model parameters
    dr is the spacing between relative radial positions to group by, as a fraction
    '''
    rz = survivalCalc(folder, time=time, a=a, b=b, c=c, dr=dr, **kwargs)
    return survivalRateRZ(rz, dr)

def survivalzPlot(rz:dict, xvar:str, axs:np.array, zunits:str, cm):
    '''plot the survival metrics as a function of z or t position. 
    rz is a dctionary of survival objects created by survivalCalc. 
    xvar is a string, z or t
    axs is a list of matplotlib axes
    zunits is a string, usually mm or nozzle_inner_width
    cm is a colormap'''
    for rbar in rz:
        if xvar=='z':
            xlist = rz[rbar].zlist
        else:
            xlist = rz[rbar].tlist
        xlist = xlist[1:]
        if len(xlist)>0: 
            for i,yl in enumerate(['vlist', 'taulist', 'Ilist']):
                ylist = getattr(rz[rbar], yl)[1:]
                axs[i].plot(xlist, ylist, c=cm(rbar), linewidth=0.75)
                xi = int(len(xlist)/2)
                axs[i].text(xlist[xi], ylist[xi], rbar, color=cm(rbar), horizontalalignment='center', verticalalignment='top') 

                if xvar=='z':
                    zunname = unitName(zunits)
                    axs[i].set_xlabel(f'$z$ position ({zunname})')
                elif xvar=='t':
                    axs[i].set_xlabel('Time (s)')
                if yl=='Ilist':
                    axs[i].set_ylabel('Surviving cells/initial cells')
                elif yl=='taulist':
                    axs[i].set_ylabel('Shear stress (Pa)')
                elif yl=='dtlist':
                    axs[i].set_ylabel('Time at step (s)')
                elif yl=='vlist':
                    axs[i].set_ylabel('Velocity (mm/s)')
        
    for ax in [axs[0], axs[1]]:
        ax.set_ylim(bottom=0)
        
    for ax in axs:
        setSquare(ax)
        
    
        
def survivalrPlot(folder:str, ax, time:float=2.5, a:float=10**-3, b:float=0.5, c:float=0.5, dr:float=0.05, xlabel:bool=True, ylabel:bool=True, fontsize:int=8, **kwargs):
    '''plot cell survival as a function of normalized radius within the nozzle
    ax is the axis to plot on
    time is the time at which to collect stress data
    a,b,c are survival model parameters
    dr is the spacing between relative radial positions, as a fraction
    xlabel = True to label the x axis
    ylabel = True to label the y axis
    '''
    xlist = []
    ylist = []
    rz = survivalCalc(folder, time=time, a=a, b=b, c=c, dr=dr, **kwargs)
    for rbar in rz:
        xlist.append(rbar)
        ylist.append(rz[rbar].totalI)
    if not 'color' in kwargs:
        color='black'
    else:
        color = kwargs['color']
    ax.plot(xlist, ylist, color=color, linewidth=0.75)
    if 'label' in kwargs and len(xlist)>2:
        xi = int(len(xlist)/2)
        ax.text(xlist[xi], ylist[xi], kwargs['label'], color=color, horizontalalignment='center', verticalalignment='top', fontsize=fontsize) 
    if xlabel:
        ax.set_xlabel('Radius/nozzle radius', fontsize=fontsize)
    if ylabel:
        ax.set_ylabel('Surviving cells/initial cells', fontsize=fontsize)
    setSquare(ax)
    return rz


    
def survivalRMulti(topFolder:str, axs, cvar:str='nozzle_angle', time:float=2.5, a:float=10**-3, b:float=0.5, c:float=0.5, dr:float=0.05, xlabel:bool=True, ylabel:bool=True, fontsize:int=8, **kwargs):
    '''plot cell survival as a function of normalized radius within the nozzle, for multiple sims. axs should be an array of 2 axes
    topFolder holds multiple simulations
    axs is the list of axes to plot on
    cvar is the variable to color by
    time is the time at which to collect stress data
    a,b,c are survival model parameters
    dr is the spacing between relative radial positions, as a fraction
    xlabel = True to label the x axis
    ylabel = True to label the y axis
    '''
    
    plt.rc('font', size=fontsize)
    
    folders = fp.caseFolders(topFolder)       # all sims in folder
    folders, _ = listTPvalues(folders, **kwargs) # remove any values that don't match
    
    cm = sns.color_palette('viridis', n_colors=len(folders)) # uses viridis color scheme
    
    flist = []
    for i,folder in enumerate(folders):
        le, u = extractTP(folder, units=True)
        nz = round(float(le[cvar]), 10)
        flist.append({'folder':folder, 'cvar':nz})
        
    if len(flist)==0:
        return
        
    flist = pd.DataFrame(flist)
    flist.sort_values(by='cvar', inplace=True)
    flist.reset_index(drop=True, inplace=True)
    flist['cvar'] = expFormatList(list(flist['cvar']))
    
    
    for i,row in flist.iterrows():
        nz = row['cvar']
        if cvar=='nozzle_angle':
            nz = f'{int(nz)}$\degree$'
        rz = survivalrPlot(row['folder'], axs[0], a=a, dr=dr, b=b, c=c, color=cm[i], label=nz, ylabel=ylabel, xlabel=xlabel, **kwargs)
        rate = survivalRateRZ(rz, dr)
        axs[1].scatter([nz], [rate], color=cm[i])
       
    if xlabel:
        axs[1].set_xlabel(f'{unitName(cvar)} ({u[cvar]})', fontsize=fontsize)
    if ylabel:
        axs[1].set_ylabel('Surviving cells/initial cells', fontsize=fontsize)
    setSquare(axs[1])

    axs[0].set_title(survivalEqLabel(a,b,c), fontsize=fontsize)
    
    

    
def survivalRMultiRow(topFolder:str, exportFolder:str, fontsize:int=8, export:bool=True, overwrite:bool=False, **kwargs):
    '''plot cell survival as a function of radius, at three weights of the equation
    topFolder holds multiple simulations
    exportFolder is the folder to export figures to
    export=True to export images to file
    overwrite=True to overwrite existing files
    '''
    
    labels = ['survivalMulti']
    fn = intm.imFn(exportFolder, labels, topFolder, **kwargs) # output file name
    if not overwrite and os.path.exists(f'{fn}.png'):
        return
    
    fig,axs = plt.subplots(2,3, figsize=(6.5, 4.5))
    if 'logx' in kwargs and kwargs['logx']:
        for ax in axs[1]:
            ax.set_xscale('log')
    plt.rc('font', size=fontsize) 
    for i,d in enumerate(weights(**kwargs)):
        survivalRMulti(topFolder, [axs[0][i], axs[1][i]], a=d['a'], b=d['b'], c=d['c'], **kwargs)
    subFigureLabels(axs, inside=False)
    fig.tight_layout()
    
    if export:
        intm.exportIm(fn, fig) # export figure


def survivalPlot(folder:str, exportFolder:str, xvar:str, time:float=2.5, a:float=10**-3, b:float=0.5, c:float=0.5
                 , zunits:str='nozzle_inner_width', fontsize:int=8, dr:float=0.05, export:bool=True, overwrite:bool=False, **kwargs):
    '''plot survival as a function of z position, relative radius, or get a single value.
    folder is the simulation folder
    exportFolder is the folder to export results to
    xvar is the variable on the x axis. value of xvar should be 'z', 't', 'rbar', or 'scalar' 
    time is the time at which we evaluate survival
    a,b, and c are model parameters for the cell survival
    zunits = mm or nozzle_inner_width
    dr is the spacing between relative r values, as a fraction
    export True to export images
    overwrite True to export images even if the file already exists    
    '''

    labels = ['survival', xvar, os.path.basename(folder), str(a), str(b), str(c), zunits]
    fn = intm.imFn(exportFolder, labels, os.path.dirname(folder), **kwargs) # output file name
    if not overwrite and os.path.exists(f'{fn}.png'):
        return
    
    plt.rc('font', size=fontsize) 
    if not xvar=='scalar':
        
        if not 'cm' in kwargs:
            cm = sns.color_palette('viridis', as_cmap=True) # uses viridis color scheme
        else:
            cm = kwargs['cm']
        
    if xvar=='z' or xvar=='t':
        if 'figsize' in kwargs:
            fig,axs = plt.subplots(1,3, figsize=kwargs['figsize'])
        else:
            fig,axs = plt.subplots(1,3, figsize=(6.5, 4))
        rz = survivalCalc(folder, time=time, a=a, b=b, c=c, zunits=zunits, dr=dr, **kwargs)
        survivalzPlot(rz, xvar, axs, zunits, cm)
        axs[2].set_title(survivalEqLabel(a,b,c), fontsize=fontsize)
        fig.tight_layout()
            
    elif xvar=='rbar':
        fig,ax = plt.subplots(1,1)
        survivalrPlot(folder, ax, time=time, a=a, b=b, c=c, zunits=zunits, dr=dr, **kwargs)
    elif xvar=='scalar':
        return survivalRate(folder, time=time, a=a, b=b, c=c, dr=dr, **kwargs)
    
    
    if export:
        intm.exportIm(fn, fig) # export figure
    
    
    
 #----------------------------------  
    
def stressInPlane(df0:pd.DataFrame) -> float:
    '''get the average shear stress in the plane, where df0 holds points in the plane'''
    errorRet = -1
    if len(df0)==0:
        return errorRet
    
    if not 'shearstressmag' in df0:
        raise ValueError('No shearstressmag in table')
    
    weightedsum = 0
    weight = 0

    df0 = intm.removeOutliers(df0, 'shearstressmag', sigma=6) # remove outliers
    
    rbarlist = list(df0.rbar.unique())
    if len(rbarlist)<2:
        return errorRet # not enough points in ring. return
    
    rbarlist.sort()
    dr = rbarlist[1]-rbarlist[0]
    for rbar in rbarlist:
        df1 = df0[df0.rbar==rbar]              # only points within a ring
        stress = df1['shearstressmag'].mean()  # average stress in the ring
        
        area = abs((rbar+dr/2)**2 - (rbar-dr/2)**2)        # area of ring
        weightedsum = weightedsum + stress*area # stress * area of ring
        weight = weight + area                 # cumulative area
    stress = weightedsum/weight
    return stress


def shearStressCalcVolume(folder:str, time:float, zunits:str, z0:float) -> Tuple[List[float], List[float]]: # RG
    '''calculate mean shear stress across the length of the nozzle, from all points in nozzle
    time is the time at which we collect the shear stress
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    z0 is the bottom z position of the nozzle, or any value in mm to set the z position relative to
    '''
    
    df,units = intm.importPtsNoz(folder, time) # get points in nozzle
    if len(df)==0:
        return [], []
    
    vals = df.groupby(by='z').mean()['shearstressmag']
    
    if len(vals)==0:
        return [],[]
    
    xlist = list(z0 - vals.index)
    ylist = list(vals)
    if not zunits=='mm':
        le = fp.legendUnique(folder)
        xlist = [i/float(le[zunits]) for i in xlist]

    return xlist, ylist


def shearStressCalcSlice(folder:str, time:float, zunits:str) -> Tuple[List[float], List[float]]: 
    '''Calculates mean shear stress across the length of the nozzle, from slice
    folder is the folder to do calculations on
    time is the time at which to calculate the shear stress
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    '''
    
    df,units = intm.importSliceNoz(folder, time) # get points in nozzle
    if len(df)==0:
        return [], []
    if not 'shearstressmag' in df:
        logging.warning(f'No shearstressmag in {folder}')
        return [],[]
    df = intm.takePlane(df, folder, dr=0.001) # add rbar column
    df['z'] = [round(z,5) for z in df['z']]
    
    # iterate through z slices and get average stress
    zlist0 = list(df.z.unique())
    zlist0.sort()
    zlist = []
    stresslist = []
    for zi in zlist0:
        df0 = df[df.z==zi]
        if len(df0)>0:
            sip = stressInPlane(df0)
            if sip>0 and (len(stresslist)==0 or sip<10*stresslist[-1]): # filter out spikes
                stresslist.append(sip)
                zlist.append(zi)
    
        
    # normalize z values by units
    if not zunits=='mm':
        le = fp.legendUnique(folder)
        zlist = [i/float(le[zunits]) for i in zlist]

    return zlist, stresslist
#-------------------
    
def nozzleLineTrace(folder:str, time:float, zabove:float, zunits:str='mm', volume:bool=False) -> pd.DataFrame:
    '''Calculates mean shear stress across the width of the nozzle
    folder is the folder to do calculations on
    time is the time at which to collect the stress
    zabove is the z position relative to the bottom of the nozzle, in units of zunits
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    volume = True to use the whole volume of the nozzle, if there is a file. otherwise, use a slice collected from the center of the nozzle
    '''
    
    if volume:
        df,units = intm.importPtsNoz(folder, time) # get points in nozzle
    else:
        df,units = intm.importSliceNoz(folder, time) # get points in nozzle
    
    if len(df)==0:
        return []
    
    le = fp.legendUnique(folder)
    
    df = intm.takePlane(df, folder, dr=0.001, xhalf=True) # add rbar column

    if not zunits==units['z']:
        # convert z units
        if zunits in le:
            con = float(le[zunits])
            df['z'] = [i/con for i in df['z']]
            
    if zabove>0:
        zabove=-zabove
    
    zvar = intm.closest(df['z'].unique(), zabove) # get exact z val
    df = df[(df['z']==zvar)] # get points at that z val

    # convert kinematic viscosity to dynamic
    rho = float(le['ink_rho'])    
    if 'nu_ink' in df:
        df['nu_ink'] = [rho*nu for nu in df['nu_ink']]

    md = float(le['nozzle_center_x_coord'])
    df['x'] = [x-md for x in df['x']]
    
    df.sort_values(by='x', inplace=True)
        
    return df


def getListsFromTrace(row:pd.Series, zstress:pd.DataFrame, xvar:str, yvar:str, cvar:str):
    '''from a dataframe tracing values, get x and y points to plot
    row is a row from a dataframe
    zstress is the trace holding the mean shear stress across the width of the nozzle
    xvar is the variable to plot on the x axis
    yvar is the variable to plot on the y axis
    cvar is the variable to color by    
    '''
    if len(zstress)>0:
        theta = row[cvar]    
        xlist = list(zstress[xvar])
        ylist = list(zstress[yvar])
        if yvar=='magu' or xvar=='rbar':
            li = int(len(zstress)/2)
        else:
            li = 0
        return theta, li, xlist, ylist
    else:
        return '', 0, [], []
    

def getDataWithinNozzle(xvar:str, yvar:str, cvar:str, row:pd.Series, time:float, zabove:float, zunits:str, volume:bool):
    '''get point data within the nozzle for the intended metrics
    xvar is the variable to plot on the x axis
    yvar is the variable to plot on the y axis
    cvar is the variable to color by   
    row holds metadata about the simulation
    time is the time at which to collect the stress
    zabove is the z position relative to the bottom of the nozzle, in units of zunits
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    volume = True to use the whole volume of the nozzle, if there is a file. otherwise, use a slice collected from the center of the nozzle
    '''
    # get data
    if yvar=='shearstressz':
        theta = row[cvar]
        if not volume:
            xlist,ylist = shearStressCalcSlice(row['folder'], time, zunits)
        else:
            z0 = row['nozzle_bottom_coord']
            xlist,ylist = shearStressCalcVolume(row['folder'], time, zunits, z0)
        li = int(len(xlist)/2)
        zstress = []
    else:
        zstress = nozzleLineTrace(row['folder'], time, zabove, zunits=zunits, volume=volume)
        theta, li, xlist, ylist = getListsFromTrace(row, zstress, xvar, yvar, cvar)
    return theta, li, zstress, xlist, ylist


def withinNozzleFolder(i:int, row:pd.Series, time:float, zabove:float, ax, yvar:str, cvar:str, cm, legendloc:str='overlay', zunits:str='mm', xvar:str='x', volume:bool=False, zstress:List=[], **kwargs) -> None:
    '''get data from inside the nozzle and plot it for a single folder
    i is the row number, used for determining color
    row holds metadata about the simulation
    time is the time at which to collect the stress
    zabove is the z position relative to the bottom of the nozzle, in units of zunits
    ax is the axis to plot the data on
    yvar is the variable to plot on the y axis
    cvar is the variable to color by   
    legendloc= overlay, right, or center legend location
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    xvar is the variable to plot on the x axis
    volume = True to use the whole volume of the nozzle, if there is a file. otherwise, use a slice collected from the center of the nozzle
    if zstress is empty, this calculates the stress list using getDataWithinNozzle, otherwise it uses the already collected zstress
    '''
    xlist = []
    if len(zstress)==0:
        theta, li, zstress, xlist, ylist = getDataWithinNozzle(xvar, yvar, cvar, row, time, zabove, zunits, volume)
    else:
        theta, li, xlist, ylist = getListsFromTrace(row, zstress, xvar, yvar, cvar)

    if len(xlist)==0:
        folder = row['folder']
        logging.warning(f'No data collected in {folder}')
        return zstress
    
    # set color label
    if cvar=='nozzle_angle':
        clabel = f'{int(theta)}$\degree$'
    else:
        clabel = theta

    # plot data
    ax.plot(xlist, ylist, label=clabel, c=cm[i], linewidth=0.75)

    # label line
    if legendloc=='overlay':
        x0 = xlist[li]
        y0 = ylist[li]
        if li==0:
            ha = 'right'
        else:
            ha = 'center'
        ax.text(x0, y0, clabel, color=cm[i], horizontalalignment=ha, verticalalignment='top') 
        
    return zstress
    
    
def labelNozAxs(ax, yvar:str, zunits:str, xvar:str, zabove:float, time:float, legendloc:str, **kwargs):
    '''put plot labels on the axis
    ax is the axis to plot the data on
    yvar is the variable to plot on the y axis
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    xvar is the variable to plot on the x axis
    zabove is the z position relative to the bottom of the nozzle, in units of zunits
    time is the time at which to collect the stress
    legendloc= overlay, right, or center legend location
    '''
    if yvar=='shearstressmag' or yvar=='shearstressz':
        ax.set_ylabel('Shear Stress (Pa)')
        if not ('logy' in kwargs and kwargs['logy']):
            ax.set_ylim(bottom=0)
    elif yvar=='nu_ink':
        ax.set_ylabel('Viscosity (Pa*s)')
        ax.set_yscale('log')
    elif yvar=='magu':
        ax.set_ylabel('Velocity (mm/s)')
        
    zunname = unitName(zunits)
    if yvar=='shearstressz':
        ax.set_xlabel(f'$z$ position ({zunname})')
        if not ('logy' in kwargs and kwargs['logy']):
            ax.vlines([0], 0, 1, transform=ax.get_xaxis_transform(),  color='#666666', linestyle='--', linewidth=0.75)
            ax.text(-0.05,0, 'nozzle exit', horizontalalignment='right', verticalalignment='bottom')
        ax.set_title(f'{time} s')
    else:
        if xvar=='rbar':
            ax.set_xlabel('$r/r_0$')
        elif xvar=='x':
            ax.set_xlabel(f'$x$ (mm)')
        
        ax.set_title(f'{zabove} {zunname} before exit, {time} s')
    
    if not legendloc=='overlay':
        ax.legend(loc='lower left', bbox_to_anchor=(0,1))
        

def withinNozzle(folders:List[str], time:float, zabove:float, axs, cvar:str, yvars:str, legendloc:str='overlay', zunits:str='mm', xvar:str='x', volume:bool=False,  **kwargs) -> None:
    '''plot line traces of value yvars across the nozzle as a function of xvar at position zabove relative to the bottom of the nozzle in zunits and time time, on axes axs, coloring the lines by variable cvar
     zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    volume=True to use points from the whole nozzle volume, False to use only a slice at y=0
    '''
    _, u = extractTP(folders[0], units=True) # get units

    if not 'cm' in kwargs:
        cm = sns.color_palette('viridis', n_colors=len(folders)) # uses viridis color scheme
    else:
        cm = kwargs['cm']
    
    tplist = pd.DataFrame([extractTP(folder) for folder in folders])
    tplist.sort_values(by=cvar, inplace=True)
    tplist.reset_index(drop=True, inplace=True)
    tplist[cvar] = expFormatList(tplist[cvar])
    maxy = 0
    
    for i,row in tplist.iterrows():
        zstress =  []
        for j,yvar in enumerate(yvars):
            zstress = withinNozzleFolder(i, row, time, zabove, axs[j], yvar, cvar, cm, legendloc=legendloc, zunits=zunits, xvar=xvar, volume=volume, zstress = zstress, **kwargs)
            
    # add labels
    for j,yvar in enumerate(yvars):
        labelNozAxs(axs[j], yvar, zunits, xvar, zabove, time, legendloc, **kwargs)

    
def withinNozzle0(topFolder:str, exportFolder:str, time:float, zabove:float, zunits:str='mm', cvar:str='nozzle_angle'
                  , overwrite:bool=False, export:bool=True, fontsize:int=8, **kwargs):
    '''plots line traces within the nozzle at a given z position and time. 
    topfolder is the folder holding the simulations. you can filter the folder using **kwargs, as described in plot_general.listTPvalues
    exportFolder is the folder to export figures to
    time is the time at which to collect the stress
    zabove is the z position relative to the bottom of the nozzle, in units of zunits
    zunits is 'mm' or any parameter in legend, e.g. 'nozzle_inner_width'
    cvar is the variable to color by
    overwrite True to overwrite existing files
    export True to export figures
    '''
    

    labels = ['trace_across']
    fn = intm.imFn(exportFolder, labels, topFolder, **kwargs) # output file name
    if (not overwrite and export) and os.path.exists(f'{fn}.png'):
        return

    plt.rc('font', size=fontsize) 
    
    folders = fp.caseFolders(topFolder)
    folders, _ = listTPvalues(folders, **kwargs) # remove any values that don't match
    
    if len(folders)==0:
        logging.warning('No files in list')
        return

    fig, axs = plt.subplots(1,3, figsize=(6.5,3.5))
    
    if 'logy' in kwargs:
        if kwargs['logy']:
            for ax in [axs[0], axs[1]]:
                ax.set_yscale('log')
    
    withinNozzle(folders, time, zabove, axs, cvar, ['shearstressz', 'shearstressmag', 'magu'], zunits=zunits, **kwargs) # plot the values on the axis
   
    for ax in axs:
        setSquare(ax)
  
    subFigureLabels(axs, inside=False)
    fig.tight_layout()

    if export:
        intm.exportIm(fn, fig) # export figure