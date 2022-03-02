#!/usr/bin/env python
'''Functions for plotting overall metrics, such as simulation time, folder name, simulation rate, cross-sectional area...'''

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
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------


    

def setSquare(ax):
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')


def plotSquare(ax:plt.Axes, x0:float, y0:float, dx:float, caption:str, color) -> None:
    '''plotSquare plots a square
    ax is axis to plot on
    x0 is the x position
    y0 is the y position
    dx is the spacing between xlist
    caption is the label
    color is the color of the circle'''
    if len(caption)>0:
        # calculate brightness of color
        l = 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2]
        if l>0.4:
            txtcolor = 'black'
        else:
            txtcolor = 'white'
        ax.text(x0, y0, caption, horizontalalignment='center', verticalalignment='center', color=txtcolor)
    box = plt.Rectangle([x0-dx/2,y0-dx/2], dx, dx, color=color, ec=None)
    ax.add_artist(box)
    
    
def plotSS(ss:pd.DataFrame, column:str, tmin:float) -> None:
    '''plot slice summaries. produces a row of 2 plots: column as a function of x and as a function of time
    ss is a dataframe of slice summary data, e.g. as imported by importSS or created by summarize
    column is the column name, e.g. 'centery'
    tmin is the minimum time to include in these plots'''
    size=4
    ss2 = ss[ss['time']>=tmin]
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, figsize=(10, 4))
    p = ss2.plot.scatter(x='xbehind', y=column, c='time', legend=False, ax=axes[0], colormap = 'RdBu')
    axes[0].set_xlabel('Position (mm)')
    axes[0].set_ylabel(column)   
    p = ss2.plot.scatter(x='time', y=column, c='x', legend=False, ax=axes[1], colormap = 'RdBu')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel(column)
    
      

#------------------------------------------
# text plot, shows which folders correspond to which viscosities    


def txtPlot(folder:str, cp:folderPlots, dt:float) -> None:
    '''txtPlot assigns a single folder to the plot
    folder is a full path name
    cp is a comboPlot object
    dt is the spacing between text items. A good value is 0.2'''
    try:
        [color, x0, y0, sigmapos] = vvplot(folder, cp) # find the position of this plot within the big plot
    except Exception as e:
        return
    xmid = x0
#     ymid = y0+dt*(sigmapos-1) # xmid and ymid are positions in the plot
    ymid = y0+dt # xmid and ymid are positions in the plot
    b = os.path.basename(folder)
    if cp.split:
        axnum = sigmapos
    else:
        axnum = 0
    cp.axs[axnum].text(xmid, ymid, b, horizontalalignment='center', verticalalignment='center', c=color) # put the folder name on the plot
   
    
def txtPlots0(topFolder:str, exportFolder:str, overwrite:bool=False, export:bool=False, **kwargs) -> None:
    '''write names of a list of folders on one big plot
    topFolder is the folder that holds all of the folders
    exportFolder is the folder to export the figure to
    overwrite true to overwrite existing files'''
    labeli = 'names'
    fn = intm.imFn(exportFolder, labeli, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    cp = comboPlot(topFolder, [-0.6, 0.6], [-0.6, 0.6], 6.5, **kwargs)
    if len(cp.flist)==0:
        return
    if cp.split:
        dt = 0
    else:
        dt = 0.2
        
    for folder in cp.flist:
        txtPlot(folder, cp, dt)
    cp.figtitle = 'Folder names'
    cp.clean()
    if export:
        intm.exportIm(fn, cp.fig, **kwargs)
    return

#------------------------------------------   
# run time plot: how long the simulation ran in simulation seconds


def runtimePlot(folder:str, cp:folderPlots, dt:float) -> None:
    '''runtimePlot assigns a single folder to the plot
    folder is a full path name
    cp is a comboPlot object
    dt is the spacing between text items. A good value is 0.2'''
    try:
        [color, x0, y0, sigmapos] = vvplot(folder, cp) # find the position of this plot within the big plot
    except:
        return
    xmid = x0
#     ymid = y0+dt*(sigmapos-1) # xmid and ymid are positions in the plot
    ymid = y0+dt # xmid and ymid are positions in the plot
    b = fp.currentTime(folder)[0]
    if cp.split:
        axnum = sigmapos
    else:
        axnum = 0
    cp.axs[axnum].text(xmid, ymid, b, horizontalalignment='center', verticalalignment='center', c=color) # put the folder name on the plot
   
    
def runtimePlots0(topFolder:str, exportFolder:str, overwrite:bool=False, **kwargs) -> None:
    '''write names of a list of folders on one big plot
    topFolder is the folder that holds all of the folders'''
    labeli = 'runtime'
    fn = intm.imFn(exportFolder, labeli, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    cp = comboPlot(topFolder, [-0.6, 0.6], [-0.6, 0.6], 6.5, **kwargs)
    if len(cp.flist)==0:
        return
    dt = 0.2
    for folder in cp.flist:
        runtimePlot(folder, cp, dt)
    cp.figtitle = 'Run times'
    cp.clean()
    intm.exportIm(fn, cp.fig, **kwargs)
    
 #------------------------------------------   
#### generic value plots

def valueCaption(val:str, tmax:str, tmin:str) -> str:
    '''Caption that shows the value, with an appropriate number of decimals. If the val is outside the captioning region, don't generate a caption. If the val is large or tmax is large, use 0 decimals. If they are small, use 2.'''
    if math.isnan(val):
        return ''
    if val>=tmax or val<=tmin:
        if val>10 or tmax>100:
            caption='%2.0f'%(val)
        else:
            caption = '%1.2f'%(val)
    else:
        caption = ''
    return caption
        
        
def plotTableVals(t1:pd.DataFrame, cp:comboPlot, tminmode:int, timeplot:bool=False) -> Dict:
    '''Plot a list of values on a comboPlot using either circle size plot or color density plot. 
    t1 is a dataframe made from timePlot outputs
    cp is the comboPlot to plot the values on
    tminmode=0 to set the minimum to 0. tminmode=1 to set the minimum to the min value in the table.
    timeplot true if we are plotting times. This is necessary for circle scaling.'''

    
    if len(t1)<1:
        raise ValueError
    if 'tmin' in cp.kwargs:
        tmin = cp.kwargs['tmin']
    else:
        if tminmode==0:
            tmin = 0
        else:
            tmin = t1['rate'].min()
    if 'tmax' in cp.kwargs:
        tmax = cp.kwargs['tmax']
    else:
        tmax = t1['rate'].max()
    
        
    t1 = t1.sort_values(by=['rate'])
    t1 = t1.reset_index(drop=True)

    # set circle size
    if not cp.split:
        rmax = cp.dx # maximum radius of circle is the spacing between points/2
        if tmax>(100*t1['rate'].median()) or timeplot:
            rmax = rmax*np.sqrt(tmax/500)
    
    #cmap = sns.cubehelix_palette(as_cmap=True, rot=-0.4)
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    #cmap = sns.cubehelix_palette(as_cmap=True)
    # label only evenly spaced across values
    spacing = max(float(np.floor(len(t1)/20)),1)
    dummy = tmin
    for index,t in t1.iterrows():
        if float(index)%spacing==0 or cp.split:
            dummy = t['rate']
        caption = valueCaption(t['rate'], tmax, dummy)
        val = (t['rate']-tmin)/(tmax-tmin) # actual value, normalized to range
        
        if cp.split:
            if val>1:
                color = [138/256, 75/256, 60/256]
            elif val<0:
                color = [33/256, 85/256, 101/256]
            else:
                color = cmap(val)
            sp = t['sigmapos']
            if sp<len(cp.axs):
                ax = cp.axs[t['sigmapos']]
                plotSquare(ax, t['x0'], t['y0'], cp.dx, caption, color)
        else:
            ax = cp.axs[0] 
            plotCircle(ax, t['x0'], t['y0'], np.sqrt(val)*rmax, caption, t['color'], sigmapos=t['sigmapos'])

    return {'tmin':tmin, 'tmax':tmax, 'cmap':cmap}


def folderToPlotVals(folder:str, cp, rate) -> Dict:
    '''this is the output dictionary for a timePlot or metricPlot entry from one simulation
    folder is a full path name
    cp is a comboPlot object
    rate is a value to plot'''
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp) # find the position of this plot within the big plot
    except Exception as e:
        raise ValueError
    return {'color':color, 'x0':x0, 'y0':y0, 'rate':rate, 'sigmapos':sigmapos}


def plotAllFolderVals(function, cp:comboPlot, tminmode:int, timeplot:bool=False) -> Dict:
    '''Go through all of the folders in the folder list stored in cp, and construct a table of values to plot.
    function is the function to use on each folder. Usually some variation on folderToPlotVals
    tminmode is 0 to use 0 as the min value, 1 to use the min value in the table as min value for choosing colors
    timeplot=True if we are plotting simulation rates. Important for circle size scaling.'''
    t1 = []
    for f in cp.flist:
        try:
            row = function(f, cp)
        except NameError as n:
            raise n
        except Exception as e:
            pass
        else:
            t1.append(row)
    t2 = pd.DataFrame(t1)
    t2 = t2.dropna()
    return plotTableVals(t2, cp, tminmode, timeplot)
    
def valueLegend(cp:comboPlot, vpout:Dict) -> None:
    '''Put a color legend for the gradient plot on the bottom'''
    if cp.split:
        ylim = cp.axs[0].get_ylim()
        if ylim[1]-ylim[0]>2.5:
            y = -0.15
        elif ylim[1]-ylim[0]>1.3:
            y = 0
        else:
            y = 0.1
        cbaxes = cp.fig.add_axes([0.2, y, 0.6, 0.05])
        nm = plt.Normalize(vmin=vpout['tmin'], vmax=vpout['tmax'])
        sm = plt.cm.ScalarMappable(cmap=vpout['cmap'], norm=nm)
        plt.colorbar(sm, cax=cbaxes, orientation="horizontal")



#------------------------------------------
# time plots: how fast the simulation ran, in real hr per simulation s
        

def timePlot(folder:str, cp:comboPlot):
    '''timePlot determines the position and size of circle to plot
    folder is a full path name
    cp is a comboPlot object'''

    le = fp.currentRate(folder)
    try:
        rate = float(le['rate'])
    except:
        return
    else:
        return folderToPlotVals(folder, cp, rate) 


def timePlots(topFolder:str, exportFolder:str, overwrite:bool=False, **kwargs) -> None:
    '''timePlots plots computation rates as circles
    topFolder is a full path name to the folder containing all the simulations
    exportFolder is the folder to export the plot to'''
    
    labeli = 'simrate'
    fn = intm.imFn(exportFolder, labeli, topFolder, **kwargs)  # file name for the plot
    if not overwrite and os.path.exists(fn+'.png'):            # quit if this plot already exists and overwrite==False
        return
    
    cp = comboPlot(topFolder, [-0.6, 0.6], [-0.6, 0.6], 6.5, **kwargs)  # create a plot
    if len(cp.flist)==0:
        return
    lfunc = lambda folder, cp: timePlot(folder,cp)                      # we are going to run timePlot on every folder
    try:
        vpout = plotAllFolderVals(lfunc, cp, 0, timeplot=True)          # plot all the files
    except:
        return
    cp.figtitle = 'Simulation time (real hr/simulation s)'
    cp.clean()                  # clean up the plot
    valueLegend(cp, vpout)      # add a color legend
    intm.exportIm(fn, cp.fig, **kwargs)   # export figure


#------------------------------------------  
# slice summary metrics plots

def summaryRow(folder:str, time:float, xbehind:float) -> Tuple[dict,dict]:
    '''turn the folder into a single row, for a given time and position'''
    row, rowunits = metricVals(folder, time, xbehind, ['x', 'xbehind', 'time', 'centery', 'centerz', 'area','maxheight', 'maxwidth', 'centeryn', 'centerzn', 'arean', 'maxheightn', 'maxwidthn', 'vertdisp', 'vertdispn', 'aspectratio', 'speed', 'speeddecay'])
    meta, u = extractTP(folder, units=True)
    if not (u['dink']=='mm' and u['nuink']=='Pa*s' and u['rhoink']=='g/mL' and u['sigma']=='mJ/m^2' and u['vink']=='mm/s'):
        raise ValueError('Units of extractTP do not match needed units for summaryRow')

    for s in ['ink', 'sup']:
        meta['gdot_'+s] = meta['v'+s]/meta['d'+s] # 1/s
        u['gdot_'+s] = '1/s'
        if meta['k'+s]>0:
            mu = meta['k'+s]*(abs(meta['gdot_'+s])**(meta['n'+s]-1)) + meta['tau0'+s]/abs(meta['gdot_'+s])
            meta['visc0_'+s] = min(mu, meta['nu'+s])  # Pa*s
        else:
            meta['visc0_'+s] = meta['nu'+s]
        u['visc0_'+s] = 'Pa*s'
        meta['CaInv_'+s] = meta['sigma']/(meta['visc0_'+s]*meta['v'+s]) # capillary number ^-1 []
        u['CaInv_'+s] = ''
        meta['Re_'+s] = 10**3*(meta['rho'+s]*meta['v'+s]*meta['d'+s])/(meta['visc0_'+s]) # reynold's number
        u['Re_'+s] = ''
    
    meta['viscRatio'] = meta['visc0_ink']/meta['visc0_sup']
    u['viscRatio'] = ''
    meta['ReRatio'] = meta['Re_ink']/meta['Re_sup']  
    u['ReRatio'] = ''
    
    out = {**meta, **row}
    units = {**u, **rowunits}
    
    return out, units

def summaryTable(topfolders:str, time:float, xbehind:float, exportFolder:str, filename:str='summaryTable') -> Tuple[pd.DataFrame, dict]:
    '''collect summary data for each topfolder and put it all into a table'''
    t = []
    u = []
    for topfolder in topfolders:
        logging.info(topfolder)
        for cf in fp.caseFolders(topfolder):
            try:
                tt0,u0 = summaryRow(cf, time, xbehind)
            except:
                pass
            else:
                if len(tt0)>0:
                    t = t+[tt0]
                if len(u0)>len(u):
                    u = u0
    tt = pd.DataFrame(t)
    if os.path.exists(exportFolder):
        plainExp(os.path.join(exportFolder, f'{filename}_x_{xbehind}_t_{time}.csv'), tt, u)
    return tt,u


def metricVals(folder:str, time:float, xbehind:float, labels:List[str], units:bool=False, xunits:str='mm') -> Dict:
    '''Find the value of slice summary metrics for a single simulation.
    folder is the full path name to a simulation folder
    time is the time of the slice
    xbehind is the position of the slice, relative to the center of the nozzle
    labels is a list of metrics to collect, e.g. ['area', 'centery']'''
    if not os.path.exists(folder):
        raise ValueError(f"Path {folder} does not exist")

    ss, u = intm.importSS(folder)
        # get slice summaries
    if len(ss)<2:
        raise ValueError(f"Slice summaries too short: {folder}")
        
    ss = ss[ss['time']==time]
    
    if len(ss)<2:
        raise ValueError(f"Slice summaries too short: {folder}")
        
    try:
        if not xunits=='mm':
            le = fp.legendUnique(folder)
            xbehind = round(xbehind*float(le[xunits]), 5) # convert xbehind to mm
        xreal = intm.closest(list(ss['xbehind'].unique()), xbehind)
        # this is the actual x value that we measured that's 
        # closest to the one we're asking for
    except Exception as e:
        raise e
        
    if abs(xreal-xbehind)>0.2:
        # if the x value is too far away, abort
        raise ValueError(f"No valid x value: {folder}")
        
    row = ss[(ss['xbehind']==xreal) & (ss['time']==time)] 
        # select the slice summary at the position and time we asked for
    if not len(row)==1:
        raise ValueError(f"Not enough rows: {folder}, {xreal}, {xbehind}, {time}")
        
    try:
        rates = {label:row.iloc[0][label] for label in labels}
        # find the value of the metric we're looking for
    except:
        logging.debug(folder)
        raise NameError(f'Error collecting metric value: {folder}')
        
    if units:
        return rates, dict([[i, u[i]] for i in labels])
    else:
        return rates



def metricPlot(folder:str, cp:comboPlot, time:float, xbehind:float, label:str) -> Dict:
    '''metricPlot determines the position and size of circle or square to plot
    folder is a full path name
    cp is a comboPlot object
    time is the time since extrusion started in s
    xbehind is the distance behind the center of the nozzle in mm
    label is the column label, e.g. 'maxz'.
    '''
    
    try:
        rate = metricVals(folder, time, xbehind, [label])
        rate = rate[label]
    except Exception as e:
        raise e
    return folderToPlotVals(folder, cp, rate)
       
        

def metricPlots(topFolder:str, exportFolder:str, time:float, xbehind:float, label:str, overwrite:bool=False, **kwargs) -> None:
    '''# metricPlots plots slice summaries as color density plots
    topFolder is a full path name to the folder containing all the simulations
    exportFolder is the folder to export plots to
    time is the time since extrusion started in s
    xbehind is the distance behind the center of the nozzle in mm
    label is the column label, e.g. 'maxz' '''
    
    labeli =f'{label}_{xbehind}_t_{time}'
    fn = intm.imFn(exportFolder, labeli, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    cp = comboPlot(topFolder, [-0.6, 0.6], [-0.6, 0.6], 6.5, gridlines=False, **kwargs)
    if len(cp.flist)==0:
        return
    lfunc = lambda f, cp: metricPlot(f, cp, time, xbehind, label)
    try:
        vpout = plotAllFolderVals(lfunc, cp, 1) # use tminmode 1 so the min of the color map is the min of the metric
    except Exception as e:
        logging.error(str(e))
        return
    cp.figtitle = f'{label}, {xbehind} mm behind nozzle, t = {time} s'
    cp.clean()
    valueLegend(cp, vpout)
    intm.exportIm(fn, cp.fig, **kwargs)
    
#--------------------------------    
    
def qualityPlots(d:pd.DataFrame, time:float, xbehind:float, cvar:str='', xvar:str='theta', **kwargs):
    '''Plot points for all labels, one plot per label
    d is dataframe holding table of angles and metrics
    time is time in s
    xbehind is the distance behind the nozzle to take the slice summaries
    '''
    
    labdict = {'arean':'Area/intended', 'vertdispn':'Vert disp/nozzle diam', 'aspectratio':'Height/width', 'speeddecay':'Speed/intended'}
    yvars = list(d.keys())
    yvars.remove(xvar)
    if cvar in yvars:
        yvars.remove(cvar)
    nvars = len(yvars)
    fig, axs = multiPlots(nvars, sharex=True, imsize=6.5)
    fig.suptitle(f'Print Quality Metrics, {xbehind} mm behind nozzle, t = {time} s')
    plt.xticks(ticks=[0,5,10,15,20,25,30])

    axs[0,0].set_title(" ") # so y axis names do not overlap title
#     colors = ['#32964d', '#27cae6', '#335862', '#38f0ac'] # colors to plot with

    for r in range(rows): # creates plots, sets labels, and adds ideal lines
        for c in range(cols):
            ax = axs[r,c]
            yvar = yvars.pop(0)
            if not cvar in d or len(cvar)==0:
                ax.scatter(d[xvar], d[yvar], c='black', marker='D')
            else:
                zplot=ax.scatter(d[xvar], d[yvar], c=d[cvar], marker='D', cmap='coolwarm')
            ax.set_ylabel(labdict[yvar])
            if yvar=='vertdispn':
                yideal=0
            else:
                yideal=1
            ax.axhline(yideal, ls='--', c='black')
            ax.text(min(d[xvar]), yideal, 'ideal', ha='left', va='bottom', color='black') # key for ideal horizontal lines
    if cvar in d and len(cvar)>0:
        cb_ax = fig.add_axes([0.95,.124,.03,.7])
        if 'cvarlabel' in kwargs:
            cvl = kwargs['cvarlabel']
        else:
            cvl = cvar
        fig.colorbar(zplot, cax=cb_ax, orientation='vertical', label=cvl)
        plt.subplots_adjust(wspace=0.25, hspace=0.05)
    else:
        fig.tight_layout()
    for axrow in axs:
        for ax in axrow[:-1]:
            setSquare(ax)
            
    for ax in axs[-1]:
        if 'xvarlabel' in kwargs:
            ax.set_xlabel(kwargs['xvarlabel'])
        else:
            ax.set_xlabel(varNickname(xvar))
    return fig

    
def qualityPlots0(topFolder:str, exportFolder, time:float, xbehind:float, labels:List[str]=['arean', 'vertdispn', 'aspectratio', 'speeddecay'], overwrite:bool=False, xvar:str='nozzle_angle', cvar:str='', **kwargs) -> None: # RG
    '''Plots slice summaries for against nozzle angles in scatter plots
    topFolder is a full path name to the folder containing all the simulations
    exportFolder is the folder to export plots to
    time is the time since extrusion started in s
    xbehind is the distance behind the center of the nozzle in mm
    labels is the slice summaries to plot in order
    overwrite is whether to overwrite if the file exists'''
    
    fn = intm.imFn(exportFolder, labels, topFolder, **kwargs) # output file name
    if not overwrite and os.path.exists(fn+'.png'):
        return


    folders = fp.caseFolders(topFolder) # names of folders to plot
    if len(folders)==0:
        logging.warning('No folders found')
        return
    folders, _ = listTPvalues(folders, **kwargs) # list of transport property lists
    
    xlist = [] # list of nozzle angles
    for folder in folders:
        meta, u = extractTP(folder, units=True)
        d = {xvar:meta[xvar]}
        kwargs['xvarlabel'] = varNickname(s) + f' ({u[xvar]})'
        if cvar in meta:
            d[cvar] = meta[cvar]
            kwargs['cvarlabel'] = varNickname(s) + f' ({u[cvar]})'
        try:
            value = metricVals(folder, time, xbehind, labels, units=False)
        except Exception as e:
            pass
        else:
            d = {**d, **value}
            xlist.append(d)

    fig = qualityPlots(pd.DataFrame(xlist), time, xbehind, cvar=cvar, xvar=xvar, **kwargs)

    intm.exportIm(fn, fig) # export figure
    
    
#------------------------------------------------

def viscRatioPlot(folders:str, time:float, xbehind:float, xunits:str='mm', xvar:str='viscRatio', yvar:str='aspectratio', cvar:str='', fontsize:int=8, **kwargs) -> None:
    '''plot viscosity ratio vs. metrics'''
    
    dlist = []
    for f in folders:
        try:
            d = metricVals(f, time, xbehind, [yvar], xunits=xunits) # get the metric value
        except Exception as e:
            pass
        else:
            vr = intm.viscRatio(f)   # get the viscosity ratio
            d = {**d, **vr}
            if not cvar=='':
                le, u = extractTP(f, units=True)
                d['cvar'] = float(le[cvar])  # get the color
            else:
                d['cvar'] = 0
            dlist.append(d)   # save these values to list
    df = pd.DataFrame(dlist)
    
    if len(df)==0:
        return
    
    if 'figsize' in kwargs:
        figsize = kwargs['figsize']
    else:
        figsize = (3,3)
    plt.rc('font', size=fontsize) 
    fig,ax = plt.subplots(1,1, figsize=figsize)
    
    if not cvar=='':
        clist = list(df['cvar'].unique())
        cm = sns.color_palette('viridis', n_colors=len(clist))
        clist.sort()
        for j,cval in enumerate(clist):
            dfi = df[df['cvar']==cval]
            ax.scatter(dfi[xvar], dfi[yvar], color=cm[j], label=cval)
        ax.legend(title=varNickname(cvar, short=True, units=u))
    else:
        ax.scatter(df[xvar], df[yvar])
        
    setSquare(ax)
    ax.set_xlabel(varNickname(xvar, short=False, units=u))
    ax.set_ylabel(varNickname(yvar, short=False, units=u))
    return fig
    
def viscRatioPlot0(topFolder:str, exportFolder:str, time:float, xbehind:float, xunits:str='mm', xvar:str='viscRatio', yvar:str='aspectratio', overwrite:bool=False, export:bool=False, cvar:str='', **kwargs) -> None:
    
    fn = intm.imFn(exportFolder, [yvar, xvar], topFolder, **kwargs) # output file name
    if export and (not overwrite and os.path.exists(fn+'.png')):
        return
    
    folders = fp.caseFolders(topFolder) # names of folders to plot
    if len(folders)==0:
        logging.warning('No folders found')
        return
    folders, _ = listTPvalues(folders, **kwargs) # list of transport property lists
    
    fig = viscRatioPlot(folders, time, xbehind, xvar=xvar, xunits=xunits, yvar=yvar, cvar=cvar, **kwargs)
    
    if export:
        intm.exportIm(fn, fig) # export figure
