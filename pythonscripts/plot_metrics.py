#!/usr/bin/env python
'''Functions for plotting overall metrics, such as simulation time, folder name, simulation rate, cross-sectional area...'''

# external packages
import sys
import os
import numpy as np
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

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plots
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

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
   
    
def txtPlots0(topFolder:str, exportFolder:str, overwrite:bool=False, **kwargs) -> None:
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
    intm.exportIm(fn, cp.fig, **kwargs)

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
    b = fp.currentTime(folder)
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
    t1 = t1.reset_index()

    
    # set circle size
    if not cp.split:
        rmax = cp.dx/2 # maximum radius of circle is the spacing between points/2
        if tmax>(100*t1['rate'].median()) or timeplot:
            rmax = rmax*np.sqrt(tmax/500)
    #         rmax = rmax*(tmax/(100*t1['rate'].median()))

    
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
    
    le = intm.importLegend(folder)
    rate = float((le[le['title']=='simulation rate (hr/s)']).val) # find the simulation rate
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


def metricVals(folder:str, time:float, xbehind:float, labels:List[str]) -> Dict:
    '''Find the value of slice summary metrics for a single simulation.
    folder is the full path name to a simulation folder
    time is the time of the slice
    xbehind is the position of the slice, relative to the center of the nozzle
    labels is a list of metrics to collect, e.g. ['area', 'centery']'''
    if not os.path.exists(folder):
        raise ValueError(f"Path {folder} does not exist")

    le, units = intm.importSS(folder)
        # get slice summaries
    if len(le)<2:
        raise ValueError(f"Slice summaries too short: {folder}")
    try:
        xreal = intm.closest(list(le['xbehind'].unique()), xbehind) 
        # this is the actual x value that we measured that's 
        # closest to the one we're asking for
    except Exception as e:
        raise e
    if abs(xreal-xbehind)>0.2:
        # if the x value is too far away, abort
        raise ValueError(f"No valid x value: {folder}")
    row = le[(le['xbehind']==xreal) & (le['time']==time)] 
        # select the slice summary at the position and time we asked for
    if not len(row)==1:
        raise ValueError(f"Not enough rows: {folder}, {xreal}, {xbehind}, {time}")
    try:
        rates = {label:row.iloc[0][label] for label in labels}
        # find the value of the metric we're looking for
    except:
        logging.debug(folder)
        raise NameError(f'Error collecting metric value: {folder}')
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
    
    labeli = label+'_'+str(xbehind)+'_t_'+str(time)
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
    cp.figtitle = label+', '+str(xbehind)+' mm behind nozzle, t = '+str(time)+' s'
    cp.clean()
    valueLegend(cp, vpout)
    intm.exportIm(fn, cp.fig, **kwargs)
    
    
    
def qualityPlots(rows:int, time:float, xbehind:float, xlist:List[int], matrix:List[List[int]], labels:List[str]): # RG
    '''Plot points for all labels, one plot per label
    rows is the number of labels
    time is time in s
    xbehind is the distance behind the nozzle to take the slice summaries
    xlist is the angles to plot
    matrix contains all slice summaries for all angles
    labels is the slice summaries to extract'''
    
    yval = ['NA']*rows # list to put the labels in
    for i, name in enumerate(labels): # allows the labels to be passed in any order
        if name=='arean':
            yval[i] = 'Normalized Area'
        elif name=='vertdispn':
            yval[i] = 'Normalized Vertical Displacement'
        elif name=='aspectratio':
            yval[i] = 'Aspect Ratio'
        elif name=='speeddecay':
            yval[i] = 'Speed Decay'
    
    fig, axs = plt.subplots(2, 2, sharex=True, constrained_layout=True)
    fig.suptitle('Print Quality Metrics, '+str(xbehind)+' mm behind nozzle, t = '+str(time)+' s')
    plt.xticks(ticks=[0,5,10,15,20,25,30])
    
    fig.text(0.54, -0.04, 'Nozzle angle (degrees)', ha='center') # x label, did not center when using fig.set_xlabel
    fig.text(0.95, 0.95, '---- ideal', ha='center', color='#2b6cb3') # key for ideal horizontal lines
    axs[0,0].set_title(" ") # so y axis names do not overlap title
    colors = ['#32964d', '#27cae6', '#335862', '#38f0ac'] # colors to plot with
    
    nvd = yval.index('Normalized Vertical Displacement') # find location of nvd, necessary because it has an ideal of 0 not 1
    
    for r in range(2): # creates plots, sets labels, and adds ideal lines
        for c in range(2):
            axs[r,c].scatter(xlist, matrix[r*2+c][:], c=colors[r*2+c], marker='D')
            axs[r,c].set_ylabel(yval[r*2+c])
            if r*2+c==nvd:
                axs[r,c].axhline(0, ls='--', c='#2b6cb3')
            else:
                axs[r,c].axhline(1, ls='--', c='#2b6cb3')
    return fig

    
def qualityPlots0(topFolder:str, exportFolder, time:float, xbehind:float, labels:List[str], overwrite:bool=False, **kwargs) -> None: # RG
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
            
    xlist = [] # list of nozzle angles
    for theta in folders:
        thetaValue = theta.split('cn')[1] # take only the angle from the folder path
        xlist.append(thetaValue)
        
    xlist = [int(i) for i in xlist] # cast xlist to ints     
    idx = np.argsort(xlist) # create indices for sorted xlist
    xlist = [xlist[i] for i in idx] # sort xlist
    folders = [folders[i] for i in idx] # sort folders so the matrix later is already sorted

    ylist = [] # y values as a list of dicts containing all labels
    for theta in folders:
        value = metricVals(theta, time, xbehind, labels) # get slice summary values for each label
        ylist.append(value)
    
    rows = len(labels) # a row for each label
    cols = len(xlist) # a column for each nozzle angle
    matrix = [[0 for x in range(cols)] for y in range(rows)] # a 2d list to hold y values for all labels
    
    for r in range(rows):
        current = [sub[labels[r]] for sub in ylist] # get all of the values for the current label without their key
        for c in range(cols):
            matrix[r][c] = current[c] # add each value to its appropriate position in the 2d list
            
    fig = qualityPlots(rows, time, xbehind, xlist, matrix, labels) # plot values
    intm.exportIm(fn, fig) # export figure


def shearStressPlots(time:float, ss:List[float], alist:List[float]): # RG
    '''Plot average shear stress over the nozzle
    time is time in s
    ss is the average shear stress of each nozzle down the nozzle length
    alist is the angles to plot'''
    
#     xmin = 0.3015 # tip of nozzle
#     xmax = 2.05623 # 3% below top of nozzle
    xmin = 0 # 3% below top of nozzle
        # xmin is 0 so that our frame of reference is of the ink ask it flows down the nozzle
    xmax = 1.75473 # tip of nozzle, translated with reference
    dx = 0.05
    xlist = np.arange(xmin, xmax, dx) # evenly spaced values
    
    xvar = 'z position (mm)'
    yvar = 'Shear Stress (Pa)'
    
    fig, ax = plt.subplots(constrained_layout=True)
    fig.suptitle('Average shear stress along the nozzle')
    cm = sns.color_palette('viridis', n_colors=7) # uses viridis color scheme

    for i,name in enumerate(alist):
        ax.plot(xlist, ss[i], label=str(name)+' degrees', c=cm[i]) # plot each shear stress by distance from top of nozzle and label by name
      
    ax.set_xlabel(xvar)
    ax.set_ylabel(yvar)
    ax.legend()
    
    return fig
    

def shearStressCalc(folder:str): # RG
    '''Calculates mean shear stress across the length of the nozzle
    folder is the folder to do calculations on'''
    
    df,units = intm.importPtsNoz(folder, 2.5) # get points in nozzle
    vals = [] # average shear stress for current angle
    realz = round(df['z'],5) # z values in CSV, rounded to avoid floating point error
    delta = len(np.unique(realz)) # number of cross-sections
    zlist = [round(0.3015+(delta-i-1)*0.05,5) for i in range(delta)] # z values expected, rounded to avoid floating point error

    for i in range(delta): # go through each cross-section
        mask = df[realz==zlist[i]] # extract only the cells at current z value Currently extracts nothing
        ss = mask['shearstressmag'] # extract only shear stress magnitude
        av = np.mean(ss) # average shear stress
        vals.append(av) # add average to vals
        
    return vals


def shearStressPlots0(topFolder:str, exportFolder, time:float, overwrite:bool=False, **kwargs) -> None: # RG
    '''Plots average shear stress along the length of the nozzle
    topFolder is a full path name to the folder containing all the simulations
    exportFolder is the folder to export plots to
    time is the time since extrusion started in s
    overwrite is whether to overwrite if the file exists'''
    
    labels = ['shear_stress']
    fn = intm.imFn(exportFolder, labels, topFolder, **kwargs) # output file name
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    folders = fp.caseFolders(topFolder)
            
    alist = [] # a list of nozzle angles
    for theta in folders:
        thetaValue = theta.split('cn')[1]
        alist.append(thetaValue)
    
    alist = [int(i) for i in alist] # cast xlist to ints     
    idx = np.argsort(alist) # create indices for sorted xlist
    alist = [alist[i] for i in idx] # sort xlist
    folders = [folders[i] for i in idx] # sort folders to match sorted angles

    allVals = [] # average shear stresses for all angles
    for theta in folders: # for each folder
        ss = shearStressCalc(theta) # calculate average shear stresses
        allVals.append(ss)
    
    fig = shearStressPlots(time, allVals, alist) # create figure
    intm.exportIm(fn, fig) # export figure
