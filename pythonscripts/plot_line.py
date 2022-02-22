#!/usr/bin/env python
'''Plotting line traces from Paraview'''

# external packages
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot_general import *
import interfacemetrics as intm

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
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

############ line plots

  
    

def linePlot(folder:str, time:float, ax:plt.Axes, color, yvar:str='vx', label:str='', **kwargs) -> None:
    '''plot the result of a line trace collected with paraview
        colorf is the function to use to determne the color of the plotted line, e.g. sigfunc
        rang is the total list of values that you get when you evaluate colorf over the whole folder
        yvar is 0 to plot velocities, 1 to plot viscosities'''
    t1,units = intm.importLine(folder, time, **kwargs)
    if len(t1)==0:
        logging.warning(f'Line file is missing in {folder}')
        return 
    t1 = t1.sort_values(by='z')
    if yvar in t1:
        ystrink = yvar
        ystrsup = yvar
    elif yvar=='nu':
        tp = extractTP(folder)
        ystrsup = 'nu_sup'
        ystrink = 'nu_ink'
        # viscosity yvar
        if 'nu_sup' in t1:
            t1['nu_sup']=t1['nu_sup']*10**3 # this is specifically when the density is 10**3 kg/m^3
        else:
            nu_suplist = [tp['nusup'] for i in range(len(t1))]
            t1['nu_sup'] = nu_suplist
        if 'nu_ink' in t1:
            t1['nu_ink']=t1['nu_ink']*10**3
        else:
            nu_inklist = [tp['nuink'] for i in range(len(t1))]
            t1['nu_ink'] = nu_inklist
    elif yvar=='shearrate':
        ystrsup = 'shearrate'
        ystrink = 'shearrate'
        t1['shearrate'] = [np.linalg.norm(np.array([[row['shearrate0'],row['shearrate1'],row['shearrate2']],[row['shearrate3'],row['shearrate4'],row['shearrate5']],[row['shearrate6'],row['shearrate7'],row['shearrate8']]])) for i,row in t1.iterrows()]
    else:
        raise NameError(f'yvar must be column in line csv or \'nu\', given {yvar}')

    inkPts = t1[t1['alpha']>0.5]

    zname = 'z'
    minx = min(inkPts['z'])
    maxx = max(inkPts['z'])
    supPtsLeft = t1[t1['z']<minx]
    supPtsRight = t1[t1['z']>maxx]
    
    for suppts in [supPtsLeft, supPtsRight]:
        ax.plot(suppts['z'], suppts[ystrsup], color=color, linewidth=0.75)
    ax.plot(inkPts['z'], inkPts[ystrink], color=color, linestyle='--', linewidth=0.75)
    pts = t1[(t1['z']==minx) | (t1['z']==maxx)]
    ax.scatter(pts['z'], pts[ystrink], color=color, label=label)
    ax.scatter(pts['z'], pts[ystrsup], color=color)
    
    
def labDict(yvar:str) -> str:
    '''dictionary for variable headers in line traces'''
    if yvar=='vx':
        return ('$x$ velocity (mm/s)')
    elif yvar=='vz':
        return ('$z$ velocity (mm/s)')
    elif yvar=='nu':
        return ('Viscosity (Pa$\cdot$s)')
    elif yvar=='p':
        return ('Pressure (Pa)')
    else:
        return (yvar)

def linePlots(folders:List[str], cvar, time:float=2.5, imsize:float=3.25, yvar:str='vx', legend:bool=True, xlabel:bool=True, ylabel:bool=True, **kwargs) -> plt.Figure:
    '''files is a list of folders (e.g. nb16, nb17) to include in the plot
    cvar is the function to use for deciding plot colors. func should be as a function of transport properties dictionary or a string
    e.g. func could be multfunc'''
    
    if 'ax' in kwargs and 'fig' in kwargs:
        ax = kwargs['ax']
        fig = kwargs['fig']
        kwargs.pop('ax')
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False, figsize=(imsize, imsize))
    
    if type(cvar) is str: 
        # this is a column header from extractTP
        tplist = pd.DataFrame([extractTP(folder) for folder in folders])
        tplist.sort_values(by=cvar, inplace=True)
        tplist.reset_index(drop=True, inplace=True)
        tplist[cvar] = expFormatList(tplist[cvar])
        cm = sns.color_palette('viridis', as_cmap=True) # uses viridis color scheme
        for i,row in tplist.iterrows():
            linePlot(row['folder'], time, ax, cm(i/len(tplist)), yvar, cmap=cm, label=row[cvar], **kwargs)
    else:
        # this is a function operated on extractTP
        funcvals = unqListFolders(folders, cvar)
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        for f in folders:
            val = folderToFunc(folder, colorf)
            fracval = decideRatio(val, rang)
            if fracval==0.5:
                color = 'gray'
            else:
                color = cmap(fracval)
            linePlot(f, time, ax, color, yvar, label=decideFormat(val), **kwargs)
    if legend:
        ax.legend(bbox_to_anchor=(1, 1), loc='upper right', ncol=1)
    if xlabel:
        ax.set_xlabel('$z$ (mm)')
    if ylabel:
        ax.set_ylabel(labDict(yvar))
    if yvar=='nu' or yvar=='shearrate':
        ax.set_yscale('log')
    return fig

def linePlots0(topFolder:str, exportFolder:str, cvar:str, time:float, imsize:float=6.5, yvar:str='vx', overwrite:bool=False, export:bool=True, fontsize:int=8, **kwargs) -> plt.Figure:
    '''plot the value over the line trace'''

    
    labels = ['line', time, yvar, cvar]
    fn = intm.imFn(exportFolder, labels, topFolder, **kwargs) # output file name
    if not overwrite and os.path.exists(fn+'.png'):
        return

    plt.rc('font', size=fontsize) 
    
    folders = fp.caseFolders(topFolder)
    folders, _ = listTPvalues(folders, **kwargs) # remove any values that don't match
    
    if type(yvar) is list:
        fig, axs = multiPlots(len(yvar), sharex=True, imsize=imsize)
        cols = len(axs[0])
        rows = len(axs)
        axlist = axs.flatten()
        for i,m in enumerate(yvar):
            xlabel = ((i/cols+1)>=np.floor(len(yvar)/cols))
            legend = (i==(len(yvar)-1))
            linePlots(folders, cvar, time=time, yvar=m, ax=axlist[i], fig=fig, legend=legend, xlabel=xlabel, **kwargs)
        for ax in axlist:
            setSquare(ax)
        subFigureLabels(axs)
        plt.subplots_adjust(hspace=0)
        fig.tight_layout()
    else:
        fig = linePlots(folders, cvar, time, imsize, yvar, **kwargs)
    if export:
        intm.exportIm(fn, fig) # export figure



#------------------------------


def linePressure(folder:str) -> dict:
    '''get the pressure differential across the nozzle'''
    file = os.path.join(folder, 'line_t_10_z_0.5.csv')
    if not os.path.exists(file):
        return {}, {}
    data, units = intm.importPointsFile(file) 
    pUpstream = data[(data.x>-2.9)&(data.x<-2.87)].p.max() # upstream pressure
    pDownstream = data[(data.x>-1.96)&(data.x<-1.9)].p.max() # downstream pressure
    dp = pUpstream-pDownstream
    pdict = {'pU':pUpstream, 'pD':pDownstream, 'dP':dp}
    units = {'pU':units['p'], 'pD':units['p'], 'dP':units['p']}
    meta, u = extractTP(folder, units=True)
    retval = {**meta, **pdict}
    units = {**u, **units}
    return retval, units

def linePressureRecursive(folder:str) -> dict:
    '''find line pressures for all sims in folder, all the way down'''
    if not os.path.isdir(folder):
        return [], {}
    r,u = linePressure(folder)
    if len(r)>0:
        return [r],u
    else:
        rlist = []
        units = {}
        for f in os.listdir(folder):
            r,u = linePressureRecursive(os.path.join(folder, f))
            if len(r)>0:
                rlist = rlist+r
                units = u
        return rlist, units

def linePressures(topfolder:str, exportFolder:str, filename:str) -> dict:
    '''find line pressures for all sims in folder and export'''
    rlist, units = linePressureRecursive(topfolder)
    tt = pd.DataFrame(rlist)
    if os.path.exists(exportFolder):
        fn = os.path.join(exportFolder, filename)
        col = pd.MultiIndex.from_tuples([(k,v) for k, v in units.items()])
        data = np.array(tt)
        df = pd.DataFrame(data, columns=col)       
        df.to_csv(fn)
        logging.info(f'Exported {fn}')
    return tt,units

    
    