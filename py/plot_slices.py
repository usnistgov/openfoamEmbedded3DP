#!/usr/bin/env python
'''Functions for plotting cross-sections'''

# external packages
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
import interfacemetrics as intm
from plot_general import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='10.0')


#-------------------------------------------

#### plot slices


def plotSlices(data:pd.DataFrame) -> plt.Figure:
    '''use this to plot all slices within a pandas dataframe containing the points on the interface
    data is a pandas dataframe, e.g. one imported from importpoints'''
    d = data
    ticks = [-1, -0.5, 0, 0.5, 1]
    size = 8
    if len(xpts(data))>1:
        cmap = 'RdBu'
        color = 'x'
        sidelegend = True
    else:
        cmap = None
        color = 'Black'
        sidelegend = False
    p = d.plot.scatter(x='y', y='z', c=color, xlim=[-1, 1], ylim=[-1, 1], legend=False,\
               colormap = cmap, figsize=[size, size], \
               xticks=ticks, yticks=ticks, sharex=False, s=0.1)
    p.set_aspect('equal', adjustable='box')
    p.set_xlabel('y (mm)')
    p.set_ylabel('z (mm)')
    if sidelegend:
        f = plt.gcf()
        cax = f.get_axes()[1]
        cax.set_ylabel('x (mm)')
        cax.figsize=[0.5, size]
    return p



#### cross-section plots

def angle360(row:pd.Series, center:List) -> float:
    '''get the angle between the vector from center to row[y,z], relative to the x axis'''
    v = [row['y']-center[0],row['z']-center[1]]
    if v==[0,0]:
        return 0
    angle = np.arccos(np.dot(v, [1,0])/(np.linalg.norm(v)))
    if v[1]<0:
        angle = 2*np.pi-angle
    return angle

def plotXSOnAx(pts:pd.DataFrame, ax, color) -> None:
    '''sort the points by radial position and plot as line'''
    left = pts.y.min()
    right = pts.y.max()
    bottom = pts.z.min()
    top = pts.z.max()
    center = [(left+right)/2, (top+bottom)/2]
    pts['theta'] = [angle360(row, center) for i,row in pts.iterrows()]
    pts.sort_values(by='theta', inplace=True, ignore_index=True)
    ptlists = [[]]
    for i,row in pts.iterrows():
        if i==0:
            ptlists[0].append(dict(row))
        else:
            j = 0
            while j<len(ptlists) and np.sqrt((row['y']-ptlists[j][-1]['y'])**2+(row['z']-ptlists[j][-1]['z'])**2)>0.1:
                # point too far away
                j+=1
            if j<len(ptlists):
                ptlists[j].append(row)
            else:
                ptlists.append([row])
    for ptsi in ptlists:
        ptsidf = pd.DataFrame(ptsi)
        ax.plot(ptsidf['y'], ptsidf['z'], color=color, linewidth=1, marker=None)

def XSPlot(xs:pd.DataFrame, folder:str, cp:comboPlot, recenter:bool=False, **kwargs) -> None:
    '''plot cross-sections for whole folder
    xs is a pandas DataFrame holding points
    folder is a full pathname
    cp is a comboPlot object'''
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except:
        logging.info(folder)
        traceback.print_exc()
        return
    if xs.z.max()>cp.dy/2:
        # xs is going to fall outside of the plot bounds
        xind = cp.xmlist.index(x0)
        yind = cp.ymlist.index(y0)
        zout = (xs.z.max()-(cp.dy/2))/cp.dy + cp.dy/10
        cp.indicesreal = cp.indicesreal.append({'x':xind, 'y':yind+zout}, ignore_index=True)

    xs['y'] = xs['y']+x0  # shift the points relative to the point on the plot where this sim's origin should be
    xs['z'] = xs['z']+y0
    try:
        xmid,ymid = intm.pdCentroid(xs)
    except:
        return
    if cp.split:
        ax = cp.axs[sigmapos]
    else:
        ax = cp.axs[0]
    if 'color' in kwargs:
        color = kwargs['color']
    plotXSOnAx(xs, ax, color)
    
    le = fp.legendUnique(folder)
    if 'adjacent_filament_offset' in le and float(le['adjacent_filament_offset'])>0:
        # neighboring filaments
        if float(le['ink_velocity'])>0 and xs.x.max()>float(le['nozzle_center_x_coord']):
            # with extrusion
            if le['adjacent_filament_orientation']=='y':
                x0 = x0-float(le['adjacent_filament_offset'])/2
            else:
                y0 = y0-float(le['adjacent_filament_offset'])/2
        else:
            # no extrusion
            if le['adjacent_filament_orientation']=='y':
                x0 = x0-float(le['adjacent_filament_offset'])
            else:
                y0 = y0-float(le['adjacent_filament_offset'])
    if recenter:
        if le['adjacent_filament_orientation']=='y':
            x0 = x0+float(le['adjacent_filament_offset'])
        else:
            y0 = y0-float(le['adjacent_filament_offset'])
        
    ax.arrow(x0, y0, xmid-x0, ymid-y0,  head_width=0.05, head_length=0.1, fc=color, ec=color, length_includes_head=True)
    
def putAbove(cp:comboPlot) -> Tuple[int,int,float,float]:
    '''get the indices and positions for the ideal plot to put the plot above the top left corner'''
    xind = int(cp.indicesreal.x.min())
    x0 = cp.xmlist[xind]
    yind = int(cp.indicesreal[cp.indicesreal.x==xind].y.max())
    y0 = cp.ymlist[yind]+cp.dy
    yind = yind+1
    return xind, yind, x0, y0

def putLeft(cp:comboPlot) -> Tuple[int,int,float,float]:
    '''get the indices to put the plot left of the bottom left corner'''
    xind = int(cp.indicesreal.x.min())
    x0 = cp.xmlist[xind]-cp.dx
    xind = xind-1
    yind = int(cp.indicesreal.y.max())
    y0 = cp.ymlist[yind]
    return xind, yind, x0, y0
    

def XSPlotIdeal(cp:comboPlot, le:dict, **kwargs) -> None:
    '''this plots an ideal cross-section on the cross-section plot
    le is a dictionary holding metadata, from legendUnique'''
    if 'adjacent_filament_offset' in le and float(le['adjacent_filament_offset'])>0:
        # adjacent filaments
        if 'adjacent_filament_offset' in cp.xvar:
            # put above
            x = cp.xfunc(extractTP(os.path.join(cp.topFolder, le['folder']), units=False))
            xind = findPos(cp.xlist, x)
            x0 = cp.xmlist[xind]
            yind = cp.indicesreal[cp.indicesreal.x==xind].y.min()
            if not pd.isna(yind):
                yind = int(yind)
            else:
                yind = int(cp.indicesreal.y.min())
            y0 = cp.ymlist[yind]-cp.dy
            yind = yind-1
        elif 'adjacent_filament_offset' in cp.yvar:
            # put left
            xind,yind,x0,y0 = putLeft(cp)
            y = cp.yfunc(extractTP(os.path.join(cp.topFolder, le['folder']), units=False))
            yind = findPos(cp.ylist, y)
            y0 = cp.ymlist[yind]
        else:
            if len(cp.ylistreal)==1:
                # put to the left
                xind,yind,x0,y0 = putLeft(cp)
            else:
                # put above
                xind,yind,x0,y0 = putAbove(cp)
    else:
        # single filaments
        if len(cp.ylistreal)==1:
            # put to the left
            xind,yind,x0,y0 = putLeft(cp)
        else:
            # put above
            xind,yind,x0,y0 = putAbove(cp)
    color1='Black'
    color2 = 'Gray'
    
    if 'adjacent_filament_offset' in le and float(le['adjacent_filament_offset'])>0:
        # plot 2 filaments
        plotCircle(cp.axs[0], x0, y0, float(le['nozzle_inner_width'])/2, '', color2)
        if le['adjacent_filament_orientation']=='y':
            plotCircle(cp.axs[0], x0-float(le['adjacent_filament_offset']), y0, float(le['nozzle_inner_width'])/2, '', color1)
        else:
            plotCircle(cp.axs[0], x0, y0-float(le['adjacent_filament_offset']), float(le['nozzle_inner_width'])/2, '', color1)
    else:
        # plot 1 filament
        plotCircle(cp.axs[0], x0, y0, float(le['nozzle_inner_width'])/2, 'Ideal', color1)
    cp.indicesreal = cp.indicesreal.append({'x':xind, 'y':yind}, ignore_index=True)


def XSPlotf(folder:str, time:float, xbehind:float, cp:comboPlot, xunits:str='mm', ref:dict={'nozzle_inner_width':0, 'ink_velocity':0, 'bath_velocity':0}, recenter:bool=False, **kwargs) -> None:
    '''plots cross-section from one file
    time is the time of the cross-section in seconds
    xbehind is the distance behind the nozzle to take the cross-section, in xunits
    cp is the comboPlot object to plot on
    xunits is 'mm' or 'nozzle_inner_width'
    ref holds metadata about a reference simulation. leave values at 0 to not do any scaling. Otherwise, scale the size of the cross-section relative to the reference values
    '''
    ptsx = intm.importPtsSlice(folder, time, xbehind, xunits=xunits)
    if len(ptsx)==0:
        return
    if float(ref['nozzle_inner_width'])>0:
        # rescale the points so all cross-sections reference the same 
        le = fp.legendUnique(folder)
        if not ('adjacent_filament_offset' in le and float(le['adjacent_filament_offset'])>0):
            # normalize cross-section if not writing fusions
            dEst = float(le['nozzle_inner_width'])*np.sqrt(float(le['ink_velocity'])/float(le['bath_velocity'])) # ideal final diameter
            dEst0 = float(ref['nozzle_inner_width'])*np.sqrt(float(ref['ink_velocity'])/float(ref['bath_velocity'])) # ideal final diameter
            if dEst0>0:
                for s in ['y', 'z']:
                    ptsx[s] = [dEst0/dEst*i for i in ptsx[s]]
        else:
            if recenter:
                if le['adjacent_filament_orientation']=='y':
                    ptsx['y'] = [i+float(le['adjacent_filament_offset']) for i in ptsx['y']]
                else:
                    ptsx['z'] = [i-float(le['adjacent_filament_offset']) for i in ptsx['z']]
    if len(ptsx)>0:
        XSPlot(ptsx, folder, cp, recenter=recenter, **kwargs)

def xbStr(xbehind:Union[List[float], float, dict]) -> str:
    '''convert xbehind to a readable string'''
    if type(xbehind) is dict:
        xbstr =''
        for key,val in xbehind.items():
            xbstr = xbstr+os.path.basename(key)+str(val)+'_'
    else:
        xbstr = str(xbehind)
    return xbstr

def XSPlots0(topFolder:str, exportFolder:str, time:float, xbehind:Union[List[float], float], xunits:str='mm', overwrite:bool=False, dx:float=0.5, **kwargs) -> None:
    '''plot all cross-sections together
    topFolder is the folder that holds all the files
    time is the time at which to take the cross-section
    xbehind is the distance behind the center of the nozzle to take the cross-section
    xunits is 'mm' or 'nozzle_inner_width'
    overwrite True to overwrite values
    dx is the spacing between cross-sections, in mm
    
    '''
    label = f'xs_{xbStr(xbehind)}{xunits}_t_{time}'
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return

    if 'dy' in kwargs:
        dy = kwargs['dy']
    else:
        dy = dx
    if 'imsize' in kwargs:
        imsize = kwargs['imsize']
        kwargs.pop('imsize')
    else:
        imsize = 6.5
    cp = comboPlot(topFolder, [-dx, dx], [-dy, dy], imsize, **kwargs)
    
    if len(cp.flist)==0:
        logging.warning("No files in list")
        return
    
    if type(xbehind) is list or type(xbehind) is dict:
        colors = ['#940f0f', '#61bab0', '#1e5a85']
    
    (cp.flist).sort(key=lambda folder:extractTP(folder)['sigma']) # sort folders by sigma value so they are stacked in the right order
    le0 = fp.legendUnique(cp.flist[0]) # use first file as reference dimensions
    if 'adjacent_filament_offset' in le0:
        idealFolders = {le0['adjacent_filament_offset']:le0}
    else:
        idealFolders = {0:le0}
    for folder in cp.flist:
        # identify if we need to add a new ideal plot
        if 'adjacent_filament_offset' in cp.xvar or 'adjacent_filament_offset' in cp.yvar:
            le = fp.legendUnique(folder)
            if not le['adjacent_filament_offset'] in idealFolders:
                idealFolders[le['adjacent_filament_offset']] = le
                
        # plot xs
        if type(xbehind) is list:
            # plot xs for multiple x positions
            for i,xb in enumerate(xbehind):
                XSPlotf(folder, time, xb, cp, xunits=xunits, ref=le0, color=colors[i], recenter=False)
        elif type(xbehind) is dict:
            # compare xs ahead of nozzle to xs behind nozzle in reference file
            i = 0
            le = fp.legendUnique(folder)
            for topf,xb in xbehind.items():
                if topf in folder:
                    # plot folder
                    XSPlotf(folder, time, xb, cp, xunits=xunits, ref=le0, color=colors[i], recenter=True)
                else:
                    # plot reference folder
                    reff = findRef(le, topf)
                    if os.path.exists(reff):
                        XSPlotf(reff, time, xb, cp, xunits=xunits, ref=le0, color=colors[i])
                    else:
                        logging.info(f'Reference file not found for {os.path.basename(folder)}')
                i = i+1
        else:
            # plot xs for single x position
            XSPlotf(folder, time, xbehind, cp, xunits=xunits, ref=le0)
    xunname = xunits.replace('nozzle_inner_width', '$d_i$')
    cp.figtitle = f'Cross sections, {xbStr(xbehind)} {xunname} behind nozzle, t = {time} s'
    
    for adj,le in idealFolders.items():
        XSPlotIdeal(cp,le)
    cp.clean()

    for ax in cp.axs:
        ax.grid(linestyle='-', linewidth='0.25', color='#949494')
        
    intm.exportIm(fn, cp.fig, **kwargs)
    
    
def findRef(le:dict, topFolder:str):
    '''find the reference folder'''
    p = os.path.join(topFolder, le['compare_to'].strip())
    if os.path.exists(p):
        return p
    else:
        for f1 in os.listdir(topFolder):
            f1full = os.path.join(topFolder, f1)
            if os.path.isdir(f1full) and not fp.isSimFolder(f1full) and not f1 in ['mesh', 'geometry']:
                r = findRef(le, f1full)
                if len(r)>0:
                    return r
        return ''

