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
import re

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from folder_stats import folderStats
from points.folder_points import folderPoints
import points.points_tools as pto
from points.slice_points import slicePoints
import tools.strings as st
import plot.plot_tools as plto
from plot.colors import plotColors
from plot.combo_plot import comboPlot

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

class XSPlot(comboPlot):
    '''plot all cross-sections together
    topFolder is the folder that holds all the files
    time is the time at which to take the cross-section
    xbehind is the distance behind the center of the nozzle to take the cross-section. if it's a list, take multiple positions. if it's a dict, take positions from different folders, where the topfolder is the key and the position is the val
    xunits is 'mm' or 'nozzle_inner_width'
    overwrite True to overwrite plots
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , time:float
                 , xbehind:Union[List[float], float]
                 , xunits:str='mm'
                 , rescale:bool=False
                 , plotNoz:bool=False
                 , plotArrows:bool=False
                 , **kwargs):
        self.time = time
        self.xbehind = xbehind
        self.xunits = xunits
        self.rescale = rescale
        self.plotNoz = plotNoz
        self.plotArrows = plotArrows
        kwargs = self.getCvar(**kwargs)
        super().__init__(topFolder, exportFolder=exportFolder, **kwargs)
        self.getFN(**kwargs)
        if not self.checkOverwrite(export=self.export, overwrite=self.overwrite):
            return
        self.setColors()   
            # determine if we should change the color scheme from default
        self.idealFolders = [[{} for i in range(self.ncol)]  for j in range(self.nrow)] 
            # folders to use for creating ideal reference circles. one dict for each axis
        self.fPoints = {}   
            # dictionary of folderPoints objects
        
        # iterate through folders and plot data
        for i,row in self.filedf.iterrows():
            self.plotFolder(row)
            
        self.addIdeals()
        xunname = xunits.replace('nozzle_inner_width', '$d_i$')
        self.figtitle = f'Cross sections, {self.xbStr()} {xunname} behind nozzle, t = {self.time} s'
        self.clean()
            
        if self.export:
            self.exportIm(**kwargs) 
                
    #-------------------
                        
    def getCvar(self, **kwargs) -> None:
        '''determine a colorvar to feed to super.init and a color var to actually use'''
        if 'cvar' in kwargs:
            self.cvarreal = kwargs['cvar']
            spl = re.split(',', self.cvarreal)
            if 'xbehind' in spl:
                self.cxbehindPos = spl.index('xbehind')   # store where the xbehind value is in the color values
                spl.remove('xbehind')
                if len(spl)==0:
                    kwargs['cvar'] = 'bn'
                else:
                    kwargs['cvar'] = ','.join(spl)
            else:
                return kwargs
            if 'colorDict' in kwargs:
                self.colorDictReal = kwargs['colorDict'].copy()
                newColorDict = {}
                for key,val in kwargs['colorDict'].items():
                    spli = re.split(', ', key)
                    if self.cxbehindPos>0:
                        if self.cxbehindPos<len(spli):
                            newkey = ', '.join(spli[:self.cxbehindPos]+spli[self.cxbehindPos+1:])
                        else:
                            newkey = ', '.join(spli[:self.cxbehindPos])
                    else:
                        if self.cxbehindPos<len(spli):
                            newkey = ', '.join(spli[self.cxbehindPos+1:])
                        else:
                            newkey = 'bn'
                    newColorDict[newkey] = val
                kwargs['colorDict'] = newColorDict
        else:
            self.cvarreal = ''
        return kwargs
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'xs_{self.xbStr()}{self.xunits}_t_{self.time}'
        
    def xbStr(self) -> str:
        '''convert xbehind to a readable string'''
        if type(self.xbehind) is dict:
            xbstr =''
            for key,val in self.xbehind.items():
                xbstr = f'{xbstr}{os.path.basename(key)}{val}_'
        else:
            xbstr = str(self.xbehind)
        return xbstr
    
    #-------------------
    
    def newColorVal(self, cval0:str, xbehind:float) -> str:
        '''get the new color val'''
        if 'xbehind' in self.cvarreal:
            spli = re.split(', ', cval0)
            if self.cxbehindPos>0:
                if self.cxbehindPos<len(spli):
                    newval = ', '.join(spli[:self.cxbehindPos]+[str(xbehind)]+spli[self.cxbehindPos+1:])
                else:
                    newval = ', '.join(spli[:self.cxbehindPos]+[str(xbehind)])
            else:
                if self.cxbehindPos<len(spli):
                    newval = ', '.join([str(xbehind)]+spli[self.cxbehindPos+1:])
                else:
                    newval = str(xbehind)
            return newval
        else:
            return cval0
    
    def setColors(self):
        '''change the color system if we're coloring by xbehind'''
        if type(self.xbehind) is list:
            # list of x positions
            clist = self.xbehind
        elif type(self.xbehind) is dict:
            # dictionary of x positions in different folders
            clist = list(self.xbehind.values())
        else:
            # use standard coloring
            return
        if not 'xbehind' in self.cvarreal and self.cvar=='':
            self.cvar = 'xbehind'
            colors = ['Black', '#cf2213']
            self.colors = plotColors(clist, 'xbehind', f'x ({self.xunits})', colorList = colors)
            self.filedf.loc[:, 'cvar'] = self.xbehind[0]  # this just creates an empty value
        else:
            vallist = [self.newColorVal(c,b) for c in self.clist for b in clist]
            if hasattr(self, 'colorDictReal'):
                self.kwargs['colorDict'] = self.colorDictReal
            self.colors = plotColors(vallist, self.cvarreal, f'x ({self.xunits})', **self.kwargs)
        
    
    #-------------------
        
    def addToSpacingLegend(self, fs:folderStats, row:pd.Series) -> None:
        '''determine if we need a new entry in the spacing legend'''
        if hasattr(fs.geo, 'spacing'):
            spacing = fs.geo.spacing
            j = row['splityindex']
            i = row['splitxindex']
            if not spacing in self.idealFolders[j][i]:
                self.idealFolders[j][i][spacing] = fs
      
    #-------------------
    
    def plotNozzle(self, row:pd.Series, pos:dict, scale:float=1, recenter:bool=False) -> None:
        '''plot the nozzle'''
        if not self.plotNoz:
            return
        fs = self.fstats[row['folder']]
        now = fs.geo.now*scale
        niw = fs.geo.niw*scale
        bot = fs.geo.nbz*scale
        x0 = pos['x0']
        y0 = pos['y0']
        dx = 0
        dy = 0
        if recenter:
            if fs.geo.Ddir=='y':
                dy = fs.geo.dDis
            else:
                dx = fs.geo.dDis
        ybot = y0+dy+bot
        ytop = ybot+0.4
        xleft = x0+dx-now/2
        xright = x0+dx+now/2
        xileft = x0+dx-niw/2
        xiright = x0+dx+niw/2
        pos['ax'].fill([xleft, xleft, xright, xright, xiright, xiright, xileft, xileft, xleft], [ytop, ybot, ybot, ytop, ytop, ybot, ybot, ytop, ytop], color='#d1d1d1', linewidth=0)
        
    def plotArrow(self, row:pd.Series, sp:slicePoints, pos:dict, recenter:bool=False) -> None:
        '''plot the arrow from the intended position to the centroid'''
        if not self.plotArrows:
            return
        fs = self.fstats[row['folder']]
        sp.plotArrow(pos, recenter, fs)

    def plotPoints(self, row:pd.Series, sp:slicePoints, xbehind:float, recenter:bool=False, **kwargs) -> None:
        '''we have a list of points, now we need to shift them to the right position on the axis'''  
        
        pos = self.getXYRow(row)
        if 'xbehind' in self.cvarreal:
            # change the color if we're coloring by xbehind
            val = self.newColorVal(row['cvar'], xbehind)
            pos['color'] = self.colors.getColor(val)
        x0 = pos['x0']
        y0 = pos['y0']
        
        scale = 1
        dx = 0
        dy = 0
        if self.rescale:
            # rescale the points so all cross-sections reference the same 
            if not hasattr(fs.geo, 'spacing'):
                # normalize cross-section if not writing fusions
                dEst = fs.geo.dEst # ideal final diameter
                dEst0 = self.fs1.geo.dEst # ideal final diameter
                if dEst0>0:
                    scale = dEst0/dEst
            else:
                if recenter:
                    if fs.geo.Ddir=='y':
                        dy = fs.geo.dDis
                    else:
                        dx = fs.geo.dDis

        if sp.df.z.max()>self.dy/2:
            # xs is going to fall outside of the plot bounds. expand the bounds
            xind = self.xmlist.index(x0)
            yind = self.ymlist.index(y0)
            zout = (sp.df.z.max()-(self.dy/2))/self.dy + self.dy/10
            self.indicesreal = self.indicesreal.append({'x':xind, 'y':yind+zout, 'ax':pos['ax'], 'color':pos['color']}, ignore_index=True)
            
        sp.shiftPoints(scale, x0+dx, y0+dy)   # shift and scale the points relative to the point on the plot where this sim's origin should be
        self.plotNozzle(row, pos, scale, recenter=recenter)
        sp.plotSlice(pos['ax'], pos['color'])
        self.plotArrow(row, sp, pos, recenter=recenter)
                
    def plotSlice(self, row:pd.Series, xbehind:float, recenter:bool=False, **kwargs) -> None:
        '''plots cross-section from one file
        xbehind is the distance behind the nozzle to take the cross-section, in xunits
        '''
        fs = self.fstats[row['folder']]
        if fs.folder in self.fPoints:
            fp = self.fPoints[fs.folder]
        else:
            fp = folderPoints(fs)
            self.fPoints[fs.folder] = fp
        sp = fp.importPtsSlice(self.time, xbehind, xunits=self.xunits)
        if len(sp.df0)==0:
            return
        

        self.plotPoints(row, sp, xbehind, recenter=recenter, **kwargs)
        
    def findRef(self, fs:folderStats, topFolder:str) -> str:
        '''find the reference folder'''
        p = os.path.join(topFolder, fs.compare_to)
        if os.path.exists(p):
            return p
        else:
            for f1 in os.listdir(topFolder):
                f1full = os.path.join(topFolder, f1)
                if os.path.isdir(f1full) and not fp.isSimFolder(f1full) and not f1 in ['mesh', 'geometry']:
                    r = self.findRef(fs, f1full)
                    if len(r)>0:
                        return r
            return ''
        
    def plotFolder(self, row) -> None:
        '''given a row in the pandas dataframe, plot the slices'''
        # identify if we need to add a new ideal plot
        fs = self.fstats[row['folder']]
        # determine if we need a new legend
        self.addToSpacingLegend(fs, row)
        # plot xs
        if type(self.xbehind) is list:
            # plot xs for multiple x positions
            for i,xb in enumerate(self.xbehind):
                self.plotSlice(row, xb, recenter=False)
        elif type(self.xbehind) is dict:
            # compare xs ahead of nozzle to xs behind nozzle in reference file
            i = 0
            for topf,xb in self.xbehind.items():
                if topf in folder:
                    # plot folder
                    self.plotSlice(row, xb, recenter=True)
                else:
                    # plot reference folder
                    reff = self.findRef(fs, topf)
                    if os.path.exists(reff):
                        row2 = row.copy()
                        row2['folder'] = reff
                        self.plotSlice(row2, xb)
                    else:
                        logging.info(f'Reference file not found for {os.path.basename(folder)}')
                i = i+1
        else:
            # plot xs for single x position
            self.plotSlice(row, xb)  

    
    #----------------------------------------------------------
    
    def idealPlotSingle(self, fs:folderStats, i:int, j:int) -> Tuple[int, int]:
        '''plot a single ideal filament cross-section'''
        # single filaments
        if len(self.ylistreal)==1:
            # put to the left
            xind,yind,x0,y0 = self.putLeft(i, axrow=j)
        else:
            # put above
            xind,yind,x0,y0 = self.putAbove(i, axrow=j)
        color1='Black'
        plto.plotCircle(self.axs[j][i], x0, y0, fs.geo.niw/2, 'Ideal', color1)
        return xind, yind, x0, y0
    
    def idealPlotAdj(self, fs:folderStats, i:int, j:int) -> Tuple[int, int]:
        '''plot two ideal filament cross-sections, overlapped appropriately'''
        ir = self.indicesreal[(self.indicesreal.axx==i)&(self.indicesreal.axy==j)]
        if self.vn.hasSpacing(self.xvar):
            # put above
            pos = self.getPos(fs)
            xind = pos['xindex']
            x0 = self.xmlist[xind]
            yind = ir[ir.x==xind].y.min()
            if not pd.isna(yind):
                yind = int(yind)
            else:
                yind = int(ir.y.min())
            y0 = self.ymlist[yind]-self.dy
            yind = yind-1
        elif self.vn.hasSpacing(self.yvar):
            # put left
            xind,yind,x0,y0 = self.putLeft()
            pos = self.getPos(fs)
            yind = pos['yindex']
            y0 = self.ymlist[yind]
        else:
            if len(self.ylistreal)==1:
                # put to the left
                xind,yind,x0,y0 = self.putLeft(i)
            else:
                # put above
                xind,yind,x0,y0 = self.putAbove(i)
                
        color1 = 'Black'
        color2 = 'Gray'
        
        plto.plotCircle(self.axs[j][i], x0, y0, fs.geo.niw/2, '', color2)
        if fs.geo.Ddir=='y':
            plto.plotCircle(self.axs[j][i], x0-fs.geo.Ddis, y0, fs.geo.niw/2, '', color1)
        else:
            plto.plotCircle(self.axs[j][i], x0, y0-fs.geo.Ddis, fs.geo.niw/2, '', color1)
            ybot = -fs.geo.Ddis-fs.geo.niw/2
            if ybot<-self.dy/2:
                self.indicesreal = self.indicesreal.append({'x':xind, 'y':ybot/(self.dy/2) - self.dy/10, 'axx':i, 'axy':j}, ignore_index=True)
                
        
        return xind, yind, x0, y0     
        
    def XSPlotIdeal(self, fs:folderStats, i:int, j:int, **kwargs) -> None:
        '''this plots an ideal cross-section on the cross-section plot
        fs is a folderStats object to use as a reference to find ideal position
        i is the axis column
        j is the axis row
        '''
        if hasattr(fs.geo, 'spacing'):
            xind, yind, x0, y0 = self.idealPlotAdj(fs, i, j)
        else:
            xind, yind, x0, y0 = self.idealPlotSingle(fs, i, j)
        if not x0 in self.xmlist:
            self.xmlist.append(x0)
            self.xlist.append('Ideal')
        if not y0 in self.ymlist:
            self.ymlist.append(y0)
            self.ylist.append('Ideal')
        self.indicesreal = self.indicesreal.append({'x':xind, 'y':yind, 'axx':i, 'axy':j}, ignore_index=True)
        
        
    def addIdeals(self):
        '''add illustrations for scale and position of ideal cross-sections'''
        for i in range(self.ncol):
            for j in range(self.nrow):
                if len(self.idealFolders[j][i])==0:
                    self.idealFolders[j][i] = {0:self.firstFS()}
                for fs in self.idealFolders[j][i].values():
                    self.XSPlotIdeal(fs, i, j)