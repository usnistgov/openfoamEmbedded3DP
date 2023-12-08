#!/usr/bin/env python
'''Collecting points from folders'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import file.plainIm as pi
import tools.strings as st
from folder_stats import geoStats, folderStats
import points.points_tools as pto
from add_units import addUnits
from points.slice_points import slicePoints

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#################################################################

class folderPoints:
    '''for extracting points from a folder'''
    
    def __init__(self, fs:Union[str, folderStats]):
        if type(fs) is str:
            fs = folderStats(fs)
        self.fs = fs
        self.folder = fs.folder
        self.geo = self.fs.geo
        self.interfacePointsFiles = {}
        self.nozzlePointsFiles = {}
        self.nozzleSlicePointsFiles = {}
        self.interfacePointsSliceFiles = {}
        
    #--------------------------------------
    
    def importPointsFile(self, file:str) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
        '''This is useful for importing interfacePoints.csv files. '''
        d,units = pi.plainIm(file, False)
        if len(d)==0:
            return d,units
        try:
            d = d[pd.to_numeric(d.time, errors='coerce').notnull()] # remove non-numeric times
        except:
            pass
        mdict = {'m':'mm', 'm/s':'mm/s'}
        for s in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'arc_length', 'magu']:
            if s in d:
                try:
                    d[s] = d[s].astype(float)
                except:
                    raise Exception(f'Non-numeric values for {s} in {file}')
                else:
                    d[s]*=1000
                units[s]=mdict.get(units[s], units[s]+'*10^3')
        d = d.sort_values(by='x')
        return d,units
    
    def getExistingPoints(self, name:str, time:float):
        '''check if we already imported this file, and if not, import and store the values'''
        if time in getattr(self, f'{name}Files'):
            d, units = getattr(self, f'{name}Files')[time]
            return d.copy(), units.copy()
        d,units = self.importPointsFile(self.fn(name, time))
        getattr(self, f'{name}Files')[time] = (d,units)
        return d,units


    def fn(self, name:str, time:float) -> str:
        '''get a file name from the name of the folder and the time in seconds'''
        bn = f'{name}_t_{int(round(time*10))}.csv'
        return os.path.join(self.folder, name, bn)
    
    def importInterfacePoints(self, time:float) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
        '''import all points in the interface at a given time'''
        return self.getExistingPoints('interfacePoints', time)

    def importNozzlePoints(self, time:float) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]: # RG
        '''import all points in the center y slice of the nozzle at a given time'''
        return self.getExistingPoints('nozzlePoints', time)

    def importNozzleSlicesPoints(self, time:float) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
        '''import all points in the center y slice of the nozzle at a given time'''
        return self.getExistingPoints('nozzleSlicePoints', time)
    
    
    #--------------------------------
    
    def tryImportPtsSlice(self, time:float, xbehind:float, xunits:str='mm') -> Union[pd.DataFrame, List[Any]]:
        '''try to import the points from an existing file'''
        # try to import the points
        xu = xunits.replace('nozzle_inner_width', 'niw')
        bn = f'interfacePoints_t_{int(round(time*10))}_x_{xbehind}_{xu}.csv'
        ipxfolder = os.path.join(self.folder, 'interfacePointsX')
        fn = os.path.join(ipxfolder, bn)
        if os.path.exists(fn):
            pts, units= pi.plainIm(fn, ic=0, checkUnits=True)
            return pts, fn
        if not os.path.exists(ipxfolder):
            os.mkdir(ipxfolder)
        return [], fn
    
    def makeRelativeX(self, pts:pd.DataFrame) -> pd.DataFrame:
        '''make the x position relative to the nozzle center'''
        xc = self.geo.nozzle_center_x_coord
        pts['x'] =  [i-xc for i in pts['x']]
        return pts
    
    def convertXunits(self, pts:pd.DataFrame, units:dict, xunits:str, xvar:str='x') -> Tuple[pd.DataFrame, dict]:
        '''convert the x points to the units given by xunits, but leave the rest in mm'''
        if not xunits=='mm' and units[xvar]=='mm':
            # convert the x units
            if xunits=='niw' or xunits=='nozzle_inner_width':
                pts[xvar] = pts[xvar]/self.geo.niw
                units[xvar] = 'niw'
            else:
                val = self.fs.getVal(xunits)
                if type(val) is float and not np.isnan(val):
                    pts[xvar] = pts[xvar]/val
                    units[xvar] = xunits
        return pts, units
    
    def selectSlice(self, pts:pd.DataFrame, xbehind:float) -> pd.DataFrame:
        '''select just the points closest to the xbehind position'''
        xlist = pto.xpts(pts)
        xreal = pto.closest(xlist, xbehind)
        
        if abs(xreal-xbehind)>0.2:
            return []
        ptsx = pts[pts['x']==xreal]
        return ptsx.copy()

    def importPtsSlice(self, time:float, xbehind:float, xunits:str='mm') -> slicePoints:
        '''import points just from one slice. 
        time is in s
        xbehind is distance behind nozzle in xunits. 
        returns a slicePoints object, which holds points, representations as polygons, and functions to find centroids etc.
        Finds the closest position to the requested position and gives up if there is no x value within 0.2 xunits'''
        pts, fn = self.tryImportPtsSlice(time, xbehind, xunits)
        if len(pts)>0:
            return slicePoints(pts)
        
        # import all points at that time
        pts,units = self.importInterfacePoints(time)
        if len(pts)==0:
            return slicePoints([])

        # use relative coordinates
        pts = self.makeRelativeX(pts)
        pts, units = self.convertXunits(pts, units, xunits)
        ptsx = self.selectSlice(pts, xbehind)
        if len(ptsx)==0:
            return slicePoints([])
        sp = slicePoints(ptsx)
        units['theta'] = 'rad'
        units['segment'] = ''
        pi.plainExp(fn, sp.combinedf(), units)    # export the points for next time
        return sp
    
    #------------------------------------------

    def importSummarySlice(self, time:float, xbehind:float, xunits:str='mm') -> Tuple[pd.DataFrame, dict]:
        '''import measurements just from one slice. 
        time is in s
        xbehind is distance behind nozzle in xunits. 
        returns a slicePoints object, which holds points, representations as polygons, and functions to find centroids etc.
        Finds the closest position to the requested position and gives up if there is no x value within 0.2 xunits'''
        pts,u = self.importSummary()
        if len(pts)==0:
            return pts,u
        pts = pts[pts.time==time]
        if len(pts)==0:
            return pts,u
        
        # use relative coordinates
        pts = self.makeRelativeX(pts)
        pts, u = self.convertXunits(pts, u, xunits)
        ptsx = self.selectSlice(pts, xbehind)
        return ptsx, u
    
    def importSummary(self) -> Tuple[pd.DataFrame, dict]:
        '''get the summary file'''
        fn = os.path.join(self.folder, 'sliceSummaries.csv')
        if not os.path.exists(fn):
            return pd.DataFrame([]), {}
        pts, u = pi.plainIm(fn)
        return pts,u
        
    
    #--------------------------------------
    
    def importFilemm(file:str, slist:List[str]) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
        '''Import a file and convert m to mm. File is a full path name. slist is the list of column names for the columns that need mm conversion'''
        d,units = plainIm(file, False)
        mdict = {'m':'mm', 'm/s':'mm/s'}
        if len(d)==0:
            return d,units
        else:
            for s in slist:
                d[s.lower()]*=1000
                units[s]=mdict.get(units[s], units[s]+'*10^3')
            return d,units


    def importLine(self, time:float, x:float=1.4, xunits:str='mm', **kwargs) -> pd.DataFrame:
        '''import a csv of a line trace pulled from ParaView. 
        time and x are the time and position the line trace are taken at. 
        xunits can be mm or nozzle_inner_width'''
        if xunits=='mm':
            file = os.path.join(self.folder, f'line_t_{int(round(time*10))}_x_{x}.csv')
        else:
            xform = '{:.1f}'.format(x)
            file = os.path.join(self.folder, f'line_t_{int(round(time*10))}_x_{xform}_di.csv')
        try:
            data, units = self.importPointsFile(file) 
        except:
            addUnits(file)
            data, units = self.importPointsFile(file) 
        return  data, units
    
    #---------------------------------------
    
    def posSlice(self, loc:float, xunits:str='mm') -> pd.DataFrame: # RG
        '''go through an interface file and take the positions of the interface points for a given x behind the nozzle
        outputs a pandas DataFrame'''
        ipfolder = os.path.join(self.folder, 'interfacePoints')
        if not os.path.exists(ipfolder):
            raise Exception('No interface points')
        sp = self.importPtsSlice(2.5, loc, xunits=xunits)
        d = sp.df
        d = d.sort_values(by='z')
        pos = d[['x', 'y', 'z']]
        return pos
