#!/usr/bin/env python
'''Collecting points from folders'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
import cv2 as cv
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import file.plainIm as pi
import tools.strings as st
from folder_stats import geoStats, folderStats
from add_units import addUnits

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#---------------------------------------------------------

class segmentPoints:
    '''for holding points from a single segmented object in a slice'''
    
    def __init__(self, df:pd.DataFrame, smoothing:float=0.03):
        self.df = df.copy()
        self.smoothing = smoothing
        
    def angle360(self, row:pd.Series, cx:float, cy:float) -> float:
        '''get the angle between the vector from center to row[y,z], relative to the x axis'''
        v = [row['y']-cx,row['z']-cy]
        if v==[0,0]:
            return 0
        angle = np.arccos(np.dot(v, [1,0])/(np.linalg.norm(v)))
        if v[1]<0:
            angle = 2*np.pi-angle
        return angle*360/(2*np.pi)
    
    def sortRadially(self, cx:float='', cy:float='', closeGaps:bool=True) -> None:
        '''sort the list of points radially so we can plot this as a line instead of points (saves space in vector images)'''
        if cx=='' or cy=='':
            # find center
            left = self.df.y.min()
            right = self.df.y.max()
            bottom = self.df.z.min()
            top = self.df.z.max()
            cx = (left+right)/2
            cy = (top+bottom)/2
        df2 = self.df.copy()
        df2['theta'] = [self.angle360(row, cx, cy) for _,row in df2.iterrows()]  # calculate angle relative to center
        df2.sort_values(by='theta', inplace=True, ignore_index=True)             # sort by angle relative to center
        
        if closeGaps:
            df2['dtheta'] = df2['theta'].diff()
            m = df2.dtheta.max()
            if m>5:
                # rearrange the points so they go from one endpoint to another endpoint
                row = df2[df2.dtheta==m]
                i = row.iloc[0].name
                df2.loc[:i, 'theta'] = df2.loc[:i, 'theta']+360
                df2.sort_values(by='theta', inplace=True, ignore_index=True)
            df2.drop(columns=['dtheta'], inplace=True)  # don't need this column anymore
        self.df = df2
        
    def sortVal(self, var:str, inc:bool) -> None:
        '''sort the points by y value. inc=True to sort increasing, otherwise decreasing y'''
        self.df.sort_values(by=var, inplace=True, ignore_index=True, ascending=inc)
        
    def polygon(self):
        '''get the shapely polygon for this segment'''
        if not hasattr(self, 'poly'):
            self.resetPoly()
        return self.poly
        
    def resetPoly(self):
        '''reset the polygon object'''
        self.poly = Polygon(self.df[['y','z']].values.tolist())
        if self.smoothing>0:
            self.poly = self.poly.buffer(self.smoothing, join_style=1).buffer(-self.smoothing, join_style=1)
        
    
    def shiftPoints(self, scale:float=1, dx:float=0, dy:float=0):
        '''scale the points and shift the points in the combined dataframe'''
        if not scale==1:
            self.df.loc[:,'y'] = self.df.y*scale
            self.df.loc[:,'z'] = self.df.z*scale
        if not dx==0:
            self.df.loc[:,'y'] = self.df.y + dx
        if not dy==0:
            self.df.loc[:,'z'] = self.df.z + dy
            
    def plotSlice(self, ax, color):
        '''plot the slice on an existing axis'''
        #x,y = self.poly.exterior.xy
        x = self.df.y
        y = self.df.z
        ax.plot(x,y, color=color, linewidth=0.75, marker=None)


class slicePoints:
    '''for holding points that belong to slices. dcrit is the distance between points when they can be in the same group in the first step, and ccrit is the distance between endpoints that will allow two segments to be combined in the second step'''
    
    def __init__(self, df:pd.DataFrame, dcrit:float=0.05, ccrit:float=0.05, smoothing:float=0.01):
        self.df0 = df
        if len(self.df0)==0:
            return
        self.dcrit = dcrit
        self.ccrit = ccrit
        self.smoothing = smoothing
        
        if 'segment' in df:
            # this dataframe is pre-sorted
            self.df = self.df0.copy()
            sortRadial=False
            self.smoothing = 0
        else:
            # we need to label the segments
            self.sortPoints()
            sortRadial=True
        self.splitPoints(sortRadial)
        # self.plot()
            
    def sortPoints(self) -> None:
        '''sort the dataframe into segments'''
        # create an initial combined segmentsPoints to sort the points radially
        self.labelSegments()
        
    def splitPoints(self, sortRadial:bool, **kwargs) -> None:
        ''' create a unique segmentPoints object for each  segment'''
        self.splist = []
        for segment in self.df.segment.unique():
            df1 = self.df[self.df.segment==segment]
            if len(df1)>10:
                sp = segmentPoints(df1, smoothing=self.smoothing)
                if sortRadial:
                    sp.sortRadially()   # sort the points radially for correct polygon rendering, calculations
                self.splist.append(sp)
                
    def manualSplit(self, refPts:dict) -> None:
        '''split the points based on their proximity to the reference points'''
        for i,row in self.df0.iterrows():
            ds = {}
            for segment, points in refPts.items():
                dist = min([(row['y']-p2[0])**2+(row['z']-p2[1])**2 for p2 in points])
                ds[dist] = segment
            s = ds[min(ds.keys())]
            self.df.loc[i, 'segment'] = s
            
    def labelNozzle(self, xL:float, xR:float, yB:float, margin:float) -> None:
        '''create new labels for the segments that are attached to the nozzle'''
        n = self.df.segment.max()+1
        self.df.loc[(self.df.y>xL-margin)&(self.df.y<xL+margin)&(self.df.z>yB-margin), 'segment'] = n  # label left
        self.df.loc[(self.df.y>xR-margin)&(self.df.y<xR+margin)&(self.df.z>yB-margin), 'segment'] = n+1  # label right
        self.df.loc[(self.df.y>xL-margin)&(self.df.y<xR+margin)&(self.df.z>yB-margin)&(self.df.z<yB+margin), 'segment'] = n+2  # label bottom
        
            
    def flatten(self) -> None:
        '''relabel all segments as the same segment using the slicePoints dataframes'''
        self.df = pd.concat([sp.df.loc[:len(sp.df)-2] for sp in self.splist])
        self.selfFlatten()
        
    def selfFlatten(self) -> None:
        '''flatten all points onto the same segment using this dataframe'''
        self.df.loc[:, 'segment'] = 1
        self.df.dropna(inplace=True)
        
    def addToCombiningList(self, l:list, currEp:list, endpoints:dict) -> Tuple[list, list, dict]:
        '''add the next splitPoints to the list of points'''
        last = currEp[1]  # get the last point in the list
        distances = dict([[self.ppdist(last, ep[0]), key] for key,ep in endpoints.items()]) # calculate distance to beginning points
        mindist = min(distances.keys())
        key = distances[mindist]
        l.append(self.splist[key].df.iloc[:len(self.splist[key].df)-2])    # add points to list
        currEp = endpoints.pop(key)
        return l,currEp,endpoints
        
    def combinedf(self) -> pd.DataFrame:
        '''get a dataframe that combines all of the sorted components'''
        endpoints = dict([[i,[sp.df.iloc[0], sp.df.iloc[-1]]] for i,sp in enumerate(self.splist)]) # get a dictionary of endpoints for each list
        l = [self.splist[0].df.iloc[:len(self.splist[0].df)-2]]  # start the list
        currEp = endpoints.pop(0)  # current endpoint
        while len(endpoints)>0:
            # add the next closest list
            l, currEp, endpoints = self.addToCombiningList(l, currEp, endpoints)
        self.df = pd.concat(l)
        return self.df
        
        
    def labelSegments(self) -> None:
        '''go through the points and sort them into separate segments by close contact'''
        if len(self.df0)==0:
            return
        sp = segmentPoints(self.df0)
        sp.sortRadially()
        self.df = sp.df.copy()
        self.df['segment'] = [-1 for i in range(len(self.df))]
        # go through the points and sort them into groups based on how far they are to each other
        self.lasts = {}
        self.firsts = {}
        self.nsegments = 0
        for i in range(len(self.df)):
            self.addToList(i)
        
        self.recombine()
        self.df.sort_values(by='segment')
            
    def ppdist(self, p1, p2) -> float:
        return np.sqrt((p1['y']-p2['y'])**2+(p1['z']-p2['z'])**2)
    
    def recombine(self):
        '''combine close groups of points'''
        if len(self.firsts)==0:
            return
        dists = pd.DataFrame([{'f':f, 'l':l, 'd':self.ppdist(self.firsts[f], self.lasts[l])} for f in self.firsts for l in self.lasts])  # find distances between firsts and lasts
        dists2 = pd.DataFrame([{'f':f, 'l':l, 'd':self.ppdist(self.lasts[f], self.lasts[l])} for f in self.firsts for l in self.lasts])   # find distances between lasts
        dists3 = pd.DataFrame([{'f':f, 'l':l, 'd':self.ppdist(self.firsts[f], self.firsts[l])} for f in self.firsts for l in self.lasts])   # find distances between firsts
        dists = pd.concat([dists, dists2, dists3])
        dists = dists[dists.d<self.ccrit]
        for i,row in dists.iterrows():
            l = row['l']
            f = row['f']
            if not l==f:
                self.df.loc[(self.df.segment==l), 'segment'] = f  # adopt the first segment label
    
    def createNewSegment(self, i:int, pt:dict):
        j = self.nsegments+1
        self.lasts[j] = pt
        self.firsts[j] = pt
        self.nsegments = j
        self.df.loc[i,'segment'] = j
    
    def addToList(self, i:int) -> None:
        '''add this point to one of the existing segments'''
        row = self.df.loc[i]
        pt = {'y':row['y'], 'z':row['z']}
        for j in self.lasts:
            ld = self.ppdist(pt, self.lasts[j])
            if ld<self.dcrit:
                # close to an endpoint
                self.lasts[j] = pt
                self.df.loc[i,'segment'] = j
                if not j in self.firsts:
                    self.firsts[j] = pt
                    self.nsegments = j
                return
            fd = self.ppdist(pt, self.firsts[j])
            if fd<self.dcrit:
                # close to a starting point
                self.firsts[j] = pt
                self.df.loc[i, 'segment'] = j
                return
        self.createNewSegment(i, pt)
        
    #-------------------------------------------------
        
    def polygon(self):
        '''get a polygon that is a combination of all of the polygons'''
        if not hasattr(self, 'poly'):
            self.resetPoly()
        return self.poly
        
    def resetPoly(self):
        '''create a new polygon object'''
        for sp in self.splist:
            sp.resetPoly()
        self.poly = unary_union([o.polygon() for o in self.splist])
        return self.poly

    def convexHull(self):
        if hasattr(self, 'ch'):
            return self.ch
        self.ch = self.polygon().convex_hull
        return self.ch
        
    def centroid(self)-> Tuple[float, float]:
        '''get the coordinates of a polygon centroid'''
        if hasattr(self, 'cx'):
            return self.cx, self.cy
        p = self.polygon()
        c = p.centroid.wkt
        slist = re.split('\(|\)| ', c)
        self.cx = float(slist[2])
        self.cy = float(slist[3])
        return self.cx,self.cy
    
    def area(self) -> float:
        '''get the area of the polygon'''
        if hasattr(self, 'a'):
            return self.a
        p = self.polygon()
        self.a = p.area
        return self.a
    
    def centroidAndArea(self) -> Tuple[float, float, float]:
        '''get the centroid and area of a slice, given as a list of points or a dataframe'''
        x,y=self.centroid()
        a = self.area()
        return x,y,a

    def roughness(self) -> float:
        '''get the excess perimeter relative to the perimeter of the convex hull'''
        p = self.polygon()
        ch = self.convexHull()
        return p.length/ch.length - 1
    
    def emptiness(self) -> float:
        '''how much of the convex hull is filled'''
        p = self.polygon()
        ch = self.convexHull()
        return 1-p.area/ch.area
    
    #----------------------------------------------
    
    def plot(self):
        fig,ax = plt.subplots(1,1)
        sl = self.df.segment.unique()
        cmap = plt.cm.get_cmap('jet', len(sl))
        for i,segment in enumerate(sl):
            color = cmap(i)
            self.df[self.df.segment==segment].plot(x='y', y='z', c=color, ax=ax, label=segment, linewidth=0.75)
        if hasattr(self, 'firsts'):
            for i in self.firsts:
                ax.plot([self.firsts[i]['y'], self.lasts[i]['y']], [self.firsts[i]['z'], self.lasts[i]['z']],color='Black')
         #   ax.scatter([self.lasts[i]['y']], [self.lasts[i]['z']], s=10, color='Yellow')
     
    #-------------------------------------------------
    
    def shiftPoints(self, scale:float=1, dx:float=0, dy:float=0):
        '''scale the points and shift the points in the combined dataframe'''
        # shift the points in each segmentPoints object
        for sp in self.splist:
            sp.shiftPoints(scale, dx, dy)
            
        # recombine the points
        self.combinedf()
        self.resetPoly()   # create new polygons with the new locations
        
    def plotSlice(self, ax, color) -> None:
        '''plot the slice on an existing axis'''
        for sp in self.splist:
            sp.plotSlice(ax, color)
            
    def plotArrow(self, pos:dict, recenter:bool, fs:folderStats) -> None:
        '''plot an arrow from the ideal position x0, y0 to the centroid of the slice'''
        x0 = pos['x0']
        y0 = pos['y0']
        color = pos['color']
        try:
            xmid,ymid = self.centroid()
        except:
            logging.info(f'Could not find centroid for {fs.folder}')
            traceback.print_exc()
            return
        if hasattr(fs.geo, 'spacing'):
            # neighboring filaments
            if float(fs.getVal('vink'))>0 and self.df.x.max()>fs.geo.ncx:
                # with extrusion
                if fs.geo.Ddir=='y':
                    x0 = x0-fs.geo.Ddis/2
                else:
                    y0 = y0-fs.geo.Ddis/2
            else:
                # no extrusion
                if fs.geo.Ddir=='y':
                    x0 = x0-fs.geo.Ddis
                else:
                    y0 = y0-fs.geo.Ddis
        if recenter:
            if fs.geo.Ddir=='y':
                x0 = x0+fs.geo.Ddis
            else:
                y0 = y0-fs.geo.Ddis

        pos['ax'].arrow(x0, y0, xmid-x0, ymid-y0
                        , head_width=0.05, head_length=0.1
                        , fc=color, ec=color, linewidth=0.75
                        , length_includes_head=True)

    #---------------------------------------------------
    
    def removeNozzleBottom(self, xL:float=-0.45, xR:float=0.45, yB:float=0.3, margin:float=0.02, plot:bool=True) -> None:
        '''remove points around the nozzle'''
        self.df0 = self.df0[~((self.df0.y>xL+margin)&(self.df0.y<xR-margin)&(self.df0.z>yB-margin))]
        self.sortPoints()
        self.splitPoints(True)
        if plot:
            self.plot()
    
    def sortSpreadAroundNozzle(self, xL:float=-0.45, xR:float=0.45, yB:float=0.3, margin:float=0.02, plot:bool=True) -> None:
        '''re-sort the points when the filament has spread up along the sides of the nozzle'''
        self.df.loc[:,'segment'] = 1
        self.labelNozzle(xL, xR, yB, margin)   # make the sides of the nozzle different objects
        self.splitPoints(True)
        if plot:
            self.plot()
        self.splist[1].sortVal('z', False)   # sort the right side descending
        self.splist[2].sortVal('z', True)    # sort the left side ascending
        if len(self.splist)>3:
            self.splist[3].sortVal('y', False)  # sort the bottom edge
        self.combinedf()
        self.selfFlatten()
        self.splitPoints(False)
        
        if plot:
            self.plot()
            
    def sort2Segments(self, l1:list, l2:list
                             , pts1:list=[], pts2:list=[]
                             , plot:bool=True
                             , **kwargs):
        '''we have two separated segments printed next to each other. separate them'''
        p1 = l1+pts1
        p2 = l2+pts2
        self.manualSplit({1:p1, 2:p2})
        self.splitPoints(True)
        if 'cx1' in kwargs and 'cy1' in kwargs:
            self.splist[0].sortRadially(cx=kwargs['cx1'], cy=kwargs['cy1'])
        if 'cx2' in kwargs and 'cy2' in kwargs:
            self.splist[1].sortRadially(cx=kwargs['cx2'], cy=kwargs['cy2']) 
        self.combinedf()
        if plot:
            self.plot()
            for p in [p1, p2]:
                plt.scatter(np.array(p)[:,0], np.array(p)[:,1], s=2, marker='X')
                
    def sort3Segments(self, l1:list, l2:list, l3:list
                             , pts1:list=[], pts2:list=[], pts3:list=[]
                             , plot:bool=True
                             , **kwargs):
        '''we have two separated segments printed next to each other. separate them'''
        p1 = l1+pts1
        p2 = l2+pts2
        p3 = l3+pts3
        self.manualSplit({1:p1, 2:p2, 3:p3})
        self.splitPoints(True)
        if 'cx1' in kwargs and 'cy1' in kwargs:
            self.splist[0].sortRadially(cx=kwargs['cx1'], cy=kwargs['cy1'])
        if 'cx2' in kwargs and 'cy2' in kwargs:
            self.splist[1].sortRadially(cx=kwargs['cx2'], cy=kwargs['cy2']) 
        if 'cx3' in kwargs and 'cy3' in kwargs:
            self.splist[3].sortRadially(cx=kwargs['cx3'], cy=kwargs['cy3']) 
        self.combinedf()
        if plot:
            self.plot()
            for p in [p1, p2, p3]:
                plt.scatter(np.array(p)[:,0], np.array(p)[:,1], s=2, marker='X')
                
                
    def sort2SegmentsInPlane(self, a1:float=0.8, x1:float=-0.46, y1:float=-0.2
                             , a2:float=1.8, x2:float=-0.44, y2:float=-0.18
                             , y0:float=-0.2, yf:float=0.5, dy:float=0.1
                             , **kwargs):
        '''we have two separated segments printed next to each other. separate them'''
        l1 = [[a1*(y-y1)**2+x1, y] for y in np.arange(y0, yf, dy)]
        l2 = [[a2*(y-y2)**2+x2, y] for y in np.arange(y0, yf, dy)]
        self.sort2Segments(l1, l2, **kwargs)
        
    def sort2SegmentsOutOfPlane(self, a1:float=0.8, x1:float=0, y1:float=-0.2
                             , a2:float=0.7, x2:float=0, y2:float=-0.18
                             , x0:float=-0.7, xf:float=-0.1, dx:float=0.1
                             , **kwargs):
        '''we have two separated segments printed on top of each other. separate them'''
        l1 = [[x,a1*(x-x1)**2+y1] for x in np.arange(x0, xf, dx)]
        l2 = [[x,a2*(x-x2)**2+y2] for x in np.arange(x0, xf, dx)]
        self.sort2Segments(l1, l2, **kwargs) 
        
    def sort3SegmentsOutOfPlane(self, a1:float=0.8, x1:float=0, y1:float=-0.2
                             , a2:float=0.7, x2:float=0, y2:float=-0.18
                             , x0:float=-0.7, xf:float=-0.1, dx:float=0.1
                             , **kwargs):
        '''we have two separated segments printed on top of each other. separate them and split the 2nd one in half horizontally'''
        l1 = [[x,a1*(x-x1)**2+y1] for x in np.arange(x0, xf, dx)]
        l2 = [[x,a2*(x-x2)**2+y2] for x in np.arange(x0, (x0+xf)/2, dx)]
        l3 = [[x,a2*(x-x2)**2+y2] for x in np.arange((x0+xf)/2, xf, dx)]
        self.sort3Segments(l1, l2, l3, **kwargs) 
        
    def simpleCombine(self, cx, cy, plot:bool=True):
        '''combine everything into one object and sort radially around cx, cy'''
        self.selfFlatten()
        self.splitPoints(True)
        self.splist[0].sortRadially(cx=cx, cy=cy, closeGaps=True)
        self.combinedf()
        if plot:
            self.plot()
 