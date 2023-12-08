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

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#################################################################

def xpts(data:pd.DataFrame) -> pd.Series:
    '''List of all of the x positions in the dataframe'''
    xl = np.sort(data.x.unique())
    return xl

def closest(lst:List[float], K:float) -> float: # RG
    '''find the closest value in a list to the float K'''
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]


def pdCentroid(xs:pd.DataFrame) -> Tuple[float, float, float]:
    '''find the centroid of a slice represented as a pandas dataframe'''
    return centroid(np.array(xs[['y', 'z']]))

def setX(pts2d:np.array, x:float) -> np.array: # RG
    '''given an array of 2d points pts2d, create an array of 3d points at position x'''
    out = np.zeros([len(pts2d), 3])
    out[:, 1:3] = pts2d
    out[:, 0] = x
    return out


def setZ(pts2d:np.array, z:float) -> np.array:
    '''given an array of 2d points pts2d, create an array of 3d points at position z'''
    out = np.zeros([len(pts2d), 3])
    out[:, 0:2] = pts2d
    out[:, 2] = z
    return out

def sortPolar(pts:np.array) -> List[float]: # RG
    '''sort a list of points from 0 to 2 pi
    pts is a 3d list of unarranged points lying in the x, y, or z plane'''
    # maps x,y,z so z is the plane the xs lies in
    for i in range(3):
        if np.all(pts[:,i] == pts[0,i]):
            x = pts[:,(i+1)%3]
            y = pts[:,(i+2)%3]
            z = pts[:,i]
            j = i
    # organizes points by polar coordinates with the center as the origin
    x0 = np.mean(x)
    y0 = np.mean(y)
    r = np.sqrt((x-x0)**2+(y-y0)**2)
    theta = np.where(y>y0, np.arccos((x-x0)/r), 2*np.pi-np.arccos((x-x0)/r))
    mask = np.argsort(theta)
    xsort = x[mask]
    ysort = y[mask]
    xsort = np.append(xsort,xsort[0])
    ysort = np.append(ysort,ysort[0])
    z = np.append(z,z[0])
    # maps x,y,z back to original
    ptstemp = np.asarray(list(zip(xsort,ysort,z)))
    ptssort = np.zeros((len(xsort),3))
    for k in range(3):
        ptssort[:,(j+k+1)%3]=ptstemp[:,k]
    return ptssort, theta

def denoise(pts:np.ndarray) -> np.ndarray: # RG
    '''takes the interface which has inner and outer points, making jagged edges,
    and roughly turns it into only outer points
    outputs the outer points'''
    x = pts[:,0]
    y = pts[:,1]
    x0 = np.mean(x)
    y0 = np.mean(y)
    points = np.zeros((90,2))
    phi = np.degrees(sortPolar(setZ(np.column_stack((x,y)),0))[1])/4 # quarter-angle values for each point
    miss = 0 # handles situations where there are no points in the slice
    for theta in range(90): # slice for effectively every 4 degrees
        eligx = []
        eligy = []
        r = []
        for i,val in enumerate(phi):
            if val>=theta and val<theta+1:
                eligx.append(x[i]) # x of point in slice
                eligy.append(y[i]) # y of point in slice
                r.append(np.sqrt((x[i]-x0)**2+(y[i]-y0)**2)) # distance of point from center
        if len(r)>0:
            d = r.index(max(r))
            outpt = [eligx[d],eligy[d]] # coordinates of point furthest away
            points[theta-miss] = outpt # add point to points
        else:
            miss+=1
    if miss==0:
        return points
    return points[:-miss]

def xspoints(p:pd.DataFrame, dist:float, ore:str) -> Union[np.ndarray, List[float]]: # RG
    '''take interface points and shift them the desired offset from the nozzle
    outputs the points and the centerpoint of the points'''
    y = p['y']
    z = p['z']
    d = dist
    if ore=='y':
        y = [a-d for a in y] # shift cross section by d
    elif ore=='z':
        z = [a-d for a in z]
    else:
        raise Exception('Valid orientations are y and z')
    vertices = list(zip(y,z)) # list of tuples
    points = np.asarray(vertices) # numpy array of coordinate pairs
    centery = np.mean(y)
    centerz = np.mean(z)
    center = [centery, centerz]
    points = denoise(points)
    return points, center