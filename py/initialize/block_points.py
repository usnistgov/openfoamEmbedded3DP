#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import numpy as np

# local packages

# logging
logging.basicConfig(level=logging.INFO)

#-------------------------------------------------------------------------------------------------  

def ptFromPts(pts:np.array, x:float, y:float, z:float) -> List[float]:
    '''select the point (x,y,z) from an array of points pts
    Input: pts, x, y, z'''
    pt = pts[(pts[:,0]==x) & (pts[:,1]==y) & (pts[:,2]==z), :]
    return pt[0]


def blockPts(pl:float, corner:List[float]) -> np.array:
    '''select the points in this block based on the first corner, and put the points in the order that openFOAM likes
    pl is an array containing points
    corner is one point which marks the first corner of the block'''
    xlist = np.unique(pl[:,0])
    ylist = np.unique(pl[:,1])
    zlist = np.unique(pl[:,2])
    cx = np.where(xlist==corner[0])
    cy = np.where(ylist==corner[1])
    cz = np.where(zlist==corner[2])
    cx = cx[0][0]
    cy = cy[0][0]
    cz = cz[0][0]
    pts = np.zeros([8, pl[1,:].size])
    for k in [0,1]:
        pts[0+4*k, :] = ptFromPts(pl, xlist[cx], ylist[cy], zlist[cz+k])
        pts[1+4*k, :] = ptFromPts(pl, xlist[cx+1], ylist[cy], zlist[cz+k])
        pts[2+4*k, :] = ptFromPts(pl, xlist[cx+1], ylist[cy+1], zlist[cz+k])
        pts[3+4*k, :] = ptFromPts(pl, xlist[cx], ylist[cy+1], zlist[cz+k])
    return pts

    