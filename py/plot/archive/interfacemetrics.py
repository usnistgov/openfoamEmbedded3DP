#!/usr/bin/env python
'''Analyzing simulated single filaments'''

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
import folder_parser as fp
from pvCleanup import addUnits
from plainIm import *
import ncreate3d as nc
from file_meta import folderStats
import points_tools as pto

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

    
########### FILE HANDLING ##############

def takePlane(df:pd.DataFrame, folder:str, dr:float=0.05, xhalf:bool=True) -> pd.DataFrame:
    '''add an rbar column to the dataframe and only take the center plane'''
    fs = folderStats(folder)
    xc = fs.geo.ncx
    yc = fs.geo.ncy
    di = fs.geo.niw
    na = fs.geo.nozzle_angle
    tana = np.tan(np.deg2rad(na))        # tangent of nozzle angle
    zbot = fs.geo.nbz
    df['z'] =[zbot - i for i in df['z']] # put everything relative to the bottom of the nozzle
    ztop = fs.geo.nozzle_length
    df = df[df['z']>-ztop*0.9]           # cut off the top 10% of the nozzle
    
    dy = di/3                           # take middle y portion of nozzle
    df = df[(df['y']>-1*(dy))&(df['y']<dy)] # only select y portion
    if xhalf:
        df = df[(df['x']>xc)] # back half of nozzle
    
    df['rbar'] = [np.sqrt((row['x']-xc)**2+(row['y']-yc)**2)/(di/2+abs(row['z'])*tana) for i,row in df.iterrows()]
    df['rbar'] = [round(int(rbar/dr)*dr,5) for rbar in df['rbar']] # round to the closest dr
        # radius as fraction of total radius at that z
    df = df[df['rbar']<0.95]
        
    return df


def averageRings(df:pd.DataFrame) -> pd.DataFrame:
    '''given a list of points, group them by z and rbar, and take the average values by group'''
    vals = df.groupby(by=['z', 'rbar'])
    
def removeOutliers(df:pd.DataFrame, col:str, sigma:int=3) -> pd.DataFrame:
    '''remove outliers from the list based on column col'''
    med = df[col].median()
    stdev = df[col].std()
    return df[(df[col]>med-sigma*stdev)&(df[col]<med+sigma*stdev)]
