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
from points.slice_points import slicePoints
from summarize.summarizer import summarizer


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#------------------------------------------------------

class summarizerAdjacent(summarizer):
    '''for adjacent lines for a single simulation'''
    
    def __init__(self, folder:str, overwrite:bool=False):
        super().__init__(folder, overwrite)
        
    def defaultHeader(self) -> list:
        return list(self.sliceUnits().keys())
    
    def sliceUnits(self) -> dict:
        '''Find the units for a slice summary dictionary based on the units of the interfacePoints csv'''
        if not hasattr(self, 'units'):
            xu = 'mm'
            tu = 's'
            su = 'mm/s'
        else:
            xu = self.units['x']
            tu = self.units['time']
            su = self.units['vx']
        return {'x':xu, 'xbehind':xu, 'time':tu, 'centery':xu, 'centerz':xu, 'centeryr':xu, 'centerzr':xu
                , 'area':xu+'^2', 'maxheight':xu, 'maxwidth':xu, 'centeryn':'', 'centerzn':''
                , 'arean':'', 'maxheightn':'', 'maxwidthn':'', 'vertdisp':xu, 'vertdispn':''
                , 'horizdisp':xu, 'horizdispn':''
                , 'aspectratio':'', 'aspectration':'', 'speed':su, 'speeddecay':''
                , 'roughness':'', 'emptiness':'', 'asymmetryh':'w', 'asymmetryv':'h'}
    
        
        
    def sliceSummary(self, sli:pd.DataFrame) -> Dict[str, float]:
        '''sliceSummary collects important stats from a slice of a filament at a certain x and time and returns as a dictionary
        sli is a subset of points as a pandas dataframe
        fs is a folderStats object'''
        ev = -100 # error value
        rv = dict([[x, ev] for x in self.defaultHeader()])
            # time is in s
            # centery, centerz, minz, vertdisp, maxz, maxheight, maxwidth are in mm
            # speed is in mm/s
            # centeryn, centerzn, vertdispn, maxheightn, maxwidthn are normalized by nozzle inner diameter
            # aspectratio, topbotratio are dimensionless
            # speeddecay is normalized by the bath speed
        rv['x'] = sli.iloc[0]['x']
        rv['xbehind'] = rv['x']-self.fs.geo.ncx
        rv['time'] = np.float64(sli.iloc[0]['time'])

        if len(sli)<10:
            #logging.error('Not enough points')
            raise ValueError('Not enough points')
            
        sm = slicePoints(sli)
        rv['centery'], rv['centerz'], rv['area'] = sm.centroidAndArea()
        rv['centeryr'] = rv['centery'] - self.fs.geo.intentycenter
        rv['centerzr'] = rv['centerz'] - self.fs.geo.intentzcenter
        rv['maxheight'] = sli['z'].max()-sli['z'].min()
        rv['maxwidth'] = sli['y'].max()-sli['y'].min()
        if rv['maxheight']==0 or rv['maxwidth']==0:
            #logging.error('Cross-section is too small')
            raise ValueError('Cross-section is too small')

        rv['centerzn'] = rv['centerzr']/self.fs.geo.niw  # this is the center of mass
        rv['centeryn'] = rv['centeryr']/self.fs.geo.niw
        rv['arean'] = rv['area']/self.fs.geo.intentarea # normalize area by nozzle area
        rv['maxheightn'] = rv['maxheight']/self.fs.geo.intenth # normalize height by nozzle diameter
        rv['maxwidthn'] = rv['maxwidth']/self.fs.geo.intentw # normalized width by intended width

        rv['vertdisp'] = (sli['z'].min() - self.fs.geo.intentzbot)
        rv['vertdispn'] = rv['vertdisp']/self.fs.geo.niw
        rv['horizdisp'] = (sli['y'].min() - self.fs.geo.intentyleft)
        rv['horizdispn'] = rv['horizdisp']/self.fs.geo.niw
        rv['aspectratio'] = rv['maxheight']/rv['maxwidth']
        rv['aspectration'] = rv['aspectratio']/(self.fs.geo.intenth/self.fs.geo.intentw)
        rv['speed'] = sli['vx'].mean() # speed of interface points in x
        rv['speeddecay'] = rv['speed']/self.fs.sup.dynamic['v'] # relative speed of interface relative to the bath speed
        rv['roughness'] = sm.roughness()
        rv['emptiness'] = sm.emptiness()
        rv['asymmetryh'] = 0.5-(rv['centery']-sli.y.min())/rv['maxwidth']
        rv['asymmetryv'] = 0.5-(rv['centerz']-sli.z.min())/rv['maxheight']

        return rv