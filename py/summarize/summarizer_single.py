#!/usr/bin/env python
'''Analyzing simulated single filaments'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import points.points_tools as pto
from summarize.summarizer import summarizer


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#------------------------------------------------------

class summarizerSingle(summarizer):
    '''for single lines'''
    
    def __init__(self, folder:str, overwrite:bool=False):
        super().__init__(folder, overwrite)
        
    def defaultHeader(self) -> list:
        return ['x', 'xbehind', 'time', 'centery', 'centerz', 'area', 'maxheight', 'maxwidth', 'centeryn', 'centerzn', 'arean', 'maxheightn', 'maxwidthn','vertdisp', 'vertdispn', 'aspectratio', 'speed', 'speeddecay']
    
    def sliceUnits(self) -> dict:
        '''Find the units for a slice summary dictionary based on the units of the interfacePoints csv'''
        xu = self.units['x']
        return {'x':xu, 'xbehind':xu, 'time':self.units['time'], \
              'centery':xu, 'centerz':xu, 'area':xu+'^2', 'maxheight':xu, 'maxwidth':xu, \
               'centeryn':'', 'centerzn':'', 'arean':'', 'maxheightn':'', 'maxwidthn':'',\
              'vertdisp':xu, 'vertdispn':'', 'aspectratio':'', 'speed':self.units['vx'], 'speeddecay':''}
    
    def xlist(self, data) -> list:
        '''get a list of x positions to probe'''
        xlist = pto.xpts(data)
        xlist = xlist[xlist>self.fs.geo.behind]
        return xlist
        
        
    def sliceSummary(self, sli:pd.DataFrame) -> Dict[str, float]:
        '''sliceSummary collects important stats from a slice of a filament at a certain x and time and returns as a dictionary
        sli is a subset of points as a pandas dataframe
       '''
        ev = -100 # error value
        rv = {'x':ev, 'xbehind':ev, 'time':ev, \
              'centery':ev, 'centerz':ev, 'area':ev, 'maxheight':ev, 'maxwidth':ev, \
               'centeryn':ev, 'centerzn':ev, 'arean':ev, 'maxheightn':ev, 'maxwidthn':ev,\
              'vertdisp':ev, 'vertdispn':ev, 'aspectratio':ev, 'speed':ev, 'speeddecay':ev} # return value
            # x is in mm
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

            
        sm = pto.shapeMeasure(sli)
        try:
            rv['centery'], rv['centerz'], rv['area'] = sm.centroidAndArea()
        except KeyboardInterrupt as e:
            raise e
        except:
            #logging.error('centroid error')
            raise ValueError('centroid error')
        rv['maxheight'] = sli['z'].max()-sli['z'].min()
        rv['maxwidth'] = sli['y'].max()-sli['y'].min()
        if rv['maxheight']==0 or rv['maxwidth']==0:
            #logging.error('Cross-section is too small')
            raise ValueError('Cross-section is too small')

        rv['centerzn'] = rv['centerz']/self.fs.geo.niw
        rv['centeryn'] = rv['centery']/self.fs.geo.niw
        rv['arean'] = rv['area']/(np.pi*(self.fs.geo.niw/2)**2) # normalize area by nozzle area
        rv['maxheightn'] = rv['maxheight']/self.fs.geo.niw # normalize height by nozzle diameter
        rv['maxwidthn'] = rv['maxwidth']/self.fs.geo.niw # normalized width by intended width

        rv['vertdisp'] = (sli['z'].min() - self.fs.geo.intentzbot)
        rv['vertdispn'] = rv['vertdisp']/self.fs.geo.niw
        rv['aspectratio'] = rv['maxheight']/rv['maxwidth']
        rv['speed'] = sli['vx'].mean() # speed of interface points in x
        rv['speeddecay'] = rv['speed']/self.fs.geo.bv # relative speed of interface relative to the bath speed

        return rv