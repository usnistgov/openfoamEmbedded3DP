#!/usr/bin/env python
'''Ideal values of measured variables'''

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

class ideals:
    '''stores the ideal values for measurements'''
    
    def __init__(self):
        self.adj = {'centeryr':0, 'centerzr':0
                    , 'centeryn':0, 'centerzn':0
                    , 'arean':1
                    , 'maxheightn':1, 'maxwidthn':1
                    , 'vertdisp':0, 'vertdispn':0
                    , 'horizdisp':0, 'horizdispn':0
                    , 'aspectration':1,  'speeddecay':1
                    , 'roughness':0, 'emptiness':0
                    , 'asymmetryh':0, 'asymmetryv':0}
        self.single = {'centery':0, 'centerz':0
                       , 'centeryn':0, 'centerzn':0
                       , 'arean':1
                       , 'maxheightn':1, 'maxwidthn':1
                       , 'vertdisp':0, 'vertdispn':0
                       , 'aspectratio':1, 'speeddecay':1}
        
    def getValue(self, var:str, printType:str=''):
        if printType=='adjacent':
            if var in self.adj:
                return self.adj[var]
            else:
                return ''
        elif printType=='single':
            if var in self.single:
                return self.single[var]
            else:
                return ''
        else:
            if var in self.adj:
                return self.adj[var]
            elif var in self.single:
                return self.single[var]
            else:
                return ''
        