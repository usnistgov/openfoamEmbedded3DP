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
from summarize.steady import steadyMetrics
from summarize.summarizer_single import summarizerSingle
from summarize.summarizer_adjacent import summarizerAdjacent


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#-------------------------------------------

class sumAndSteadySingle:
    '''for a single line, summarize and get steady metrics'''
    
    def __init__(self, folder:str, overwrite:bool=False):
        self.sa = summarizerSingle(folder, overwrite=overwrite)
        if self.sa.success:
            for s in ['time', 'xbehind']:
                steadyMetrics(folder, overwrite=overwrite, mode=s, dother=1, vdcrit=0.01, col='vertdispn')
                
class sumAndSteadyAdjacent:
    '''for 2 adjacent lines, summarize and get steady metrics'''
    
    def __init__(self, folder:str, overwrite:bool=False):
        print(folder)
        self.sa = summarizerAdjacent(folder, overwrite=overwrite)
        if self.sa.success:
            for s in ['time', 'xbehind']:
                steadyMetrics(folder, overwrite=overwrite, mode=s, dother=1, vdcrit=0.01, col='vertdispn')
    
            