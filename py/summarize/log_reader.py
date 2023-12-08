#!/usr/bin/env python
'''Functions for plotting convergence for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import numpy as np
import os
import re
import pandas as pd
import csv
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging, platform, socket, sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from file.file_handling import folderHandler
from file.plainIm import plainIm, plainExp

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------------------------------------------------------------   

class logReader:
    '''scrapes info from logs'''
    
    def __init__(self, folder:str, overwrite:bool=False) -> None:
        self.folder = folder
        self.fh = folderHandler(folder)
        self.importFile(overwrite=overwrite)
        
    def fn(self):
        '''get the name of the interfoam file'''
        cf = self.fh.caseFolder()
        fn = os.path.join(cf, 'log_interFoam')
        if not os.path.exists(fn):
            fn = os.path.join(os.path.dirname(cf), 'log_interFoam')
        return fn
    
    def exportFN(self):
        '''get the name of the file to export to'''
        cf = self.fh.caseFolder()
        return os.path.join(cf, 'log_read.csv')
    
    def importFile(self, overwrite:bool=False) -> int:
        '''import values to the dataframe'''
        fn = self.exportFN()
        if os.path.exists(fn) and not overwrite:
            self.df,self.u = plainIm(fn, ic=0)
        else:
            self.readLog()
            self.exportFile()
        
    def exportFile(self):
        '''export values to a csv'''
        plainExp(self.exportFN(), self.df, self.u)
        
    def logEntry(self) -> dict:
        '''the logEntry dictionary is used to store information scraped from log files'''
        return {'courantmin': 0, 'courantmax': 0, 'deltaT': 0, 'simTime': 0, 'ralpha': 0, 'rprgh': 0, 'realtime': 0, 'cells': 0}
    
    def selectIf(self, strs:List[str], i:int) -> float:
        '''Get the float out of the list of strings, if the list of strings is long enough. otherwise, raise an error'''
        if len(strs)>i:
            try:
                f = float(strs[i])
            except Exception as e:
                print(e)
                raise NameError
            else:
                return f
        else:
            raise NameError
    
    def readLog(self):
        '''scrape values from the log'''
        file = self.fn()
        if not os.path.exists(file):
            return
        li = []
        with open(file, 'r') as f:
            for i in range(50): # skip all the headers and startup output
                line = f.readline()
            while line:
                try:
                    if line.startswith('Courant'): # we've hit a new time step
                        newEntry = self.logEntry()
    #                         li.append(logEntry()) # create a new object and store it in the list
    #                         lectr+=1
                        if len(li)>0:
                            newEntry['cells'] = li[-1]['cells'] 
                            # copy the number of cells from the last step and only adjust if the log says the number changed
                        strs = re.split('Courant Number mean: | max: |\n', line)
                        newEntry['courantmin'] = self.selectIf(strs, 1)
                        newEntry['courantmax'] = self.selectIf(strs, 2)
                    elif line.startswith('deltaT'):
                        strs = re.split('deltaT = |\n', line)
                        newEntry['deltaT'] = self.selectIf(strs, 1)
                    elif line.startswith('Time = '):
                        strs = re.split('Time = |\n', line)
                        newEntry['simTime'] = self.selectIf(strs, 1)
                        if len(li)>0:
                            lastTime = li[-1]['simTime']
                            while lastTime>newEntry['simTime'] and len(li)>1:
                                li = li[:-1]
                                lastTime = li[-1]['simTime']
                    elif line.startswith('Unrefined from '):
                        strs = re.split('Unrefined from | to | cells.\n', line)
                        newEntry['cells'] = self.selectIf(strs, 2)
                        if len(li)>0 and li[-1]['cells']==0:
                            for i in range(len(li)):
                                li[i]['cells'] = float(strs[1])
                                # the log never says the initial number of cells, 
                                # but it says the previous number if it changes the number of cells, 
                                # so to get the initial value, look for the first time the mesh is refined
                    elif line.startswith('smoothSolver'):
                        strs = re.split('Final residual = |, No Iterations', line)
                        newEntry['ralpha'] = self.selectIf(strs, 1)
                    elif line.startswith('DICPCG:  Solving for p_rgh,'):
                        strs = re.split('Final residual = |, No Iterations', line)
                        rprgh = self.selectIf(strs, 1)
                        newEntry['rprgh'] = rprgh
                    elif line.startswith('ExecutionTime'):
                        strs = re.split('ExecutionTime = | s', line)
                        newEntry['realtime'] = self.selectIf(strs, 1)
                        if newEntry['ralpha']>0 and newEntry['rprgh']>0:
                            while len(li)>0 and newEntry['simTime']<li[-1]['simTime']:
                                li = li[:-1]
                            li.append(newEntry)
                except NameError:
                    pass
                line = f.readline()
            li = li[1:]    
        self.df = pd.DataFrame(li)
        self.u = {'cells':'', 'courantmin':'', 'courantmax':'', 'deltaT':'s', 'simTime':'s', 'cells':'', 'ralpha':'', 'rprgh':'', 'realtime':'s'}
