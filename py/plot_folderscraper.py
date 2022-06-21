#!/usr/bin/env python
'''Functions for plotting convergence for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket, sys
from backwardsRead import fileReadBackwards

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from folderparser import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------------------------------------------------------------       
    
    
    
###################################################
###################################################
# These functions are used to scrape log files for times and fitting metrics
###################################################


def logEntry() -> Dict:
    '''the logEntry dictionary is used to store information scraped from log files'''
    return {'courantmin': 0, 'courantmax': 0, 'deltaT': 0, 'simTime': 0, 'ralpha': 0, 'rprgh': 0, 'realtime': 0, 'cells': 0}

def interFile(folder:str) -> str:
    '''interFile finds the interFoam log
    folder can be a case folder or its parent'''
    cf = caseFolder(folder)
    fn = os.path.join(cf, 'log_interFoam')
    if not os.path.exists(fn):
        fn = os.path.join(os.path.dirname(cf), 'log_interFoam')
    return fn

def selectIf(strs:List[str], i:int) -> float:
    '''Get the float out of the list of strings, if the list of strings is long enough. otherwise, raise an error'''
    if len(strs)>i:
        try:
            f = float(strs[i])
        except Exception as e:
            raise NameError
        else:
            return f
    else:
        raise NameError


def logRead(folder:str) -> List[logEntry]:
    '''logRead extracts values from log files and stores them in a list of logEntry objects
    folder can be a case folder or its parent'''
    intf = interFile(folder)
    li = []   # this will be a list of logEntry objects
    if os.path.exists(intf):
        with open(intf, 'r') as f:
            for i in range(50): # skip all the headers and startup output
                line = f.readline()
            while line:
                try:
                    if line.startswith('Courant'): # we've hit a new time step
                        newEntry = logEntry()
#                         li.append(logEntry()) # create a new object and store it in the list
#                         lectr+=1
                        if len(li)>0:
                            newEntry['cells'] = li[-1]['cells'] 
                            # copy the number of cells from the last step and only adjust if the log says the number changed
                        strs = re.split('Courant Number mean: | max: |\n', line)
                        newEntry['courantmin'] = selectIf(strs, 1)
                        newEntry['courantmax'] = selectIf(strs, 2)
                    if line.startswith('deltaT'):
                        strs = re.split('deltaT = |\n', line)
                        newEntry['deltaT'] = selectIf(strs, 1)
                    elif line.startswith('Time = '):
                        strs = re.split('Time = |\n', line)
                        newEntry['simTime'] = selectIf(strs, 1)
#                         if len(li)>0 and newEntry['simTime']<li[-1]['simTime']: # if the time is less than the previous time, restart the table
#                             print('restarting, ', newEntry['simTime'], li[-1]['simTime'])
# #                             li = []
                    elif line.startswith('Unrefined from '):
                        strs = re.split('Unrefined from | to | cells.\n', line)
                        newEntry['cells'] = selectIf(strs, 2)
                        if len(li)>0 and li[-1]['cells']==0:
                            for i in range(len(li)):
                                li[i]['cells'] = float(strs[1])
                                # the log never says the initial number of cells, 
                                # but it says the previous number if it changes the number of cells, 
                                # so to get the initial value, look for the first time the mesh is refined
                    elif line.startswith('smoothSolver'):
                        strs = re.split('Final residual = |, No Iterations', line)
                        newEntry['ralpha'] = selectIf(strs, 1)
                    elif line.startswith('DICPCG:  Solving for p_rgh,'):
                        strs = re.split('Final residual = |, No Iterations', line)
                        newEntry['rprgh'] = selectIf(strs, 1)
                    elif line.startswith('ExecutionTime'):
                        strs = re.split('ExecutionTime = | s', line)
                        newEntry['realtime'] = selectIf(strs, 1)
                        li.append(newEntry)
                except Exception as e:
                    # if we hit an error, skip this line
                    print(f'error hit: {e}')
                    pass
                line = f.readline()
            li = li[1:]    
        #### plot 
        # plotAll(folder, li) 
        # rAlphaPlot(folder, li)
    return pd.DataFrame(li)


def la(li:List[Any], at:str) -> List[Any]:
    '''list an attribute for each object in a list
    li is a list of objects
    at is the attribute that we want to get from each object'''
    return list(getattr(o, at) for o in li)

#--------------------------------------------------
# plots


def plotConvergence(folder:str, li:Union[List, pd.DataFrame], export:bool=False):
    '''plot 4 plots that show the simulation metrics over time
    folder can be a case folder or its parent
    li is a dataframe holding log values. Give empty list to automatically generate a list'''
    fn = os.path.join(folder, 'images', 'convergence.png')
    if export and os.path.exists(fn): # if the file already exists, return
        return
    if len(li)==0:
        li = logRead(folder)
    if len(li)==0: # if there is no log, return
        return
    fig, axs = plt.subplots(4, sharex=True)
    fig.suptitle = os.path.basename(folder)
    fig.set_size_inches(3, 9)
    fs = 12
    axs[0].plot(li['simTime'], li['ralpha'], 'maroon', label='Alpha') 
    axs[0].plot(li['simTime'], li['rprgh'], label='p_rgh') 
    axs[0].set_yscale("log")
#     axs[0].set_xlabel('Time in simulation (s)', fontsize=fs) 
    axs[0].set_ylabel('Residual', fontsize=fs)     
    axs[0].legend()
    axs[0].set_title(os.path.basename(folder), fontsize=fs)

    axs[1].plot(li['simTime'], li['cells']) 
#     axs[1].set_xlabel('Time in simulation (s)', fontsize=fs) 
    axs[1].set_ylabel('Cells in mesh', fontsize=fs) 

    axs[2].plot(li['simTime'], li['courantmin'], 'maroon', label='Co min') 
    axs[2].plot(li['simTime'], li['courantmax'], label='Co max') 
#     axs[2].set_xlabel('Time in simulation (s)', fontsize=fs) 
    axs[2].set_ylabel('Courant number', fontsize=fs)     
    axs[2].legend()

    axs[3].plot(li['simTime'], li['deltaT']) 
    axs[3].set_xlabel('Time in simulation (s)', fontsize=fs) 
    axs[3].set_ylabel('Delta t', fontsize=fs)
    
    plt.close()
    if export:
        fig.savefig(fn, bbox_inches='tight')
        logging.info(f'{shortName(folder)}\images\convergence.png')
    return fig



def rAlphaPlot(folders:List[str], lis:List[Any]):
    '''plot just the alpha residual over time
    folders ia a list of folders
    lis is a list of lists of logEntry objects'''
    if len(lis)==0:
        lis = [logRead(folder) for folder in folders]
    fs = 12
    fig, axs = plt.subplots(1)
    fig.set_size_inches(8, 8)
    colors = ['firebrick', 'darkblue', 'burlywood', 'dodgerblue']
    for i, li in enumerate(lis):
        axs.plot(la(li, 'simTime'), la(li, 'ralpha'), colors[i], label=folders[i])
    axs.set_xlabel('Time in simulation (s)', fontsize=fs) 
    axs.set_ylabel('Alpha residual', fontsize=fs) 
    axs.legend()
    plt.close()
    return fig