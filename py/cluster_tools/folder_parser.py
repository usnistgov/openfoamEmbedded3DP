#!/usr/bin/env python
'''Functions for handling files and folders used in OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. folderparser identifies log files from interFoam, etc. and collects information into csv tables
'''

# external packages
import os, sys
import numpy as np
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentdir)
from tools.config import cfg
from folder_stats import folderStats


#-------------------------------------------------------------------------------------------------  

    
def mkdirif(path:str) -> int:
    '''make a directory if it doesn't exist
    path is a full path name
    returns 1 if error, returns 0 for no error'''
    try:
        os.mkdir(path)
    except OSError as e:
        return 1
    else:
        logging.info ("Created directory %s" % path)
    return 0


def copy(src:str, dest:str) -> None:
    '''Copy directory src to directory dest. Both should be full folder paths'''
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            logging.error('Directory not copied. Error: %s' % e)
    return
 
    
######################################
#-------------------------------------------------------------------------------------------------  

    
def modifyControlDict(folder:str, tend:float) -> int:
    '''change the endtime in the controlDict
    returns 0 if the dictionary was changed, 1 if not'''
    cfi = os.path.join(folder, 'case')
    if os.path.exists(cfi):
        cf = cfi
    else:
        cf = folder
    cdfile = os.path.join(cf, 'system', 'controlDict')
    cdfilenew = os.path.join(cf, 'system', 'controlDict2')
    if not os.path.exists(cdfile):
        return 1
    retval = 0
    # if the endtime is already at the value, abort this loop and delete the new file
    with open(cdfile, 'r') as fold:
        with open(cdfilenew, 'w') as fnew:
            for line in fold:
                if line.startswith('endTime'):
                    linenew = 'endTime\t'+str(tend)+';\n'
                    if linenew==line:
                        retval = 1
                        break
                    else:
                        line = linenew
                fnew.write(line)
                
    # if we changed the endtime, overwrite the old file
    if retval==0:            
        os.remove(cdfile)
        os.rename(cdfilenew, cdfile)
        print('Set end time to '+str(tend)+' in '+cdfile)
    else:
        os.remove(cdfilenew)
    return retval

    
    
    