#!/usr/bin/env python
'''Functions for generating legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import os
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket, sys
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import file.file_handling as fh
from file.backwards_read import fileReadBackwards
from file.file_export import *
from scrape import scrape

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)



#-------------------------------------------------------------------------------------------------  

def legendFN(folder:str) -> str:
    return os.path.join(folder, 'legend.csv')

def importIf(folder:str) -> List[List[str]]:
    '''importIf imports an existing legend file if it exists 
    folder is a full path, e.g. "C:\\...\\nb64"'''
    fn = legendFN(folder)
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            return(list(csv.reader(f)))
    else:
        print('no legend')
        return []

def legendTableToDict(t:np.array, units:bool=False) -> Union[Tuple[dict,dict], dict]:
    '''convert the table into a dictionary with unique headers'''
    headers = [i[0] for i in t]
    headers2 = headers
    section = ''
    for i in range(len(t)):
        if t[i][0]=='':
            section = t[i+1][0]
        elif t[i][0] in ['sup', 'ink', 'controlDict', 'dynamicMeshDict']:
            section = t[i][0]
        else:
            if headers.count(t[i][0])>1:
                newname = t[i][0].replace(' ', '_')
                if len(section)>0:
                    newname = f'{section}_{newname}'
                if headers2.count(newname)>0 or len(section)==0:
                    if len(t[i])>2:
                        newname = f'{newname}_{t[i][2]}'
                    else:
                        newname = f'{newname}_0'
                t[i][0] = newname
                headers2 = [i[0] for i in t]
            else:
                t[i][0] = t[i][0].replace(' ', '_')
        if len(t[i])==2:
            t[i] = t[i] +['']
    values = {a[0]:a[1] for a in t}        

    if units:
        if len(t[4])==3:
            u = {a[0]:a[2] for a in t}
        else:
            u = legendUnique(cfg.path.legend_units, units=False)
        return values, u
    else:
        return values

def legendUnique(folder:str, units:bool=False) -> Union[Tuple[dict,dict], dict]:
    '''legendUnique imports a legend, rewrites the variable names so they are all unique and ready to compile into a pandas dataframe, and outputs a dictionary
    if units=True, also imports a dictionary of units'''
    t = importIf(folder)
    if len(t)==0:
        return {}

    if units:
        values, u = legendTableToDict(t, units=units)
    else:
        values = legendTableToDict(t, units=units)

    # values['compare_to'] = os.path.basename(os.path.dirname(folder)) # name of the folder it's in RG
    if units:
        return values,u
    else:
        return values

def populate(folder:str, *varargin, readLogs:bool=True, overwrite:bool=False) -> List[List[str]]:
    '''populate scrapes all the data from the folder and exports it to a table called legend.csv
    folder is a full path name
    can also add strings that evaluate specific functions, e.g. scrapeTP, so that you can just scrape the transportproperties if there is already a legend'''
    if not fh.folderHandler(folder).isSimFolder():
        raise Exception("Not a simulation folder")
    s = scrape(folder)   # create an object to store variables
    fn = legendFN(folder)     # export file name
    # s.compareto[1] = os.path.basename(os.path.dirname(folder)) # RG
    if readLogs:
        s.scrapeLogs()   # scrape the logs
        s.scrapeCD() # scrape the control dictionary
    if not overwrite and os.path.exists(fn):
        leOld, uOld = legendUnique(folder, units=True) # import old legend and units to dictionary
        t = s.table()
        leNew, uNew = legendTableToDict(t, units=True) # convert new legend and units to dictionary
        leNew['compare_to'] = leOld['compare_to'] # RG
        for key in leNew:
            if not key=='compare_to':
                # keep old values if not same, except for scraped timings
                if key in leOld and not leOld[key]==leNew[key]:
                    if not (('time' in key or 'rate' in key) and len(leNew[key])>0 and float(leNew[key])>0):
                        if len(leOld[key])>0:
                            if leOld[key][0]==' ':
                                # remove space from beginning of unit
                                leOld[key] = leOld[key][1:]
                            leNew[key] = leOld[key]
                if key in uOld and not uOld[key]==uNew[key]:
                    if uOld[key][0]==' ':
                        # remove space from beginning of unit
                        uOld[key] = uOld[key][1:]
                    uNew[key] = uOld[key]
        t = [[key, leNew[key], uNew[key]] for key in leNew]                        
    else:
        # if there is no legend file, we have to go through all the files and scrape data
        s.scrapeAll()
        t = s.table() # generate a table from all the data we scraped
    exportCSV(fn, t) # export the updated table
    return t



def populateToTable(folder:str, repopulate:bool = False) -> List[List[str]]:
    '''take one folder and either import the existing legend or create a new one
    folder is a full path name
    repopulate is true if you want to overwrite existing legend.csv files'''
    if repopulate:
        # overwrite times in legend.csv, or if there is no file, create a new legend file
        t2 = populate(folder)
    else:
        # don't overwrite anything, just import the existing file
        t2 = importIf(folder)
    return t2


def populateList(liInit:List[str], exportFilename:str, repopulate:bool = False) -> None:
    '''populateList scrapes all data for all files in a list and creates a combined table
    liInit is a list of folders to scrape
    exportFilename is the destination to export the summary file
    repopulate is true if you want to overwrite existing legend.csv files'''
    li = []
    for folder in liInit:
        if os.path.exists(folder):
            li.append(folder)
    t1 = populateToTable(li[0], repopulate)   
        # first collect one table so you have the format
    tbig = [['' for i in range(len(li)+1)] for j in range(len(t1))] 
        # tbig combines all the legend files into one table
    for i,row in enumerate(t1):
        tbig[i][0] = row[0] 
            # import the variable names and the values 
            # for the first file into the big table
        tbig[i][1] = row[1]
    for j,l in enumerate(li[1:]): 
        # go through the rest of the files and put the values in the big table
        t2 = populateToTable(l, repopulate)
        for i,row in enumerate(t2):
            tbig[i][j+2] = row[1]
    exportCSV(exportFilename, tbig) 
        # export the combined table
    return




#------------------------------

# RETROFITTING FUNCTIONS


def updateGeoFile(folder:str) -> None:
    '''update the geometry file to put units in the third column'''
    fn = os.path.join(folder, 'geometry.csv')
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            geo = list(csv.reader(f))
    else:
        # no geometry file
        return
    if len(geo[0])==3:
        # geometry file already has units
        return
    new = []
    for row in geo:
        spl = re.split(' \(|\)', row[0])
        if len(spl)>1:
            newrow = [spl[0], row[1], spl[1]]
        else:
            newrow = [spl[0], row[1], '']
        new.append(newrow)
    exportCSV(fn, new)
    logging.info(f'Exported {fn}')
    
def addUnitsToLegend(folder:str) -> None:
    '''add units to the legend file'''
    updateGeoFile(folder)
    populate(folder)
        

