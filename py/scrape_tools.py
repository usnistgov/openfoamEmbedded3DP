#!/usr/bin/env python
'''Functions for generating legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import os,sys
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import re


# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)



#-------------------------------------------------------------------------------------------------  

def ca(col:List[List[str]], hlist:List[str], li:List[List[str]]) -> None:
        '''ca is used by scrape.table() to compile all of the data in the object into a table
        col is a table with 2 columns. 
        hlist is a list of headers that describe this chunk of data. they are added as rows with an empty value 
        li is a list of [variable name, value] to add to col'''
        for i in hlist:
            col.append([i, ''])
        for i in li:
            col.append(i)

def placeInList(l:List[List[str]], s:str, v:str) -> None: 
    '''placeInList puts a value v into a list with 2 columns, where s is the name of the variable
        this function will only place the value into the list if the variable name s is in the list and there isn't already a value there
        l is a list
        s is a string
        v is a value'''
    
    # find the index in the list where this value should go
    i = -1
    found=False
    while i<len(l)-1 and not found:
        i+=1
        if len(l[i])>0 and l[i][0]==s:
            found=True    
            
    # put the value in the list
    if not found:
        return
    else:
        l[i][1] = v # this puts the value in the list
    return


def cancelUnits(s:str) -> str:
    '''cancelUnits removes the units list (e.g. [ 0 2 -1 0 0 0 0]) from a value
    s is a string'''
    strs = re.split(' ', s)
    return strs[-1]


def listLevel(startString:str, endString:str, line:str, f:TextIO, l:List[List[str]]) -> str:
    '''listLevel is a tool that looks for sections of interest within files and scrapes values out of them
    startString is a string that tells us that we've reached the section of interest. It should be at the beginning of the line
    endString tells us that we've reached the end of the section of interest. It should be at the beginning of the line.
    line is the starting line
    f is a file stream created by open()
    l is a list of variable names and values that we're going to scrape values into
    returns the line we just read'''
    
    while not line.startswith(startString):
        line = f.readline()
    
    while not line.startswith(endString):
        strs = re.split(';|\t', line) # split the line at tabs and semicolons
        ii = 0
        si = 0
        while ii<len(strs) and len(strs[ii])==0:
            ii+=1 # find the first non-empty string in the list
        if ii+1<=len(strs)-1:
            s0 = strs[ii] # if there are enough entries in the split line to contain a name and value, place this entry into l
            s1 = cancelUnits(strs[ii+1])
            placeInList(l, s0, s1)
        line = f.readline()
    return line  


def readLevel0(s:str, line:str, f:TextIO, obj:List[List[str]]) -> str:
    '''readLevel0 is a simpler version of listLevel, where instead of placing many values in a list,
    we're looking for a single value. This function targets lines in files that have no tabs at the beginning, just "name\tvalue"
    s is a trigger string that tells us we've found the value we're looking for
    line is a starting line
    f is a file stream created by open()
    obj is a [1x2] list into which we'll store our value
    returns the line we just read'''
    while not line.startswith(s):
        line = f.readline()
    strs = re.split(';|\t', line) # split the line at ; and tabs
    obj[1] = strs[1] # the value will always be the second item
    return line