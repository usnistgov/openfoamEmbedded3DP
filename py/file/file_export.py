#!/usr/bin/env python
'''Functions for importing files and csvs. Do not add pandas to this file, it will break paraview'''

# external packages
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import numpy as np
import csv

# local packages

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#----------------------------------------------

def plainImDict(fn:str, unitCol:int=-1, valCol:Union[int,list]=1) -> Tuple[dict,dict]:
    '''import values from a csv into a dictionary'''
    if type(valCol) is list:
        d = dict([[i,{}] for i in valCol])
    else:
        d = {}
    u = {}
    with open(fn, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            # save all rows as class attributes
            title = row[0]
            if title in d:
                title = f'{title}_1'
            if unitCol>0:
                u[row[0]] = row[unitCol]
            if type(valCol) is int:
                if valCol==1:
                    val = row[valCol]
                else:
                    val = (','.join(row[valCol:])).replace('\"', '')
                d[row[0]]=tryfloat(val)
            elif type(valCol) is list:
                for i in valCol:
                    d[i][row[0]]=tryfloat(row[i])
    return d,u

def plainExpDict(fn:str, vals:dict, units:dict={}, diag:bool=True, quotechar:str='|') -> None:
    '''export the dictionary to file'''
    with open(fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar=quotechar, quoting=csv.QUOTE_MINIMAL)
        for st,val in vals.items():
            if st in units:
                row = [st, units[st], val]
            else:
                if len(units)>0:
                    row = [st, '', val]
                else:
                    row = [st, val]
            writer.writerow(row)
    if diag:
        logging.info(f'Exported {fn}')
        
def exportFile(folder:str, file:str, text:str) -> None:
    '''exportFile exports a text file
    folder is the folder to save the file to (e.g. 'C:\\...\\myfolder')
    file is a file base name (e.g. 'myfile.txt') within that folder
    text is the text to write to the file'''
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    logging.info("Exported file %s" % fn)


def exportCSV(fn:str, table:List[Any]) -> None:
    '''exportCSV exports a csv file
    fn is the full path of the file to export
    table is the table to export, as a list'''
    with open(fn, mode='w', newline='') as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in table:
            w.writerow(i)
    logging.info('Exported file %s' % fn)
    