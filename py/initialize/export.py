#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
from stl import mesh
import numpy as np
import os,sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from boundary_input import BoundaryInput
from initialize_tools import scale

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

def mkdirif(path:str) -> int:
    '''makes a directory if it doesn't already exist
    path is the directory to create'''
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        logging.info("Created directory %s" % path)

def saveStls(folder:str, bl:List[BoundaryInput]) -> None:
    '''exports stls
    folder is a full path name
    bl is a list of boundaryInput objects '''
    for b in bl:
        me = b.meshi
        me2 = mesh.Mesh(np.zeros(len(me), dtype=mesh.Mesh.dtype))
        for i,m in enumerate(me['vectors']):
            me2.vectors[i] = m*scale()
        title = os.path.join(folder, f'{b.label}.stl')
        me2.save(title)
                             
def replaceCR(file:str) -> None:
    '''replace windows carriage return with linux carriage return. Otherwise bash scripts will throw error. From https://stackoverflow.com/questions/36422107/how-to-convert-crlf-to-lf-on-a-windows-machine-in-python'''
    WINDOWS_LINE_ENDING = b'\r\n'
    UNIX_LINE_ENDING = b'\n'

    with open(file, 'rb') as open_file:
        content = open_file.read()
    
    content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

    with open(file, 'wb') as open_file:
        open_file.write(content)                    

def exportFile(folder:str, file:str, text:str, linux:bool=False) -> None:
    '''exports text files
    folder is a full path name
    file is a basename
    text is the text to export'''
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    logging.debug("Exported %s" % fn)
    if linux:
        replaceCR(fn)
    
