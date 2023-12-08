#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''

# external packages
import os
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging

# local packages
from file_creator import fileCreator
from fluid import Fluid

# logging
logging.basicConfig(level=logging.INFO)

#-------------------------------------------------------------------------------------------------  

def genericExport(ii:Union[int,str], sup:Fluid, ink:Fluid, sigma:float, topFolder:str, exportMesh:bool=False, **kwargs) -> None:
    ''' Export a folder, given a support fluid, ink fluid, and surface tension. 
        ii is for the folder label. If you want it to be labeled nb#, input a number. Otherwise, input a string.
        sup is a fluid object that holds info about the support transport properties
        ink is a fluid object that holds info about the ink transport properties
        sigma is in J/m^2 (e.g. 0.04)
        topFolder is the folder to save this new folder in
        exportMesh true to export geometry folders into this folder'''
    out = fileCreator(ii, topFolder, exportMesh=exportMesh, **kwargs)
    out.addTransportProperties(ink, sup, sigma)
    out.export()

def genericMesh(parentFolder:str, **kwargs) -> fileCreator:
    '''This generates a folder with mesh and geometry files, but nothing else important. If you want to customize the nozzle, input keyword variables. e.g. genericMesh(myFolder, bathWidth=16)'''
    out = fileCreator('temp', parentFolder, exportMesh=True, onlyMesh=True, **kwargs) # generate mesh files
    out.export() # export all of the mesh files
    out.plot() # make a plot of the boundaries
    return out
