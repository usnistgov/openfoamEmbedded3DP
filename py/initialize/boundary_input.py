#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
# local packages


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  
        

class BoundaryInput:
    '''# stores information about a boundary '''
    
    def __init__(self, labelin:str, typin:str):
        '''Input: labelin is the name of the boundary e.g. outerWall. typin is the boundary type, e.g. patch'''
        self.label = labelin # name of the boundary e.g. outerWall
        self.typ = typin # name of boundary type e.g. patch for blockMeshDict
        self.flist = [] # list of faces, which are lists of point indices
        self.alphalist = [] # list of alpha.ink properties
        self.Ulist = [] # list of U properties
        self.plist = [] # list of p_rgh properties
        self.meshi = [] # mesh for exporting to stl
        self.reflev = 0