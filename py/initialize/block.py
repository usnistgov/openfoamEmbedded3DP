#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging

# local packages


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  

class Block:
    '''stores information about a block during blockMesh '''
    
    def __init__(self):
        self.vertices = []
        self.meshi = [1,1,1]
        self.grading = [1,1,1]