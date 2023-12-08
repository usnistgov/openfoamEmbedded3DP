#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import os,sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_tools import OpenFOAMFile

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class compileMeshQualityDict(OpenFOAMFile):
    '''compile meshQualityDict'''
    
    def __init__(self):
        super().__init__()
        s = self.header("dictionary", "meshQualityDict")
        s = s + "#includeEtc \"caseDicts/mesh/generation/meshQualityDict\"\n\n" # RG
        s = s + self.closeLine()
        self.s = s