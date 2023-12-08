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
from boundary_input import BoundaryInput
from dict_list import DictList
from initialize_tools import OpenFOAMFile

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class compileSurfaceFeaturesDict(OpenFOAMFile):
    '''compile surfaceFeaturesDict'''
    
    def __init__(self, bl:List[BoundaryInput]):
        super().__init__()
    
        bnames = [o.label for o in bl]
        s = self.header("dictionary", "surfaceFeaturesDict")
        s = s + DictList("surfaces", 1, [f'"{b}.stl"' for b in bnames]).prnt(0)
        s = s + DictList("", 0, [["includedAngle", 180]]).prnt(-1)
        s = s + DictList("", 0, [["writeObj", "yes"]]).prnt(-1)
        s = s + self.closeLine()
        self.s = s