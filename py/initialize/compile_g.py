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
from dict_list import DictList
from initialize_tools import OpenFOAMFile

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class compileG(OpenFOAMFile):
    '''compile g'''
    
    def __init__(self):
        super().__init__()
        s = self.header("uniformDimensionedVectorField", "g")
        s = s + DictList("", 0, [["dimensions", "[0 1 -2 0 0 0 0]"], ["value", "(0 0 -9.81)"]]).prnt(-1)
        s = s + self.closeLine()
        self.s = s