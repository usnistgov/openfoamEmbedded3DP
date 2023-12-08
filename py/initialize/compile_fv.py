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
from fv_vars import FVVars
from dict_list import DictList
from initialize_tools import OpenFOAMFile
from file_group import FileGroup

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class compileFvSchemes(OpenFOAMFile):
    '''gets the text for the fvSchemes file. Input empty string to get file for mesh folder'''
    
    def __init__(self, fvv:Union[FVVars, str]):
        super().__init__()
        s = self.header("dictionary", "fvSchemes")
        if isinstance(fvv, FVVars):
            s = s + DictList("", 0, fvv.fvSchemeList()).prnt(-1)
        else:
            for st in ["gradSchemes", "divSchemes", "laplacianSchemes"]:
                s = s + st + "\n{\n}\n"
        s = s + self.closeLine()
        self.s = s

class compileFvSolution(OpenFOAMFile):
    '''gets the text for fvSolution. Input empty string to get empty file for mesh folder'''
    
    def __init__(self, fvv:Union[FVVars, str]):
        super().__init__()
        s = self.header("dictionary", "fvSolution")
        if isinstance(fvv, FVVars):
            s = s + fvv.solverlist().prnt(0)
            s = s + fvv.pimple().prnt(0)
            s = s + DictList("relaxationFactors", 0, [DictList("equations", 0, [["\".*\" 1"]])]).prnt(0)
        s = s + self.closeLine()
        self.s = s