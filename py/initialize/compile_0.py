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

class compiler0(OpenFOAMFile):
    '''puts files in the 0 folder into the correct format
    classn is the class that goes into the header, e.g. "volScalarField" or "pointScalarField"
    obj is the object that goes into the header, e.g. "U"
    dims is the dimensions of the field, e.g. "[1 -1 -2 0 0 0 0]", which means kg/(m*s^2)
    intfield is the internal field, e.g. "0" or "(0 0 0)"
    # f is the function to run on each boundaryInput object to get a list of properties of interest'''
    
    def __init__(self, bl:List[BoundaryInput], classn:str, obj:str, dims:str, intfield:str, f:Callable):
        super().__init__()    
        s = self.header(classn, obj) # header
        simplelist = DictList("", 0, [["dimensions", dims], ["internalField", "uniform " + intfield]]) # list of simple variables to define
        s = s + simplelist.prnt(-1) # format the simple list 
        bflist = DictList("boundaryField", 0, []) # list of boundary dictionary entries
        for b in bl:
            bflist.proplist.append(f(b)) # add the appropriate list for each boundary, e.g. alpha properties for inkFlow
        s = s + bflist.prnt(0) # format the dictionary list
        s = s + self.closeLine()
        self.s = s
        
class compileAlphaOrig(compiler0):
    
    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__(bl, "volScalarField", "alpha.ink", "[0 0 0 0 0 0 0]", "0", lambda bi:bi.alphalist)
    
class compileU(compiler0):
    
    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__(bl, "volVectorField", "U", "[0 1 -1 0 0 0 0]", "(0.01 0 0)", lambda bi:bi.Ulist)
    
class compileP(compiler0):
    
    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__(bl, "volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]", "0", lambda bi:bi.plist)

class compileCellLevel(compiler0):
    
    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__(bl, "volScalarField", "cellLevel", "[0 0 0 0 0 0 0]", "0", lambda bi:[["type", "zeroGradient"]])

class compilePointLevel(compiler0):
    
    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__(bl, "pointScalarField", "pointLevel", "[0 0 0 0 0 0 0]", "0", lambda bi:[["type", "calculated"]])