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


class compileSurfaceFeatureExtractDict(OpenFOAMFile):
    '''compile surfaceFeatureExtractDict'''

    def __init__(self, bl:List[BoundaryInput]) -> str:
        super().__init__()
        bnames = [o.label for o in bl]
        s = self.header("dictionary", "surfaceFeatureExtractDict")
        s = s + DictList("", 0, \
                         list(map(lambda x: DictList(x+".stl", 0, \
                                                     [["extractionMethod", "extractFromSurface"],\
                                                      DictList("extractFromSurfaceCoeffs", 0, [["includedAngle", 180]]),\
                                                      ["writeObj", "yes"]]), bnames))).prnt(-1)
        s = s + self.closeLine()
        self.s = s