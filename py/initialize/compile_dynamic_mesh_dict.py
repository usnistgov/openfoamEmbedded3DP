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
from mesh_vars import MeshVars
from dict_list import DictList
from initialize_tools import OpenFOAMFile

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class compileDynamicMeshDict(OpenFOAMFile):
    
    def __init__(self, mv:MeshVars) -> str:
        super().__init__()
        s = self.header("dictionary", "dynamicMeshDict")
        simplelist = DictList("", 0, [["dynamicFvMesh", "dynamicRefineFvMesh"]]) # list to hold simple variables for dynamic meshing 
        s = s + simplelist.prnt(-1)
        correctfluxes = DictList("correctFluxes", 1, \
                                 [["phi", "none"], ["nHatf", "none"], ["rhoPhi", "none"],\
                                  ["alphaPhi", "none"], ["ghf", "none"], ["flux(alpha.ink)", "none"],\
                                 ["alphaPhi0.ink", "none"], ["alphaPhiUn", "none"], ["dVf_", "none"]])
        dmd = mv.dmd()
        dmd.append(correctfluxes)
                    # list to hold mesh coefficient dictionary entries for dynamic meshing
        s = s + DictList("dynamicRefineFvMeshCoeffs", 0, dmd).prnt(0)
        s = s + self.closeLine()
        self.s = s