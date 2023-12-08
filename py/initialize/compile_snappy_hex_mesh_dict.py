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
from boundary_input import BoundaryInput
from mesh_vars import MeshVars
from dict_list import DictList

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class compileSnappyHexMeshDict(OpenFOAMFile):
    '''compile SnappyHexMeshDict'''
    
    def __init__(self, bl:List[BoundaryInput], mv:MeshVars) -> str:
        super().__init__()
        bnames = [o.label for o in bl]
        s = self.header("dictionary", "snappyHexMeshDict")
        s = s + DictList("", 0, [\
                                 ["castellatedMesh", mv.castellatedMesh],\
                                 ["snap", mv.snap],\
                                 ["addLayers", mv.addLayers]]).prnt(-1)
        s = s + DictList("geometry", 0,\
                         map(lambda x: DictList(x+".stl", 0, [["type",  "triSurfaceMesh"], ["name ", x]]), bnames)).prnt(0)

        # castellated mesh controls
        cmc = mv.cmc()
        featuredictlist = DictList("features", 1, \
                list(map(lambda b: DictList("", 0, [["file", f"\"{b.label}.eMesh\""], ["level", b.reflev]]), bl)))
        refinementsurfaces = DictList("refinementSurfaces", 0,\
                list(map(lambda b: DictList(b.label, 0, [["level", f"({b.reflev} {b.reflev})"]]), bl)))
        cmc.append(featuredictlist)
        cmc.append(refinementsurfaces)
        cmc.append(DictList("refinementRegions", 0, [[]]))   
        s = s + DictList("castellatedMeshControls", 0, cmc).prnt(0)

        s = s + DictList("snapControls", 0,  mv.sc()).prnt(0)

        # add layers controls
        alc = mv.alc()
        layerslist = []
        for o in bl:
            if o.reflev>0:
                layerslist.append(DictList(o.label, 0, [["nSurfaceLayers", mv.nSurfaceLayers]]))
        alc.append(DictList("layers", 0, layerslist))
        s = s + DictList("addLayersControls", 0, alc).prnt(0)

        # mesh quality controls
        mqc = mv.mqc()
        mqc.append(["#include \"meshQualityDict\""])
        mqc.append(DictList("relaxed", 0, [["maxNonOrtho", 65]]))
        s = s + DictList("meshQualityControls", 0, mqc).prnt(0)

        s = s + DictList("writeFlags", 1, ["scalarLevels", "layerSets", "layerFields"]).prnt(0)
        s = s + DictList("", 0, [["mergeTolerance", "1E-6"]]).prnt(-1)
        s = s + self.closeLine()
        self.s = s