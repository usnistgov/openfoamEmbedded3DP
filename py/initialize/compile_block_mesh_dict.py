#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import numpy as np
import os,sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from block import Block
from boundary_input import BoundaryInput
from dict_list import DictList
from initialize_tools import OpenFOAMFile, scale

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class compileBlockMeshDict(OpenFOAMFile):
    '''get the blockMeshDict text
    pl is a point list
    blocks is a list of Block objects
    bl is a list of BoundaryInput objects'''
    
    def __init__(self, pl:np.array, blocks:List[Block], bl:List[BoundaryInput]) -> str:
        super().__init__()
        s = self.header("dictionary", "blockMeshDict")
        s = s + f"scale {scale()};\n\n"
        s = s + DictList("vertices", 1, pl[:, 0:3]).prnt(0)
        s = s + self.blocks2txt(blocks)
        s = s + "edges\n(\n);\n\n"
        s = s + self.boundarycpp(bl)
        s = s + "mergePatchPairs\n(\n);\n\n" 
        s = s + self.closeLine()
        self.s = s
        
    def vec2cpp(self, v:List[float]) -> str:
        '''convert a vector to the format it needs to be in for OpenFOAM to read it
        v is a list'''

        s = "("
        for vi in v:
            s = f'{s}{vi} '
        s = s + ")"
        return s

    def block2txt(self, block:Block) -> str:
        '''convert block to openfoam string'''
        s = "hex " + self.vec2cpp(block.vertices[:, 3].astype(int)) + " " + self.vec2cpp(block.meshi)
        s = s + " simpleGrading " + self.vec2cpp(block.grading)
        return s

    def blocks2txt(self, blocks:List[Block]) -> str:
        '''convert list of blocks to openfoam string'''
        pl = DictList("blocks", 1, [])
        for b in blocks:
            pl.proplist.append(self.block2txt(b))
        return pl.prnt(0)

    def faceSelector(self, block:Block, st:str) -> np.array:
        '''gets a list of vertices for a list
        block is a Block object
        st is a string indicating which face to use'''
        if st == "x-":
            li = [0, 4, 7, 3]
        elif st == "x+":
            li = [1, 2, 6, 5]
        elif st == "y-": 
            li = [1, 5, 4, 0]
        elif st == "y+":
            li = [3, 7, 6, 2]
        elif st == "z-": 
            li = [0, 3, 2, 1]
        elif st == "z+": 
            li = [4, 5, 6, 7]
        return block.vertices[li, :]

    def boundaryList(self, blocks:List[Block]) -> List[BoundaryInput]:
        '''compiles a list of all the boundaries in the system for blockMesh
        because we are using snappyHexMesh, we only need one boundary in blockMesh'''
        allb = BoundaryInput("allBoundary", "patch")
        for st in ["x-", "x+", "y-", "y+", "z-", "z+"]:
            allb.flist.append(self.faceSelector(blocks[0], st)) # front and back faces
        return [allb]

    def boundarycpp(self, bl:List[BoundaryInput]) -> str: 
        '''get the whole cpp section for list of boundaries'''
        bb = DictList("boundary", 1, [])
        for b in bl:
            # for each boundary in the boundary list
            thisboundary = DictList(b.label, 0, []) # label of the boundary
            thisboundary.proplist.append(["type patch"]) # boundary type
            thesefaces = DictList("faces", 1, []) # list of faces
            for f in b.flist: # extract list of faces, which are lists of points
                thesefaces.proplist.append(f[:,3].astype(int)) # for each list of points, take the last column
            thisboundary.proplist.append(thesefaces) # save faces to this boundary
            bb.proplist.append(thisboundary) # save this boundary to list of boundaries
        return bb.prnt(0)