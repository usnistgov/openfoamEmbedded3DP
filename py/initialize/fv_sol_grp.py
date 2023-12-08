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

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class FVSolGrp:
    '''fvsolution variables for a given solve variable (e.g. interFoam)'''

    def __init__(self, st:str, solver:str):
        '''Inputs: st, solver
        st is a string that tells us what variable is being solved for, e.g. 'alphaink', 'pcorr', 'prgh', 'prghfinal', or 'U'
        solver is the type of solver being used: 'interFoam' or 'interIsoFoam'. '''
        
        self.badcharlist = []
                     
        if st=="alphaink":
            self.dicttitle = "\"alpha.ink.*\""
            self.nAlphaSubCycles = 1
            self.cAlpha = 1
            if solver=="interIsoFoam":
                # these variables are specific to isoAdvector
                self.isofaceTol = "1e-6" # Error tolerance on alpha when cutting surface cells into sub-cells
                self.surfCellTol = "1e-6" # Only cells with surfCellTol < alpha < 1-surfCellTol 
                                            # are treated as surface cells
                self.nAlphaBounds = 3  # Number of times the ad-hoc bounding step should
                                        # try to correct unboundedness. Strictly volume
                                        # conserving (provided that sum(phi) = 0 for a cell).
                self.snapTol = "1e-12" # Optional: cells with alpha < snapAlphaTol are
                                        # snapped to 0 and cells with 1 - alpha <
                                        # snapAlphaTol are snapped to 1
                self.clip = "true" # Optional: clip remaining unboundedness
            elif solver=="interFoam":
                # these variables are specific to MULES
                self.nAlphaCorr = 2
                self.MULESCorr = "yes"
                self.nLimiterIter = 5
                self.solver = "smoothSolver"
                self.smoother = "symGaussSeidel"                 
                self.tolerance = "1e-8"
                self.relTol = 0
        elif st=="pcorr":
            self.dicttitle = "\"pcorr.*\""
            self.solver = "PCG"
            self.preconditioner = "DIC"
            self.tolerance = "1e-5"
            self.relTol = 0
        elif st=="prgh":
            self.dicttitle = "p_rgh"
            self.solver = "PCG"
            self.preconditioner = "DIC"
            if solver=="interIsoFoam":
                self.tolerance = "1e-8"
            else:
                self.tolerance = "1e-7"
            self.relTol = 0.05
        elif st == "prghfinal":
            self.dicttitle = "p_rghFinal"
            self.relTol = 0
            self.badcharlist = [["$p_rgh"]]
        elif st == "U":
            self.dicttitle = "U"
            self.solver = "smoothSolver"
            self.smoother = "symGaussSeidel"                 
            self.tolerance = "1e-6"
            self.relTol = 0
    
    def dl(self) -> DictList:
        '''Gives a DictList object which is used for printing variables to file'''
        l = self.badcharlist
        for attr, value in self.__dict__.items():
            if attr!="dicttitle" and attr!="badcharlist":
                l.append([attr, value])
        return DictList(self.dicttitle, 0, l)