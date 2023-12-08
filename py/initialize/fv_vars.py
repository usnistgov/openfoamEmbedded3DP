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
from fv_sol_grp import FVSolGrp
from dict_list import DictList

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class FVVars: 
    '''Holds all fvsolution and fvschemes variables'''
    

    def __init__(self, solver:str) -> None:
        '''Input: solver is the type of solver being used: 'interFoam' or 'interIsoFoam' '''
        self.solver = solver
        
        # FVSchemes
        self.ddtSchemesDefault = "Euler"
        self.gradSchemesDefault = "Gauss linear"
        self.laplacianSchemesDefault = "Gauss linear corrected"
        self.interpolationSchemesDefault = "linear"
        self.snGradSchemesDefault = "corrected"
        self.divSchemes = [["div(rhoPhi,U) Gauss linearUpwind grad(U)"], ["div(phi,alpha) Gauss vanLeer"], \
                      ["div(phirb,alpha) Gauss linear"], ["div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear"]]
        self.slist = ["ddtSchemes", "gradSchemes", "laplacianSchemes", "interpolationSchemes", "snGradSchemes"]
        self.postlist = ["divSchemes"]
        if solver=="interIsoFoam":
            self.postlist.append("fluxRequired")
            self.fluxRequired = [["default", "no"], ["p_rgh"], ["pcorr"], ["alpha.ink"]]

        # FVSolution
        self.alphaink = FVSolGrp("alphaink", solver)
        self.pcorr = FVSolGrp("pcorr", solver)
        self.prgh = FVSolGrp("prgh", solver)
        self.prghfinal = FVSolGrp("prghfinal", solver)
        self.U = FVSolGrp("U", solver)
        
        # PIMPLE
        self.momentumPredictor = "no"
        self.nOuterCorrectors = 1
        self.nCorrectors = 3
        self.nNonOrthogonalCorrectors = 0
    

    def fvSchemeList(self) -> List[DictList]:
        '''format the fvSchemes variables into a list of DictLists so they are ready to print to file'''
        l = []       
        for v in self.slist:
            l.append(DictList(v, 0, [["default", getattr(self, v+"Default")]]))
        for v in self.postlist:
            l.append(DictList(v, 0, getattr(self, v)))
        return l


    def solverlist(self) -> DictList:
        '''for fvsolutions, get a list of DictLists created by the FVSolGrp objects so they are ready to print to file'''
        l = []
        for o in [self.alphaink, self.pcorr, self.prgh, self.prghfinal, self.U]:
            l.append(o.dl())
        return DictList("solvers",  0, l)

    
    def pimple(self) -> DictList:
        '''get a DictList for the PIMPLE variables'''
        l = []
        for o in ["momentumPredictor", "nOuterCorrectors", "nCorrectors", "nNonOrthogonalCorrectors"]:
            l.append([o, getattr(self, o)])
        return DictList("PIMPLE", 0, l)