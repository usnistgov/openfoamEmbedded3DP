#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import os,sys
import logging
import numpy as np

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from block import Block
from block_points import blockPts
from boundary_input import BoundaryInput
from cd_vars import CDVars
from compile_0 import compileAlphaOrig, compileU, compileP, compileCellLevel, compilePointLevel
from compile_all_run import compileSlurm, compileAllAllRun, compileAllClean, compileAllAllClean, compileAllRunMesh, compileAllRun
from compile_block_mesh_dict import compileBlockMeshDict
from compile_control_dict import compileControlDict
from compile_dynamic_mesh_dict import compileDynamicMeshDict
from compile_fv import compileFvSchemes, compileFvSolution
from compile_g import compileG
from compile_mesh_quality_dict import compileMeshQualityDict
from compile_set_fields_dict import compileSetFieldsDict
from compile_snappy_hex_mesh_dict import compileSnappyHexMeshDict
from compile_surface_features_dict import compileSurfaceFeaturesDict
from compile_surface_feature_extract_dict import compileSurfaceFeatureExtractDict
from compile_transport_properties import compileTransportProperties
from compile_turbulence_properties import compileTurbulenceProperties
from file_group import FileGroup
from fluid import Fluid
from fv_vars import FVVars
from geometry_file import geometryFile
from mesh_vars import MeshVars
from noz_vars import NozVars
from real_boundaries import realBoundaries



# logging
logging.basicConfig(level=logging.INFO)

#-------------------------------------------------------------------------------------------------  

class fileCreator:
    '''Generate all files but transport files.
    ii is either a folder name (e.g. 'folderA') or a number (e.g. 25) to be appended to the folder name
    exportMesh true to export a mesh folder inside of this simulation folder
    onlyMesh true to only export a mesh folder
    topFolder is the parent of this simulation folder
    folderBase is the first part of the folder name, if ii is not a string
    startTime is the simulation start time in s
    endTime is the simulation end time in s
    dt is the initial time step in s
    writeDt is the time step for saving results to file in s
    solver could be interFoam or interIsoFoam
    Additional keyword vars are passed into NozVars. Examples are bathWidth=16, vink=10
    '''
    
    def __init__(self, ii:Union[str, float]
                 , topFolder:str
                 , exportMesh:bool=False
                 , onlyMesh:bool=False
                 , folderBase:str="nb"
                 , startTime:float=0
                 , endTime:float=2.5
                 , dt:float=0.001
                 , writeDt:float=0.1
                 , solver:float="interFoam"
                 , **kwargs):
        self.ii = ii
        self.topFolder = topFolder
        self.exportMesh = exportMesh
        self.onlyMesh = onlyMesh
        self.folderBase = folderBase
        self.startTime = startTime
        self.endTime = endTime
        self.dt = dt
        self.writeDt = writeDt
        self.solver = solver
    
        self.folderName()
        self.initialize(**kwargs)
        self.createNozzleBlockFile(**kwargs)

        self.solverObjects()
        self.compileSolverFiles()
    
    def labels(self, ink:Fluid, sup:Fluid) -> str:
        '''get a file that stores the ink and support labels'''
        return f'ink,{ink.label}\nsup,{sup.label}'
    
    def addTransportProperties(self, ink:Fluid, sup:Fluid, sigma:float) -> None:
        ''''add the fluid properties'''
        self.fg.labels = self.labels(ink, sup)
        self.fg.transportProperties = compileTransportProperties(ink.transportGroup('ink'), sup.transportGroup('sup'), sigma).s
        
    def export(self):
        self.fg.exportAllFiles()
        
    def plot(self):
        self.fg.makePlot()
    
    def folderName(self) -> None:
        '''determine the name of the new folder to be created'''
        if not os.path.exists(self.topFolder):
            raise FileNotFoundError(f'Folder {self.topFolder} does not exist')

        if self.onlyMesh:
            self.folder = self.topFolder
        else:
            if type(self.ii) is str:
                self.folder = os.path.join(self.topFolder, self.ii)
            else:
                self.folder = os.path.join(self.topFolder, f'{self.folderBase}{self.ii}')
                
    def initialize(self, **kwargs) -> None:
        '''initialize the NozVars and MeshVars objects'''
        self.geo = NozVars(**kwargs)
        if 'meshSize' in kwargs:
            self.mv = MeshVars(meshSize=kwargs['meshSize'])
        else:
            self.mv = MeshVars(meshSize=self.geo.niw/3) # make the initial mesh size 1/3 of the inner nozzle diameter
    #     mv.locationInMesh = f'({geo.ncx} {geo.ncy} {geo.nbo})'
        

    def createNozzleBlockFile(self, **kwargs) -> None:
        '''gets the text for most of the files in the folder
        exportMesh is true if we want to create a mesh folder'''
        self.fg = FileGroup(self.folder, exportMesh=self.exportMesh, onlyMesh=self.onlyMesh, **kwargs)
        self.fg.geofile = geometryFile(self.geo).s
        self.br = realBoundaries(self.geo, self.exportMesh, **kwargs).boundaryList # this is a real boundary list for the imported stl boundaries
        self.createSystem()
        self.createConstant()
        self.create0()
        
        if self.exportMesh:
            self.addMesh()
            #plot
            self.fg.geo = self.geo
            self.fg.pl = self.pl
            self.fg.blocks = self.blocks
            self.fg.br = self.br 
            
            
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
    
    def addSystemMesh(self) -> None:
        '''create mesh files'''
        self.pl = self.geo.createNozzle() # point list
        bcs = self.geo.blockCornerList() # list of block corners
        self.blocks = []
        for bc in bcs:
            blocki = Block()
            blocki.vertices = blockPts(self.pl, bc)
            msz = self.mv.meshsize # mesh size
            blocki.meshi = [round(self.geo.bw/msz), round(self.geo.bh/msz), round(self.geo.bd/msz)]
            self.blocks.append(blocki) 
        bl = self.boundaryList(self.blocks) # this is a dummy boundary list for the original blockMeshDict
        self.fg.blockMeshDict = compileBlockMeshDict(self.pl, self.blocks, bl).s
        self.fg.fvSchemesMesh = compileFvSchemes("").s
        self.fg.fvSolutionMesh = compileFvSolution("").s
        
    def createSystem(self) -> None:
        '''create files in system folder'''
        self.fg.setFieldsDict = compileSetFieldsDict(self.geo).s
        if self.exportMesh:
            self.addSystemMesh()
            
    def createConstant(self) -> None:
        '''create files in the constant folder'''
        self.fg.dynamicMeshDict = compileDynamicMeshDict(self.mv).s
        self.fg.g = compileG().s
        self.fg.turbulenceProperties = compileTurbulenceProperties().s
        
    def addMesh(self) -> None:
        '''create mesh files for the 0 folder'''
        self.fg.snappyHexMeshDict = compileSnappyHexMeshDict(self.br, self.mv).s
        self.fg.surfaceFeatureExtractDict = compileSurfaceFeatureExtractDict(self.br).s
        self.fg.surfaceFeaturesDict = compileSurfaceFeaturesDict(self.br).s
        self.fg.meshQualityDict = compileMeshQualityDict().s
        self.fg.cellLevel = compileCellLevel(self.br).s
        self.fg.pointLevel = compilePointLevel(self.br).s
        
        
    def create0(self) -> None:
        '''create files in the 0 folder'''
        self.fg.meshes = self.br # keep the real boundaries for exporting
        self.fg.alphainkorig = compileAlphaOrig(self.br).s
        self.fg.U = compileU(self.br).s
        self.fg.prgh = compileP(self.br).s

    def solverObjects(self) -> None:
        '''gets the control dictionary variable object and fvsolution and fvschemes dictionary variable object
        starttime is the start time for the simulation in s
        endtime is the end time for the simulation in s
        dt is the initial solve time step in s
        writedt is the write time step in s
        solver is interFoam or interIsoFoam'''
        self.cdv = CDVars(self.startTime, self.endTime, self.dt, self.writeDt)
        self.cdv.application = self.solver
        self.fvv = FVVars(self.solver)
        
    def compileSolverFiles(self):
        '''generates the text for the controlDict, FVsolution, FVschemes, and bash scripts'''
        self.fg.controlDict = compileControlDict(self.cdv).s
        self.fg.fvSchemes = compileFvSchemes(self.fvv).s
        self.fg.fvSolution = compileFvSolution(self.fvv).s
        self.fg.slurm = compileSlurm(self.fg.folder, self.fg.slurmFolder).s
        self.fg.allallrun = compileAllAllRun().s
        self.fg.allclean = compileAllClean(self.cdv.endTime, self.cdv.writeInterval).s
        self.fg.allallclean = compileAllAllClean().s
        self.fg.allrunmesh = compileAllRunMesh(self.fg.folder).s
        self.fg.allrun = compileAllRun(self.fg.folder, self.cdv.application).s
    #     self.fg.cont = compileContinue(self.fg.folder, cdv.application)
        self.cdv.application = "icoFoam"
        self.fg.controlDictMesh = compileControlDict(self.cdv).s
