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
from tools.config import cfg
from block import Block
from boundary_input import BoundaryInput
import export as ex
import folder_scraper as fs
from noz_vars import NozVars
from file_plotter import filePlotter


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  

class FileGroup:
    '''This holds all of the strings that get outputted to text files and the meshes used to generate stls'''

    def __init__(self, folder:str, exportMesh:bool=False, onlyMesh:bool=False, **kwargs):
        '''Input is the folder that all of these files will go into'''
        
        self.exportMesh = exportMesh
        self.onlyMesh = onlyMesh
        self.folder = folder
        self.geofile = ""
        self.allclean = "" # for both mesh and case folders
        self.allallrun = "" # for the folder above mesh and case
        self.allrun = ""
        self.allrunmesh = ""
        self.cont = ""

        self.alphainkorig = ""
        self.prgh = ""
        self.U = ""

        self.g = ""
        self.transportProperties = ""
        self.turbulenceProperties = ""
        self.dynamicMeshDict = ""

        self.blockMeshDict = ""
        self.controlDict = ""
        self.controlDictMesh = ""
        self.fvSchemes = ""
        self.fvSchemesMesh = ""
        self.fvSolution = ""
        self.fvSolutionMesh = ""
        self.setFieldsDict = ""

        self.cellLevel = ""
        self.pointLevel = ""
        self.meshQualityDict = ""
        self.snappyHexMeshDict = ""
        self.surfaceFeatureExtractDict = ""
        self.surfaceFeaturesDict = ""

        if 'slurmFolder' in kwargs:
            self.slurmFolder = kwargs['slurmFolder']
        else:
            self.slurmFolder = os.path.join(cfg.path.slurmFolder, os.path.basename(folder))

        self.plot = ""

        self.meshes = []

    def exportAllFiles(self) -> str:
        '''exports all of the files that we generated. exportMesh is true to export mesh files in this folder'''
        f = self.folder
        folderList = [f]
        
        if not self.onlyMesh:
            casef = os.path.join(f, "case")
            f0 = os.path.join(casef, "0")
            fconst = os.path.join(casef,"constant")
            fsyst = os.path.join(casef,"system")
            folderList = folderList + [casef, f0, fconst, fsyst]
             
        if self.exportMesh:
            fgeom = os.path.join(f,"geometry")
            fmesh = os.path.join(f,"mesh")
            fmeshconst = os.path.join(fmesh,"constant")
            fmeshconsttri = os.path.join(fmeshconst,"triSurface")
            fmeshsys = os.path.join(fmesh,"system")
            fmesh0 = os.path.join(fmesh,"0")
            folderList = folderList + [fgeom, fmesh, fmeshconst, fmeshconsttri, fmeshsys, fmesh0]

        list(map(ex.mkdirif, folderList)) # create folders

        
        if not self.onlyMesh:
#             ex.exportFile(f, 'labels.csv', self.labels)
            ex.exportFile(f, 'run.slurm', self.slurm, linux=True)
            ex.exportFile(f, "Allclean.sh", self.allallclean, linux=True) 
            ex.exportFile(casef, "Allclean.sh", self.allclean, linux=True) 
#             ex.exportFile(casef, "Allrun", self.allrun)
            ex.exportFile(casef, 'Allrun.sh', self.allrun, linux=True)
#             ex.exportFile(casef, "Continue", self.cont) 
            ex.exportFile(f0, "alpha.ink.orig", self.alphainkorig) 
            ex.exportFile(f0, "alpha.ink", self.alphainkorig) 
            ex.exportFile(f0, "p_rgh", self.prgh) 
            ex.exportFile(f0, "U", self.U) 
            ex.exportFile(fconst, "g", self.g) 
            ex.exportFile(fconst, "transportProperties", self.transportProperties) 
            ex.exportFile(fconst, "turbulenceProperties", self.turbulenceProperties) 
            ex.exportFile(fconst, "dynamicMeshDict", self.dynamicMeshDict) 
            ex.exportFile(fsyst, "controlDict", self.controlDict)     
            ex.exportFile(fsyst, "fvSchemes", self.fvSchemes)     
            ex.exportFile(fsyst, "fvSolution", self.fvSolution)    
            ex.exportFile(fsyst, "setFieldsDict", self.setFieldsDict) 

        if self.exportMesh:
            if not self.onlyMesh:
                ex.exportFile(f, "Allrun.sh", self.allallrun, linux=True)
                ex.exportFile(f, "geometry.csv", self.geofile)
#             ex.exportFile(fmesh, "Allclean.sh", self.allclean, linux=True) 
            ex.exportFile(fmesh, "Allrun.sh", self.allrunmesh, linux=True)
            ex.exportFile(fmesh0, "pointLevel", self.pointLevel)
            ex.exportFile(fmesh0, "cellLevel", self.cellLevel)
            ex.exportFile(fmeshsys, "blockMeshDict", self.blockMeshDict) 
            ex.exportFile(fmeshsys, "controlDict", self.controlDictMesh) 
            ex.exportFile(fmeshsys, "fvSchemes", self.fvSchemesMesh) 
            ex.exportFile(fmeshsys, "fvSolution", self.fvSolutionMesh)
            ex.exportFile(fmeshsys, "meshQualityDict", self.meshQualityDict)
            ex.exportFile(fmeshsys, "snappyHexMeshDict", self.snappyHexMeshDict)
            ex.exportFile(fmeshsys, "surfaceFeatureExtractDict", self.surfaceFeatureExtractDict)
            ex.exportFile(fmeshsys, "surfaceFeaturesDict", self.surfaceFeaturesDict)
            ex.saveStls(fgeom, self.meshes)
            ex.saveStls(fmeshconsttri, self.meshes)

        if not self.onlyMesh:
            fs.populate(self.folder)
                 
                           
    def makePlot(self, elev:int=20, azim:int=145, showMesh:bool=False
                 , vertexNumbers:bool=False, blockNumbers:bool=False
                 , fs:int=10, figsize:float=10, **kwargs):
        '''plot the requested geometry
        geo is a nozzleVars object
        pl is a point list used for blockMesh
        blocks is a list of blocks used for blockMesh
        br is real boundaries used for generating stls and setting boundary conditions'''
        self.plotter = filePlotter(self.geo, self.pl, self.blocks, self.br)
        self.plotter.makePlot(**kwargs)
        self.plot = self.plotter.fig
