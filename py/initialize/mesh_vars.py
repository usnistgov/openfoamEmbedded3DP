#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging

# local packages

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class MeshVars:
    '''mesh variables'''
    
    def __init__(self, meshSize:float=0.2):
    
        self.fwreflev = 4
        self.xareflev = 2 # RG
        self.nCellsBetweenLevels = 5
        self.meshsize = meshSize

        ### snappyHexMesh
        # https://www.openfoam.com/documentation/user-guide/snappyHexMesh.php#x14-500004.4
        # castellatedMeshControls
        self.castellatedMesh = "true" 
            # Create the castellated mesh?
        self.allowFreeStandingZoneFaces = "true" 
            # Allow the generation of free-standing zone faces
        self.locationInMesh = "(0 0 0)" 
            # point which must be meshed, as opposed to inside walls
        self.maxLocalCells = 100000
            # Maximum number of cells per processor during refinement
        self.maxGlobalCells = 2000000 # *10 RG
            # Overall cell limit during refinement (i.e. before removal)
        self.minRefinementCells = 0 
            # If ≥ number of cells to be refined, surface refinement stops
        self.resolveFeatureAngle = 30
            # Applies maximum level of refinement to cells
            # that can see intersections whose angle exceeds this
            
                # snapControls
        self.snap = "true" 
            # Do the surface snapping stage?
        self.explicitFeatureSnap = "true" 
            # Use castellatedMeshControls features
        self.implicitFeatureSnap = "false"
            # Detect (geometric only) features by sampling the surface
        self.multiRegionFeatureSnap = "true" 
            # Detect features between multiple surfaces when using the explicitFeatureSnap
        self.nFeatureSnapIter = 10
            # Number of feature edge snapping iterations
        self.nRelaxIter = 5
            # Maximum number of snapping relaxation iterations
        self.nSmoothPatch = 3
            # Number of patch smoothing iterations before finding correspondence to surface
        self.nSolveIter = 100
            # Number of mesh displacement relaxation iterations
        self.tolerance = 2.0    
            # Ratio of distance for points to be attracted by 
            # surface feature point or edge, to local maximum edge length
            
        # addLayersControls
        self.addLayers = "true" 
            # Add surface layers?
        self.expansionRatio = 1
            # Expansion factor for layer mesh
        self.featureAngle = 30
            # Angle above which surface is not extruded
        self.finalLayerThickness = 0.3
            # Thickness of layer furthest from the wall, 
            # either relative or absolute according to the relativeSizes entry
        self.maxFaceThicknessRatio = 0.5
            # Face thickness ratio above which surface is not extruded,
            # useful for warped cells
        self.maxThicknessToMedialRatio = 0.3 
            # Reduce layer growth where ratio thickness to medial distance is large
        self.minMedianAxisAngle = 90
            # Angle used to pick up medial axis points
        self.minThickness = 0.25
            # Minimum overall thickness of all layers, below which surface is not extruded
        self.nBufferCellsNoExtrude = 0
            # Create buffer region for new layer terminations
        self.nGrow = 0
            # Number of layers of connected faces that are not grown 
            # if points are not extruded; helps convergence of 
            # layer addition close to features
        self.nLayerIter = 50
            # Overall max number of layer addition iterations
        self.nRelaxIter = 5
            # Maximum number of snapping relaxation iterations
        self.nRelaxedIter = 20
            # Max number of iterations after which the controls in the 
            # relaxed sub dictionary of meshQuality are used
        self.nSmoothNormals = 15
            # Number of smoothing iterations of interior mesh movement direction
        self.nSmoothSurfaceNormals = 10
            # Number of smoothing iterations of surface normals
        self.nSmoothThickness = 10
            # Smooth layer thickness over surface patches
        self.nSurfaceLayers = 10
            # number of surface layers for each patch
        self.relativeSizes = "true"
            # Are layer thicknesses relative to undistorted cell size 
            # outside layer or absolute?



        # meshQualityControls
        self.errorReduction = 0.75
            # Amount to scale back displacement at error points
        self.maxNonOrtho = 65
            # Maximum non-orthogonality allowed; 180 disables
        self.nSmoothScale = 4
            # Number of error distribution iterations

        self.mergeTolerance = "1E-6" # Merge tolerance as fraction of bounding box of initial mesh



        ### dynamicMeshDict
        # https://openfoamwiki.net/index.php/Parameter_Definitions_-_dynamicRefineFvMesh
        self.refineInterval = 5
            # how many iterations between mesh refinements
        self.field = "alpha.ink"
            # field to refine on
        self.lowerRefineLevel = 0.001
            # values below trigger refinement
        self.upperRefineLevel = 0.999
            # values above trigger unrefinement
        self.unrefineLevel = 5
            # number of times cells can be coarsened
        self.nBufferLayers = 1
            # number of layers around a refined cell
        self.maxRefinement = 4
            # max number of times cells can be refined
        self.maxCells = self.maxGlobalCells*10 # RG
            # total number of cells in the mesh
        self.dumpLevel = "false"
            # writes the refinement level for each cell as a volScalarField
            
    def varList(self, li) -> List[List[str]]:
        '''Given a list of variable names li, construct a table with those variables.'''
        
        out = [["", ""] for x in range(len(li))]
        for i, l in enumerate(li):
            out[i] = [l, getattr(self, l)]
        return out

        
    def cmc(self) -> List[List[str]]:
        '''generates a list of variables for the castellatedMeshControls file'''
        
        return(self.varList(["maxLocalCells", "maxGlobalCells", "minRefinementCells", "nCellsBetweenLevels", "resolveFeatureAngle", "locationInMesh", "allowFreeStandingZoneFaces"]))
       
    
    def sc(self) -> List[List[str]]:
        '''generates a list of variables for the snapControls file'''
        
        return(self.varList(["nSmoothPatch", "tolerance", "nSolveIter", "nRelaxIter", "nFeatureSnapIter", "implicitFeatureSnap", "explicitFeatureSnap", "multiRegionFeatureSnap"]))
    
    

    def alc(self) -> List[List[str]]:
        '''generates a list of variables for the addLayersControls table in dynamicMeshDict'''
        
        return(self.varList(["relativeSizes", "expansionRatio", "finalLayerThickness", "minThickness", "nGrow",\
                             "featureAngle", "nRelaxIter", "nSmoothSurfaceNormals", "nSmoothNormals",\
                             "nSmoothThickness", "maxFaceThicknessRatio", "maxThicknessToMedialRatio",\
                             "minMedianAxisAngle", "nBufferCellsNoExtrude", "nLayerIter", "nRelaxedIter"]))
    
    def mqc(self) -> List[List[str]]:
        '''generates a list of variables for the meshQualityControls table in dynamicMeshDict'''
        
        return(self.varList(["nSmoothScale", "errorReduction"]))
    
    def dmd(self) -> List[List[str]]:
        '''generates a list of variables for the dynamicMeshDictionary table in dynamicMeshDict'''
        
        return(self.varList(["refineInterval", "field", "lowerRefineLevel", "upperRefineLevel", "unrefineLevel",\
                             "nBufferLayers", "maxRefinement", "maxCells", "dumpLevel"]))
    
    