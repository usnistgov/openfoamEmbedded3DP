#### ncreate3d.py ####
# these functions are used to generate OpenFOAM input files for a nozzle in a 3D bath


import math
import numpy as np
from stl import mesh
import matplotlib.pyplot as plt
import statistics as st
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import functools
import os
import sys
sys.path.append(r'C:\Users\lmf1\Documents\OpenFOAM')
import folderparser as fp
import folderscraper as fs
from typing import List, Dict, Tuple, Union, Any, TextIO

#-------------------------------------------------
#---------------GLOBAL VARIABLES------------------
#-------------------------------------------------

SCALE = 0.001 # scale all units by this amount * m

#-------------------------------------------------
#--------------------CLASSES----------------------
#-------------------------------------------------

# DictList is an object that ensures uniform formatting of variables into files
class DictList:
        # title is the title of the variable group
        # form is 0 for bracketed entries e.g. group definitions, 1 for parentheses e.g. point lists
        # proplist is a list of properties that we're storing within this group
    def __init__(self, title:str, form:int, proplist:List[str]):
        self.title = title
        self.proplist = proplist
        if form==0: 
            self.pl = "{" # parens (group) left
            self.pr = "}" # parens (group) right
            self.sp = "\t" # spacing
            self.rl = "" # row left
            self.rr = ";" # row right
        else:
            self.pl = "("
            self.pr = ");"
            self.sp = " "
            self.rl = "("
            self.rr = ")"
    
    # format the title and variables into the format that OpenFOAM likes
        # level is the number of tabs before the title of the group
    def prnt(self, level:int):
        if level<0:
            s = ""
        else:
            s = tabs(level) + self.title + "\n" + tabs(level) + self.pl + "\n"
        for p in self.proplist:
            if isinstance(p, DictList):  # if this property is a DictList, recurse
                s = s + p.prnt(level + 1)
            else: # row level
                s = s + tabs(level + 1)
                if isinstance(p, str): # single item
                    s = s + p + "\n"
                else: # list of items
                    if len(p)>0:
                        # one row
                        s = s + self.rl
                        for pi in p[0:-1]:
                            s = s + str(pi) + self.sp
                        s = s + str(p[-1]) + self.rr + "\n"
                    else:
                        s = s + "\n"
                if level<0:
                    s = s + "\n"
        if level>=0:
            s = s + tabs(level) + self.pr + "\n\n"
        return s

#--------------------------------------------------------------------------------------------------------
# mesh variables
class MeshVars:
    fwreflev = 2
    nCellsBetweenLevels = 10
    meshsize = 0.2
    
    ### snappyHexMesh
    # https://www.openfoam.com/documentation/user-guide/snappyHexMesh.php#x14-500004.4
    # castellatedMeshControls
    castellatedMesh = "true" 
        # Create the castellated mesh?
    allowFreeStandingZoneFaces = "true" 
        # Allow the generation of free-standing zone faces
    locationInMesh = "(0 0 0)" 
        # Merge tolerance as fraction of bounding box of initial mesh
    maxLocalCells = 100000 
        # Maximum number of cells per processor during refinement
    maxGlobalCells = 2000000 
        # Overall cell limit during refinement (i.e. before removal)
    minRefinementCells = 0 
        # If ≥ number of cells to be refined, surface refinement stops
    resolveFeatureAngle = 30 
        # Applies maximum level of refinement to cells
        # that can see intersections whose angle exceeds this
        
    def cmc(self):
        return(self.varlist(["maxLocalCells", "maxGlobalCells", "minRefinementCells", "nCellsBetweenLevels", "resolveFeatureAngle", "locationInMesh", "allowFreeStandingZoneFaces"]))
       
    # snapControls
    snap = "true" 
        # Do the surface snapping stage?
    explicitFeatureSnap = "true" 
        # Use castellatedMeshControls features
    implicitFeatureSnap = "false" 
        # Detect (geometric only) features by sampling the surface
    multiRegionFeatureSnap = "true" 
        # Detect features between multiple surfaces when using the explicitFeatureSnap
    nFeatureSnapIter = 10 
        # Number of feature edge snapping iterations
    nRelaxIter = 5
        # Maximum number of snapping relaxation iterations
    nSmoothPatch = 3 
        # Number of patch smoothing iterations before finding correspondence to surface
    nSolveIter = 100 
        # Number of mesh displacement relaxation iterations
    tolerance = 2.0    
        # Ratio of distance for points to be attracted by 
        # surface feature point or edge, to local maximum edge length
        
    def sc(self):
        return(self.varlist(["nSmoothPatch", "tolerance", "nSolveIter", "nRelaxIter", "nFeatureSnapIter", "implicitFeatureSnap", "explicitFeatureSnap", "multiRegionFeatureSnap"]))
    
    # addLayersControls
    addLayers = "false" 
        # Add surface layers?
    expansionRatio = 1 
        # Expansion factor for layer mesh
    featureAngle = 30 
        # Angle above which surface is not extruded
    finalLayerThickness = 0.3 
        # Thickness of layer furthest from the wall, 
        # either relative or absolute according to the relativeSizes entry
    maxFaceThicknessRatio = 0.5 
        # Face thickness ratio above which surface is not extruded,
        # useful for warped cells
    maxThicknessToMedialRatio = 0.3 
        # Reduce layer growth where ratio thickness to medial distance is large
    minMedianAxisAngle = 90 
        # Angle used to pick up medial axis points
    minThickness = 0.25 
        # Minimum overall thickness of all layers, below which surface is not extruded
    nBufferCellsNoExtrude = 0 
        # Create buffer region for new layer terminations
    nGrow = 0 
        # Number of layers of connected faces that are not grown 
        # if points are not extruded; helps convergence of 
        # layer addition close to features
    nLayerIter = 50 
        # Overall max number of layer addition iterations
    nRelaxIter = 5 
        # Maximum number of snapping relaxation iterations
    nRelaxedIter = 20 
        # Max number of iterations after which the controls in the 
        # relaxed sub dictionary of meshQuality are used
    nSmoothNormals = 3 
        # Number of smoothing iterations of interior mesh movement direction
    nSmoothSurfaceNormals = 1 
        # Number of smoothing iterations of surface normals
    nSmoothThickness = 10
        # Smooth layer thickness over surface patches
    nSurfaceLayers = 10
        # number of surface layers for each patch
    relativeSizes = "true"
        # Are layer thicknesses relative to undistorted cell size 
        # outside layer or absolute?
    
    
    
    # meshQualityControls
    errorReduction = 0.75
        # Amount to scale back displacement at error points
    maxNonOrtho = 65
        # Maximum non-orthogonality allowed; 180 disables
    nSmoothScale = 4
        # Number of error distribution iterations

    mergeTolerance = "1E-6" # Merge tolerance as fraction of bounding box of initial mesh
    
    
    
    ### dynamicMeshDict
    # https://openfoamwiki.net/index.php/Parameter_Definitions_-_dynamicRefineFvMesh
    refineInterval = 5
        # how many iterations between mesh refinements
    field = "alpha.ink"
        # field to refine on
    lowerRefineLevel = 0.001
        # values below trigger refinement
    upperRefineLevel = 0.999
        # values above trigger unrefinement
    unrefineLevel = 5
        # number of times cells can be coarsened
    nBufferLayers = 1
        # number of layers around a refined cell
    maxRefinement = 4
        # max number of times cells can be refined
    maxCells = maxGlobalCells
        # total number of cells in the mesh
    dumpLevel = "false"
        # writes the refinement level for each cell as a volScalarField
    
    ########### functions
    # add layers table
    def alc(self):
        return(self.varlist(["relativeSizes", "expansionRatio", "finalLayerThickness", "minThickness", "nGrow",\
                             "featureAngle", "nRelaxIter", "nSmoothSurfaceNormals", "nSmoothNormals",\
                             "nSmoothThickness", "maxFaceThicknessRatio", "maxThicknessToMedialRatio",\
                             "minMedianAxisAngle", "nBufferCellsNoExtrude", "nLayerIter", "nRelaxedIter"]))
    
    # mesh quality table
    def mqc(self):
        return(self.varlist(["nSmoothScale", "errorReduction"]))
    
    # dynamic mesh table
    def dmd(self):
        return(self.varlist(["refineInterval", "field", "lowerRefineLevel", "upperRefineLevel", "unrefineLevel",\
                             "nBufferLayers", "maxRefinement", "maxCells", "dumpLevel"]))
    
    # given a list of variable names, construct a table with those variables
        # li is a list of variable names (strings)
    def varlist(self, li):
        out = [["", ""] for x in range(len(li))]
        for i, l in enumerate(li):
            out[i] = [l, getattr(self, l)]
        return out

    #--------------------------------------------------------------------------------------------------------
# controlDict variables    
class CDVars:
   
        # starttime is in seconds
        # endtime is in seconds
        # dt is the initial time step of the simulation in seconds
        # writedt is the time step at which to write results to file
    def __init__(self, starttime, endtime, dt, writedt):
        # solver
        self.application = "interFoam"
        
        # time control
        self.startFrom = "latestTime"
        self.startTime = starttime
        self.stopAt = "endTime"
        self.endTime = endtime
        self.deltaT = dt # Time step of the simulation
        
        # time step control
        self.adjustTimeStep = "no" # yes/no† to adjust time step according to maximum Courant number in transient simulation.
        self.maxCo = 1 # Maximum Courant number allowed.
        self.maxAlphaCo = 1
        self.maxDeltaT = 1 # Maximum time step allowed in transient simulation
        
        # data writing
        #self.writeControl = "adjustable" # Writes data every writeInterval seconds of simulated time.
        self.writeControl = "adjustableRunTime"
        self.writeInterval = writedt # Scalar used in conjunction with writeControl described above.
        self.purgeWrite = 0 # Integer representing a limit on the number of time directories that are stored by overwriting time directories on a cyclic basis. Example of t0 = 5s, Δt = 1s and purgeWrite 2;: data written into 2 directories, 6 and 7, before returning to write the data at 8 s in 6, data at 9 s into 7, etc. To disable the time directory limit, specify purgeWrite 0;† For steady-state solutions, results from previous iterations can be continuously overwritten by specifying purgeWrite 1;
        self.writeFormat = "ascii" # ASCII format, written to writePrecision significant figures.
        self.writePrecision = 6 # Integer used in conjunction with writeFormat described above, 6† by default
        self.writeCompression = "uncompressed" # Specifies the compression of the data files.
        self.timeFormat = "general" # Specifies scientific format if the exponent is less than -4 or greater than or equal to that specified by timePrecision.
        self.timePrecision = 6 # Integer used in conjunction with timeFormat described above, 6† by default
        
        # data reading
        self.runTimeModifiable = "yes" # yes†/no switch for whether dictionaries, e.g.controlDict, are re-read by OpenFOAM at the beginning of each time step.

    # vlist constructs a table with all of the variables held in this object
    # the table has two columns: the variable name and the value
    def vlist(self):
        l = []        
        for attr, value in self.__dict__.items():
            l.append([attr, value])
        return l

    #--------------------------------------------------------------------------------------------------------
# fvsolution variables for a given solve variable
class FVSolGrp:
        # st is a string that tells us what variable is being solved for, e.g. 'alphaink', 'pcorr', 'prgh', 'prghfinal', or 'U'
        # solver is the type of solver being used: 'interFoam' or 'interIsoFoam'
    def __init__(self, st:str, solver:str):
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
    
    # dl gives a DictList object which is used for printing variables to file
    def dl(self):
        l = self.badcharlist
        for attr, value in self.__dict__.items():
            if attr!="dicttitle" and attr!="badcharlist":
                l.append([attr, value])
        return DictList(self.dicttitle, 0, l)

    #--------------------------------------------------------------------------------------------------------
# all fvsolution and fvschemes variables
class FVVars: 
        # solver is the type of solver being used: 'interFoam' or 'interIsoFoam'
    def __init__(self, solver:str):
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
    
    # format the fvSchemes variables into a list of DictLists so they are ready to print to file
    def fvSchemeList(self):
        l = []       
        for v in self.slist:
            l.append(DictList(v, 0, [["default", getattr(self, v+"Default")]]))
        for v in self.postlist:
            l.append(DictList(v, 0, getattr(self, v)))
        return l

    # for fvsolutions, get a list of DictLists created by the FVSolGrp objects so they are ready to print to file
    def solverlist(self):
        l = []
        for o in [self.alphaink, self.pcorr, self.prgh, self.prghfinal, self.U]:
            l.append(o.dl())
        return DictList("solvers",  0, l)

    # get a DictList for the PIMPLE variables
    def pimple(self):
        l = []
        for o in ["momentumPredictor", "nOuterCorrectors", "nCorrectors", "nNonOrthogonalCorrectors"]:
            l.append([o, getattr(self, o)])
        return DictList("PIMPLE", 0, l)

    #--------------------------------------------------------------------------------------------------------
# geometry variables
# this is where the geometry of the nozzle is defined!

class NozVars:
    # ALL VALUES IN MM  
        # wm bath width multiplier
        # hm bath height multiplier
        # dm bath depth multiplier
        # fm front of nozzle bath width multiplier
        # v print speed
        # npts number of poitns in the circle used to define the nozzle
    def __init__(self, wm:float, hm:float, dm:float, fm:float, v:float, npts:int):
        self.niw = 0.603 # nozzle inner width
        self.nt = 0.152 # nozzle thickness
        self.bw = wm*self.niw # bath width (x)
        self.bh = hm*self.niw # bath height (y)
        self.bd = dm*self.niw # bath depth (z)
        self.nl = self.bd/2-self.niw/2
        
        self.ble = -self.bw/2 # bath left coord
        self.bri = self.bw/2 # bath right
        self.bfr = -self.bh/2 # bath front
        self.bba = self.bh/2 # bath back
        self.bbo = -self.bd/2 # bath bottom
        self.bto = self.bd/2 # bath top
        self.nbo = self.bto - self.nl # nozzle bottom coord
        self.ncx = self.ble + fm*self.niw # nozzle center bottom x coord
        self.ncy = self.bfr + self.bh/2 # nozzle center bottom y coord
        
        # here, we define the geometry of the nozzle using four circles: 
        # the inner and outer edges of the inlet (top) and outlet (bottom) of the nozzle
        # if the nozzle is a cylinder, the inner edge of the top is the same as the inner edge of the bottom, and the outer edge at the top is the same as the outer edge of the bottom
        # if the nozzle is a cone, the points at the top might cover a bigger circle!
        
        # the way that snappyHexMesh works, you give it stl files to define the geometry of the boundaries
        # that means that our cylinder actually needs to be made of triangles. 
        # You can do that by defining your circles at the top and bottom using a discrete number of points.
        # We use npts=50 in noz3dscript.ipynb
        # we also add the center point of the nozzle, which you need to create triangles for the inkFlow boundary
        # all of these lists of points are in x,y coords
        
        
        bottomRadius = self.niw/2
        topRadius = bottomRadius
        
        self.inptsb = circlePoints(bottomRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the bottom of the nozzle
        self.outptsb = circlePoints(bottomRadius + self.nt, npts)+[self.ncx, self.ncy] 
            # outer radius points at the bottom
        self.inptst = circlePoints(topRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the top of the nozzle (the inlet)
        self.outptst = circlePoints(topRadius + self.nt, npts)+[self.ncx, self.ncy]
            # outer radius points at the top of the nozzle

        self.bv = v*SCALE # bath velocity
        self.iv = v*SCALE*bottomRadius/topRadius # ink velocity

        #--------------------------------------------------------------------------------------------------------
# this holds all of the strings that get outputted to text files, 
# as well as the meshes used to generate stls
class FileGroup:
    def __init__(self, folder:str):
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

        self.g = compileg()
        self.transportProperties = ""
        self.turbulenceProperties = compileTurbulenceProperties()
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
        self.meshQualityDict = compileMeshQualityDict()
        self.snappyHexMeshDict = ""
        self.surfaceFeatureExtractDict = ""

        self.plot = plt.figure

        self.meshes = []
        
    def makePlot(self):
        self.plot = plotGeo(self.geo, self.pl, self.blocks, self.br)
        
    
#--------------------------------------------------------------------------------------------------------        
######## blocks and boundaries
        
# stores information about a block during blockMesh 
class Block:
    def __init__(self):
        self.vertices = []
        self.meshi = [1,1,1]
        self.grading = [1,1,1]
        
# stores information about a boundary 
class BoundaryInput:
    def __init__(self, labelin:str, typin:str):
        self.label = labelin # name of the boundary e.g. outerWall
        self.typ = typin # name of boundary type e.g. patch for blockMeshDict
        self.flist = [] # list of faces, which are lists of point indices
        self.alphalist = [] # list of alpha.ink properties
        self.Ulist = [] # list of U properties
        self.plist = [] # list of p_rgh properties
        self.meshi = [] # mesh for exporting to stl
        self.reflev = 0
      
  #--------------------------------------------------------------------------------------------------------  
#-------------------------------------------------
#--------------------FUNCTIONS--------------------
#-------------------------------------------------
    # bi is a BoundaryInput object
def allalpha(bi:BoundaryInput) -> List:
    return bi.alphalist

    # bi is a BoundaryInput object
def allU(bi:BoundaryInput) -> List:
    return bi.Ulist
 
    # bi is a BoundaryInput object
def allp(bi:BoundaryInput) -> List:
    return bi.plist  

# return n tabs in a row
    # n is the number of tabs
def tabs(n:int) -> str:
    s = ""
    for i in range(n):
        s = s + "\t"
    return s

# iterate over all items stored in an object
    # obj an object
    # f is a function with 2 inputs: the attribute name and the value
def classIterate(obj, f):
     for attr, value in obj.__dict__.items():
        f(attr, value)
        
# calculate the grading for the mesh within a block
    # W is the total width of the block, 
    # R is the ratio of largest to smallest cell
    # t is the smallest cell size        
def gradingcalc(self, W:float, R:float, t:float) -> float:
    er = (W-t)/(W-R*t)
    n = round(math.log(1-W*(1-er)/t, er))
    return n


#-------------------------------------------------
#----------------CREATING GEOMETRIES--------------
# get a list of points in a circle
    # r is the radius of the circle
    # npts is the number of points
def circlePoints(r:float, npts:int) -> List:
    pts = np.zeros([npts+1, 2])
    tstep = 2*math.pi/(npts)
    theta = 0
    for i in range(npts):
        pts[i] = ([math.cos(theta)*r, math.sin(theta)*r])
        theta = theta + tstep
    pts[npts] = pts[0]
    return pts
    
# get list of vertices in the block mesh
    # nv is a NozVars object
def createnozzle(nv:NozVars) -> List:
    xlist = np.array([nv.ble-0.01,  nv.bri+0.01])
    ylist = np.array([nv.bfr-0.01, nv.bba+0.01])
    zlist = np.array([nv.bbo-0.01, nv.bto+0.01])
    pts = np.zeros([xlist.size*ylist.size*zlist.size,4])
    ptsi = 0
    for z in zlist:
        for y in ylist:
            for x in xlist:
                pts[ptsi, :] = [x, y, z, ptsi]
                ptsi = ptsi+1
    return pts

# select the points that are the bottom corner of the blockcornerlist
    # pl is a point list
def blockcornerlist(pl:np.array) -> List:
    xmax = max(pl[:, 0])
    ymax = max(pl[:, 1])
    zmax = max(pl[:, 2])
    corners = pl[(pl[:,0]<xmax) & (pl[:,1]<ymax) & (pl[:,2]<zmax), :]
    return corners
 
# select the point (x,y,z) from an array of points pts
def ptfrompts(pts:np.array, x:float, y:float, z:float) -> List[float]:
    pt = pts[(pts[:,0]==x) & (pts[:,1]==y) & (pts[:,2]==z), :]
    return pt[0]

# select the points in this block based on the first corner, and put 
# the points in the order that openFOAM likes
    # pl is an array containing points
    # corner is one point which marks the first corner of the block
def blockpts(pl:float, corner:List[float]) -> List:
    xlist = np.unique(pl[:,0])
    ylist = np.unique(pl[:,1])
    zlist = np.unique(pl[:,2])
    cx = np.where(xlist==corner[0])
    cy = np.where(ylist==corner[1])
    cz = np.where(zlist==corner[2])
    cx = cx[0][0]
    cy = cy[0][0]
    cz = cz[0][0]
    pts = np.zeros([8, pl[1,:].size])
    for k in [0,1]:
        pts[0+4*k, :] = ptfrompts(pl, xlist[cx], ylist[cy], zlist[cz+k])
        pts[1+4*k, :] = ptfrompts(pl, xlist[cx+1], ylist[cy], zlist[cz+k])
        pts[2+4*k, :] = ptfrompts(pl, xlist[cx+1], ylist[cy+1], zlist[cz+k])
        pts[3+4*k, :] = ptfrompts(pl, xlist[cx], ylist[cy+1], zlist[cz+k])
    return pts

#--------------------------------------------------
#----------------------STLs------------------------

# given two two-element lists and one one-element list, determine the two triangles that constitute that face
# for example, two x values, two y values, and a z value
def axisface(xyz:List[float]):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    data = np.zeros(2, dtype=mesh.Mesh.dtype)
    if len(x)==1:
        data['vectors'][0] = np.array([[x[0], y[0], z[0]], [x[0], y[1], z[0]], [x[0], y[0], z[1]]])
        data['vectors'][1] = np.array([[x[0], y[1], z[1]], [x[0], y[1], z[0]], [x[0], y[0], z[1]]])
    if len(y)==1:
        data['vectors'][0] = np.array([[x[0], y[0], z[0]], [x[1], y[0], z[0]], [x[0], y[0], z[1]]])
        data['vectors'][1] = np.array([[x[1], y[0], z[1]], [x[1], y[0], z[0]], [x[0], y[0], z[1]]])
    if len(z)==1:
        data['vectors'][0] = np.array([[x[0], y[0], z[0]], [x[0], y[1], z[0]], [x[1], y[0], z[0]]])
        data['vectors'][1] = np.array([[x[1], y[1], z[0]], [x[0], y[1], z[0]], [x[1], y[0], z[0]]])
    return data

# given an array of 4 points pts, construct a face
# the line from point 1 to point 2 must cross through the center of the face
def ptface(pts):
    data = np.zeros(2, dtype=mesh.Mesh.dtype)
    data['vectors'][0] = np.array([pts[0], pts[1], pts[2]])
    data['vectors'][1] = np.array([pts[1], pts[2], pts[3]])
    return data

# given an array of 2d points pts2d, create an array of 3d points at position z
def setz(pts2d, z):
    out = np.zeros([len(pts2d), 3])
    out[:, 0:2] = pts2d
    out[:, 2] = z
    return out

# an arc face could be a donut on a plane (e.g. the bottom of the nozzle) if both lists have the same x,y,or z
# it could be a circle or cone (e.g. the nozzle inlet) if ina has length 1
# it could be a cylinder or frustrum if the two lists have different x,y,and z
    # ina is a list of points on the inner radius of the arc
    # outa is a list of points on the outer radius of the arc  
def arcface(ina, outa):
    data = np.zeros(2*len(ina), dtype = mesh.Mesh.dtype)
    if len(ina)==1:
        ina2 = np.ones([len(outa), 3])
        ina2[:,:] = ina
        ina = ina2
    for i in range(len(outa)-1):
        d = ptface([ina[i], outa[i], ina[i+1], outa[i+1]])
        data['vectors'][2*i] = d['vectors'][0]
        data['vectors'][2*i+1] = d['vectors'][1]
    return data

# get a mesh array that describes a hole in a plane
    # cpts is a list of circle points. cpts should be in order from 0 to 2 pi
    # x is a list of 2 x values for the plane
    # y is a list of 2 y values for the plane
    # z is a scalar that the plane lies on
def holeinplane(cpts, x, y, z):
    n = len(cpts)-1 # number of points on the circle
    nchunks = math.floor(n/4) # split the circle into four chunks, number of points in a chunk
    data = np.zeros(n+4, dtype = mesh.Mesh.dtype)
    x.sort()
    y.sort()
    corners = [[x[1], y[0], z], [x[1], y[1], z], [x[0], y[1], z], [x[0], y[0], z],  [x[1], y[0], z]]
    for i in range(4):
        data['vectors'][i*(nchunks+1)] = np.array([cpts[i*nchunks], corners[i], corners[i+1]])
        for j in range(nchunks):
            data['vectors'][i*(nchunks+1)+j+1] = np.array([cpts[i*nchunks+j],cpts[i*nchunks+j+1], corners[i+1]])
    premain = cpts[4*nchunks:]
    for j in range(len(premain)-1):
        data['vectors'][4*nchunks+4+j] = np.array([premain[j], premain[j+1], corners[-1]])
    return data

# combine all meshes into one list
    # meshlist is a list of lists of triangles
def combinemeshes(meshlist):
    n = functools.reduce(lambda a,b : a+len(b), meshlist, 0) # total number of mesh triangles
    data = np.zeros(n, dtype=mesh.Mesh.dtype)
    di = 0
    for m in meshlist:
        for mi in m:
            data[di] = mi
            di+=1
    return data


    
#--------------------------------------------------    
#---------------TRANSLATING TO openFOAM------------

# header for openfoam files
    # obj is the object name (string)
    # cl is the class name (string)
def header(cl:str, obj:str):
    s = ("/*--------------------------------*- C++ -*----------------------------------*\n"
        +"| =========                 |                                                 |\n"
        +"| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        +"|  \\\\    /   O peration     | Version:  v1912                                 |\n"
        +"|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
        +"|    \\\\/     M anipulation  |                                                 |\n"
        +"*---------------------------------------------------------------------------*/\n"
        +"FoamFile\n"
        +"{\n"
        +"    version     2.0;\n"
        +"    format      ascii;\n"
        +"    class       " + cl + ";\n"
        +"    object      " + obj + ";\n"
        +"}\n"
        +"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //;\n\n"
        )
    return s

# this removes residual data from previous runs
def compileAllClean():
    s = ("cd \"${0%/*}\" || exit                                # Run from this directory"
         + ". $WM_PROJECT_DIR/bin/tools/CleanFunctions; "
         + "cleanCase"
        )
    return s

# for folders that contain both a mesh and case folder, create a function that runs both allrun functions
def compileAllAllRun():
    s = ("./mesh/Allrun; ./case/Allrun")
    return s

def flistloop(s:str, functionlist:List[str], folder:str) -> str:
    for f in functionlist:
        s = s + 'echo \"running ' + f + ' in ' + folder + '\"; '
        s = s + f + ">>log_" + f + "; "
    return s

# this is the allrun bash script for the case folder
def compileAllRun(folder, solver):
    s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " \
         + "echo \""+folder+"\"; "
          + "cd \"${0%/*}\" || exit; "
         + "cp -r ../../mesh/constant/polyMesh constant; " 
         + "cp 0/alpha.ink.orig 0/alpha.ink; " )
    functionlist = ["setFields", solver, "foamToVTK"]
    s = flistloop(s, functionlist, folder)
    return s

# this script picks up in the middle of a solve and keeps solving
def compileContinue(folder, solver):
    s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " 
          + "cd \"${0%/*}\" || exit; "
         )
    functionlist = [solver, "foamToVTK"]
    s = flistloop(s, functionlist, folder)
    return s

# this script runs the meshing functions in the mesh folder
def compileAllRunMesh(folder):
    s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " 
        + "cd \"${0%/*}\" || exit; ")
    functionlist = ["surfaceFeatureExtract", "blockMesh", "snappyHexMesh -overwrite", "foamToVTK"]
    s = flistloop(s, functionlist, folder)
    return s

CLOSELINE = "// ************************************************************************* //";


    
#-------------------------------------------------   
#-----------------SYSTEM--------------------------


#-------------------------------------------------
#---------------blockMeshDict---------------------

# convert a vector to the format it needs to be in for OpenFOAM to read it
    # v is a list
def vec2cpp(v):
    # convert 1D array v to openfoam string
    s = "("
    for vi in v:
        s = s + str(vi) + " "
    s = s + ")"
    return s

# convert block to openfoam string
    # block is a Block object
def block2txt(block):
    s = "hex " + vec2cpp(block.vertices[:, 3].astype(int)) + " " + vec2cpp(block.meshi)
    s = s + " simpleGrading " + vec2cpp(block.grading)
    return s

# convert list of blocks to openfoam string
    # block is a list of Block objects
def blocks2txt(blocks):
    pl = DictList("blocks", 1, [])
    for b in blocks:
        pl.proplist.append(block2txt(b))
    return pl.prnt(0)
    
EDGES = "edges\n(\n);"

# faceselector gets a list of vertices for a list
    # block is a block object
    # st is a string indicating which face to use
def faceselector(block, st):
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

# boundarylist compiles a list of all the boundaries in the system for blockMesh
# because we are using snappyHexMesh, we only need one boundary in blockMesh
    # blocks is a list of Block objects
def boundarylist(blocks):
    all = BoundaryInput("allBoundary", "patch")
    for st in ["x-", "x+", "y-", "y+", "z-", "z+"]:
        all.flist.append(faceselector(blocks[0], st)) # front and back faces
    return [all]

# select the boundaries of a wall
    # geo is a nozzleVars object
    # stri is a string indicating which face to take
def walsel(geo, st):
    xli = [geo.ble, geo.bri]
    yli = [geo.bfr, geo.bba]
    zli = [geo.bbo, geo.bto]
    if st == "x-":
        xli = [geo.ble]
    elif st == "x+":
        xli = [geo.bri]
    elif st == "y-": 
        yli = [geo.bfr]
    elif st == "y+":
        yli = [geo.bba]
    elif st == "z-": 
        zli = [geo.bbo]
    elif st == "z+": 
        zli = [geo.bto]
    return [xli, yli, zli]

# get a list of boundaries for snappyHexMesh, setFields
    # geo is a NozVars object geo.bv and geo.iv are bath velocities and ink velocities, scaled
    # mv is a MeshVars object
    # exportmesh is true if we want to create a mesh folder
def realboundaries(geo, mv, exportmesh):
    bf = BoundaryInput("bathFlow", "")
    bf.alphalist = DictList(bf.label, 0, [["type", "fixedValue"], ["value", "uniform 0"]])
    bf.Ulist = DictList(bf.label, 0, [["type", "fixedValue"], ["value", "uniform (" + str(geo.bv) + " 0 0)"]])
    bf.plist = DictList(bf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    if exportmesh:
        bf.meshi = combinemeshes(list(map(lambda s: axisface(walsel(geo, s)), ["x-", "y-", "y+", "z-"])))
    #print(bf.meshi)
    
    inkf = BoundaryInput("inkFlow", "")
    inkf.alphalist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform 1"]])
    inkf.Ulist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(geo.iv) + ")"]])
    inkf.plist = DictList(inkf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    if exportmesh:
        inkf.meshi = arcface(setz(np.zeros([len(geo.inptst), 2])+[geo.ncx, geo.ncy], geo.bto), setz(geo.inptst, geo.bto))
    #print(inkf.meshi)
    
    at = BoundaryInput("atmosphere", "")
    at.alphalist = DictList(at.label, 0, [["type", "inletOutlet"], ["value", "uniform 0"], ["inletValue", "uniform 0"]])
    at.Ulist = DictList(at.label, 0, [["type", "pressureInletOutletVelocity"], ["value", "uniform (0 0 0)"]])
    at.plist = DictList(at.label, 0, [["type", "totalPressure"], ["p0", "uniform 0"]])
    if exportmesh:
        at.meshi = combinemeshes([holeinplane(setz(geo.outptst, geo.bto), [geo.ble, geo.bri], [geo.bfr, geo.bba], geo.bto), \
                           axisface(walsel(geo, "x+"))])
    #print(at.meshi)
    
    fw = BoundaryInput("fixedWalls", "")
    fw.alphalist = DictList(fw.label, 0, [["type", "zeroGradient"]])
    fw.Ulist = DictList(fw.label, 0, [["type", "noSlip"]])  
    fw.plist = DictList(fw.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    if exportmesh:
        fw.meshi = combinemeshes([arcface(setz(geo.inptsb, geo.nbo), setz(geo.outptsb, geo.nbo)),\
                            arcface(setz(geo.inptsb, geo.nbo), setz(geo.inptst, geo.bto)),\
                            arcface(setz(geo.outptsb, geo.nbo), setz(geo.outptst, geo.bto))])
    fw.reflev = mv.fwreflev
    #print(fw.meshi)
    
    return [bf, inkf, at, fw]
    
# get the whole cpp section for list of boundaries 
    # bl is a list of BoundaryInput objects
def boundarycpp(bl): 
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

# get the blockMeshDict text
    # pl is a point list
    # blocks is a list of Block objects
    # bl is a list of BoundaryInput objects
def compileBlockMeshDict(pl, blocks, bl):
    s = header("dictionary", "blockMeshDict")
    s = s + "scale " + str(SCALE) + ";\n\n"
    s = s + DictList("vertices", 1, pl[:, 0:3]).prnt(0)
    s = s + blocks2txt(blocks)
    s = s + "edges\n(\n);\n\n"
    s = s + boundarycpp(bl)
    s = s + "mergePatchPairs\n(\n);\n\n" + CLOSELINE
    return s


#-------------------------------------------------
#---------------controlDict-----------------------

# compileControlDict gets the text for controlDict
    # cdv is a CDVars object
def compileControlDict(cdv):
    s = header("dictionary", "controlDict")
    v = cdv.vlist()
    #v.append(["#sinclude", "sampling"])
    v.append(['#includeIfPresent', '\"sampling\"'])
    s = s + DictList("", 0,v).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#---------------fvSchemes-----------------------  

def compileFvSchemes(fvv):
    s = header("dictionary", "fvSchemes")
    if isinstance(fvv, FVVars):
        s = s + DictList("", 0, fvv.fvSchemeList()).prnt(-1)
    else:
        for st in ["gradSchemes", "divSchemes", "laplacianSchemes"]:
            s = s + st + "\n{\n}\n"
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#---------------fvSolution-----------------------
def compileFvSolution(fvv):
    s = header("dictionary", "fvSolution")
    if isinstance(fvv, FVVars):
        s = s + fvv.solverlist().prnt(0)
        s = s + fvv.pimple().prnt(0)
        s = s + DictList("relaxationFactors", 0, [DictList("equations", 0, [["\".*\" 1"]])]).prnt(0)
    s = s + CLOSELINE
    return s

def compileFV(fvv, out):
    out.fvSchemes = compileFvSchemes(fvv)
    out.fvSolution = compileFvSolution(fvv)
    return out

#-------------------------------------------------
#---------------all solver files-----------------------

# solverobjects gets the control dictionary variables and fvsolution and fvschemes dictionary variables
    # starttime is the start time for the simulation in s
    # endtime is the end time for the simulation in s
    # dt is the initial solve time step in s
    # writedt is the write time step in s
    # solver is "interFoam" or "interIsoFoam"
def solverobjects(starttime, endtime, dt, writedt, solver):
    cdv = CDVars(starttime, endtime, dt, writedt)
    cdv.application = solver
    fvv = FVVars(solver)
    return [cdv, fvv]

# compileSolverfiles gets the text for the controlDict, FVsolution, FVschemes, and bash scripts
    # cdv is a CDVars object
    # fvv is a FVVars object
    # out is a FileGroup object
def compileSolverFiles(cdv, fvv, out):
    out.controlDict = compileControlDict(cdv)
    out = compileFV(fvv, out)
    out.allallrun = compileAllAllRun()
    out.allclean = compileAllClean()
    out.allrunmesh = compileAllRunMesh(out.folder)
    out.allrun = compileAllRun(out.folder, cdv.application)
    out.cont = compileContinue(out.folder, cdv.application)
    cdv.application = "icoFoam"
    out.controlDictMesh = compileControlDict(cdv)
    return out

#-------------------------------------------------
#------------------setFieldsDict----------------------

# compile setFieldsDict
    # geo is a NozVars object
def compileSetFieldsDict(geo):
    s = header("dictionary", "setFieldsDict")
    s = s + DictList("defaultFieldValues", 1, ["volScalarFieldValue alpha.ink 0"]).prnt(0)
    r = DictList("regions", 1, [])                 
        # r is the region where the ink originally is
    c2c = DictList("cylinderToCell", 0, [])         
        # c2c is a cylinderToCell dictionary list
    for i,z in enumerate([geo.nbo, geo.bto]):
        c2c.proplist.append(["p" + str(i+1) + " ( " + str(geo.ncx*SCALE) + " " + str(geo.ncy*SCALE) + " " + str(z*SCALE) + " )"])
                # top and bottom points of central cylinder axis
    c2c.proplist.append(["radius " + str(str((geo.niw+geo.nt)/2*SCALE))])        
                # radius of cylinder
    c2c.proplist.append(DictList("fieldValues", 1, ["volScalarFieldValue alpha.ink 1"])) 
                # value of alpha inside the cylinder is 1 because it is full of ink
    r.proplist.append(c2c)
    s = s + r.prnt(0)
    s = s + CLOSELINE
    return s

#---------------------------------------------
#-------------------meshQualityDict-----------

# compile meshQualityDict
def compileMeshQualityDict():
    s = header("dictionary", "meshQualityDict")
    s = s + "#includeEtc \"caseDicts/meshQualityDict\"\n\n"
    s = s + CLOSELINE
    return s

#---------------------------------------------
#----------------snappyHexMeshDict-------------

# compileSnappyHexMeshDict
    # bl is a list of BoundaryInput objects
    # mv is a MeshVars object
def compileSnappyHexMeshDict(bl, mv):
    bnames = [o.label for o in bl]
    s = header("dictionary", "snappyHexMeshDict")
    s = s + DictList("", 0, [\
                             ["castellatedMesh", mv.castellatedMesh],\
                             ["snap", mv.snap],\
                             ["addLayers", mv.addLayers]]).prnt(-1)
    s = s + DictList("geometry", 0,\
                     map(lambda x: DictList(x+".stl", 0, [["type",  "triSurfaceMesh"], ["name ", x]]), bnames)).prnt(0)
    
    # castellated mesh controls
    cmc = mv.cmc()
    featuredictlist = DictList("features", 1, \
            list(map(lambda b: DictList("", 0, [["file", "\""+b.label+".eMesh\""], ["level", b.reflev]]), bl)))
    refinementsurfaces = DictList("refinementSurfaces", 0,\
            list(map(lambda b: DictList(b.label, 0, [["level", "("+str(b.reflev)+" "+str(b.reflev)+")"]]), bl)))
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
    s = s + CLOSELINE
    return s

#---------------------------------------------------
#-----------------surfaceFeatureExtractDict--------

# compile surfaceFeatureExtractDict
    # bl is a list of BoundaryInput objects
def compileSurfaceFeatureExtractDict(bl):
    bnames = [o.label for o in bl]
    s = header("dictionary", "surfaceFeatureExtractDict")
    s = s + DictList("", 0, \
                     list(map(lambda x: DictList(x+".stl", 0, \
                                                 [["extractionMethod", "extractFromSurface"],\
                                                  DictList("extractFromSurfaceCoeffs", 0, [["includedAngle", 180]]),\
                                                  ["writeObj", "yes"]]), bnames))).prnt(-1)
    s = s + CLOSELINE
    return s
    
    

#-------------------------------------------------   
#-----------------CONSTANT--------------------------
#-------------------------------------------------  

#-------------------------------------------------
#------------------------g------------------------

# compile g
def compileg():
    s = header("uniformDimensionedVectorField", "g")
    s = s + DictList("", 0, [["dimensions", "[0 1 -2 0 0 0 0]"], ["value", "(0 0 -9.81)"]]).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#------------------transportProperties------------

# transportGroupNewt gets a DictList object that describes the transport properties
    # title is the name of the phase, e.g. 'ink'
    # nu is the viscosity in m^2/s
    # rho is the density in kg/m^3
def transportGroupNewt(title, nu, rho):
    return DictList(title, 0, [["transportModel", "Newtonian"], ["nu", str(nu)], ["rho", str(rho)]])

# transportGroupHB gets a DictList that describes Herschel-Bulkley transport properties
    # title is the name of the phase
    # the HB model follows nu = min(nu0, tau0/gammadot + k*gammadot^(n-1))
    # inputs to the model are nu0 in m^2/s, tau0 in m^2/s^2, k in m^2/s, n unitless
    # rho is the density in kg/m^3
def transportGroupHB(title, nu0, tau0, k, n, rho):
    return DictList(title, 0, \
                  [["transportModel", "HerschelBulkley"], \
                  DictList("HerschelBulkleyCoeffs", 0, \
                           [["nu0", "[ 0 2 -1 0 0 0 0 ] " + str(nu0)],\
                            ["tau0", "[ 0 2 -2 0 0 0 0 ] " + str(tau0)], \
                            ["k", "[ 0 2 -1 0 0 0 0 ] " + str(k)],\
                            ["n", "[ 0 0 0 0 0 0 0 ] " + str(n)]]\
                           ),\
                  ["rho", str(rho)]])

# compile transportProperties
    # ink and sup are DictLists created by transportGroupNewt
    # sigma is the surface tension in J/m^2
def compileTransportProperties(ink, sup, sigma):
    s = header("dictionary", "transportProperties")
    s = s + DictList("", 0, [["phases (ink sup)"],ink, sup, ["sigma", str(sigma)]]).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#------------turbulenceProperties-----------------
def compileTurbulenceProperties():
    s = header("dictionary", "turbulenceProperties")
    s = s + DictList("", 0, [["simulationType", "laminar"]]).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#------------dynamicMeshDict-----------------
    # mv is a meshVars object
def compileDynamicMeshDict(mv):
    s = header("dictionary", "dynamicMeshDict")
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
    s = s + CLOSELINE
    return s
        

#-------------------------------------------------   
#----------------------0--------------------------
#-------------------------------------------------  

# format0 puts files in the 0 folder into the correct format
    # geo is a NozVars object
    # bl is the list of boundaryInput objects
    # classn is the class that goes into the header, e.g. "volScalarField" or "pointScalarField"
    # obj is the object that goes into the header, e.g. "U"
    # dims is the dimensions of the field, e.g. "[1 -1 -2 0 0 0 0]", which means kg/(m*s^2)
    # intfield is the internal field, e.g. "0" or "(0 0 0)"
    # f is the function to run on each boundaryInput object to get a list of properties of interest
def format0(geo, bl, classn, obj, dims, intfield, f):
    # geo is a NozVars object
    # bl is a boundary list (list of BoundaryInput objects)
    s = header(classn, obj) # header
    simplelist = DictList("", 0, [["dimensions", dims], ["internalField", "uniform " + intfield]]) # list of simple variables to define
    s = s + simplelist.prnt(-1) # format the simple list 
    bflist = DictList("boundaryField", 0, []) # list of boundary dictionary entries
    for b in bl:
        bflist.proplist.append(f(b)) # add the appropriate list for each boundary, e.g. alpha properties for inkFlow
    s = s + bflist.prnt(0) # format the dictionary list
    s = s + CLOSELINE
    return s

def compilealphaorig(geo, bl):
    return format0(geo, bl, "volScalarField", "alpha.ink", "[0 0 0 0 0 0 0]", "0", allalpha)
    
def compileU(geo, bl):
    return format0(geo, bl, "volVectorField", "U", "[0 1 -1 0 0 0 0]", "(0 0 0)", allU)
    
def compilep(geo, bl):
    return format0(geo, bl, "volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]", "0", allp)

def compileCellLevel(geo, bl):
    return format0(geo, bl, "volScalarField", "cellLevel", "[0 0 0 0 0 0 0]", "0", lambda x:[["type", "zeroGradient"]])

def compilePointLevel(geo, bl):
    return format0(geo, bl, "pointScalarField", "pointLevel", "[0 0 0 0 0 0 0]", "0", lambda x:[["type", "calculated"]])

        
#-------------------------------------------------        
#------------compile all the files--------------
#-------------------------------------------------

# createNozzleBlockFile gets the text for most of the files in the folder
    # geo is a NozVars object
    # mv is a meshVars object
    # exportmesh is true if we want to create a mesh folder
def createNozzleBlockFile(geo, mv, exportmesh, folder):

    fg = FileGroup(folder)
    fg.geofile = geometryFile(geo)
    
    #system
    fg.setFieldsDict = compileSetFieldsDict(geo)
    if exportmesh:
        pl = createnozzle(geo) # point list
        bcs = blockcornerlist(pl) # list of block corners
        blocks = []
        for bc in bcs:
            blocki = Block()
            blocki.vertices = blockpts(pl, bc)
            msz = mv.meshsize # mesh size
            blocki.meshi = [round(geo.bw/msz), round(geo.bh/msz), round(geo.bd/msz)]
            blocks.append(blocki) 
        bl = boundarylist(blocks) # this is a dummy boundary list for the original blockMeshDict
        fg.blockMeshDict = compileBlockMeshDict(pl, blocks, bl)
    
    
    #constant
    fg.dynamicMeshDict = compileDynamicMeshDict(mv)
    
    #0
    br = realboundaries(geo, mv, exportmesh) # this is a real boundary list for the imported stl boundaries
    fg.meshes = br # keep the real boundaries for exporting
    fg.alphainkorig = compilealphaorig(geo, br)
    fg.U = compileU(geo, br)
    fg.prgh = compilep(geo, br)
    fg.cellLevel = compileCellLevel(geo, br)
    fg.pointLevel = compilePointLevel(geo, br)
    if exportmesh:
        fg.snappyHexMeshDict = compileSnappyHexMeshDict(br, mv)
        fg.surfaceFeatureExtractDict = compileSurfaceFeatureExtractDict(br)
    #plot
    if exportmesh:
        fg.geo = geo
        fg.pl = pl
        fg.blocks = blocks
        fg.br = br
    #fg.plot = plotGeo(geo, pl, blocks, br)
    return fg

# geometryFile gets a csv string of all of the geometry variables we care about
    # geo is a NozVars object
def geometryFile(geo):
    l = [['nozzle inner width (mm)', geo.niw],\
         ['nozzle thickness (mm)', geo.nt], \
         ['bath width (mm)', geo.bw], \
         ['bath depth (mm)', geo.bd], \
         ['nozzle length (mm)', geo.nl],\
         ['bath left coord (mm)', geo.ble], \
         ['bath right coord (mm)', geo.bri],\
         ['bath front coord (mm)', geo.bfr], \
         ['bath back coord (mm)', geo.bba],\
         ['bath bottom coord (mm)', geo.bbo], \
         ['bath top coord (mm)', geo.bto], \
         ['nozzle bottom coord (mm)', geo.nbo],\
         ['nozzle center x coord (mm)', geo.ncx],\
         ['nozzle center y coord (mm)', geo.ncy], \
         ['bath velocity (m/s)', geo.bv], \
         ['ink velocity (m/s)', geo.iv]]
    s = ""
    for li in l:
        s = s + li[0] + ', ' + str(li[1]) + '\n'
    return s
         


#-------------------------------------------------        
#---------------------PLOTS-----------------------
#-------------------------------------------------

# plot the requested geometry
    # geo is a nozzleVars object
    # pl is a point list used for blockMesh
    # blocks is a list of blocks used for blockMesh
    # br is real boundaries used for generating stls and setting boundary conditions
def plotGeo(geo, pl, blocks, br):
    
    fig = plt.figure(figsize=[14,14])
    ax = fig.add_subplot(111, projection='3d')
    maxdim = max([geo.bw, geo.bh, geo.bd])
    ax.set_xlim3d(geo.ble, geo.ble+maxdim)
    ax.set_ylim3d(geo.bfr, geo.bfr+maxdim)
    ax.set_zlim3d(geo.bbo, geo.bbo+maxdim)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.scatter(pl[:,0], pl[:,1], pl[:,2], s=20, c='k')
    fs = 20;
    for p in pl:
        ax.text(p[0], p[1], p[2], "{:.0f}".format(p[3]), color='k', fontsize=fs)
    for i in range(len(blocks)):
        vs = blocks[i].vertices
        ax.text(st.mean(vs[:,0]), st.mean(vs[:,1]), st.mean(vs[:,2]), i, color='r', fontsize=fs)
    clist = ['tomato', 'navajowhite', 'mediumaquamarine', 'deepskyblue', 'white']
    plist = []
    for bi in br:
        # each boundary bi in list bl
        col = clist.pop(0) # pop a color from the front of the list
        plist.append(mpatches.Patch(color=col, label=bi.label))
        for m in bi.meshi['vectors']:
            #print(m)
            p3c = Poly3DCollection([list(zip(m[:,0],m[:,1],m[:,2]))], alpha=0.35)
            p3c.set_facecolor(col)
            p3c.set_edgecolor('gray')
            ax.add_collection3d(p3c)      
    ax.legend(handles=plist, fontsize=fs)
    return fig

#------------------------------------------------
#--------------------EXPORTING-------------------

# exports stls
    # folder is a full path name
    # bl is a list of boundaryInput objects 
def savestls(folder, bl):
    for b in bl:
        me = b.meshi
        me2 = mesh.Mesh(np.zeros(len(me), dtype=mesh.Mesh.dtype))
        for i,m in enumerate(me['vectors']):
            me2.vectors[i] = m*SCALE
        title = os.path.join(folder, b.label+".stl")
        me2.save(title)

# exports text files
    # folder is a full path name
    # file is a basename
    # text is the text to export
def exportfile(folder, file, text):
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    print("Exported file %s" % fn)
    
# makes a directory if it doesn't already exist
    # path is the directory to create
def mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)

# exports all of the files that we generated
    # f is a folder string
    # fg is a FileGroup object
def exportallfiles(fg, exportmesh):
    f = fg.folder
    casef = os.path.join(f, "case")
    f0 = os.path.join(casef, "0")
    fconst = os.path.join(casef,"constant")
    fsyst = os.path.join(casef,"system")
    folderlist = [f, casef, f0, fconst, fsyst]
    if exportmesh:
        fgeom = os.path.join(f,"geometry")
        fmesh = os.path.join(f,"mesh")
        fmeshconst = os.path.join(fmesh,"constant")
        fmeshconsttri = os.path.join(fmeshconst,"triSurface")
        fmeshsys = os.path.join(fmesh,"system")
        fmesh0 = os.path.join(fmesh,"0")
        for f in [fgeom, fmesh, fmeshconst, fmeshconsttri, fmeshsys, fmesh0]:
            folderlist.append(f)
    
    list(map(mkdirif, folderlist)) # create folders
    
    exportfile(f, "geometry.csv", fg.geofile) 
    exportfile(casef, "Allclean", fg.allclean) 
    exportfile(casef, "Allrun", fg.allrun) 
    exportfile(casef, "Continue", fg.cont) 
    exportfile(f0, "alpha.ink.orig", fg.alphainkorig) 
    exportfile(f0, "alpha.ink", fg.alphainkorig) 
    exportfile(f0, "p_rgh", fg.prgh) 
    exportfile(f0, "U", fg.U) 
    exportfile(fconst, "g", fg.g) 
    exportfile(fconst, "transportProperties", fg.transportProperties) 
    exportfile(fconst, "turbulenceProperties", fg.turbulenceProperties) 
    exportfile(fconst, "dynamicMeshDict", fg.dynamicMeshDict) 
    exportfile(fsyst, "controlDict", fg.controlDict)     
    exportfile(fsyst, "fvSchemes", fg.fvSchemes)     
    exportfile(fsyst, "fvSolution", fg.fvSolution)    
    exportfile(fsyst, "setFieldsDict", fg.setFieldsDict) 
        
    if exportmesh:
        exportfile(f, "Allrun", fg.allallrun) 
        exportfile(fmesh, "Allclean", fg.allclean) 
        exportfile(fmesh, "Allrun", fg.allrunmesh)
        exportfile(fmesh0, "pointLevel", fg.pointLevel)
        exportfile(fmesh0, "cellLevel", fg.cellLevel)
        exportfile(fmeshsys, "blockMeshDict", fg.blockMeshDict) 
        exportfile(fmeshsys, "controlDict", fg.controlDictMesh) 
        exportfile(fmeshsys, "fvSchemes", compileFvSchemes("")) 
        exportfile(fmeshsys, "fvSolution", compileFvSolution(""))
        exportfile(fmeshsys, "meshQualityDict", fg.meshQualityDict)
        exportfile(fmeshsys, "snappyHexMeshDict", fg.snappyHexMeshDict)
        exportfile(fmeshsys, "surfaceFeatureExtractDict", fg.surfaceFeatureExtractDict)
        savestls(fgeom, fg.meshes)
        savestls(fmeshconsttri, fg.meshes)
    
    fs.populate(f)
    



def allbuttransport(ii, exportmesh, topfolder):
    folder = os.path.join(topfolder, "nb"+str(ii))

    bathwidth = 16 # times nozzle inner diameter
    bathheight = 7
    bathdepth = 7
    frontwidth = 4
    printspeed = 10 # mm/s
    npts = 50 # number of points in the nozzle circle

    geo = NozVars(bathwidth, bathheight, bathdepth, frontwidth, printspeed, npts)
    mv = MeshVars()
    mv.maxRefinement = 4
    out = createNozzleBlockFile(geo, mv, exportmesh, folder)

    starttime = 0
    endtime = 2.5
    dt = 0.001
    writedt = 0.1
    solver = "interFoam" # options: interFoam, interIsoFoam
    [cdv, fvv] = solverobjects(starttime, endtime, dt, writedt, solver)
    cdv.startFrom = "latestTime"
    cdv.adjustTimeStep = "yes"
    #cdv.runTimeModifiable = "yes"
    out = compileSolverFiles(cdv, fvv, out)
    return out

####### for newtonian ink and support
def newtnewtexport(ii, ivisc, svisc, sigma, topfolder, exportmesh):
    out = allbuttransport(ii, exportmesh, topfolder)
    ink = transportGroupNewt("ink", ivisc, 1000)
    sup = transportGroupNewt("sup", svisc, 1000)
    out.transportProperties = compileTransportProperties(ink, sup, sigma)
    exportallfiles(out, exportmesh)

####### for Herschel Bulkley support and newtonian ink
def HBnewtexport(ii, tau0, k, n, ivisc, nu0, sigma, topfolder, exportmesh):
    out = allbuttransport(ii, exportmesh, topfolder)
    ink = transportGroupNewt("ink", ivisc, 1000)
    sup = transportGroupHB("sup", nu0, tau0, k, n, 1000)
    out.transportProperties = compileTransportProperties(ink, sup, sigma)
    exportallfiles(out, exportmesh)
    
####### for newtonian support and HerschelBulkley ink
def newtHBexport(ii, tau0, k, n, svisc, nu0, sigma, topfolder, exportmesh):
    out = allbuttransport(ii, exportmesh, topfolder)
    sup = transportGroupNewt("sup", svisc, 1000)
    ink = transportGroupHB("ink", nu0, tau0, k, n, 1000)
    out.transportProperties = compileTransportProperties(ink, sup, sigma)
    exportallfiles(out, exportmesh)
    
####### for newtonian support and HerschelBulkley ink
def HBHBexport(ii, tau0, k, n, nu0sup, nu0ink, sigma, topfolder, exportmesh):
    out = allbuttransport(ii, exportmesh, topfolder)
    sup = transportGroupHB("sup", nu0sup, tau0, k, n, 1000)
    ink = transportGroupHB("ink", nu0ink, tau0, k, n, 1000)
    out.transportProperties = compileTransportProperties(ink, sup, sigma)
    exportallfiles(out, exportmesh)
    
####### for newtonian support and HerschelBulkley ink
def HBHByielded(ii, tau0sup, tau0ink, ksup, kink, nsup, nink, nu0sup, nu0ink, sigma, topfolder, exportmesh):
    out = allbuttransport(ii, exportmesh, topfolder)
    sup = transportGroupHB("sup", nu0sup, tau0sup, ksup, nsup, 1000)
    ink = transportGroupHB("ink", nu0ink, tau0ink, kink, nink, 1000)
    out.transportProperties = compileTransportProperties(ink, sup, sigma)
    exportallfiles(out, exportmesh)
 
        
        
    