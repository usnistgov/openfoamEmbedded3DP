#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
import os
from stl import mesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
import matplotlib.patches as mpatches
import functools
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import numpy as np
import logging
import sys
import pandas as pd
from descartes import PolygonPatch

# local packages
import folderscraper as fs
from config import cfg
import interfacemetrics as intm # RG

# logging
logging.basicConfig(level=logging.INFO)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  


#-------------------------------------------------
#---------------GLOBAL VARIABLES------------------
#-------------------------------------------------

SCALE = 0.001 # scale all units by this amount * m
 
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## GENERAL
    
def tabs(n:int) -> str:
    '''return n tabs in a row'''
    s = ""
    for i in range(n):
        s = s + "\t"
    return s


class DictList:
    '''DictList is an object that ensures uniform formatting of variables into files'''
        
    def __init__(self, title:str, form:int, proplist:List[str]):
        '''Inputs: title, form, proplist
        title is the title of the variable group
        form is 0 for bracketed entries e.g. group definitions, 1 for parentheses e.g. point lists
        proplist is a list of properties that we're storing within this group'''
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
    
    
    def prnt(self, level:int) -> str:
        '''Format the title and variables into the format that OpenFOAM likes. 
        Input: the number of tabs before the title of the group. If you want no title, use level -1
        Output: formatted string'''
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
            self.slurmFolder = cfg.path.slurmFolder

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

        list(map(mkdirif, folderList)) # create folders

        
        if not self.onlyMesh:
            exportFile(f, 'labels.csv', self.labels)
            exportFile(casef, "Allclean", self.allclean, linux=True) 
#             exportFile(casef, "Allrun", self.allrun)
            exportFile(casef, 'Allrun.sh', self.allrun, linux=True)
            exportFile(casef, 'run.slurm', self.slurm, linux=True)
#             exportFile(casef, "Continue", self.cont) 
            exportFile(f0, "alpha.ink.orig", self.alphainkorig) 
            exportFile(f0, "alpha.ink", self.alphainkorig) 
            exportFile(f0, "p_rgh", self.prgh) 
            exportFile(f0, "U", self.U) 
            exportFile(fconst, "g", self.g) 
            exportFile(fconst, "transportProperties", self.transportProperties) 
            exportFile(fconst, "turbulenceProperties", self.turbulenceProperties) 
            exportFile(fconst, "dynamicMeshDict", self.dynamicMeshDict) 
            exportFile(fsyst, "controlDict", self.controlDict)     
            exportFile(fsyst, "fvSchemes", self.fvSchemes)     
            exportFile(fsyst, "fvSolution", self.fvSolution)    
            exportFile(fsyst, "setFieldsDict", self.setFieldsDict) 

        if self.exportMesh:
            if not self.onlyMesh:
                exportFile(f, "Allrun", self.allallrun, linux=True)
                exportFile(f, "geometry.csv", self.geofile)
            exportFile(fmesh, "Allclean", self.allclean, linux=True) 
            exportFile(fmesh, "Allrun", self.allrunmesh, linux=True)
            exportFile(fmesh0, "pointLevel", self.pointLevel)
            exportFile(fmesh0, "cellLevel", self.cellLevel)
            exportFile(fmeshsys, "blockMeshDict", self.blockMeshDict) 
            exportFile(fmeshsys, "controlDict", self.controlDictMesh) 
            exportFile(fmeshsys, "fvSchemes", self.fvSchemesMesh) 
            exportFile(fmeshsys, "fvSolution", self.fvSolutionMesh)
            exportFile(fmeshsys, "meshQualityDict", self.meshQualityDict)
            exportFile(fmeshsys, "snappyHexMeshDict", self.snappyHexMeshDict)
            exportFile(fmeshsys, "surfaceFeatureExtractDict", self.surfaceFeatureExtractDict)
            exportFile(fmeshsys, "surfaceFeaturesDict", self.surfaceFeaturesDict)
            saveStls(fgeom, self.meshes)
            saveStls(fmeshconsttri, self.meshes)

            
        try:
            fs.populate(f)
        except:
            pass
        
    
    def makePlot(self) -> None:
        '''Makes a plot of the simulation geometry'''
        self.plot = plotGeo(self.geo, self.pl, self.blocks, self.br)
        
        
        
        
    
#--------------------------------------------------    
#---------------TRANSLATING TO openFOAM------------


def header(cl:str, obj:str) -> str:
    '''header for openfoam files
    cl is the class name (string)
    obj is the object name (string)
    '''
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


def compileAllClean() -> str:
    '''this removes residual data from previous runs'''
    s = ("cd \"${0%/*}\" || exit                                # Run from this directory"
         + ". $WM_PROJECT_DIR/bin/tools/CleanFunctions; "
         + "cleanCase"
        )
    return s

def compileAllAllRun() -> str:
    '''for folders that contain both a mesh and case folder, create a function that runs both Allrun functions'''
    s = ("chmod +x ./mesh/Allrun; chmod +x ./case/run.slurm; chmod +x ./case/Allrun.sh;\n") # RG
    s = s + ("./mesh/Allrun; sbatch ./case/run.slurm") # RG
    return s

def fListLoop(s:str, functionlist:List[str], folder:str, ifstarted:bool=False) -> str:
    '''Write a bash script to go through multiple functions. ifstarted True if you do not want to run this function if the sim has started'''
    for f in functionlist:
        s1 = 'echo \"running ' + f + ' in ' + folder + '\";\n' + f + ">>log_" + f
        if ifstarted:
            s = s + '[ ! d \"0.1\"] && ('+s1+'); '
        else:
            s = s + s1 + ';\n'
    return s


# def compileAllRun(folder:str, solver:str) -> str:
#     '''this is the allrun bash script for the case folder'''
#     s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " \
#          + "echo \""+folder+"\"; "
#           + "cd \"${0%/*}\" || exit; ")
#     s = (s + '[ ! d \"0.1\"] && ('\
#          + "cp -r ../../mesh/constant/polyMesh constant; " 
#          + "cp 0/alpha.ink.orig 0/alpha.ink); " )
#     s = fListLoop(s, ["setFields"], folder, ifstarted=True)
#     s = s = fListLoop(s, [solver, "foamToVTK"], folder, ifstarted=False)
#     return s

def compileAllRun(folder:str, solver:str) -> str: # RG
    '''this is the allrun bash script for the case folder'''
    f = os.path.basename(folder)
    s = '#!/bin/bash\n\n'
    s = s + '. $WM_PROJECT_DIR/bin/tools/RunFunctions;\n'
    s = s + 'cd "${0%/*}" || exit;\n' # RG
    s = s + 'echo '+f+'\n'
    s = s + 'if [ ! -d "0.1" ]; then\n'
    s = s + '\tcp -r ../mesh/constant/polyMesh constant;\n'
    s = s + '\tcp 0/alpha.ink.orig 0/alpha.ink;\n'
    s = s + '\techo \"running setFields in '+f+'\";\n'
    s = s + '\tsetFields>>log_setFields;\n'
    s = s + 'fi \n'
    s = fListLoop(s, [solver, "foamToVTK"], f, ifstarted=False)
    return s

def compileSlurm(folder:str, parentdir:str) -> str:
    '''this is the slurm script for the case folder'''
    workdir = (os.path.join(parentdir, os.path.basename(folder), 'case')).replace("\\","/")
    s = f'#!/bin/bash\n#SBATCH -p local\n#SBATCH --time=14-00:00:00\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n#SBATCH --job-name={os.path.basename(folder)}\n#SBATCH --workdir={workdir}\n\n'
    s = s + 'export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n\nsrun bash '+workdir+'/Allrun.sh' # RG
    return s


# def compileContinue(folder:str, solver:str) -> str:
#     '''this script picks up in the middle of a solve and keeps solving'''
#     s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " 
#           + "cd \"${0%/*}\" || exit; "
#          )
#     functionlist = [solver, "foamToVTK"]
#     s = fListLoop(s, functionlist, folder)
#     return s


def compileAllRunMesh(folder:str) -> str:
    '''this script runs the meshing functions in the mesh folder'''
    s = (". $WM_PROJECT_DIR/bin/tools/RunFunctions; " 
        + "cd \"${0%/*}\" || exit;\n")
    functionlist = ["surfaceFeatures", "blockMesh", "snappyHexMesh -overwrite", "foamToVTK"] # RG
    s = fListLoop(s, functionlist, folder)
    return s

CLOSELINE = "// ************************************************************************* //";
    
    
    #----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## GEOMETRY
    
   #--------------------------------------------------------------------------------------------------------        
######## blocks and boundaries
        

class Block:
    '''stores information about a block during blockMesh '''
    
    def __init__(self):
        self.vertices = []
        self.meshi = [1,1,1]
        self.grading = [1,1,1]
        

class BoundaryInput:
    '''# stores information about a boundary '''
    
    def __init__(self, labelin:str, typin:str):
        '''Input: labelin is the name of the boundary e.g. outerWall. typin is the boundary type, e.g. patch'''
        self.label = labelin # name of the boundary e.g. outerWall
        self.typ = typin # name of boundary type e.g. patch for blockMeshDict
        self.flist = [] # list of faces, which are lists of point indices
        self.alphalist = [] # list of alpha.ink properties
        self.Ulist = [] # list of U properties
        self.plist = [] # list of p_rgh properties
        self.meshi = [] # mesh for exporting to stl
        self.reflev = 0
        
        
#-------------------------------------------------        
#---------------------PLOTS-----------------------
#-------------------------------------------------




class NozVars:
    '''this is where the geometry of the nozzle is defined'''
    
    def __init__(self, bathWidth:float=16, bathHeight:float=7, bathDepth:float=7, frontWidth:float=4, vink:float=10, vbath:float=10, npts:int=50, nozzleInnerWidth:float=0.603, nozzleThickness:float=0.152, nozzleAngle:float=0, horizontal:bool=False, adjacent:str='None', distance:float=0, **kwargs):
        ''' Allowed input variables:
            bathWidth: (default=16) bath width in nozzle inner diameters
            bathHeight: (default=7) bath height in nozzle inner diameters
            bathDepth: (default=7) bath depth in nozzle inner diameters
            frontWidth: (default=4) front of nozzle bath width in nozzle inner diameters
            vink: (default=10) ink extrusion speed in mm/s
            vbath: (default=10) bath translation speed in mm/s
            npts: (default=50) number of points in the circle used to define the nozzle
            nozzleInnerWidth: (default=0.603) inner diameter of the nozzle in mm
            nozzleThickness: (default=0.152) nozzle wall thickness in mm
            nozzleAngle: (default=0) nozzle angle in degrees
            horizontal: (default=False) whether to have the nozzle horizontal
            adjacent: (default=None) adjacent filament orientation
            distance: (default=-1) offset of adjacent filament in nozzle inner diameters
        '''
        self.niw = nozzleInnerWidth # nozzle inner width
        self.nt = nozzleThickness # nozzle thickness
        self.bw = bathWidth*self.niw # bath width (x)
        self.bh = bathHeight*self.niw # bath height (y)
        self.bd = bathDepth*self.niw # bath depth (z)
        self.nl = self.bd/2-self.niw/2 # nozzle length
        self.na = nozzleAngle # nozzle angle RG
        self.hor = horizontal # vertical or horizontal RG
        self.adj = adjacent # adjacent filament orientation RG
        self.dst = distance*self.niw # adjacent filament offset RG
        try:
            reffolder = kwargs.get('reffolder') # RG
            self.cor = reffolder[reffolder.find('nb'):] # RG
        except AttributeError:
            self.cor = 'None'
        
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
        
        bottomRadius = self.niw/2 # nozzle radius outlet RG
        topRadius = bottomRadius+(self.nl*np.tan(np.deg2rad(self.na))) # nozzle radius inlet RG
        
        if 2*topRadius>self.bh-4*self.niw: # resizes support bath for large nozzle angles RG
            oldbh = self.bh
            self.bh = 2*topRadius+4*self.niw
            scale = self.bh/oldbh
            self.bw = self.bw*scale
            frontWidth = frontWidth*scale
        
        if self.hor:
            temp = self.bw # flip x and z dimensions
            self.bw = self.bd
            self.bd = temp
            
        self.ble = -self.bw/2 # bath left coord
        self.bri = self.bw/2 # bath right
        self.bfr = -self.bh/2 # bath front
        self.bba = self.bh/2 # bath back
        self.bbo = -self.bd/2 # bath bottom
        self.bto = self.bd/2 # bath top
        self.nbo = self.bto - self.nl # nozzle bottom z coord
        self.ncx = self.ble + frontWidth*self.niw # nozzle center bottom x coord
        self.ncy = self.bfr + self.bh/2 # nozzle center bottom y coord
        
        if self.hor:
            self.ncx = 0 # center nozzle
        
        npts = int(np.ceil(npts*topRadius/bottomRadius)) # RG
        
        self.inptsb = circlePoints(bottomRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the bottom of the nozzle
        self.outptsb = circlePoints(bottomRadius + self.nt, npts)+[self.ncx, self.ncy] 
            # outer radius points at the bottom
        self.inptst = circlePoints(topRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the top of the nozzle (the inlet)
        self.outptst = circlePoints(topRadius + self.nt, npts)+[self.ncx, self.ncy]
            # outer radius points at the top of the nozzle

        self.bv = vbath*SCALE # bath velocity
        self.iv = vink*SCALE*bottomRadius**2/topRadius**2 # initial ink velocity RG

        #--------------------------------------------------------------------------------------------------------

        


def plotGeo(geo:NozVars, pl:np.array, blocks:List[Block], br:List[BoundaryInput]) -> plt.Figure:
    '''plot the requested geometry
    geo is a nozzleVars object
    pl is a point list used for blockMesh
    blocks is a list of blocks used for blockMesh
    br is real boundaries used for generating stls and setting boundary conditions'''
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
        ax.text(np.mean(vs[:,0]), np.mean(vs[:,1]), np.mean(vs[:,2]), i, color='r', fontsize=fs)
    clist = ['tomato', 'navajowhite', 'mediumaquamarine', 'deepskyblue', 'white']
    plist = []
    for bi in br:
        # each boundary bi in list bl
        col = clist.pop(0) # pop a color from the front of the list
        plist.append(mpatches.Patch(color=col, label=bi.label))
        for m in bi.meshi['vectors']:
            p3c = Poly3DCollection([list(zip(m[:,0],m[:,1],m[:,2]))], alpha=0.35)
            p3c.set_facecolor(col)
            p3c.set_edgecolor('gray')
            ax.add_collection3d(p3c)      
    ax.legend(handles=plist, fontsize=fs)
    return fig
   #--------------------------------------------------------------------------------------------------------
# geometry variables        
        
 
      

  #--------------------------------------------------------------------------------------------------------  
#-------------------------------------------------
#--------------------FUNCTIONS--------------------
#-------------------------------------------------




#-------------------------------------------------
#----------------CREATING GEOMETRIES--------------

def circlePoints(r:float, npts:int) -> List:
    '''get a list of points in a circle. Inputs: r, npts
    r is the radius of the circle
    npts is the number of points'''
    pts = np.zeros([npts+1, 2])
    tstep = 2*np.pi/(npts)
    theta = 0
    for i in range(npts):
        pts[i] = ([np.cos(theta)*r, np.sin(theta)*r])
        theta = theta + tstep
    pts[npts] = pts[0]
    return pts
    

def createNozzle(nv:NozVars) -> np.array:
    '''get list of vertices in the block mesh
    Input: nv is a NozVars object'''
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


def blockCornerList(pl:np.array) -> List:
    '''select the points that are the bottom corner of the blockCornerList
    Input: pl is a point list'''
    xmax = max(pl[:, 0])
    ymax = max(pl[:, 1])
    zmax = max(pl[:, 2])
    corners = pl[(pl[:,0]<xmax) & (pl[:,1]<ymax) & (pl[:,2]<zmax), :]
    return corners
 

def ptFromPts(pts:np.array, x:float, y:float, z:float) -> List[float]:
    '''select the point (x,y,z) from an array of points pts
    Input: pts, x, y, z'''
    pt = pts[(pts[:,0]==x) & (pts[:,1]==y) & (pts[:,2]==z), :]
    return pt[0]


def blockPts(pl:float, corner:List[float]) -> np.array:
    '''select the points in this block based on the first corner, and put the points in the order that openFOAM likes
    pl is an array containing points
    corner is one point which marks the first corner of the block'''
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
        pts[0+4*k, :] = ptFromPts(pl, xlist[cx], ylist[cy], zlist[cz+k])
        pts[1+4*k, :] = ptFromPts(pl, xlist[cx+1], ylist[cy], zlist[cz+k])
        pts[2+4*k, :] = ptFromPts(pl, xlist[cx+1], ylist[cy+1], zlist[cz+k])
        pts[3+4*k, :] = ptFromPts(pl, xlist[cx], ylist[cy+1], zlist[cz+k])
    return pts

#--------------------------------------------------
#----------------------STLs------------------------


def axisFace(xyz:List[float]) -> np.array:
    '''given two two-element lists and one one-element list, determine the two triangles that constitute that face. For example, two x values, two y values, and a z value'''
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


def ptFace(pts:np.array) -> np.array:
    '''given an array of 4 points pts, construct a face. The line from point 1 to point 2 must cross through the center of the face'''
    data = np.zeros(2, dtype=mesh.Mesh.dtype)
    data['vectors'][0] = np.array([pts[0], pts[1], pts[2]])
    data['vectors'][1] = np.array([pts[1], pts[2], pts[3]])
    return data


def setZ(pts2d:np.array, z:float) -> np.array:
    '''given an array of 2d points pts2d, create an array of 3d points at position z'''
    out = np.zeros([len(pts2d), 3])
    out[:, 0:2] = pts2d
    out[:, 2] = z
    return out

def setX(pts2d:np.array, x:float) -> np.array: # RG
    '''given an array of 2d points pts2d, create an array of 3d points at position x'''
    out = np.zeros([len(pts2d), 3])
    out[:, 1:3] = pts2d
    out[:, 0] = x
    return out


def arcFace(ina:np.array, outa:np.array) -> np.array:
    '''an arc face could be a donut on a plane (e.g. the bottom of the nozzle) if both lists have the same x,y,or z
    it could be a circle or cone (e.g. the nozzle inlet) if ina has length 1
    it could be a cylinder or frustum if the two lists have different x,y,and z
    ina is a list of points on the inner radius of the arc
    outa is a list of points on the outer radius of the arc'''
    data = np.zeros(2*len(ina), dtype = mesh.Mesh.dtype)
    if len(ina)==1:
        ina2 = np.ones([len(outa), 3])
        ina2[:,:] = ina
        ina = ina2
    for i in range(len(outa)-1):
        d = ptFace([ina[i], outa[i], ina[i+1], outa[i+1]])
        data['vectors'][2*i] = d['vectors'][0]
        data['vectors'][2*i+1] = d['vectors'][1]
    return data

def sortPolar(pts:np.array) -> List[float]: # RG
    '''sort a list of points from 0 to 2 pi
    pts is a 3d list of unarranged points lying in the x, y, or z plane'''
    # maps x,y,z so z is the plane the xs lies in
    for i in range(3):
        if np.all(pts[:,i] == pts[0,i]):
            x = pts[:,(i+1)%3]
            y = pts[:,(i+2)%3]
            z = pts[:,i]
            j = i
    # organizes points by polar coordinates with the center as the origin
    x0 = np.mean(x)
    y0 = np.mean(y)
    r = np.sqrt((x-x0)**2+(y-y0)**2)
    theta = np.where(y>y0, np.arccos((x-x0)/r), 2*np.pi-np.arccos((x-x0)/r))
    mask = np.argsort(theta)
    xsort = x[mask]
    ysort = y[mask]
    xsort = np.append(xsort,xsort[0])
    ysort = np.append(ysort,ysort[0])
    z = np.append(z,z[0])
    # maps x,y,z back to original
    ptstemp = np.asarray(list(zip(xsort,ysort,z)))
    ptssort = np.zeros((len(xsort),3))
    for k in range(3):
        ptssort[:,(j+k+1)%3]=ptstemp[:,k]
    return ptssort, theta


def holeInPlane(cpts:np.array, x:List[float], y:List[float], z:List[float]) -> np.array:
    '''get a mesh array that describes a hole in a plane
    cpts is a list of circle points. cpts should be in order from 0 to 2 pi
    x is a list of 2 x values for the plane
    y is a list of 2 y values for the plane
    z is a scalar that the plane lies on'''
    n = len(cpts)-1 # number of points on the circle
    nchunks = int(np.floor(n/4)) # split the circle into four chunks, number of points in a chunk
    data = np.zeros(n+4, dtype = mesh.Mesh.dtype)
    x.sort()
    y.sort()
    z.sort() # RG
    if np.all(cpts[:,0] == cpts[0,0]): # lies in x plane RG
        corners = [[x[0], y[1], z[0]], [x[0], y[1], z[1]], [x[0], y[0], z[1]], [x[0], y[0], z[0]],  [x[0], y[1], z[0]]]
    elif np.all(cpts[:,1] == cpts[0,1]): # lies in y plane
        corners = [[x[1], y[0], z[0]], [x[1], y[0], z[1]], [x[0], y[0], z[1]], [x[0], y[0], z[0]],  [x[1], y[0], z[0]]]
    elif np.all(cpts[:,2] == cpts[0,2]): # lies in z plane
        corners = [[x[1], y[0], z[0]], [x[1], y[1], z[0]], [x[0], y[1], z[0]], [x[0], y[0], z[0]],  [x[1], y[0], z[0]]]
    for i in range(4):
        data['vectors'][i*(nchunks+1)] = np.array([cpts[i*nchunks], corners[i], corners[i+1]])
        for j in range(nchunks):
            data['vectors'][i*(nchunks+1)+j+1] = np.array([cpts[i*nchunks+j], cpts[i*nchunks+j+1], corners[i+1]])
    premain = cpts[4*nchunks:]
    for j in range(len(premain)-1):
        data['vectors'][4*nchunks+4+j] = np.array([premain[j], premain[j+1], corners[-1]])
    return data


def combineMeshes(meshList:List[List]) -> np.array:
    '''combine all meshes into one list
    meshlist is a list of lists of triangles'''
    n = functools.reduce(lambda a,b : a+len(b), meshList, 0) # total number of mesh triangles
    data = np.zeros(n, dtype=mesh.Mesh.dtype)
    di = 0
    for m in meshList:
        for mi in m:
            data[di] = mi
            di+=1
    return data
    

    
#-------------------------------------------------   
#-----------------SYSTEM--------------------------


#-------------------------------------------------
#---------------blockMeshDict---------------------


def vec2cpp(v:List[float]) -> str:
    '''convert a vector to the format it needs to be in for OpenFOAM to read it
    v is a list'''

    s = "("
    for vi in v:
        s = s + str(vi) + " "
    s = s + ")"
    return s


def block2txt(block:Block) -> str:
    '''convert block to openfoam string'''
    s = "hex " + vec2cpp(block.vertices[:, 3].astype(int)) + " " + vec2cpp(block.meshi)
    s = s + " simpleGrading " + vec2cpp(block.grading)
    return s

def blocks2txt(blocks:List[Block]) -> str:
    '''convert list of blocks to openfoam string'''
    pl = DictList("blocks", 1, [])
    for b in blocks:
        pl.proplist.append(block2txt(b))
    return pl.prnt(0)
    
EDGES = "edges\n(\n);"


def faceSelector(block:Block, st:str) -> np.array:
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


def boundaryList(blocks:List[Block]) -> List[BoundaryInput]:
    '''compiles a list of all the boundaries in the system for blockMesh
    because we are using snappyHexMesh, we only need one boundary in blockMesh'''
    allb = BoundaryInput("allBoundary", "patch")
    for st in ["x-", "x+", "y-", "y+", "z-", "z+"]:
        allb.flist.append(faceSelector(blocks[0], st)) # front and back faces
    return [allb]


def walsel(geo:NozVars, st:str):
    '''select the boundaries of a wall
    geo is a NozVars object
    st is a string indicating which face to take'''
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


def realBoundaries(geo:NozVars, exportMesh:bool, **kwargs) -> List[BoundaryInput]:
    '''get a list of boundaries for snappyHexMesh, setFields
    geo is a NozVars object geo.bv and geo.iv are bath velocities and ink velocities, scaled to m/s
    mv is a MeshVars object
    exportMesh is true if we want to create a mesh folder'''
    bf = BoundaryInput("bathFlow", "")
    bf.alphalist = DictList(bf.label, 0, [["type", "fixedValue"], ["value", "uniform 0"]])
    if geo.hor:
        bf.Ulist = DictList(bf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(geo.bv) + ")"]]) # flow in -z RG
    else:
        bf.Ulist = DictList(bf.label, 0, [["type", "fixedValue"], ["value", "uniform (" + str(geo.bv) + " 0 0)"]])
    bf.plist = DictList(bf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    if exportMesh:
        if geo.hor:
            faces = list(map(lambda s: axisFace(walsel(geo, s)), ["x+", "y-", "y+"])) # side walls RG
        else:
            faces = list(map(lambda s: axisFace(walsel(geo, s)), ["y-", "y+", "z-"]))
        if geo.adj=='y': # RG
            reffolder = kwargs.get('reffolder')
            x = geo.ble+(geo.ncx-geo.ble)*geo.niw+8*geo.niw
            p = intm.posSlice(reffolder, x)
            xspts, cent = intm.xspoints(geo, p, geo.dst)
            xspts3d = setX(xspts,geo.ble)
            xspts = sortPolar(xspts3d)[0]
            faces.insert(0, holeInPlane(xspts, [geo.ble, geo.ble], [geo.bfr, geo.bba], [geo.bto, geo.bbo]))
        else:
            faces.append(axisFace(walsel(geo, "x-")))
        bf.meshi = combineMeshes(faces)
        
    inkf = BoundaryInput("inkFlow", "")
    inkf.alphalist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform 1"]])
    inkf.Ulist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(geo.iv) + ")"]])
    inkf.plist = DictList(inkf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    if exportMesh:
        inkf.meshi = arcFace(setZ(np.zeros([len(geo.inptst), 2])+[geo.ncx, geo.ncy], geo.bto), setZ(geo.inptst, geo.bto))
        
    inkf.reflev = 0
    
    at = BoundaryInput("atmosphere", "")
    at.alphalist = DictList(at.label, 0, [["type", "inletOutlet"], ["value", "uniform 0"], ["inletValue", "uniform 0"]])
    at.Ulist = DictList(at.label, 0, [["type", "pressureInletOutletVelocity"], ["value", "uniform (0 0 0)"]])
    at.plist = DictList(at.label, 0, [["type", "totalPressure"], ["p0", "uniform 0"]])
    if exportMesh:
        if geo.hor:
            at.meshi = combineMeshes([holeInPlane(setZ(geo.outptst, geo.bto), [geo.ble, geo.bri], [geo.bfr, geo.bba], [geo.bto, geo.bto]), \
                           axisFace(walsel(geo, "z-"))]) # top and bottom walls RG
        else:
            at.meshi = combineMeshes([holeInPlane(setZ(geo.outptst, geo.bto), [geo.ble, geo.bri], [geo.bfr, geo.bba], [geo.bto, geo.bto]), \
                           axisFace(walsel(geo, "x+"))])
#         at.meshi = combineMeshes([axisFace(walsel(geo, "z+")), \
#                            axisFace(walsel(geo, "x+"))])
    
    fw = BoundaryInput("fixedWalls", "")
    fw.alphalist = DictList(fw.label, 0, [["type", "zeroGradient"]])
    fw.Ulist = DictList(fw.label, 0, [["type", "noSlip"]])  
    fw.plist = DictList(fw.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
    
    if exportMesh:
        fw.meshi = combineMeshes([arcFace(setZ(geo.inptsb, geo.nbo), setZ(geo.outptsb, geo.nbo)),\
                            arcFace(setZ(geo.inptsb, geo.nbo), setZ(geo.inptst, geo.bto)),\
                            arcFace(setZ(geo.outptsb, geo.nbo), setZ(geo.outptst, geo.bto))])
    fw.reflev = 2
    
    if geo.adj=='y': # RG
        inkx = BoundaryInput("inkxs", "")
        inkx.alphalist = DictList(inkx.label, 0, [["type", "fixedValue"], ["value", "uniform 1"]])
        inkx.Ulist = DictList(inkx.label, 0, [["type", "fixedValue"], ["value", "uniform (" + str(geo.bv) + " 0 0)"]])
        inkx.plist = DictList(inkx.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
        if exportMesh:
            xspts2 = np.copy(xspts)
            xspts2[:,0]+=geo.niw/3 # gives xs depth so snappyHexMesh can snap to its surface
            inkx.meshi = combineMeshes([arcFace(xspts, xspts2), arcFace(setX(np.zeros([len(xspts),2])+cent, geo.ble), xspts)])
        
        inkx.reflev = 2
        return [bf, inkf, inkx, at, fw]
 
    return [bf, inkf, at, fw]
    

def boundarycpp(bl:List[BoundaryInput]) -> str: 
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


def compileBlockMeshDict(pl:np.array, blocks:List[Block], bl:List[BoundaryInput]) -> str:
    '''get the blockMeshDict text
    pl is a point list
    blocks is a list of Block objects
    bl is a list of BoundaryInput objects'''
    s = header("dictionary", "blockMeshDict")
    s = s + "scale " + str(SCALE) + ";\n\n"
    s = s + DictList("vertices", 1, pl[:, 0:3]).prnt(0)
    s = s + blocks2txt(blocks)
    s = s + "edges\n(\n);\n\n"
    s = s + boundarycpp(bl)
    s = s + "mergePatchPairs\n(\n);\n\n" + CLOSELINE
    return s   
 
    

#---------------------------------------------------
#-----------------surfaceFeatureExtractDict--------

def compileSurfaceFeatureExtractDict(bl:List[BoundaryInput]) -> str:
    '''compile surfaceFeatureExtractDict'''
    bnames = [o.label for o in bl]
    s = header("dictionary", "surfaceFeatureExtractDict")
    s = s + DictList("", 0, \
                     list(map(lambda x: DictList(x+".stl", 0, \
                                                 [["extractionMethod", "extractFromSurface"],\
                                                  DictList("extractFromSurfaceCoeffs", 0, [["includedAngle", 180]]),\
                                                  ["writeObj", "yes"]]), bnames))).prnt(-1)
    s = s + CLOSELINE
    return s



#---------------------------------------------------
#-----------------surfaceFeaturesDict--------------- RG

def compileSurfaceFeaturesDict(bl:List[BoundaryInput]) -> str:
    '''compile surfaceFeaturesDict'''
    bnames = [o.label for o in bl]
    s = header("dictionary", "surfaceFeaturesDict")
    s = s + DictList("surfaces", 1, ['"' + b + '.stl"' for b in bnames]).prnt(0)
    s = s + DictList("", 0, [["includedAngle", 180]]).prnt(-1)
    s = s + DictList("", 0, [["writeObj", "yes"]]).prnt(-1)
    s = s + CLOSELINE
    return s


#    bnames = [o.label for o in bl]
#    s = header("dictionary", "surfaceFeaturesDict")
#    s = s + "surfaces",\
#                     "("
#    s = s + DictList("", 0, \
#                     list(map(lambda x: DictList(x+".stl", 0, \))).prnt(-1)
#                     ");", \
#                     DictList(["includedAngle", 180], \
#                     ["writeObj", "yes"]),

#    s = s + CLOSELINE
#    return s



#-------------------------------------------------   
#----------------------0--------------------------
#-------------------------------------------------  


def format0(bl:List[BoundaryInput], classn:str, obj:str, dims:str, intfield:str, f:Callable) -> str:
    '''puts files in the 0 folder into the correct format
    classn is the class that goes into the header, e.g. "volScalarField" or "pointScalarField"
    obj is the object that goes into the header, e.g. "U"
    dims is the dimensions of the field, e.g. "[1 -1 -2 0 0 0 0]", which means kg/(m*s^2)
    intfield is the internal field, e.g. "0" or "(0 0 0)"
    # f is the function to run on each boundaryInput object to get a list of properties of interest'''
    s = header(classn, obj) # header
    simplelist = DictList("", 0, [["dimensions", dims], ["internalField", "uniform " + intfield]]) # list of simple variables to define
    s = s + simplelist.prnt(-1) # format the simple list 
    bflist = DictList("boundaryField", 0, []) # list of boundary dictionary entries
    for b in bl:
        bflist.proplist.append(f(b)) # add the appropriate list for each boundary, e.g. alpha properties for inkFlow
    s = s + bflist.prnt(0) # format the dictionary list
    s = s + CLOSELINE
    return s


def compileAlphaOrig(bl:List[BoundaryInput]) -> str:
    return format0(bl, "volScalarField", "alpha.ink", "[0 0 0 0 0 0 0]", "0", lambda bi:bi.alphalist)
    
def compileU(bl:List[BoundaryInput]) -> str:
    return format0(bl, "volVectorField", "U", "[0 1 -1 0 0 0 0]", "(0.01 0 0)", lambda bi:bi.Ulist) # RG
    
def compileP(bl:List[BoundaryInput]) -> str:
    return format0(bl, "volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]", "0", lambda bi:bi.plist)

def compileCellLevel(bl:List[BoundaryInput]) -> str:
    return format0(bl, "volScalarField", "cellLevel", "[0 0 0 0 0 0 0]", "0", lambda bi:[["type", "zeroGradient"]])

def compilePointLevel(bl:List[BoundaryInput]) -> str:
    return format0(bl, "pointScalarField", "pointLevel", "[0 0 0 0 0 0 0]", "0", lambda bi:[["type", "calculated"]])


#-------------------------------------------------
#------------------setFieldsDict----------------------

def compileSetFieldsDict(geo:NozVars) -> str: # fills the nozzle with ink at t=0 RG
    '''compile setFieldsDict'''
    s = header("dictionary", "setFieldsDict")
    s = s + DictList("defaultFieldValues", 1, ["volScalarFieldValue alpha.ink 0"]).prnt(0)
    r = DictList("regions", 1, [])                 
        # r is the region where the ink originally is
        
    if geo.na != 0:
        delta = 0.9*geo.nt/np.tan(np.deg2rad(geo.na)) # length of cylinders
    else:
        delta = 1 # avoids divide by 0 error
        
    steps = max(np.ceil((geo.bto-geo.nbo)/delta).astype('int'),1) # number of cylinders to create
    for x in range(steps): # creates thin enough cylinders of ink so that inside of the nozzle becomes fully ink but none of the bath becomes ink
        c2c = DictList("cylinderToCell", 0, [])         
            # c2c is a cylinderToCell dictionary list
        c2c.proplist.append(["p1" + " (" + str(geo.ncx*SCALE) + " " + str(geo.ncy*SCALE) + " " + str(round((geo.nbo+(geo.bto-geo.nbo)*x/steps)*SCALE,8)) + ")"])
        c2c.proplist.append(["p2" + " (" + str(geo.ncx*SCALE) + " " + str(geo.ncy*SCALE) + " " + str(round((geo.bto-(geo.bto-geo.nbo)*(steps-1-x)/steps)*SCALE,8)) + ")"])
                    # top and bottom points of central cylinder axis
        c2c.proplist.append(["radius " + str(str((geo.niw/2+(x+1)*geo.nl/steps*np.tan(np.deg2rad(geo.na)))*SCALE))])
                    #radius of cylinder
        c2c.proplist.append(DictList("fieldValues", 1, ["volScalarFieldValue alpha.ink 1"]))
                    # value of alpha inside the cylinder is 1 because it is full of ink
        r.proplist.append(c2c)
        
    # create an initial line of ink RG
    # if geo.adj=='y':
    #     s2c = DictList("surfaceToCell", 0, [])
    #     s2c.proplist.append(['file "inkLine.stl"'])
    #     s2c.proplist.append(["outsidePoints ((0 0 -0.0015))"])
    #     s2c.proplist.append(["includeCut true"])
    #     s2c.proplist.append(["includeInside true"])
    #     s2c.proplist.append(["includeOutside false"])
    #     s2c.proplist.append(["nearDistance -1"])
    #     s2c.proplist.append(["curvature 1"])
    #     s2c.proplist.append(["fileType stl"])
    #     s2c.proplist.append(DictList("fieldValues", 1, ["volScalarFieldValue alpha.ink 1"]))
    #     r.proplist.append(s2c)
                         
    s = s + r.prnt(0)
    s = s + CLOSELINE
    return s
   

def geometryFile(geo:NozVars) -> str:
    '''geometryFile gets a csv string of all of the geometry variables we care about'''
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
         ['nozzle angle (degrees)', geo.na], \
         ['horizontal', geo.hor], \
         ['adjacent filament oriantation', geo.adj], \
         ['adjacent filament offset (mm)', geo.dst], \
         ['corresponding simulation', geo.cor], \
         ['bath velocity (m/s)', geo.bv], \
         ['ink velocity (m/s)', geo.iv]]
    s = ""
    for li in l:
        s = s + li[0] + ', ' + str(li[1]) + '\n'
    return s
         
    
    

    
    
    
    #----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## MESH
 
 #--------------------------------------------------------------------------------------------------------

class MeshVars:
    '''mesh variables'''
    
    fwreflev = 2 # RG
    inkfreflev = 0
    inkxreflev = 2 # RG
    nCellsBetweenLevels = 5
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
        
    def cmc(self) -> List[List[str]]:
        '''generates a list of variables for the castellatedMeshControls file'''
        
        return(self.varList(["maxLocalCells", "maxGlobalCells", "minRefinementCells", "nCellsBetweenLevels", "resolveFeatureAngle", "locationInMesh", "allowFreeStandingZoneFaces"]))
       
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
        
    def sc(self) -> List[List[str]]:
        '''generates a list of variables for the snapControls file'''
        
        return(self.varList(["nSmoothPatch", "tolerance", "nSolveIter", "nRelaxIter", "nFeatureSnapIter", "implicitFeatureSnap", "explicitFeatureSnap", "multiRegionFeatureSnap"]))
    
    # addLayersControls
    addLayers = "false" # RG
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
    nSmoothNormals = 15
        # Number of smoothing iterations of interior mesh movement direction
    nSmoothSurfaceNormals = 10
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
    unrefineLevel = 4 # RG
        # number of times cells can be coarsened
    nBufferLayers = 1
        # number of layers around a refined cell
    maxRefinement = 4
        # max number of times cells can be refined
    maxCells = maxGlobalCells*10 # RG
        # total number of cells in the mesh
    dumpLevel = "false"
        # writes the refinement level for each cell as a volScalarField

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
    
    
    def varList(self, li) -> List[List[str]]:
        '''Given a list of variable names li, construct a table with those variables.'''
        
        out = [["", ""] for x in range(len(li))]
        for i, l in enumerate(li):
            out[i] = [l, getattr(self, l)]
        return out
    
    

 



 #---------------------------------------------
#-------------------meshQualityDict-----------


def compileMeshQualityDict() -> str:
    '''compile meshQualityDict'''
    s = header("dictionary", "meshQualityDict")
    s = s + "#includeEtc \"caseDicts/mesh/generation/meshQualityDict\"\n\n" # RG
    s = s + CLOSELINE
    return s

#---------------------------------------------
#----------------snappyHexMeshDict-------------


def compileSnappyHexMeshDict(bl:List[BoundaryInput], mv:MeshVars) -> str:
    '''compile SnappyHexMeshDict'''
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


#-------------------------------------------------
#------------dynamicMeshDict-----------------

def compileDynamicMeshDict(mv:MeshVars) -> str:
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
        

    
    
    #----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## SOLVERS
        
       

    #--------------------------------------------------------------------------------------------------------
    
class CDVars:
    '''controlDict variables '''
        
    def __init__(self, startTime:float, endTime:float, dt:float, writeDt:float):
        '''Inputs: startTime, endTime, dt, writeDt
        startTime is in seconds
        endTime is in seconds
        dt is the initial time step of the simulation in seconds
        writeDt is the time step at which to write results to file'''
        # solver
        self.application = "interFoam"
        
        # time control
        self.startFrom = "latestTime"
        self.startTime = startTime
        self.stopAt = "endTime"
        self.endTime = endTime
        self.deltaT = dt # Time step of the simulation
        
        # time step control
        self.adjustTimeStep = "yes" # yes/no† to adjust time step according to maximum Courant number in transient simulation.
        self.maxCo = 1 # Maximum Courant number allowed.
        self.maxAlphaCo = 1
        self.maxDeltaT = 1 # Maximum time step allowed in transient simulation
        
        # data writing
        #self.writeControl = "adjustable" # Writes data every writeInterval seconds of simulated time.
        self.writeControl = "adjustableRunTime"
        self.writeInterval = writeDt # Scalar used in conjunction with writeControl described above.
        self.purgeWrite = 0 # Integer representing a limit on the number of time directories that are stored by overwriting time directories on a cyclic basis. Example of t0 = 5s, Δt = 1s and purgeWrite 2;: data written into 2 directories, 6 and 7, before returning to write the data at 8 s in 6, data at 9 s into 7, etc. To disable the time directory limit, specify purgeWrite 0;† For steady-state solutions, results from previous iterations can be continuously overwritten by specifying purgeWrite 1;
        self.writeFormat = "ascii" # ASCII format, written to writePrecision significant figures.
        self.writePrecision = 6 # Integer used in conjunction with writeFormat described above, 6† by default
        self.writeCompression = "uncompressed" # Specifies the compression of the data files.
        self.timeFormat = "general" # Specifies scientific format if the exponent is less than -4 or greater than or equal to that specified by timePrecision.
        self.timePrecision = 6 # Integer used in conjunction with timeFormat described above, 6† by default
        
        # data reading
        self.runTimeModifiable = "yes" # yes†/no switch for whether dictionaries, e.g.controlDict, are re-read by OpenFOAM at the beginning of each time step.

    
    def varList(self) -> List[List[str]]:
        '''Constructs a table with all of the variables held in this object. The table has two columns: the variable name and the value'''
        l = []        
        for attr, value in self.__dict__.items():
            l.append([attr, value])
        return l

    #--------------------------------------------------------------------------------------------------------

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

    #--------------------------------------------------------------------------------------------------------

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
    
    
    

#-------------------------------------------------
#---------------controlDict-----------------------


def compileControlDict(cdv:CDVars) -> str:
    '''gets the text for controlDict'''
    s = header("dictionary", "controlDict")
    v = cdv.varList()
    #v.append(["#sinclude", "sampling"])
    v.append(['#includeIfPresent', '\"sampling\"'])
    s = s + DictList("", 0,v).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#---------------fvSchemes-----------------------  

def compileFvSchemes(fvv:Union[FVVars, str]) -> str:
    '''gets the text for the fvSchemes file. Input empty string to get file for mesh folder'''
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
def compileFvSolution(fvv:Union[FVVars, str]) -> str:
    '''gets the text for fvSolution. Input empty string to get empty file for mesh folder'''
    s = header("dictionary", "fvSolution")
    if isinstance(fvv, FVVars):
        s = s + fvv.solverlist().prnt(0)
        s = s + fvv.pimple().prnt(0)
        s = s + DictList("relaxationFactors", 0, [DictList("equations", 0, [["\".*\" 1"]])]).prnt(0)
    s = s + CLOSELINE
    return s

def compileFV(fvv:FVVars, out:FileGroup) -> FileGroup:
    '''compiles fvsolution and fvschemes files'''
    out.fvSchemes = compileFvSchemes(fvv)
    out.fvSolution = compileFvSolution(fvv)
    return out

#-------------------------------------------------
#---------------all solver files-----------------------


def solverObjects(startTime:float, endTime:float, dt:float, writeDt:float, solver:str) -> Tuple[CDVars, FVVars]:
    '''gets the control dictionary variable object and fvsolution and fvschemes dictionary variable object
    starttime is the start time for the simulation in s
    endtime is the end time for the simulation in s
    dt is the initial solve time step in s
    writedt is the write time step in s
    solver is interFoam or interIsoFoam'''
    cdv = CDVars(startTime, endTime, dt, writeDt)
    cdv.application = solver
    fvv = FVVars(solver)
    return cdv, fvv


def compileSolverFiles(cdv:CDVars, fvv:FVVars, out:FileGroup) -> FileGroup:
    '''generates the text for the controlDict, FVsolution, FVschemes, and bash scripts'''
    out.controlDict = compileControlDict(cdv)
    out = compileFV(fvv, out)
    out.slurm = compileSlurm(out.folder, out.slurmFolder)
    out.allallrun = compileAllAllRun()
    out.allclean = compileAllClean()
    out.allrunmesh = compileAllRunMesh(out.folder)
    out.allrun = compileAllRun(out.folder, cdv.application)
#     out.cont = compileContinue(out.folder, cdv.application)
    cdv.application = "icoFoam"
    out.controlDictMesh = compileControlDict(cdv)
    return out





#-------------------------------------------------   
#-----------------CONSTANT--------------------------
#-------------------------------------------------  

#-------------------------------------------------
#-----------------gravity-------------------------


def compileG() -> str:
    '''compile g'''
    s = header("uniformDimensionedVectorField", "g")
    s = s + DictList("", 0, [["dimensions", "[0 1 -2 0 0 0 0]"], ["value", "(0 0 -9.81)"]]).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#------------------transportProperties------------


def transportGroupNewt(title:str, nu:Union[float, str], rho:Union[float, str]) -> DictList:
    '''gets a DictList object that describes the transport properties
    title is the name of the phase, e.g. 'ink'
    nu is the kinematic viscosity in m^2/s
    rho is the density in kg/m^3'''
    return DictList(title, 0, [["transportModel", "Newtonian"], ["nu", str(nu)], ["rho", str(rho)]])


def transportGroupHB(title:str, nu0:Union[float, str], tau0:Union[float, str], k:Union[float, str], n:Union[float, str], rho:Union[float, str]) -> DictList:
    '''transportGroupHB gets a DictList that describes Herschel-Bulkley transport properties
    title is the name of the phase
    the HB model follows nu = min(nu0, tau0/gammadot + k*gammadot^(n-1))
    inputs to the model are nu0 in m^2/s, tau0 in m^2/s^2, k in m^2/s, n unitless
    rho is the density in kg/m^3'''
    return DictList(title, 0, \
                  [["transportModel", "HerschelBulkley"], \
                  DictList("HerschelBulkleyCoeffs", 0, \
                           [["nu0", "[ 0 2 -1 0 0 0 0 ] " + str(nu0)],\
                            ["tau0", "[ 0 2 -2 0 0 0 0 ] " + str(tau0)], \
                            ["k", "[ 0 2 -1 0 0 0 0 ] " + str(k)],\
                            ["n", "[ 0 0 0 0 0 0 0 ] " + str(n)]]\
                           ),\
                  ["rho", str(rho)]])


def compileTransportProperties(ink:DictList, sup:DictList, sigma:Union[float, str]) -> str:
    '''compile transportProperties
    ink and sup are DictLists created by transportGroupNewt and transportGroupHB
    sigma is the surface tension in J/m^2'''
    s = header("dictionary", "transportProperties")
    s = s + DictList("", 0, [["phases (ink sup)"],ink, sup, ["sigma", str(sigma)]]).prnt(-1)
    s = s + CLOSELINE
    return s

#-------------------------------------------------
#------------turbulenceProperties-----------------

def compileTurbulenceProperties() ->str:
    s = header("dictionary", "turbulenceProperties")
    s = s + DictList("", 0, [["simulationType", "laminar"]]).prnt(-1)
    s = s + CLOSELINE
    return s
 
     
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## PUTTING IT ALL TOGETHER
    
#------------------------------------------------
#--------------------EXPORTING-------------------



def saveStls(folder:str, bl:List[BoundaryInput]) -> None:
    '''exports stls
    folder is a full path name
    bl is a list of boundaryInput objects '''
    for b in bl:
        me = b.meshi
        me2 = mesh.Mesh(np.zeros(len(me), dtype=mesh.Mesh.dtype))
        for i,m in enumerate(me['vectors']):
            me2.vectors[i] = m*SCALE
        title = os.path.join(folder, b.label+".stl")
        me2.save(title)


def exportFile(folder:str, file:str, text:str, linux:bool=False) -> None:
    '''exports text files
    folder is a full path name
    file is a basename
    text is the text to export'''
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    logging.debug("Exported %s" % fn)
    if linux:
        replaceCR(fn)
    
def replaceCR(file:str) -> None:
    '''replace windows carriage return with linux carriage return. Otherwise bash scripts will throw error. From https://stackoverflow.com/questions/36422107/how-to-convert-crlf-to-lf-on-a-windows-machine-in-python'''
    WINDOWS_LINE_ENDING = b'\r\n'
    UNIX_LINE_ENDING = b'\n'

    with open(file, 'rb') as open_file:
        content = open_file.read()
    
    content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

    with open(file, 'wb') as open_file:
        open_file.write(content)
    

def mkdirif(path:str) -> int:
    '''makes a directory if it doesn't already exist
    path is the directory to create'''
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        logging.info("Created directory %s" % path)

#-------------------------------------------------        
#------------compile all the files--------------
#-------------------------------------------------


def createNozzleBlockFile(geo:NozVars, mv:MeshVars, folder:str, exportMesh:bool=False, onlyMesh:bool=False, **kwargs) -> FileGroup:
    '''gets the text for most of the files in the folder
    exportMesh is true if we want to create a mesh folder'''
    
    fg = FileGroup(folder, exportMesh=exportMesh, onlyMesh=onlyMesh, **kwargs)
    fg.geofile = geometryFile(geo)
    
    #system
    fg.setFieldsDict = compileSetFieldsDict(geo)
    if exportMesh:
        pl = createNozzle(geo) # point list
        bcs = blockCornerList(pl) # list of block corners
        blocks = []
        for bc in bcs:
            blocki = Block()
            blocki.vertices = blockPts(pl, bc)
            msz = mv.meshsize # mesh size
            blocki.meshi = [round(geo.bw/msz), round(geo.bh/msz), round(geo.bd/msz)]
            blocks.append(blocki) 
        bl = boundaryList(blocks) # this is a dummy boundary list for the original blockMeshDict
        fg.blockMeshDict = compileBlockMeshDict(pl, blocks, bl)
        fg.fvSchemesMesh = compileFvSchemes("")
        fg.fvSolutionMesh = compileFvSolution("")
    # if geo.dst: # RG
    #     if not 'reffolder' in kwargs:
    #         raise Exception('No reference folder given')
    #     reffolder = kwargs.get('reffolder')
    #     x = geo.ble+(geo.ncx-geo.ble)*geo.niw+8*geo.niw
    #     p = intm.posSlice(reffolder,x)
    #     intm.adjacentStl(folder,geo,p)
    
    
    #constant
    fg.dynamicMeshDict = compileDynamicMeshDict(mv)
    fg.g = compileG()
    fg.turbulenceProperties = compileTurbulenceProperties()
     
    #0
    br = realBoundaries(geo, exportMesh, **kwargs) # this is a real boundary list for the imported stl boundaries
    fg.meshes = br # keep the real boundaries for exporting
    fg.alphainkorig = compileAlphaOrig(br)
    fg.U = compileU(br)
    fg.prgh = compileP(br)

    if exportMesh:
        fg.snappyHexMeshDict = compileSnappyHexMeshDict(br, mv)
        fg.surfaceFeatureExtractDict = compileSurfaceFeatureExtractDict(br)
        fg.surfaceFeaturesDict = compileSurfaceFeaturesDict(br)
        fg.meshQualityDict = compileMeshQualityDict()
        fg.cellLevel = compileCellLevel(br)
        fg.pointLevel = compilePointLevel(br)
        
    #plot
    if exportMesh:
        fg.geo = geo
        fg.pl = pl
        fg.blocks = blocks
        fg.br = br
    return fg


def allButTransport(ii:Union[str, float], topFolder:str,\
                    exportMesh:bool=False, \
                    onlyMesh:bool=False, \
                    folderBase:str="nb", \
                    startTime:float=0, \
                    endTime:float=2.5, \
                    dt:float=0.001, \
                    writeDt:float=0.1, \
                    solver:float="interFoam",\
                    **kwargs) -> FileGroup:
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
    if not os.path.exists(topFolder):
        raise Exception(f'Folder {topFolder} does not exist')
    
    if onlyMesh:
        folder = topFolder
    else:
        if type(ii) is str:
            folder = os.path.join(topFolder, ii)
        else:
            folder = os.path.join(topFolder, folderBase+str(ii))
            
    geo = NozVars(**kwargs)
    mv = MeshVars()
    out = createNozzleBlockFile(geo, mv, folder, exportMesh, onlyMesh, **kwargs)

    cdv, fvv = solverObjects(startTime, endTime, dt, writeDt, solver)
    out = compileSolverFiles(cdv, fvv, out)
    return out



class Fluid:
    '''OpenFOAM needs to use kinematic units (see: kinematic vs. dynamic viscosity). If viscosities are given in dynamic units (e.g. Pa*s for viscosity, Pa for stress), then they need to be normalized by the density. We indicate whether units are kinematic or dynamic using 'units'. Units that mention Pa or dynamic will be considered dynamic. If no units are given, kinematic units are assumed.'''
    def __init__(self, units:str='kinematic', **kwargs):
        if 'rho' in kwargs:
            self.rho = kwargs['rho']
            if self.rho<500:
                logging.warning('Density is very low. Density should be given in kg/m^3')
        else:
            self.rho = 1000
            
        if units in ['Pa', 'dynamic', 'Pa*s', 'Pas']:
            div = self.rho
        else:
            div = 1

        if 'tau0' in kwargs and 'n' in kwargs and 'k' in kwargs and 'nu0' in kwargs:
            self.model='HerschelBulkley'
            self.tau0 = kwargs['tau0']/div
            self.n = kwargs['n']
            self.k = kwargs['k']/div
            self.nu0 = kwargs['nu0']/div
        elif 'nu' in kwargs:
            self.model='Newtonian'
            self.nu=kwargs['nu']/div
        else:
            raise ValueError('Invalid inputs. Required input for Newtonian is (nu). Required inputs for Herschel-Bulkley are (tau0, n, k, nu0).')
            
        if 'label' in kwargs:
            self.label = kwargs['label']
        else:
            self.label = ''

    def transportGroup(self, name:str):
        if self.model=='Newtonian':
            return transportGroupNewt(name, self.nu, self.rho)
        elif self.model=='HerschelBulkley':
            return transportGroupHB(name, self.nu0, self.tau0, self.k, self.n, self.rho)

def labels(ink:Fluid, sup:Fluid) -> str:
    '''get a file that stores the ink and support labels'''
    return f'ink,{ink.label}\nsup,{sup.label}'

def genericExport(ii:Union[int,str], sup:Fluid, ink:Fluid, sigma:float, topFolder:str, exportMesh:bool=False, **kwargs) -> None: # RG
    ''' Export a folder, given a support fluid, ink fluid, and surface tension. 
        ii is for the folder label. If you want it to be labeled nb#, input a number. Otherwise, input a string.
        sup is a fluid object that holds info about the support transport properties
        ink is a fluid object that holds info about the ink transport properties
        sigma is in J/m^2 (e.g. 0.04)
        topFolder is the folder to save this new folder in
        exportMesh true to export geometry folders into this folder'''
    out = allButTransport(ii, topFolder, exportMesh=exportMesh, **kwargs)
    out.labels = labels(ink, sup)
    out.transportProperties = compileTransportProperties(ink.transportGroup('ink'), sup.transportGroup('sup'), sigma)
    out.exportAllFiles() 

def genericMesh(parentFolder:str, **kwargs) -> FileGroup:
    '''This generates a folder with mesh and geometry files, but nothing else important. If you want to customize the nozzle, input keyword variables. e.g. genericMesh(myFolder, bathWidth=16)'''
    out = allButTransport('temp', parentFolder, exportMesh=True, onlyMesh=True, **kwargs) # generate mesh files
    out.exportAllFiles() # export all of the mesh files
    out.makePlot() # make a plot of the boundaries
    return out


#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
#----------------------------------------------------------------------------------  
 ############## ARCHIVE
    
    
    
    

# def classIterate(obj, f):
#     '''iterate over all items stored in an object. Inputs: obj, f
#     obj is an object
#     f is a function with 2 inputs: the attribute name and the value'''
#      for attr, value in obj.__dict__.items():
#         f(attr, value)



        
# # calculate the grading for the mesh within a block
#     # W is the total width of the block, 
#     # R is the ratio of largest to smallest cell
#     # t is the smallest cell size        
# def gradingcalc(self, W:float, R:float, t:float) -> float:
#     er = (W-t)/(W-R*t)
#     n = round(math.log(1-W*(1-er)/t, er))
#     return n



 ####### for newtonian ink and support
# def newtNewtExport(ii, ivisc, svisc, sigma, topFolder, exportMesh):
#     out = allButTransport(ii, exportMesh, topFolder)
#     ink = transportGroupNewt("ink", ivisc, 1000)
#     sup = transportGroupNewt("sup", svisc, 1000)
#     out.transportProperties = compileTransportProperties(ink, sup, sigma)
#     out.exportAllFiles(exportMesh)

# ####### for Herschel Bulkley support and newtonian ink
# def HBnewtexport(ii, tau0, k, n, ivisc, nu0, sigma, topFolder, exportMesh):
#     out = allButTransport(ii, exportMesh, topFolder)
#     ink = transportGroupNewt("ink", ivisc, 1000)
#     sup = transportGroupHB("sup", nu0, tau0, k, n, 1000)
#     out.transportProperties = compileTransportProperties(ink, sup, sigma)
#     out.exportAllFiles(exportMesh)
    
# ####### for newtonian support and HerschelBulkley ink
# def newtHBexport(ii, tau0, k, n, svisc, nu0, sigma, topFolder, exportMesh):
#     out = allButTransport(ii, exportMesh, topFolder)
#     sup = transportGroupNewt("sup", svisc, 1000)
#     ink = transportGroupHB("ink", nu0, tau0, k, n, 1000)
#     out.transportProperties = compileTransportProperties(ink, sup, sigma)
#     out.exportAllFiles(exportMesh)
    
# ####### for newtonian support and HerschelBulkley ink
# def HBHBexport(ii, tau0, k, n, nu0sup, nu0ink, sigma, topFolder, exportMesh):
#     out = allButTransport(ii, exportMesh, topFolder)
#     sup = transportGroupHB("sup", nu0sup, tau0, k, n, 1000)
#     ink = transportGroupHB("ink", nu0ink, tau0, k, n, 1000)
#     out.transportProperties = compileTransportProperties(ink, sup, sigma)
#     out.exportAllFiles(exportMesh)
    
# ####### for newtonian support and HerschelBulkley ink
# def HBHByielded(ii, tau0sup, tau0ink, ksup, kink, nsup, nink, nu0sup, nu0ink, sigma, topFolder, exportMesh):
#     out = allButTransport(ii, exportMesh, topFolder)
#     sup = transportGroupHB("sup", nu0sup, tau0sup, ksup, nsup, 1000)
#     ink = transportGroupHB("ink", nu0ink, tau0ink, kink, nink, 1000)
#     out.transportProperties = compileTransportProperties(ink, sup, sigma)
#     out.exportAllFiles(exportMesh)