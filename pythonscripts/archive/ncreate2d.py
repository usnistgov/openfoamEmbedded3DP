import math
import numpy as np
import matplotlib.pyplot as plt
import statistics as st
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

#-------------------------------------------------
#---------------GLOBAL VARIABLES------------------
#-------------------------------------------------

DATAFOLDER=""
SCALE = 0.001 # scale all units by this amount * m

#-------------------------------------------------
#--------------------CLASSES----------------------
#-------------------------------------------------

class NozVars:
    # ALL VALUES IN MM 
    
    def gradingcalc(self, W, R, t):
    # calculate the grading for the mesh within a block
    # W is the total width of the block, R is the ratio of largest to smallest
    # cell, and t is the smallest cell size
        er = (W-t)/(W-R*t)
        n = round(math.log(1-W*(1-er)/t, er))
        return n
              
    def __init__(self, bw, t, v, gnx, gny):
        self.niw = 0.603 # nozzle inner width
        self.nt = 0.152 # nozzle thickness
        self.bw = bw # bath width
        self.bh = bw # bath height
        self.nl = self.bh/2
        self.t = t # under nozzle mesh size
        self.bv = v*SCALE # bath velocity
        self.iv = v*SCALE # ink velocity
        self.gnx = gnx # grading to nozzle x
        self.gny = gny # grading to nozzle y
        self.mxinnoz = round(self.niw/t) # mesh cells in nozzle x
        self.mxnt = round(self.nt/t) # mesh cells under nozzle walls x
        self.mxout = self.gradingcalc((self.bw - 2*self.nt - self.niw)/2, gnx, t)
        self.my = self.gradingcalc(self.nl, gny, t)
        
class Block:
    def __init__(self):
        self.vertices = []
        self.mesh = [1,1,1]
        self.grading = [1,1,1]
        
class BoundaryInput:
    def __init__(self, labelin, typin):
        self.label = labelin # name of the boundary e.g. outerWall
        self.typ = typin # name of boundary type e.g. patch
        self.flist = [] # list of faces of list of point indices

    
class ControlDictVars:
    def __init__(self, start, end, dt, writeint):        
        self.startTime = str(start)
        self.endTime = str(end)
        self.deltaT = str(dt)
        self.writeInterval = str(writeint)
        application = "interFoam"
        self.startFrom = "startTime"
        self.stopAt = "endTime"
        self.writeControl = "adjustable"
        self.purgeWrite = "0"
        self.writeFormat = "ascii"
        self.writePrecision = "6"        
        self.writeCompression = "off"
        self.timeFormat = "general"
        self.timePrecision = "6"
        self.runTimeModifiable = "no"
        self.adjustTimeStep = "no"
        self.maxCo = "1"
        self.maxAlphaCo = "1"
        self.maxDeltaT = "1"
        
    
class FileGroup:
    blockMeshDict = ""
    controlDict = ""
    setFieldsDict = ""
    g = ""
    transportProperties = ""
    turbulenceProperties = ""
    alphainkorig = ""
    prgh = ""
    U = ""
    allclean = ""
    allrun = ""
    plot = plt.figure
    
    
    
#-------------------------------------------------
#--------------------FUNCTIONS--------------------
#-------------------------------------------------

#-------------------------------------------------
#----------------CREATING GEOMETRIES--------------
def createnozzle(nv):
    # get list of vertices in the block mesh
    # nv is a NozVars object
    nl1 = nv.bw/2 - nv.niw/2 - nv.nt
    nl2 = nv.bw/2 - nv.niw/2
    nr1 = nv.bw/2 + nv.niw/2
    nr2 = nv.bw/2 + nv.niw/2 + nv.nt
    xlist = np.array([0, nl1, nl2, nr1, nr2, nv.bw])
    ylist = np.array([0, nv.nl, nv.bh])
    zlist = np.array([0, 1])
    pts = np.zeros([xlist.size*ylist.size*zlist.size,4])
    ptsi = 0
    for z in zlist:
        for y in ylist:
            for x in xlist:
                pts[ptsi, :] = [x, y, z, ptsi]
                ptsi = ptsi+1
    return pts

def blockcornerlist(pl):
    # select the points that are the bottom corner of the blockcornerlist
    # pl is a point list
    xmax = max(pl[:, 0])
    ymax = max(pl[:, 1])
    zmax = max(pl[:, 2])
    corners = pl[(pl[:,0]<xmax) & (pl[:,1]<ymax) & (pl[:,2]<zmax), :]
    return corners
 

def ptfrompts(pts, x, y, z):
    # select the point from a list
    pt = pts[(pts[:,0]==x) & (pts[:,1]==y) & (pts[:,2]==z), :]
    return pt[0]
    
def blockpts(pl, corner):
    # select the points in this block based on the first corner, and put 
    # the points in the order that openFOAM likes
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
#---------------TRANSLATING TO openFOAM------------

def header(cl, obj):
    # header for openfoam files
    # obj is the object name (string)
    # cl is the class name (string)
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
        +"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //;\n"
        )
    return s

CLOSELINE = "// ************************************************************************* //";

def vec2cpp(v):
    # convert 1D array v to openfoam string
    s = "("
    for vi in v:
        s = s + str(vi) + " "
    s = s + ")"
    return s
    
#-------------------------------------------------   
#-----------------SYSTEM--------------------------

#-------------------------------------------------
#---------------controlDict-----------------------
def compileControlDict(cv):
    s = header("dictionary", "controlDict")
    for attr, value in cv.__dict__.items():
        s = s + attr + "\t" + value + ";\n\n"
    s = s + "#sinclude\t\"sampling\"\n\n"
    s = s + CLOSELINE
    return s
    
def compileCDVars(starttime, endtime, dt, writedt):
    cv = ControlDictVars(starttime, endtime, dt, writedt)
    return compileControlDict(cv)

#-------------------------------------------------
#---------------blockMeshDict---------------------
def vertices2txt(pl):
    # convert list of vertices to openfoam string
    s = "vertices\n(\n"
    for i in range(pl[:,0].size):
        s = s + "\t" + vec2cpp(pl[i, 0:3]) + "\n"
    s = s + ");"
    return s

def block2txt(block):
    # convert block to openfoam string
    # block is a Block object
    s = "\n\thex " + vec2cpp(block.vertices[:, 3].astype(int)) + " " + vec2cpp(block.mesh)
    s = s + " simpleGrading " + vec2cpp(block.grading)
    return s
    
def blocks2txt(blocks):
    # convert list of blocks to openfoam string
    s = "blocks\n("
    for b in blocks:
        s = s + block2txt(b)
    s = s + "\n);"
    return s
    
EDGES = "edges\n(\n);"

def face2str(face):
    # face is a list of points in a face
    s = "\t\t\t" + vec2cpp(face[:,3].astype(int)) + "\n"
    return s
    
def faceselector(block, st):
    # block is a block object, st is a string indicating which face to users
    # gets a list of vertices for a list
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
    
    
def boundary(bi):
    # establish boundary conditions in openfoam
    # bi is a BoundaryInput objet
    s = "\t" + bi.label + "\n\t{\n\t\ttype " + bi.typ + ";\n\t\tfaces\n\t\t(\n";
    for f in bi.flist:
        s = s + face2str(f)
    s = s + "\t\t);\n\t}\n";
    return s

def boundarylist(blocks):
    bathflow = BoundaryInput("bathFlow", "patch")
    atmos = BoundaryInput("atmosphere", "patch")
    inkflow = BoundaryInput("inkFlow", "patch")
    walls = BoundaryInput("fixedWalls", "wall")
    fb = BoundaryInput("frontAndBack", "empty")
    corners = np.zeros([len(blocks), 3])
    for i in range(len(blocks)):
        corners[i] = blocks[i].vertices[0, 0:3] # corners is a list of block 1st corners
    xlist = np.sort(np.unique(corners[:,0])) # list of corner indices
    ylist = np.sort(np.unique(corners[:,1]))
    for i in range(len(blocks)):
        b = blocks[i]
        bv = b.vertices[0]
        if bv[0]==xlist[0]: 
            bathflow.flist.append(faceselector(blocks[i], "x-")) # left wall
        elif bv[0]==xlist[-1]: 
            atmos.flist.append(faceselector(blocks[i], "x+")) # right wall
            
        if bv[1]==ylist[0]: 
            bathflow.flist.append(faceselector(blocks[i], "y-")) # bottom wall
            if (bv[0]==xlist[1] or bv[0]==xlist[3]): 
                walls.flist.append(faceselector(blocks[i], "y+")) # under nozzle wall
        elif bv[0]==xlist[2]: 
            inkflow.flist.append(faceselector(blocks[i], "y+")) # nozzle inlet
        else:
            atmos.flist.append(faceselector(blocks[i], "y+")) # bath top
            
        if bv[1]==ylist[1] and (bv[0]==xlist[2] or bv[0]==xlist[4]): 
            walls.flist.append(faceselector(blocks[i], "x-")) # right of nozzle wall
        if bv[1]==ylist[1] and (bv[0]==xlist[0] or bv[0]==xlist[2]): 
            walls.flist.append(faceselector(blocks[i], "x+")) # left of nozzle wall
            
        for st in ["z-", "z+"]:
            fb.flist.append(faceselector(blocks[i], st)) # front and back faces
    return [bathflow, atmos, inkflow, walls, fb]
    
def boundarycpp(bl):
    # get the whole cpp section for list of boundaries
    s = "boundary\n(\n"    
    for b in bl:
        s = s + boundary(b) # write boundaries to file
    return s

MERGEPATCHPAIRS = "mergePatchPairs\n(\n);\n" + CLOSELINE;

def compileBlockMeshDict(pl, blocks, bl):
    s = header("dictionary", "blockMeshDict") + "\n\n"
    s = s + "scale " + str(SCALE) + ";\n\n"
    s = s + vertices2txt(pl) + "\n\n"
    s = s + blocks2txt(blocks) + "\n\n"
    s = s + EDGES + "\n\n"
    s = s + boundarycpp(bl) + "\n\n"
    s = s + MERGEPATCHPAIRS + "\n\n"
    return s

#-------------------------------------------------
#------------------SETFIELDS----------------------
def compileSetFieldsDict(blocks):
    s = header("dictionary", "setFieldsDict") + "\n\n"
    s = s + "defaultFieldValues\n(\n\tvolScalarFieldValue alpha.ink 0\n);\n\n"
    s = s + "regions\n(\n\tboxToCell\n\t{\n\t\tbox "
    for i in [0,7]:
        s = s + vec2cpp(blocks[4].vertices[i, 0:3]*SCALE) + " "
    s = s + ";\n\t\tfieldValues\n\t\t(\n\t\t\tvolScalarFieldValue alpha.ink 1\n\t\t);\n\t}\n);\n\n"
    s = s + CLOSELINE
    return s
        
#-------------------------------------------------        
#---------------------SUMMARY---------------------
#-------------------------------------------------
def createNozzleBlockFile(geo):
    pl = createnozzle(geo) # point list
    bcs = blockcornerlist(pl) # list of block corners
    xl = np.unique(pl[:,0])
    yl = np.unique(pl[:,1])
    bcs = bcs[~(((bcs[:,0]==xl[1])+(bcs[:,0]==xl[3]))*(bcs[:,1]==yl[1]))];
    blocks = []
    for bc in bcs:
        blocki = Block()
        blocki.vertices = blockpts(pl, bc)
        
        if bc[0] in [xl[0], xl[-2]]:
            blocki.mesh[0] = geo.mxout
        elif bc[0] in [xl[1], xl[-3]]:
            blocki.mesh[0] = geo.mxnt
        else:
            blocki.mesh[0] = geo.mxinnoz
        blocki.mesh[1] = geo.my
        
        if bc[0]==xl[0]:
            blocki.grading[0] = 1/geo.gnx
        elif bc[0]==xl[-2]:
            blocki.grading[0] = geo.gnx
        else:
            blocki.grading[0] = 1
        if bc[1]==yl[0]:
            blocki.grading[1] = 1/geo.gny
        else:
            blocki.grading[1] = geo.gny 
        blocks.append(blocki)
    
    fg = FileGroup()
    bl = boundarylist(blocks)
    fg.blockMeshDict = compileBlockMeshDict(pl, blocks, bl)
    fg.setFieldsDict = compileSetFieldsDict(blocks)
    fg.plot = plotGeo(geo, pl, blocks, bl)
    return fg


#-------------------------------------------------        
#---------------------PLOTS-----------------------
#-------------------------------------------------
def plotGeo(geo, pl, blocks, bl):
    fig = plt.figure(figsize=[14,14])
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(0, geo.bw)
    ax.set_ylim3d(0, geo.bw)
    ax.set_zlim3d(0, geo.bw)
    ax.scatter(pl[:,0], pl[:,1], pl[:,2], s=20, c='k')
    fs = 20;
    for p in pl:
        ax.text(p[0], p[1], p[2], "{:.0f}".format(p[3]), color='k', fontsize=fs)
    for i in range(len(blocks)):
        vs = blocks[i].vertices
        ax.text(st.mean(vs[:,0]), st.mean(vs[:,1]), st.mean(vs[:,2]), i, color='r', fontsize=fs)
    clist = ['tomato', 'navajowhite', 'mediumaquamarine', 'deepskyblue', 'white']
    plist = []
    for bi in bl:
        ff = bi.flist
        col = clist.pop(0) # pop a color from the front of the list
        plist.append(mpatches.Patch(color=col, label=bi.label))
        for f in ff:
            p3c = Poly3DCollection([list(zip(f[:,0],f[:,1],f[:,2]))], alpha=0.35)
            p3c.set_facecolor(col)
            ax.add_collection3d(p3c)
    ax.legend(handles=plist, fontsize=fs)
    return fig
    
print("ncreate imported")


