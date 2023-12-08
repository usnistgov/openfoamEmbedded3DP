#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
from stl import mesh
import os,sys
import numpy as np
import functools

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from noz_vars import NozVars
from boundary_input import BoundaryInput
from dict_list import DictList
from points.folder_points import folderPoints
from points.points_tools import xspoints, setZ, setX, sortPolar

# logging
logging.basicConfig(level=logging.INFO)

#-------------------------------------------------------------------------------------------------  


class realBoundaries:
    '''get a list of boundaries for snappyHexMesh, setFields
    geo is a NozVars object geo.bv and geo.iv are bath velocities and ink velocities, scaled to m/s
    mv is a MeshVars object
    exportMesh is true if we want to create a mesh folder'''
    
    
    def __init__(self, geo:NozVars, exportMesh:bool, reffolder:str='', **kwargs):
        self.geo = geo
        self.exportMesh = exportMesh
        self.reffolder = reffolder
        self.createBathFlow()
        self.createInkFlow()
        self.createAtmosphere()
        self.createFixedWalls()

        if self.geo.adj!='None': # RG
            self.createXSAtmosphere()
            self.createXSFlow()
            self.boundaryList = [self.bf, self.inkf, self.xf, self.xa, self.at, self.fw]
        else:
            self.boundaryList = [self.bf, self.inkf, self.at, self.fw]  
    
    #--------------------------------------------------------------------
    
    def addAdjBathHoriz(self, faces:list) -> list:
        '''horizontal displacement'''
        if self.geo.adj=='y':
            self.xspts[1] = [i-cent[1] for i in self.xspts[1]]
        else:
            self.xspts[0] = [i-cent[0] for i in self.xspts[0]]
        self.xspts3d = setZ(self.xspts, self.geo.bto)
        self.xspts = sortPolar(self.xspts3d)[0]
        faces.insert(0, self.holeInPlane(self.xspts, [self.geo.ble, self.geo.bri], [self.geo.bfr, self.geo.bba], [self.geo.bto, self.geo.bto]))
        faces.append(self.axisFace(self.walsel("x-")))
        return faces
    
    def addAdjBathVert(self, faces:list) -> list:
        '''vertical displacement'''
        self.xspts3d = setX(self.xspts, self.geo.ble)
        self.xspts = sortPolar(self.xspts3d)[0]
        faces.insert(0, self.holeInPlane(self.xspts, [self.geo.ble, self.geo.ble], [self.geo.bfr, self.geo.bba], [self.geo.bto, self.geo.bbo]))
        return faces
    
    def addAdjBath(self, faces:list) -> list:
        '''add a face for the initial filament'''
        fp = folderPoints(self.reffolder)
        p = fp.posSlice(8, xunits='nozzle_inner_width')
        self.xspts, self.cent = xspoints(p, self.geo.dst, self.geo.adj)
        if self.geo.hor:
            faces = self.addAdjBathHoriz(faces)
        else:
            faces = self.addAdjBathVert(faces)
        return faces
    
    def createBathFlow(self):
        '''create the bathFlow face'''
        self.bf = BoundaryInput("bathFlow", "")
        self.bf.alphalist = DictList(self.bf.label, 0, [["type", "fixedValue"], ["value", "uniform 0"]])
        if self.geo.hor:
            self.bf.Ulist = DictList(self.bf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(self.geo.bv) + ")"]]) # flow in -z RG
        else:
            self.bf.Ulist = DictList(self.bf.label, 0, [["type", "fixedValue"], ["value", "uniform (" + str(self.geo.bv) + " 0 0)"]])
        self.bf.plist = DictList(self.bf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
        if self.exportMesh:
            if self.geo.hor:
                faces = list(map(lambda s: self.axisFace(self.walsel(s)), ["x+", "y-", "y+"])) # side walls RG
            else:
                faces = list(map(lambda s: self.axisFace(self.walsel(s)), ["y-", "y+", "z-"]))
            if self.geo.adj!='None':
                faces = self.addAdjBath(faces)
            else:
                faces.append(self.axisFace(self.walsel("x-")))
            self.bf.meshi = self.combineMeshes(faces)
            
            
    #--------------------------------------------------------------------       

    def createInkFlow(self):
        '''create the inkFlow boundary'''
        geo = self.geo
        inkf = BoundaryInput("inkFlow", "")
        inkf.alphalist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform 1"]])
        inkf.Ulist = DictList(inkf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(geo.iv) + ")"]])
        inkf.plist = DictList(inkf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
        if self.exportMesh:
            inkf.meshi = self.arcFace(setZ(np.zeros([len(geo.inptst), 2])+[geo.ncx, geo.ncy], geo.bto), setZ(geo.inptst, geo.bto))
        self.inkf = inkf
        
        
    def createAtmosphere(self):
        '''create the atmosphere boundary'''
        geo = self.geo
        at = BoundaryInput("atmosphere", "")
        at.alphalist = DictList(at.label, 0, [["type", "inletOutlet"], ["value", "uniform 0"], ["inletValue", "uniform 0"]])
        at.Ulist = DictList(at.label, 0, [["type", "pressureInletOutletVelocity"], ["value", "uniform (0 0 0)"]])
        at.plist = DictList(at.label, 0, [["type", "totalPressure"], ["p0", "uniform 0"]])
        if self.exportMesh:
            if geo.hor:
                at.meshi = self.combineMeshes([self.holeInPlane(setZ(geo.outptst, geo.bto), [geo.ble, geo.bri], [geo.bfr, geo.bba], [geo.bto, geo.bto]), \
                               self.axisFace(self.walsel("z-"))]) # top and bottom walls RG
            else:
                at.meshi = self.combineMeshes([self.holeInPlane(setZ(geo.outptst, geo.bto), [geo.ble, geo.bri], [geo.bfr, geo.bba], [geo.bto, geo.bto]), \
                               self.axisFace(self.walsel("x+"))])
    #         at.meshi = self.combineMeshes([self.axisFace(self.walsel("z+")), \
    #                            self.axisFace(self.walsel("x+"))])
        self.at = at
        
    def createFixedWalls(self):
        '''create the fixedWalls boundary'''
        geo = self.geo
        fw = BoundaryInput("fixedWalls", "")
        fw.alphalist = DictList(fw.label, 0, [["type", "zeroGradient"]])
        fw.Ulist = DictList(fw.label, 0, [["type", "noSlip"]])  
        fw.plist = DictList(fw.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])

        if self.exportMesh:
            fw.meshi = self.combineMeshes([self.arcFace(setZ(geo.inptsb, geo.nbo), setZ(geo.outptsb, geo.nbo)),\
                                self.arcFace(setZ(geo.inptsb, geo.nbo), setZ(geo.inptst, geo.bto)),\
                                self.arcFace(setZ(geo.outptsb, geo.nbo), setZ(geo.outptst, geo.bto))])
        fw.reflev = 2
        self.fw = fw
        
    #--------------------------------------------------------------------
    
    def createXSAtmosphere(self):
        '''create an atmosphere boundary for the exit of the adj filament'''
        geo = self.geo
        xa = BoundaryInput("xsAtmosphere", "")
        xa.alphalist = DictList(xa.label, 0, [["type", "inletOutlet"], ["value", "uniform 0"], ["inletValue", "uniform 0"]])
        xa.Ulist = DictList(xa.label, 0, [["type", "pressureInletOutletVelocity"], ["value", "uniform (0 0 0)"]])
        xa.plist = DictList(xa.label, 0, [["type", "totalPressure"], ["p0", "uniform 0"]])
        if self.exportMesh:
            xspts2 = np.copy(self.xspts)
            if geo.hor:
                xspts2[:,2]-=geo.niw/3
                xa.meshi = self.combineMeshes([self.arcFace(self.xspts, xspts2), self.arcFace(setZ(np.zeros([len(self.xspts),2])+self.cent, geo.bto), self.xspts)])
            else:
                xspts2[:,0]+=geo.niw/3
                xa.meshi = self.combineMeshes([self.arcFace(self.xspts, xspts2), self.arcFace(setX(np.zeros([len(self.xspts),2])+self.cent, geo.ble), self.xspts)])
        xa.reflev = 2  # two levels of mesh refinement
        self.xa = xa
        self.xspts2 = xspts2
        
    def createXSFlow(self):
        '''create an atmosphere flow for flow of the adj filament'''
        geo = self.geo
        xf = BoundaryInput("xsFlow", "")
        xf.alphalist = DictList(xf.label, 0, [["type", "fixedValue"], ["value", "uniform 1"]])
        if geo.hor:
            xf.Ulist = DictList(xf.label, 0, [["type", "fixedValue"], ["value", "uniform (0 0 -" + str(geo.bv) + ")"]])
        else:
            xf.Ulist = DictList(xf.label, 0, [["type", "fixedValue"], ["value", "uniform (" + str(geo.bv) + " 0 0)"]])
        xf.plist = DictList(xf.label, 0, [["type", "fixedFluxPressure"], ["value", "uniform 0"]])
        if self.exportMesh:
            if geo.hor:
                xf.meshi = self.arcFace(setZ(np.zeros([len(self.xspts),2])+self.cent, geo.bto-geo.niw/3), self.xspts2)
            else:
                xf.meshi = self.arcFace(setX(np.zeros([len(self.xspts),2])+self.cent, geo.ble+geo.niw/3), self.xspts2)
        self.xf = xf
    
    #--------------------------------------------------------------------
    
        
    def axisFace(self, xyz:List[float]) -> np.array:
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
        
    def holeInPlane(self, cpts:np.array, x:List[float], y:List[float], z:List[float]) -> np.array:
        '''get a mesh array that describes a hole in a plane
        cpts is a list of circle points. cpts should be in order from 0 to 2 pi
        x is a list of 2 x values for the plane
        y is a list of 2 y values for the plane
        z is a list of 2 z values for the plane
        repeat one value twice for the plane which the hole lies on'''
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
    
    def combineMeshes(self, meshList:List[List]) -> np.array:
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
    
    def ptFace(self, pts:np.array) -> np.array:
        '''given an array of 4 points pts, construct a face. The line from point 1 to point 2 must cross through the center of the face'''
        data = np.zeros(2, dtype=mesh.Mesh.dtype)
        data['vectors'][0] = np.array([pts[0], pts[1], pts[2]])
        data['vectors'][1] = np.array([pts[1], pts[2], pts[3]])
        return data
    
    
    def arcFace(self, ina:np.array, outa:np.array) -> np.array:
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
            d = self.ptFace([ina[i], outa[i], ina[i+1], outa[i+1]])
            data['vectors'][2*i] = d['vectors'][0]
            data['vectors'][2*i+1] = d['vectors'][1]
        return data
    


    def walsel(self, st:str):
        '''select the boundaries of a wall
        geo is a NozVars object
        st is a string indicating which face to take'''
        xli = [self.geo.ble, self.geo.bri]
        yli = [self.geo.bfr, self.geo.bba]
        zli = [self.geo.bbo, self.geo.bto]
        if st == "x-":
            xli = [self.geo.ble]
        elif st == "x+":
            xli = [self.geo.bri]
        elif st == "y-": 
            yli = [self.geo.bfr]
        elif st == "y+":
            yli = [self.geo.bba]
        elif st == "z-": 
            zli = [self.geo.bbo]
        elif st == "z+": 
            zli = [self.geo.bto]
        return [xli, yli, zli]