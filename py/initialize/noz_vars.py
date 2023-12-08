#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
import os,sys
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import numpy as np

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentdir)
from initialize_tools import scale

# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  



class NozVars:
    '''this is where the geometry of the nozzle is defined
    Allowed input variables:
            bathWidth: (default=16) bath width in nozzle inner diameters (x)
            bathHeight: (default=7) bath height in nozzle inner diameters (y)
            bathDepth: (default=7) bath depth in nozzle inner diameters (z)
            frontWidth: (default=4) front of nozzle bath width in nozzle inner diameters
            vink: (default=10) ink extrusion speed in mm/s
            vbath: (default=10) bath translation speed in mm/s
            npts: (default=50) number of points in the circle used to define the nozzle
            nozzleInnerWidth: (default=0.603) inner diameter of the nozzle in mm
            nozzleThickness: (default=0.152) nozzle wall thickness in mm
            nozzleAngle: (default=0) nozzle angle in degrees
            horizontal: (default=False) whether to have the nozzle horizontal
            adjacent: (default=None) adjacent filament orientation
            distance: (default=-1) offset of adjacent filament in nozzle inner diameters'''
    
    def __init__(self
                 , bathWidth:float=16
                 , bathHeight:float=7
                 , bathDepth:float=7
                 , frontWidth:float=4
                 , vink:float=10
                 , vbath:float=10
                 , npts:int=50
                 , nozzleInnerWidth:float=0.603
                 , nozzleThickness:float=0.152
                 , nozzleAngle:float=0
                 , horizontal:bool=False
                 , adjacent:str='None'
                 , distance:float=0
                 , reffolder:str = None
                 , **kwargs):
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
        self.cor = reffolder
        
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
        
        self.inptsb = self.circlePoints(bottomRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the bottom of the nozzle
        self.outptsb = self.circlePoints(bottomRadius + self.nt, npts)+[self.ncx, self.ncy] 
            # outer radius points at the bottom
        self.inptst = self.circlePoints(topRadius, npts)+[self.ncx, self.ncy] 
            # inner radius points at the top of the nozzle (the inlet)
        self.outptst = self.circlePoints(topRadius + self.nt, npts)+[self.ncx, self.ncy]
            # outer radius points at the top of the nozzle

        self.bv = self.scale(vbath) # bath velocity
        self.iv = self.scale(vink*bottomRadius**2/topRadius**2) # initial ink velocity RG
        
    def scale(self, val:Union[float, str]) -> float:
        scl = scale()
        if type(val) is str:
            return scl*getattr(self, val)
        return scl*val
    
    def circlePoints(self, r:float, npts:int) -> list:
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
    
    def createNozzle(self) -> np.array:
        '''get list of vertices in the block mesh'''
        xlist = np.array([self.ble-0.01,  self.bri+0.01])
        ylist = np.array([self.bfr-0.01, self.bba+0.01])
        zlist = np.array([self.bbo-0.01, self.bto+0.01])
        pts = np.zeros([xlist.size*ylist.size*zlist.size,4])
        ptsi = 0
        for z in zlist:
            for y in ylist:
                for x in xlist:
                    pts[ptsi, :] = [x, y, z, ptsi]
                    ptsi = ptsi+1
        self.pl = pts
        return pts
    
    def blockCornerList(self) -> list:
        '''select the points that are the bottom corner of the blockCornerList'''
        if not hasattr(self, 'pl'):
            self.createNozzle()
        xmax = max(self.pl[:, 0])
        ymax = max(self.pl[:, 1])
        zmax = max(self.pl[:, 2])
        corners = self.pl[(self.pl[:,0]<xmax) & (self.pl[:,1]<ymax) & (self.pl[:,2]<zmax), :]
        return corners
    