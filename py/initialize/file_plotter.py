#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import os,sys
import logging
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
import matplotlib.patches as mpatches

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from block import Block
from boundary_input import BoundaryInput
from noz_vars import NozVars


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  


class filePlotter:
    
    def __init__(self, geo, pl, blocks, br, elev:int=10, azim:int=190, fs:int=10,  **kwargs):
        self.geo = geo
        self.pl = pl
        self.blocks = blocks
        self.br = br
        self.elev = elev
        self.azim = azim
        self.fs = fs
        self.makePlot(**kwargs)

    def makePlot(self, showMesh:bool=False
                 , vertexNumbers:bool=False, blockNumbers:bool=False
                 , figsize:float=10, **kwargs):
        '''plot the requested geometry
        geo is a nozzleVars object
        pl is a point list used for blockMesh
        blocks is a list of blocks used for blockMesh
        br is real boundaries used for generating stls and setting boundary conditions'''
        self.fig = plt.figure(figsize=[figsize,figsize*1.5])
        plt.rc('font',family='Arial')
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.view_init(elev=self.elev, azim=self.azim)
        self.ax.set_clip_on(False)
        maxdim = max([self.geo.bw, self.geo.bh, self.geo.bd])
        self.ax.set_xlim3d(self.geo.ble, self.geo.ble+maxdim)
        self.ax.set_ylim3d(self.geo.bfr, self.geo.bfr+maxdim)
        self.ax.set_zlim3d(self.geo.bbo, self.geo.bbo+maxdim)
        self.labelX()
        self.labelY()
        self.labelZ()
        self.ax.scatter(self.pl[:,0], self.pl[:,1], self.pl[:,2], s=self.fs, c='k')
        if blockNumbers:
            for p in self.pl:
                self.ax.text(p[0], p[1], p[2], "{:.0f}".format(p[3]), color='k', fontsize=self.fs)
        if vertexNumbers:
            for i, block in enumerate(self.blocks):
                vs = block.vertices
                self.ax.text(np.mean(vs[:,0]), np.mean(vs[:,1]), np.mean(vs[:,2]), i, color='r', fontsize=self.fs)
        if len(self.br)==4:
            clist = ['tomato', 'navajowhite', 'mediumaquamarine', 'deepskyblue', 'white']
        else:
            clist = ['tomato', 'navajowhite', '#d1b721', 'gray', 'mediumaquamarine', 'deepskyblue']
        plist = []
        for i, bi in enumerate(self.br):
            # each boundary bi in list bl
            col = clist.pop(0) # pop a color from the front of the list
            plist.append(mpatches.Patch(color=col, label=bi.label))
            for m in bi.meshi['vectors']:
                p3c = Poly3DCollection([list(zip(m[:,0],m[:,1],m[:,2]))], alpha=0.35)
                p3c.set_facecolor(col)
                if showMesh:
                    p3c.set_edgecolor('gray')
                else:
                    p3c.set_edgecolor(None)
                self.ax.add_collection3d(p3c)      
            self.ax.text(self.geo.ble+maxdim, self.geo.bfr+maxdim, self.geo.bbo+i/len(self.br)*self.geo.bd, bi.label, color=col, fontfamily='Arial', fontsize=self.fs)   # legend
        self.ax.grid(False)
        self.ax.axis(False)
        self.ax.set_facecolor('White')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_zticks([])
        self.ax.get_proj = lambda: np.dot(Axes3D.get_proj(self.ax), np.diag([1,1,1.3, 1]))
  
        
    def labelX(self):
        x = [self.geo.ble, self.geo.ble+self.geo.bw]
        if self.elev<0:
            # label top
            z = self.geo.bbo+self.geo.bd
            va = 'bottom'
        else:
            z = self.geo.bbo
            va = 'top'
        if self.azim>0 and self.azim<180:
            # label back
            y = self.geo.bfr+self.geo.bh + 0.5
            ha = 'right'
        else:
            y = self.geo.bfr - 0.5
            ha = 'left'
        self.ax.quiver(x[0], y, z, self.geo.bw, 0,0, color='Black', arrow_length_ratio=0.1)
        self.ax.text(np.mean(x), y, z, f'{self.geo.bw} mm', color='Black', fontfamily='Arial', fontsize=self.fs, verticalalignment=va, horizontalalignment=ha)
        self.ax.text(x[1], y, z, 'x', color='Black', fontfamily='Arial', fontsize=self.fs, verticalalignment=va, horizontalalignment=ha)
        
    def labelY(self):
        y = [self.geo.bfr, self.geo.bfr+self.geo.bh]
        if self.elev<0:
            # label top
            z = self.geo.bbo+self.geo.bd
            va = 'bottom'
        else:
            z = self.geo.bbo
            va = 'top'
        if self.azim>90 and self.azim<270:
            # label back
            x = self.geo.ble - 0.5
            ha = 'left'
        else:
            x = self.geo.ble+self.geo.bw + 0.5
            ha = 'right'
        self.ax.quiver(x, y[0], z, 0, self.geo.bh, 0, color='Black', arrow_length_ratio=0.1)
        self.ax.text(x, np.mean(y), z, f'{self.geo.bh} mm', color='Black', fontfamily='Arial', fontsize=self.fs, verticalalignment=va, horizontalalignment=ha)
        self.ax.text(x, y[1], z, 'y', color='Black', fontfamily='Arial', fontsize=self.fs, verticalalignment=va, horizontalalignment=ha)
        
    def labelZ(self):
        '''put a label on the z axis'''
        z = [self.geo.bbo, self.geo.bbo+self.geo.bd]
        x0 = self.geo.ble - 0.5
        xf = self.geo.ble+self.geo.bw + 0.5
        y0 = self.geo.bfr - 0.5
        yf = self.geo.bfr+self.geo.bh + 0.5
        ha = 'right'
        if self.azim>=0 and self.azim<90:
            y = y0
            x = xf
        elif self.azim>=90 and self.azim<180:
            y = yf
            x = xf
        elif self.azim>=180 and self.azim<270:
            y = yf
            x = x0
        else:
            y = y0
            x = x0
        self.ax.quiver(x, y, z[0], 0, 0, self.geo.bd, color='Black', arrow_length_ratio=0.1)
        self.ax.text(x, y, np.mean(z), f'{self.geo.bd} mm', color='Black', fontfamily='Arial', fontsize=self.fs, horizontalalignment=ha)
        self.ax.text(x, y, z[1], 'z', color='Black', fontfamily='Arial', fontsize=self.fs, horizontalalignment=ha)