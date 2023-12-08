#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import os,sys
import pandas as pd

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
sys.path.append(currentdir)
from fluid import Fluid
from tools.config import cfg
from creator import genericExport, genericMesh
from file.plainIm import plainIm

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class adjacentFluid:
    '''a class for describing fluids for printing adjacent lines'''
    
    def __init__(self, model:str):
        if model=='HerschelBulkley':
            self.rho = '1000'; 
            self.m = 'HerschelBulkley'; 
            self.nu0 = '10'
            self.tau0 = '0.01'
            self.k = '0.00375'
            self.n = 0.45
            rho = float(self.rho)
            tau0 = float(self.tau0)*10**3
            k = float(self.k)*10**3
            n = float(self.n)
            nu0 = float(self.nu0)*10**3
            self.obj = Fluid(units = "Pa", tau0=tau0, k=k, n=n, nu0=nu0, rho=rho)
        elif model=='Newtonian':
            self.rho = '1000'
            self.m = 'Newtonian'
            self.nu = '0.001'
            rho = float(self.rho)
            nu = float(self.nu)*10**3
            self.obj = Fluid(units = "Pa", label = '', nu = nu , rho = rho)
            
            
class adjacentCreator:
    '''for creating simulations for adjacent paper
    orientation must be in ['y', 'z']
    vink should be in [0,10]
    sigma should be in ['0.04', '0']
    inkModel and supportModel should be in ['HerschelBulkley', 'Newtonian']
    dist is the spacing between filaments in nozzle inner diameters
    ii is the simulation number
    '''
    
    def __init__(self, orientation:str, vink:float, sigma:str, dist:float, ii:int, inkModel:str='Newtonian', supModel:str='Newtonian'):
        self.fluids = {}
        self.fixSigma(sigma)
        self.orientation = orientation
        self.vink = vink
        self.ink = adjacentFluid(inkModel)
        self.sup = adjacentFluid(supModel)
        self.dist = dist
        self.iiname = f'aj{ii}'
        self.findModel()
        self.topfolder = os.path.join(cfg.path.c, 'adjacent') # this is the parent folder that holds everything
        
    def fixSigma(self, s):
        sigma = float(s)
        if sigma==0:
            self.sigma = '0'
            return
        if sigma>0.1:
            # wrong units
            sigma = sigma/1000
        self.sigma = str(sigma)
        
    def export(self):
        genericExport(self.iiname, self.sup.obj, self.ink.obj
                         , self.sigma, self.topfolder, reffolder=self.rf
                         , adjacent=self.orientation, distance=self.dist, vink=self.vink
                                 , exportMesh=True, slurmFolder=f'/working/lmf1/adjacent/{self.iiname}')
        
    def exportMesh(self, dist:float=0.875, adj:str='y'):
        fc = genericMesh(cfg.path.c, reffolder=self.rf, adjacent=self.orientation
                         , distance=self.dist, vink=self.vink) # generate mesh files
        fn = os.path.join(cfg.path.fig, 'adjacent', 'nozzleDiagramAdj.svg')
        fc.fg.plot.savefig(fn) # export the diagram of the nozzle
        logging.info(f'Exported {fn}')
        
        
    def extractSim(self, **kwargs) -> str: # RG
        '''Extract simulation number from the general legend based on ink and support rheology
        outputs a simulation string'''
        file = os.path.join(cfg.path.server, 'viscositysweep\legend_general.csv')
        d,u = plainIm(file, ic=0)
        d = d[(d['ink_transportmodel']==self.ink.m) & (d['sup_transportmodel']==self.sup.m)]
        d = d[(d['ink_rho']==float(self.ink.rho)) & (d['sup_rho']==float(self.sup.rho)) & (d['sigma']==float(self.sigma))]

        if self.ink.m=='Newtonian':
            d = d[(d['ink_nu']==float(self.ink.nu))]
        elif self.ink.m=='HerschelBulkley':
            d = d[(d['ink_tau0']==float(self.ink.tau0)) & (d['ink_k']==float(self.ink.k)) & (d['ink_n']==float(self.ink.n))]
        else:
            raise ValueError(f'Unexpected ink model {self.ink.m}')
        if self.sup.m=='Newtonian':
            d = d[(d['sup_nu']==float(self.sup.nu))]    
        elif self.sup.m=='HerschelBulkley':
            d = d[(d['sup_tau0']==float(self.sup.tau0)) & (d['sup_k']==float(self.sup.k)) & (d['sup_n']==float(self.sup.n))]
        else:
            raise ValueError(f'Unexpected sup model {self.sup.m}')

        if d.empty:
            raise Exception('No simulations match that rheology')
        sim = d['folder'].values[0]
        return sim
        
    def findModel(self):
        '''find the reference simulation path'''
        ref = self.extractSim()
        if self.ink.m=='Newtonian':
            if self.sup.m=='Newtonian':
                folder = 'newtnewtsweep'
            elif self.sup.m=='HerschelBulkley':
                folder = 'HBnewtsweep'
        elif self.ink.m=='HerschelBulkley':
            if self.sup.m=='Newtonian':
                folder = 'newtHBsweep'
            elif self.sup.m=='HerschelBulkley':
                folder = 'HBHBsweep'
            else:
                raise ValueError(f'Unexpected support model {self.sup.m}')
        else:
            raise ValueError(f'Unexpected ink model {self.ink.m}')
        
        self.rf = os.path.join(cfg.path.server, 'viscositysweep', folder, ref)       
