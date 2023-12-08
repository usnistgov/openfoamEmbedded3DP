#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import os,sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from transport_group import transportGroupNewt, transportGroupHB

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

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
        
        
