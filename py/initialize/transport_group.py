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
from dict_list import DictList

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class transportGroupNewt(DictList):
    '''gets a DictList object that describes the transport properties
    title is the name of the phase, e.g. 'ink'
    nu is the kinematic viscosity in m^2/s
    rho is the density in kg/m^3'''
    
    def __init__(self, title:str
                 , nu:Union[float, str]
                 , rho:Union[float, str]):
        super().__init__(title, 0, [["transportModel", "Newtonian"], ["nu", str(nu)], ["rho", str(rho)]])

        
class transportGroupHB(DictList):
    ''''transportGroupHB gets a DictList that describes Herschel-Bulkley transport properties
    title is the name of the phase
    the HB model follows nu = min(nu0, tau0/gammadot + k*gammadot^(n-1))
    inputs to the model are nu0 in m^2/s, tau0 in m^2/s^2, k in m^2/s, n unitless
    rho is the density in kg/m^3'''
    
    def __init__(self, title:str
                 , nu0:Union[float, str]
                 , tau0:Union[float, str]
                 , k:Union[float, str]
                 , n:Union[float, str]
                 , rho:Union[float, str]):
        super().__init__(title, 0, \
                  [["transportModel", "HerschelBulkley"], \
                  DictList("HerschelBulkleyCoeffs", 0, \
                           [["nu0", "[ 0 2 -1 0 0 0 0 ] " + str(nu0)],\
                            ["tau0", "[ 0 2 -2 0 0 0 0 ] " + str(tau0)], \
                            ["k", "[ 0 2 -1 0 0 0 0 ] " + str(k)],\
                            ["n", "[ 0 0 0 0 0 0 0 ] " + str(n)]]\
                           ),\
                  ["rho", str(rho)]])