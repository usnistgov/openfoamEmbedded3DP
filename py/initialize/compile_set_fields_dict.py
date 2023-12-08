#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import os,sys
import numpy as np

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from noz_vars import NozVars
from dict_list import DictList
from initialize_tools import OpenFOAMFile

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class compileSetFieldsDict(OpenFOAMFile):
    
    
    def __init__(self, geo:NozVars): # fills the nozzle with ink at t=0 RG
        '''compile setFieldsDict'''
        s = self.header("dictionary", "setFieldsDict")
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
            ncx = geo.scale('ncx')
            ncy = geo.scale('ncy')
            d1 = round(geo.scale((geo.nbo+(geo.bto-geo.nbo)*x/steps)),8)
            d2 = round(geo.scale(geo.bto-(geo.bto-geo.nbo)*(steps-1-x)/steps),8)
            rtop = geo.scale(geo.niw/2+(x+1)*geo.nl/steps*np.tan(np.deg2rad(geo.na)))

            c2c.proplist.append([f'p1 ({ncx} {ncy} {d1})'])
            c2c.proplist.append([f'p2 ({ncx} {ncy} {d2})'])
                        # top and bottom points of central cylinder axis
            c2c.proplist.append([f'radius {rtop}'])
                        #radius of cylinder
            c2c.proplist.append(DictList("fieldValues", 1, ["volScalarFieldValue alpha.ink 1"]))
                        # value of alpha inside the cylinder is 1 because it is full of ink
            r.proplist.append(c2c)

        s = s + r.prnt(0)
        s = s + self.closeLine()
        self.s = s