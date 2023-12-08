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
from noz_vars import NozVars

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 

class geometryFile:
    '''geometryFile gets a csv string of all of the geometry variables we care about'''
    
    def __init__(self, geo:NozVars):
    
        l = [['nozzle_inner_width', geo.niw, 'mm'],\
             ['nozzle_thickness', geo.nt, 'mm'], \
             ['bath_width', geo.bw, 'mm'], \
             ['bath_height', geo.bh, 'mm'], \
             ['bath_depth', geo.bd, 'mm'], \
             ['nozzle_length', geo.nl, 'mm'],\
             ['bath_left_coord', geo.ble, 'mm'], \
             ['bath_right_coord', geo.bri, 'mm'],\
             ['bath_front_coord', geo.bfr, 'mm'], \
             ['bath_back_coord', geo.bba, 'mm'],\
             ['bath_bottom_coord', geo.bbo, 'mm'], \
             ['bath_top_coord', geo.bto, 'mm'], \
             ['nozzle_bottom_coord', geo.nbo, 'mm'],\
             ['nozzle_center_x_coord', geo.ncx, 'mm'],\
             ['nozzle_center_y_coord', geo.ncy, 'mm'], \
             ['nozzle_angle', geo.na, 'degrees'], \
             ['horizontal', geo.hor, ''], \
             ['adjacent_filament_orientation', geo.adj, ''], \
             ['adjacent_filament_offset', geo.dst, 'mm'], \
             ['corresponding_simulation', geo.cor, ''], \
             ['bath_velocity', geo.bv, 'm/s'], \
             ['ink_velocity', geo.iv, 'm/s']] # RG
        s = ""
        for li in l:
            s = s + f'{li[0]}, {li[1]}, {li[2]}\n'
        self.s = s