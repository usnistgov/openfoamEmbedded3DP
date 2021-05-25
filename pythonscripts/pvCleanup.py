#!/usr/bin/env python
'''Cleaning up paraview csv files'''

import sys
import os
import csv
import logging

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import folderparser as fp

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Development"

#################################################################

def addUnits(csvfile:str):
    '''Add units to a csv file'''
    blist = []
    with open(csvfile, 'r') as b:
        canreturn = True
        csv_reader = csv.reader(b)
        try:
            line1 = next(csv_reader)
        except:
            return
        if 'interfacePoints' in csvfile or 'line' in csvfile:
            headerdict = {'alpha.ink':'alpha', 'Points:0':'x', 'Points:1':'y', 'Points:2':'z',\
                          'U:0':'vx', 'U:1':'vy', 'U:2':'vz', 'Time':'time', 'nu1':'nu_ink', 'nu2':'nu_sup',\
                          'VectorGradient:0':'shearrate0','VectorGradient:1':'shearrate1','VectorGradient:2':'shearrate2',\
                          'VectorGradient:3':'shearrate3','VectorGradient:4':'shearrate4','VectorGradient:5':'shearrate5',\
                          'VectorGradient:6':'shearrate6','VectorGradient:7':'shearrate7','VectorGradient:8':'shearrate8' }
            header = [headerdict.get(h, h) for h in line1]
            if header==line1:
                canreturn = True
            else:
                line1=header
                canreturn = False
        blist.append(line1)
        try:
            line2 = next(csv_reader)
        except:
            return
        
        if 'interfacePoints' in csvfile or 'line' in csvfile:
            unitdict = {'time':'s', 'x':'m', 'y':'m', 'z':'m', \
                        'vx':'m/s', 'vy':'m/s', 'vz':'m/s', 'alpha':'',\
                        'nu1':'m^2/s', 'nu2':'m^2/s', 'nu_ink':'m^2/s', 'nu_sup':'m^2/s',\
                        'p':'kg/(m*s^2)', 'p_rgh':'kg/(m*s^2)', 'rAU':'m^3*s/kg', 'arc_length':'m',\
                       'shearrate0':'1/s','shearrate1':'1/s','shearrate2':'1/s',\
                        'shearrate3':'1/s','shearrate4':'1/s','shearrate5':'1/s',\
                        'shearrate6':'1/s','shearrate7':'1/s','shearrate8':'1/s'
                       }
        if 'sliceSummaries' in csvfile:
            xu = 'mm'
            unitdict = {'x':xu, 'xbehind':xu, 'time':'s', 'centery':xu, 'centerz':xu, 'area':xu+'^2', 'maxheight':xu, 'maxwidth':xu, 'centeryn':'', 'centerzn':'', 'arean':'', 'maxheightn':'', 'maxwidthn':'', 'vertdisp':xu, 'vertdispn':'', 'aspectratio':'', 'speed':'mm/s', 'speeddecay':''}
        if 'steadyTimes' in csvfile:
            unitdict = {'x':'mm', 't0':'s', 'tf':'s'}
        if 'steadyPositions' in csvfile:
            unitdict = {'x0':'mm', 'xf':'mm', 't':'s'}
        unitline = [unitdict.get(s, '') for s in line1]   
        blist.append(unitline)
        if line2==unitline:
            line3 = next(csv_reader)
            if not line3==line2 and canreturn:
                return # this table already has units, so we're done
            # otherwise, this table has two rows of units, and the second one needs to be ignored
        else:
            # this table does not have units
            blist.append(line2)
        blist = blist+list(csv_reader)
    with open(csvfile, 'w', newline='', encoding='utf-8') as b:
        writer = csv.writer(b)
        for row in blist:
            writer.writerow(row)
    logging.debug(csvfile)