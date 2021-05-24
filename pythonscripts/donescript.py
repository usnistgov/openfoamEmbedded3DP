#!/usr/bin/env python
'''Functions for moving folders between computers, servers, for OpenFOAM simulations of embedded 3D printing of single filaments. '''

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import foldermover as fm
import folderparser as fp

LOGGERDEFINED = False
LOGGERDEFINED = fp.openLog('donescript.log', LOGGERDEFINED)

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  

if len(sys.argv)>1:
    loopTime = float(sys.argv[2])
else:
    loopTime = 1

fm.doneFolder(sys.argv[1], 2.5, loopTime=loopTime)
    