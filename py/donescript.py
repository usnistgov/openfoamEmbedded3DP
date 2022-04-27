#!/usr/bin/env python
'''Functions for moving folders between computers, servers, for OpenFOAM simulations of embedded 3D printing of single filaments. '''

# global packages
import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# local packages
import foldermover as fm
import folderparser as fp

# logging
LOGGERDEFINED = False
LOGGERDEFINED = fp.openLog('donescript.log', LOGGERDEFINED)


#-------------------------------------------------------------------------------------------------  

if len(sys.argv)>2:
    loopTime = float(sys.argv[2])
else:
    loopTime = 1

fm.doneFolder(sys.argv[1], 2.5, loopTime=loopTime)
    