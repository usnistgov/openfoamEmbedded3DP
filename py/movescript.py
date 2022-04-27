#!/usr/bin/env python
'''Moves folders to server. Loops continuously.'''

# external packages
import os
import time
import socket
import sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import foldermover as fm
import folderparser as fp
from config import cfg

# logging
LOGGERDEFINED = fp.openLog(os.path.realpath(__file__), LOGGERDEFINED)



#---------------------------------------


SERVERFOLDER = os.path.join(cfg.path.server, 'yieldingsweep', 'HBHByielded')
CFOLDER = os.path.join(cfg.path.c, 'HBHByielded')
EFOLDER = os.path.join(cfg.path.e, 'HBHByielded')

if len(sys.argv)>1:
    waittime = float(sys.argv[1])
else:
    waittime = 0

fm.copyCtoServerFolders(SERVERFOLDER, CFOLDER,['k', 'n', 'tau0'], waittime, EFOLDER=EFOLDER)