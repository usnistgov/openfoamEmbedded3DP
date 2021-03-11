import foldermover as fm
import folderparser as fp
import os
import time
import socket
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
LOGGERDEFINED = False
compname = socket.gethostname()
LOGGERDEFINED = fp.openLog(os.path.join(dir_path,'folderparser_'+compname+'.log'), LOGGERDEFINED)

SERVERFOLDER = r'\\cfs2e.nist.gov\642\NIST_Projects\Additive Manufacturing and Rheology\OpenFOAM\simulations\yieldingsweep\HBHByielded'
CFOLDER = r'C:\Users\lmf1\Documents\OpenFOAM\HBHByielded'
EFOLDER = r'E:\Leanne\OpenFOAM\HBHByielded'

if len(sys.argv)>1:
    waittime = float(sys.argv[1])
else:
    waittime = 0

fm.copyCtoServerFolders(SERVERFOLDER, CFOLDER,['k', 'n', 'tau0'], waittime, EFOLDER=EFOLDER)