import foldermover as fm
import folderparser as fp
import os
import time
import socket
import sys
from config import cfg

dir_path = os.path.dirname(os.path.realpath(__file__))
LOGGERDEFINED = False
compname = socket.gethostname()
LOGGERDEFINED = fp.openLog(os.path.join(dir_path,'folderparser_'+compname+'.log'), LOGGERDEFINED)

SERVERFOLDER = os.path.join(cfg.server, 'yieldingsweep', 'HBHByielded')
CFOLDER = os.path.join(cfg.c, 'HBHByielded')
EFOLDER = os.path.join(cfg.e, 'HBHByielded')

if len(sys.argv)>1:
    waittime = float(sys.argv[1])
else:
    waittime = 0

fm.copyCtoServerFolders(SERVERFOLDER, CFOLDER,['k', 'n', 'tau0'], waittime, EFOLDER=EFOLDER)