#!/usr/bin/env pvpython
'''Collecting interface points into csvs from vtk files and generating images from vtk files. Scripting for many folders and many images and tables.'''

# external packages
import os
import sys
import time
from datetime import datetime
import logging

# local packages
from paraview_general import *
import paraview_csv as pc
import paraview_screenshots as ss

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
# load the virtual environment
    # this needs to be after paraview_general, because it will otherwise break numpy
virtualEnv = os.path.join(parentdir, 'env', 'Scripts', 'activate_this.py')
if sys.version_info.major < 3:
    execfile(virtualEnv, dict(__file__=virtualEnv))
else:
    exec(open(virtualEnv).read(), {'__file__': virtualEnv})

from config import cfg
import folderparser as fp

# logging
LOGGERDEFINED = fp.openLog(os.path.realpath(__file__), False, level='DEBUG')

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Development"


#################################################################

looping = False
loopTime = 6 #hours

## csv

forceOverwrite = False
getCSVs = True


# screenshots

# for help, see ssVars in paraview_screenshots
# first entry is the type of image, second is list of times in s (leave empty to collect all). Optional keywords are described in ssVars

runList = []
# runList.append(ss.ssVars('volumes', [], volViewList=['a']))
# runList.append(ss.ssVars('volumes', [], volViewList=['y']))
# runList.append(ss.ssVars('volumes', [1.0, 2.5], volViewList=['y']))
# for s in ['viscy', 'viscx', 'uslicey', 'uslicex']:
#     runList.append(ss.ssVars(s, [1, 2.5]))
# for s in ['py', 'uslicey', 'uzslicey']:
#     runList.append(ss.ssVars(s, [1.0, 2.5]))
for s in ['shearRatex', 'shearRatey']:
    runList.append(ss.ssVars(s, [1.0, 2.5]))
# runlist.append(ss.ssVars('meshes', [2.5]))
# runlist.append(ss.ssVars('vectors', [2.5]))
# runlist.append(ss.ssVars('tubes', [2.5], tubeh=0.001, volviewlist=['a']))

# folders
folders = []

nlist = [1020]
# nlist=range(0, 1000)
# nlist = [455]


SERVERFOLDER = cfg.path.server
if not os.path.exists(SERVERFOLDER):
    logging.error('Server folder in config.yml does not exist')
    raise FileNotFoundError('Server folder in config.yml does not exist')

# topfolders = [os.path.join(SERVERFOLDER, 'viscositysweep',  s) for s in ['newtHBsweep', 'newtnewtsweep', 'HBHBsweep', 'HBnewtsweep']]

topfolders = [os.path.join(SERVERFOLDER, 'yieldingsweep', 'HBHByielded', s) for s in ['k', 'n', 'tau0']]
topfolders = topfolders + [os.path.join(SERVERFOLDER, 'yieldingsweep', 'LapRD')]
for topfolder in topfolders:
    for f in os.listdir(topfolder):
        if f.startswith('nb'):
            try:
                n1 = float(f[2:])
            except:
                n1 = f[2:]
            if n1 in nlist:
                folders.append(os.path.join(topfolder, f))


        
######################################################
####################### SCRIPT #######################

logging.info('Exporting images and csvs.')
logging.info('Images:')
for r in runList:
    logging.info(r.prnt())
logging.info(f'CSVs: Collect: {getCSVs}. Overwrite: {forceOverwrite}')
logging.info(f'Folders: {[os.path.basename(f) for f in folders]}')

while True:
    for folder in folders:
        logging.debug('Checking '+folder)
        if getCSVs:
            pc.csvFolder(folder, forceOverwrite) # create csvs
        ss.folderScript(folder, runList)
    if not looping:
        break
    now = datetime.now()
    current_time = now.strftime("%D, %H:%M:%S")
    logging.info('Waiting '+str(loopTime)+' hours for more files...')
    logging.info("------ Current Time ="+str(current_time))
    time.sleep(60*60*loopTime)
