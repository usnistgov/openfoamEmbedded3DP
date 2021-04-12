#!/usr/bin/env pvpython
'''Collecting interface points into csvs from vtk files and generating images from vtk files. Scripting for many folders and many images and tables.'''


import os
import sys
# from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime

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
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)



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


# screenshots
runList = []
runList.append(ss.ssVars('volumes', [], volViewList=['y']))
for s in ['viscy', 'viscx', 'uslicey', 'uslicex']:
    runList.append(ss.ssVars(s, [0.5, 1, 2.5]))
# runlist.append(ss.ssVars('meshes', [2.5]))
#runlist.append(ss.ssVars('vectors', [2.5]))
# runlist.append(ss.ssVars('tubes', [2.5], tubeh=0.001, volviewlist=['a']))

# folders
folders = []

#nlist = range(1000, 1100)
nlist=range(0, 1000)



SERVERFOLDER = cfg.path.server
if not os.path.exists(SERVERFOLDER):
    logging.error('Server folder in config.yml does not exist')
    raise FileNotFoundError('Server folder in config.yml does not exist')

topfolders = [os.path.join(SERVERFOLDER, 'viscositysweep',  s) for s in ['newtHBsweep', 'newtnewtsweep', 'HBHBsweep', 'HBnewtsweep']]
for topfolder in topfolders:
    for f in os.listdir(topfolder):
        if f.startswith('nb'):
            n1 = float(f[2:])
            if n1 in nlist:
                folders.append(os.path.join(topfolder, f))


        
######################################################
####################### SCRIPT #######################

logging.info(f'Exporting images and csvs.\nImages: {runList}.\nCSVs: Overwrite {forceOverwrite}.\nFolders: {[os.path.basename(f) for f in folders]}')

while True:
    for folder in folders:
        logging.debug('Checking '+folder)
        pc.csvFolder(folder, forceOverwrite) # create csvs
        ss.folderScript(folder, runList)
    if not looping:
        break
    now = datetime.now()
    current_time = now.strftime("%D, %H:%M:%S")
    logging.info('Waiting '+str(loopTime)+' hours for more files...')
    logging.info("------ Current Time ="+str(current_time))
    time.sleep(60*60*loopTime)
