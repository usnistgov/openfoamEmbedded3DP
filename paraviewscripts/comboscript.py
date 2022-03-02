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
import paraview_line as pl

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

#------
## csv

forceOverwrite = False
getCSVs = True
csvTimes = [2.5]

# modes = ['nozzle', 'interface'] # CSVs to create
modes = ['nozzleSlice', 'interface']

#------
## screenshots

getSSs = True
# for help, see ssVars in paraview_screenshots
# first entry is the type of image, second is list of times in s (leave empty to collect all). Optional keywords are described in ssVars

runList = []
runList.append(ss.ssVars('volumes', [], volViewList=['y']))
# runList.append(ss.ssVars('volumes', [1.0, 2.5], volViewList=['y']))
for s in ['viscy', 'viscx', 'py', 'uslicey', 'uslicex', 'shearStressx', 'shearStressy']:
    runList.append(ss.ssVars(s, [1, 2.5]))
# for s in ['py', 'uslicey', 'uzslicey']:
#     runList.append(ss.ssVars(s, [1.0, 2.5]))
# for s in ['shearRatex', 'shearRatey']:
#     runList.append(ss.ssVars(s, [1.0, 2.5]))
# runlist.append(ss.ssVars('meshes', [2.5]))
# runlist.append(ss.ssVars('vectors', [2.5]))
runList.append(ss.ssVars('tubes', [2.5], tubeh=0.001, volviewlist=['a']))

#------
## line traces

getLine = True

pl_tlist = [2.5]        # times at which to collect the traces
# xlist = [-0.001]      
pl_xlist = [2]       # positions at which to collect the trace, in nozzle inner diameters, relative to the center
# zlist = [0.0005]
pl_zlist = []

pl_forceOverwrite = False # True to overwrite existing files

#----------------------------------

# folders
SERVERFOLDER = cfg.path.server
if not os.path.exists(SERVERFOLDER):
    logging.error('Server folder in config.yml does not exist')
    raise FileNotFoundError('Server folder in config.yml does not exist')


# nlist = list(range(0,1000))
# nlist = [43, 53, 63, 487]
nlist = [240]
topfolders = [os.path.join(cfg.path.server, 'conicalnozzle', s) for s in ['orig', 'speed_sweep', 'diameter', 'newtonian', 'k', 'newt_diameter']]
# topfolders = [os.path.join(cfg.path.server, 'viscositysweep', s) for s in ['newtnewtsweep']]

folders = filterSimNums(topfolders, nlist)

######################################################
####################### SCRIPT #######################

logging.info('Exporting images and csvs.')
logging.info('Images:')
for r in runList:
    logging.info(r.prnt())
logging.info(f'CSVs: Collect: {getCSVs}. Overwrite: {forceOverwrite}')
logging.info(f'Folders: {[os.path.basename(f) for f in folders]}')
logging.info(f'Modes: {modes}')
logging.info(f'Line traces:\n\
            X positions: {pl_xlist} di behind nozzle.\n\
            Z positions: {pl_zlist} di above nozzle bottom.\n\
            Time list: {pl_tlist} s.')

while True:

    for folder in folders:
        logging.debug(f'Checking {folder}')
        if getCSVs:
            pc.csvFolder(folder, modes, forceOverwrite, times0=csvTimes) # create csvs RG
        if getSSs:
            ss.folderScript(folder, runList)
        if getLine:
            for xpos in pl_xlist:
                pl.csvfolder(folder, 'x', xpos, pl_tlist, forceOverwrite=pl_forceOverwrite)
            for zpos in pl_zlist:
                pl.csvfolder(folder, 'z', zpos, pl_tlist, forceOverwrite=pl_forceOverwrite)
    if not looping:
        break
    now = datetime.now()
    current_time = now.strftime("%D, %H:%M:%S")
    logging.info(f'Waiting {looptime} hours for more files...')
    logging.info(f'------ Current Time ={current_time}')
    time.sleep(60*60*loopTime)

