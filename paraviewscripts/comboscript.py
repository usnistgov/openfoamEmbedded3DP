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
sys.path.append(os.path.join(parentdir, 'py'))  # add python folder

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



#################################################################

looping = False
loopTime = 6 #hours

#------
## csv

forceOverwrite = False
getCSVs = False
getNoz = False
# csvTimes = [2.5]
csvTimes = [round(0.1*i,1) for i in range(1,26)]
nozTimes = [2.5]

#------
## screenshots

getSSs = True
ss_forceOverwrite=True
# for help, see ssVars in paraview_screenshots
# first entry is the type of image, second is list of times in s (leave empty to collect all). Optional keywords are described in ssVars

runList = []
# runList.append(ss.ssVars('volumes', [], volViewList=['y']))
# runList.append(ss.ssVars('volumes', [1.0, 2.5], volViewList=['y']))
# for s in ['viscy', 'viscx', 'py', 'uslicey', 'uslicex', 'shearStressx', 'shearStressy']:
for s in ['shearStressy']:
    runList.append(ss.ssVars(s, [2.5], yieldSurface=True))
# for s in ['py', 'uslicey', 'uzslicey']:
#     runList.append(ss.ssVars(s, [1.0, 2.5]))
# for s in ['shearRatex', 'shearRatey']:
#     runList.append(ss.ssVars(s, [1.0, 2.5]))
# runlist.append(ss.ssVars('meshes', [2.5]))
# runlist.append(ss.ssVars('vectors', [2.5]))
# runList.append(ss.ssVars('tubes', [2.5], tubeh=0.001, volviewlist=['a']))

#------
## line traces

getLine = False

pl_tlist = [2.5]        # times at which to collect the traces     
pl_xlist = [2]       # positions at which to collect the trace, in nozzle inner diameters, relative to the center
pl_zlist = []

pl_forceOverwrite = False # True to overwrite existing files

#----------------------------------

# folders
SERVERFOLDER = cfg.path.server
if not os.path.exists(SERVERFOLDER):
    logging.error('Server folder in config.yml does not exist')
    raise FileNotFoundError('Server folder in config.yml does not exist')


# nlist = list(range(0,1000))
nlist = [213]
topfolders = [os.path.join(cfg.path.server, 'conicalnozzle', s) for s in ['HB_angle', 'HB_diameter', 'HB_k', 'HB_speed', 'newt_angle', 'newt_diameter', 'newt_visc', 'visc_speed']]

folders = filterSimNums(topfolders, nlist)

######################################################
####################### SCRIPT #######################

logging.info('Exporting images and csvs.')
logging.info(f'Images. Overwrite: {ss_forceOverwrite}')
for r in runList:
    logging.info(r.prnt())
logging.info(f'Interface CSVs: Collect: {getCSVs}. Overwrite: {forceOverwrite}')
logging.info(f'Nozzle CSVs: Collect: {getNoz}. Overwrite: {forceOverwrite}')
logging.info(f'Folders: {[os.path.basename(f) for f in folders]}')
logging.info(f'Line traces:\n\
            X positions: {pl_xlist} di behind nozzle.\n\
            Z positions: {pl_zlist} di above nozzle bottom.\n\
            Time list: {pl_tlist} s.')

while True:

    for folder in folders:
        logging.debug(f'Checking {folder}')
        if getCSVs:
            pc.csvFolder(folder, ['interface'], forceOverwrite, times0=csvTimes) # create csvs RG
        if getNoz:
            pc.csvFolder(folder, ['nozzleSlice'], forceOverwrite, times0=nozTimes) # create nozzle point csv
        if getSSs:
            ss.folderScript(folder, runList, overwrite=ss_forceOverwrite)  # screenshots
        if getLine:
            for xpos in pl_xlist:
                pl.csvfolder(folder, 'x', xpos, pl_tlist, forceOverwrite=pl_forceOverwrite) # line trace at constant x
            for zpos in pl_zlist:
                pl.csvfolder(folder, 'z', zpos, pl_tlist, forceOverwrite=pl_forceOverwrite) # line trace at constant z
    if not looping:
        break
    now = datetime.now()
    current_time = now.strftime("%D, %H:%M:%S")
    logging.info(f'Waiting {looptime} hours for more files...')
    logging.info(f'------ Current Time ={current_time}')
    time.sleep(60*60*loopTime)

