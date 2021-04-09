#!/usr/bin/env pvpython
'''Collecting line traces through the bath in vtk files. Scripting for many folders and many images and tables.'''

import os
import logging

from paraview.simple import * # import the simple module from the paraview

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# load the virtual environment
virtualEnv = os.path.join(parentdir, 'env', 'Scripts', 'activate_this.py')
if sys.version_info.major < 3:
    execfile(virtualEnv, dict(__file__=virtualEnv))
else:
    exec(open(virtualEnv).read(), {'__file__': virtualEnv})
    

from config import cfg
from paraview_line import csvfolder

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"


#################################################################

#------------------------------

tlist = [1]        # times at which to collect the traces
xpos = -0.001      # position at which to collect the trace

forceoverwrite = False

folders = []

nlist = range(0, 1000)

SERVERFOLDER = cfg.path.server
topfolders = [os.path.join(SERVERFOLDER, s) for s in ['newtnewtsweep', 'HBnewtsweep', 'newtHBsweep', 'HBHBsweep']]
for topfolder in topfolders:
    for f in os.listdir(topfolder):
        if f.startswith('nb'):
            n1 = float(f[2:])
            if n1 in nlist:
                folders.append(os.path.join(topfolder, f))
                

for folder in folders:
    csvfolder(folder, xpos)
print('Done exporting csv files')
                
#-----------------------------


# class stateVars():
#     def __init__(self, folder):
#         self.folder = folder        
#         self.casevtmseries = ""
#         self.renderView1 = ""
#         self.animationScene1 = ""
#         self.timeKeeper1 = ""
#         self.timeStamp = ""
#         self.times = ""
#         self.slice = ""

# def hideAll():
#     ss = GetSources()
#     for s in ss:
#         Hide(ss[s])


# def series(folder):
#     cf = casefolder(folder)
#     return os.path.join(cf, 'VTK', os.path.basename(cf)+'.vtm.series')




# def initializeP(sv):       
#     LoadPalette(paletteName='WhiteBackground')
#         # make background white
#     paraview.simple._DisableFirstRenderCameraReset()
#         # disable automatic camera reset on 'Show'  
#     hideAll()
#         # hide all existing sources
#     sv.animationScene1 = GetAnimationScene()
#         # get animation scene   
#     sv.timeKeeper1 = GetTimeKeeper()
#         # get the time-keeper  
#     sv.animationScene1.UpdateAnimationUsingDataTimeSteps()
#         # update animation scene based on data timesteps
#     sv.renderView1 = GetActiveViewOrCreate('RenderView')
#         # get active view
#     sv.renderView1.ViewSize = [1216, 1216]
#         # set view size
#     sv.renderView1.OrientationAxesVisibility = 0
#         # hide orientation axes
#     sv.renderView1.CameraParallelProjection = 1
#         # turn off perspective
#     return sv


# def setTime(time, sv):
#     if time in sv.times:
#         sv.animationScene1.AnimationTime = time
#         sv.timeKeeper1.Time = time
#     else:
#         print(f'time {time} is not in list')

# def addToFile(x, tempfile, w, xdisp, sv):
#     sv.slice.SliceType.Origin = [x, 0, 0]
#     SaveData(tempfile, FieldAssociation="Point Data", ChooseArraysToWrite=1, AddTime=1)
#     with open(tempfile, mode='r') as f2:
#         ftemp = csv.reader(f2)
#         next(ftemp)
#         for row in ftemp:
#             r = row
#             if len(r)>1:
#                 r[1] = str(float(r[1])+xdisp)
#                 w.writerow(r)

# def casefolder(folder):
#     # if there is a folder in this folder named 'case', return that folder
#     casefold = os.path.join(folder, 'case')
#     if os.path.exists(casefold):
#         return casefold
#     # if this folder contains a folder called 'constant', this is the case folder, so return this folder
#     constfold = os.path.join(folder, 'constant')
#     if os.path.exists(constfold):
#         return folder
#     else:
#         return ''
    
# def parseVTKSeries(folder):
#     cf = casefolder(folder)
#     bn = os.path.basename(cf)
#     seriesfile = os.path.join(cf, 'VTK', bn+'.vtm.series')
#     times = []
#     if os.path.exists(seriesfile):
#         with open(seriesfile, 'r') as f:
#             for line in f:
#                 if 'name' in line:
#                     times.append(float(re.split('time\" : | }', line)[1]))
#         return times
#     else:
#         return []

# def cleansession():
#     Disconnect()
#     Connect()







