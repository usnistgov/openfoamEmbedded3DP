import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from paraview_general import *
import paraview_csv as pc
import paraview_screenshots_v2 as ss

looping = True

## csv

forceoverwrite = False


# screenshots
runlist = []
runlist.append(ss.ssVars('volumes', [], volviewlist=['y']))
for s in ['viscy', 'viscx', 'uslicey', 'uslicex']:
    runlist.append(ss.ssVars(s, [0.5, 1, 2.5]))
# runlist.append(ss.ssVars('meshes', [2.5]))
#runlist.append(ss.ssVars('vectors', [2.5]))
# runlist.append(ss.ssVars('tubes', [2.5], tubeh=0.001, volviewlist=['a']))

# folders
folders = []

#nlist = range(1000, 1100)
nlist=range(0, 1000)

#topfolders = [r'C:\Users\lmf1\Documents\OpenFOAM\HBHBsweep']
SERVERFOLDER = r'\\cfs2e.nist.gov\642\NIST_Projects\Additive Manufacturing and Rheology\OpenFOAM\simulations\viscositysweep'
# topfolders = [os.path.join(SERVERFOLDER, s) for s in ['newtnewtsweep', 'HBnewtsweep', 'newtHBsweep', 'HBHBsweep']]
#SERVERFOLDER = r'\\cfs2e.nist.gov\642\NIST_Projects\Additive Manufacturing and Rheology\OpenFOAM\simulations\yieldingsweep\HBHByielded'
#topfolders = [os.path.join(SERVERFOLDER, s) for s in ['n', 'k', 'tau0']]
topfolders = [os.path.join(SERVERFOLDER, s) for s in ['HBnewtsweep', 'HBHBsweep', 'newtnewtsweep', 'newtHBsweep']]
for topfolder in topfolders:
    for f in os.listdir(topfolder):
        if f.startswith('nb'):
            n1 = float(f[2:])
            if n1 in nlist:
                folders.append(os.path.join(topfolder, f))


        
######################################################
####################### SCRIPT #######################


while True:
    for folder in folders:
        pc.csvfolder(folder, forceoverwrite)
        ss.folderscript(folder, runlist)
    if not looping:
        break
    now = datetime.now()
    current_time = now.strftime("%D, %H:%M:%S")
    print('Waiting 6 hours for more files...')
    print("------ Current Time =", current_time)
    time.sleep(60*60*6)
