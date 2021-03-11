import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from paraview_general import *
from typing import List, Dict, Tuple, Union, Any, TextIO


############################
## CSV ##############
##################


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


#### initialize all of paraview
def initializeAll(folder):
    print('CSVs: Initializing paraview for', folder)
    sv = stateVars(folder)        # make subdirectories
    sv =  initializeP(sv)        # initialize Paraview
    sv =  initSeries(sv)        # import the vtk series files  
    return sv


def initSeries(sv):
    casevtmseries = initSeries0(sv)

    clip1 = Clip(Input=casevtmseries)
    clip1.Scalars = ['POINTS', 'alpha.ink']
    clip1.Value = 0.4
    clip1.ClipType = 'Scalar'
    clip1.Invert = 0
    clip2 = Clip(Input=clip1)
    clip2.Scalars = ['POINTS', 'alpha.ink']
    clip2.Value = 0.6
    clip2.ClipType = 'Scalar'
    clip2.Invert = 1
    slice1 = Slice(Input=clip2)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.SliceType.Origin = [0.00334, 0, 0]
    slice1.SliceType.Normal = [1, 0, 0] # look down the y axis
    Hide(casevtmseries, sv.renderView1)
    Hide(clip1, sv.renderView1)
    Hide(clip2, sv.renderView1)
    sv.casevtmseries = casevtmseries
    sv.slice = slice1
    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    return sv

################
# creating csvs

def exportcsv(folder, file, table):
    fn = os.path.join(folder, file)
    with open(fn, mode='w', newline='') as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in table:
            w.writerow(i)
    print('Exported file %s' % fn)


def addToFile(x, tempfile, w, xdisp, sv, skipheader) -> bool:
    ret = skipheader
    sv.slice.SliceType.Origin = [x, 0, 0]
    SaveData(tempfile, FieldAssociation="Point Data", ChooseArraysToWrite=1, AddTime=1)
    with open(tempfile, mode='r') as f2:
        ftemp = csv.reader(f2)
        header=next(ftemp)
        if 'Points:0' not in header:
            return False
        xpos = header.index('Points:0')
        if not skipheader:
            if len(header)>4:
                header2 = []
                for h in header:
                    h2 = h.replace('alpha.ink', 'alpha')
                    h2 = h2.replace('Points:0', 'x')
                    h2 = h2.replace('Points:1', 'y')
                    h2 = h2.replace('Points:2', 'z')
                    h2 = h2.replace('U:0', 'vx')
                    h2 = h2.replace('U:1', 'vy')
                    h2 = h2.replace('U:2', 'vz')
                    h2 = h2.replace('Time', 'time')
                    header2.append(h2)
                w.writerow(header2)
                ret = True
            else:
                ret = False
        for row in ftemp:
            r = row
            if len(r)>1:
                r[xpos] = str(float(r[xpos])+xdisp)
                w.writerow(r)
    return ret

####
# script

def csvfolder(folder:str, forceoverwrite:bool):
    try:
        if not os.path.exists(folder):
            return
        xmin = -0.004824
        xmax = -xmin
        dx = 0.0002
        xchunk = 0.0015
        initialized = False
        f1 = os.path.join(folder, 'interfacePoints')
        tempfile = os.path.join(f1, "temp.csv")
        times = readTimes(folder)
        # times = parseVTKSeries(folder)
        # if len(times)==0:
        #     generateVTKseries(folder, False)
        #     times = parseVTKSeries(folder)
        if len(times)>0:
            for time in times:
                ipfile = os.path.join(f1, "interfacePoints_t_"+str(int(round(time*10)))+".csv")
                    # this is the file that all points for this time will be saved in
                if not os.path.exists(ipfile) or forceoverwrite: 
                    # only run this if the file hasn't been created already or we're being forced to
                    if not initialized: # if paraview hasn't already been initialized, initialize it
                        sv = initializeAll(folder)
                        print('-----', times)
                        sv.times = times
                        initialized = True
                    setTime(time, sv) 
                    # error checking: sometimes we lose the ability to set the time.
                    # This seems to happen only when I walk away from the computer, 
                    # because the computer is a naughty child
                    #print('--t', time, ', file:', ipfile)
                    print(f'--t {time}, file: {ipfile}' )
                    if not sv.timeKeeper1.Time==time:
                        raise Exception('Timekeeper broken. Start again.')
                    with open(ipfile, mode='w', newline='') as f:
                        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                        #w.writerow(['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'alpha', 'p', 'rAU'])
                        skipheader = False
                        for x in (np.arange(xmin, xmax, dx)): 
                            skipheader = addToFile(x, tempfile, w, 0, sv, skipheader)
                    os.remove(tempfile)
            ResetSession()
            cleansession()
            try:
                del sv
            except:
                return
    except Exception as e:
        print(e)
        ResetSession()
        cleansession()
        return