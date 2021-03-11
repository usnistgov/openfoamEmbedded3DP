import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from typing import List, Dict, Tuple, Union, Any, TextIO
sys.path.insert(1, r'C:\Users\lmf1\OneDriveNIST\NIST\data\openfoam\pythonscripts')
from folderparser import redoVTKSeriesNoLog, parseVTKSeries, vtkfiles

class stateVars():
    def __init__(self, folder):
        mksubdirs(folder)
        self.folder = folder        
        self.casevtmseries = ""
        self.renderView1 = ""
        self.animationScene1 = ""
        self.timeKeeper1 = ""
        self.timeStamp = ""
        self.times = ""
        self.slice = ""
        self.initialized = False
        
    # get viscosity settings from legend
    def readVisc(self):
        le = importIf(self.folder)
        if len(le)==0:
            return 1
        i0 = 0
        while not le[i0][0]=='transportProperties':
            i0+=1 # find the index where the transport properties start
        self.inkmodel = le[i0+2][1] # find the row that contains the ink model.
        if self.inkmodel=='Newtonian':
            # if it's newtonian, the viscosity is stored 3 rows below 'transportProperties'
            nuinki = i0+3
            self.inkfunc = 'inknu'
        else:
            # otherwise the zero shear viscosity is 4 rows below
            nuinki = i0+4
            self.inkfunc = 'nu1'
        nuink = le[nuinki][1] # get the string viscosity in mPa s
        self.inknu = nuink

        # find the support transport properties
        # in some legends, there are no Herschel Bulkley boxes, and in some there are
        # we need to find where the support section starts
        nusupi = nuinki 

        # iterate through rows until we find the support transportmodel
        while not le[nusupi][0]=='transportModel':
            nusupi+=1

        self.supmodel = le[nusupi][1]
        if self.supmodel=='Newtonian':
            nusupi = nusupi+1
            self.supfunc = 'supnu'
        else:
            nusupi = nusupi+2
            self.supfunc = 'nu2'
        nusup = le[nusupi][1]
        self.supnu = nusup
        
        return 0

    
    
#         self.leg = importIf(self.folder)
#         if len(self.leg)>0:
            
#             self.inkmodel = self.leg[92][1]
#             if self.inkmodel=='HerschelBulkley':
#                 self.inknu = self.leg[94][1]
#                 self.inkfunc = 'nu1'
#             else:
#                 self.inknu = self.leg[93][1]
#                 self.inkfunc = 'inknu'
#             print(self.leg[95][0])
#             if self.leg[95][0]=='sup':
#                 sup0 = 96
#             else:
#                 sup0 = 100
#             self.supmodel = self.leg[sup0][1]
#             if self.supmodel=='HerschelBulkley':
#                 self.supnu = self.leg[sup0+2][1]
#                 self.supfunc = 'nu2'
#             else:
#                 self.supnu = self.leg[sup0+1][1]
#                 self.supfunc = 'supnu'
#             return 0
#         else:
#             return 1
    
        

        
# importIf imports and existing legend file if it exists 
# or creates a table if it does not exist
    # folder is a full path
def importIf(folder:str) -> List[List[str]]:
    fn = os.path.join(folder, 'legend.csv')
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            return(list(csv.reader(f)))
    else:
        return []

    # make the images folder      
def mksubdirs(folder:str) -> None:
    mkdirif(os.path.join(folder, 'images'))   
    mkdirif(os.path.join(folder, 'interfacePoints'))

# find the case folder
def casefolder(folder:str) -> str:
    # if there is a folder in this folder named 'case', return that folder
    casefold = os.path.join(folder, 'case')
    if os.path.exists(casefold):
        return casefold
    # if this folder contains a folder called 'constant', this is the case folder, so return this folder
    constfold = os.path.join(folder, 'constant')
    if os.path.exists(constfold):
        return folder
    vtkfold = os.path.join(folder, 'VTK')
    if os.path.exists(vtkfold):
        return folder
    intfold = os.path.join(folder, 'interfacePoints')
    if os.path.exists(intfold):
        return folder
    else:
        #print('no case folder in ', folder)
        return ''

# make the directory if there is no existing directory
def mkdirif(path:str) -> None:
    try:
        os.mkdir(path)
    except OSError:
        return 
    else:
        print ("Created directory %s" % path)

        
# find the path of the vtk series file
def series(folder:str) -> str:
    cf = casefolder(folder)
    vtkfolder = os.path.join(cf, 'VTK')
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if '.series' in file:
                return os.path.join(vtkfolder, file)
        # if there is a vtk folder but no series file, generate one
        redoVTKSeriesNoLog(folder)
        return parseVTKSeries(folder)
        # return generateVTKseries(folder, False)
    return ""
        
        
# # get the list of times from the vtk series file
# def parseVTKSeries(folder:str) -> List[float]:
#     seriesfile = series(folder)
#     times = []
#     if os.path.exists(seriesfile):
#         with open(seriesfile, 'r') as f:
#             for line in f:
#                 if 'name' in line:
#                     times.append(float(re.split('time\" : | }', line)[1]))
#         return times
#     else:
#         return []


# generate a new .VTK.series file if there is no series file
# searchForSeries is true if we should search for the series file
# def generateVTKseries(folder:str, searchForSeries:bool) -> str:
#     print('generating new series files')
#     if searchForSeries:
#         s = series(folder)
#         if os.path.exists(s):
#             return s
#     # check that there are vtk files
#     cf = casefolder(folder)
#     vtkfolder = os.path.join(cf, 'VTK')
#     if not os.path.exists(vtkfolder):
#         return ''
#     # find the log file
#     log = os.path.join(folder, 'log_foamToVTK')
#     if not os.path.exists(log):
#         log = os.path.join(cf, 'log_foamToVTK')
#         if not os.path.exists(log):
#             return VTKguessTimeSeries(folder)
#     fn = os.path.join(vtkfolder, 'case.vtk.series')
#     times = []
#     with open(fn,"w") as fout:
#         out = '{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n'
#         with open(log, 'r') as f:
#             for line in f:
#                 if line.startswith('Time'):
#                     time = re.split(': |\n', line)[1]
#                 if 'Internal' in line:
#                     filename = re.split('VTK/|\"', line)[2]
#                     if not time in times:
#                         out = out + '    { \"name\" : \"'+filename+'\", \"time\" : '+time+' },\n'
#                         times.append(time)
#         out = out[0:-2]+'\n  ]\n}' # eliminate the last comma
#         fout.write(out)
#     print('Exported ' + fn)
#     return fn


def readTimes(folder) -> List[float]:
    times = parseVTKSeries(folder)
    if len(times)==0 or len(times)<vtkfiles(folder):
        redoVTKSeriesNoLog(folder)
        times = parseVTKSeries(folder)
    return times

##-----------------------
# paraview

# reset the paraview session
def cleansession() -> None:
    Disconnect()
    Connect()

    
# hide all existing displays
def hideAll() -> None:
    ss = GetSources()
    for s in ss:
        Hide(ss[s])

# import the vtk series file
# this returns the datareader object that refers to the series
def initSeries0(sv:stateVars) -> Any:
    sfile = series(sv.folder)
    if not os.path.exists(sfile):
        raise NameError('vtm.series file does not exist')
    if 'vtm' in sfile:
        casevtmseries = XMLMultiBlockDataReader(FileName=sfile)
        casevtmseries.CellArrayStatus = ['U', 'alpha.ink', 'p_rgh']
        casevtmseries.PointArrayStatus = ['U', 'alpha.ink', 'p_rgh']
    else:
        bn = os.path.basename(sfile)
        casevtmseries = LegacyVTKReader(registrationName=bn, FileNames=[sfile])
    return casevtmseries
    

# go to a time
def setTime(time:float, sv:stateVars) -> None:
    if time in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        if max(sv.times)>0.1*len(sv.times):
            redoVTKSeriesNoLog(sv.folder)
            sv.times = parseVTKSeries(sv.folder)
            if time in sv.times:
                sv.animationScene1.AnimationTime = time
                sv.timeKeeper1.Time = time
            else:
                print(f'time {time} is not in list {sv.times}')
                sv.animationScene1.AnimationTime = sv.times[-1]
                sv.timeKeeper1.Time = sv.times[-1]
        
        
        
## initialize paraview
def initializeP(sv:stateVars) -> stateVars:  
    ResetSession()     
    LoadPalette(paletteName='WhiteBackground')
        # make background white
    paraview.simple._DisableFirstRenderCameraReset()
        # disable automatic camera reset on 'Show'  
    hideAll()
        # hide all existing sources
    sv.animationScene1 = GetAnimationScene()
        # get animation scene   
    sv.timeKeeper1 = GetTimeKeeper()
        # get the time-keeper  
    sv.animationScene1.UpdateAnimationUsingDataTimeSteps()
        # update animation scene based on data timesteps
    sv.renderView1 = GetActiveViewOrCreate('RenderView')
        # get active view
    sv.renderView1.ViewSize = [1216, 1216]
        # set view size
    sv.renderView1.OrientationAxesVisibility = 0
        # hide orientation axes
    sv.renderView1.CameraParallelProjection = 1
        # turn off perspective
    return sv
