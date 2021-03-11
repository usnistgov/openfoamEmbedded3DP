import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview





forceoverwrite = False
folders = []
# topfolder = 'E:\\Leanne\\OpenFOAM\\newtHBsweep'
# flist = range(344, 350)
# for f in flist:
#     folders.append(os.path.join(topfolder, 'nb'+str(f)))
topfolder = r'C:\Users\lmf1\Documents\OpenFOAM\HBHBsweep'
flist = range(347, 474)
for f in flist:
    folders.append(os.path.join(topfolder, 'nb'+str(f)))


class csv_stateVars():
    def __init__(self, folder):
        self.folder = folder        
        self.casevtmseries = ""
        self.renderView1 = ""
        self.animationScene1 = ""
        self.timeKeeper1 = ""
        self.timeStamp = ""
        self.times = ""
        self.slice = ""

def csv_hideAll():
    ss = GetSources()
    for s in ss:
        Hide(ss[s])

def csv_mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)

def csv_series(folder):
    cf = csv_casefolder(folder)
    vtkfolder = os.path.join(cf, 'VTK')
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if '.vtm.series' in file:
                return os.path.join(vtkfolder, file)
    return ""

#### initialize all of paraview
def csv_initializeAll(folder):
    print('Folder ', folder)
    print(' Initializing paraview ')
    sv = csv_stateVars(folder)        # make subdirectories
    sv =  csv_initializeP(sv)        # initialize Paraview
    sv =  csv_initSeries(sv)        # import the vtk series files  
    return sv


def csv_initializeP(sv):       
    LoadPalette(paletteName='WhiteBackground')
        # make background white
    paraview.simple._DisableFirstRenderCameraReset()
        # disable automatic camera reset on 'Show'  
    csv_hideAll()
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

def csv_initSeries(sv):
    casevtmseries = XMLMultiBlockDataReader(FileName=series(sv.folder))
    casevtmseries.CellArrayStatus = ['alpha.ink', 'U']
    casevtmseries.PointArrayStatus = ['alpha.ink', 'U']
    # sv.times = casevtmseries.TimestepValues
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

def csv_exportcsv(folder, file, table):
    fn = os.path.join(folder, file)
    with open(fn, mode='w', newline='') as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in table:
            w.writerow(i)
    print('Exported file %s' % fn)

def csv_setTime(time, sv):
    if time in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        print(f'time {time} is not in list {sv.times}')

def csv_addToFile(x, tempfile, w, xdisp, sv):
    sv.slice.SliceType.Origin = [x, 0, 0]
    SaveData(tempfile, FieldAssociation="Point Data", ChooseArraysToWrite=1, AddTime=1)
    with open(tempfile, mode='r') as f2:
        ftemp = csv.reader(f2)
        next(ftemp)
        for row in ftemp:
            r = row
            if len(r)>1:
                r[1] = str(float(r[1])+xdisp)
                w.writerow(r)

def csv_casefolder(folder):
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
    
def csv_parseVTKSeries(folder):
    seriesfile = csv_series(folder)
    times = []
    if os.path.exists(seriesfile):
        with open(seriesfile, 'r') as f:
            for line in f:
                if 'name' in line:
                    times.append(float(re.split('time\" : | }', line)[1]))
        return times
    else:
        return []

def csv_cleansession():
    Disconnect()
    Connect()

def csv_csvfolder(folder):
    try:
        if not os.path.exists(folder):
            return
        f1 = os.path.join(folder, 'interfacePoints')
        csv_mkdirif(f1)
        xmin = -0.004824
        xmax = -xmin
        dx = 0.0002
        xchunk = 0.0015
        initialized = False
        tempfile = os.path.join(f1, "temp.csv")
        times = csv_parseVTKSeries(folder) # these times are floats
        if len(times)>0:
            for time in times:
                ipfile = os.path.join(f1, "interfacePoints_t_"+str(int(round(time*10)))+".csv")
                    # this is the file that all points for this time will be saved in
                if not os.path.exists(ipfile) or forceoverwrite: 
                    # only run this if the file hasn't been created already or we're being forced to
                    if not initialized: # if paraview hasn't already been initialized, initialize it
                        sv = csv_initializeAll(folder)
                        print(times)
                        sv.times = times
                        initialized = True
                    csv_setTime(time, sv) 
                    # error checking: sometimes we lose the ability to set the time.
                    # This seems to happen only when I walk away from the computer, 
                    # because the computer is a naughty child
                    print('  time ', time, ', file: ', ipfile)
                    if not sv.timeKeeper1.Time==time:
                        raise Exception('Timekeeper broken. Start again.')
                    with open(ipfile, mode='w', newline='') as f:
                        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                        w.writerow(['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'alpha', 'p', 'rAU'])
                        for x in (np.arange(xmin, xmax, dx)): 
                            csv_addToFile(x, tempfile, w, 0, sv)                            
                    os.remove(tempfile)
            ResetSession()
            csv_cleansession()
            try:
                del sv
            except:
                return
    except Exception as e:
        print(e)
        ResetSession()
        csv_cleansession()
        return

######################################################
####################### SCRIPT #######################
def csv_mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)

for folder in folders:
    csv_csvfolder(folder)
print('Done exporting csv files')
