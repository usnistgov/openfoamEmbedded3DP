import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import logging
logger = logging.getLogger(__name__)



tlist = [1]
xpos = -0.001

forceoverwrite = False

folders = []

# nlist = [247, 241, 253, 221, 227, 709, 533, 37, 45, 53, 61]
nlist = range(0, 1000)

#topfolders = [r'C:\Users\lmf1\Documents\OpenFOAM\HBHBsweep']
SERVERFOLDER = r'\\cfs2e.nist.gov\642\NIST_Projects\Additive Manufacturing and Rheology\OpenFOAM\simulations'
topfolders = [os.path.join(SERVERFOLDER, s) for s in ['newtnewtsweep', 'HBnewtsweep', 'newtHBsweep', 'HBHBsweep']]
# topfolders = [os.path.join(SERVERFOLDER, s) for s in ['HBnewtsweep']]
for topfolder in topfolders:
    for f in os.listdir(topfolder):
        if f.startswith('nb'):
            n1 = float(f[2:])
            if n1 in nlist:
                folders.append(os.path.join(topfolder, f))


class stateVars():
    def __init__(self, folder):
        self.folder = folder        
        self.casevtmseries = ""
        self.renderView1 = ""
        self.animationScene1 = ""
        self.timeKeeper1 = ""
        self.timeStamp = ""
        self.times = ""
        self.slice = ""

def hideAll():
    ss = GetSources()
    for s in ss:
        Hide(ss[s])

def mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)

def series(folder):
    cf = casefolder(folder)
    return os.path.join(cf, 'VTK', os.path.basename(cf)+'.vtm.series')

#### initialize all of paraview
def initializeAll(folder, xpos):
    print('Folder ', folder)
    print(' Initializing paraview ')
    sv = stateVars(folder)        # make subdirectories
    sv =  initializeP(sv)        # initialize Paraview
    sv =  initSeries(sv, xpos)        # import the vtk series files  
    return sv


def initializeP(sv):       
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

def initSeries(sv, xpos):
    casevtmseries = XMLMultiBlockDataReader(FileName=series(sv.folder))
    casevtmseries.CellArrayStatus = ['alpha.ink', 'U']
    casevtmseries.PointArrayStatus = ['alpha.ink', 'U']
    # sv.times = casevtmseries.TimestepValues
    plotOverLine2 = PlotOverLine(Input=casevtmseries, Source='High Resolution Line Source')
    # plotOverLine2.Source.Point1 = [0.002588, 0, 0.0021104998886585236]
    # plotOverLine2.Source.Point2 = [0.002588, 0, -0.0021104998886585236]
    
    plotOverLine2.Source.Point1 = [xpos, 0, 0.0021104998886585236]
    plotOverLine2.Source.Point2 = [xpos, 0, -0.0021104998886585236]

    # show data in view
    plotOverLine1Display = Show(plotOverLine2, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    plotOverLine1Display.Representation = 'Surface'
    plotOverLine1Display.ColorArrayName = [None, '']
    plotOverLine1Display.OSPRayScaleArray = 'U'
    plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plotOverLine1Display.SelectOrientationVectors = 'None'
    plotOverLine1Display.ScaleFactor = 0.00042209997773170473
    plotOverLine1Display.SelectScaleArray = 'None'
    plotOverLine1Display.GlyphType = 'Arrow'
    plotOverLine1Display.GlyphTableIndexArray = 'None'
    plotOverLine1Display.GaussianRadius = 2.1104998886585235e-05
    plotOverLine1Display.SetScaleArray = ['POINTS', 'U']
    plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
    plotOverLine1Display.OpacityArray = ['POINTS', 'U']
    plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
    plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
    plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    plotOverLine1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    plotOverLine1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    lineChartView1 = CreateView('XYChartView')

    sv.plotOverLine2 = plotOverLine2
    Hide(casevtmseries, sv.renderView1)
    sv.casevtmseries = casevtmseries
    return sv


def exportcsv(folder, file, sv):

    fn = os.path.join(folder, file)
    SaveData(fn, proxy=sv.plotOverLine2, PointDataArrays=['U', 'alpha.ink', 'arc_length', 'p', 'p_rgh'], AddTime=1)

def setTime(time, sv):
    if time in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        print(f'time {time} is not in list')

def addToFile(x, tempfile, w, xdisp, sv):
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

def casefolder(folder):
    # if there is a folder in this folder named 'case', return that folder
    casefold = os.path.join(folder, 'case')
    if os.path.exists(casefold):
        return casefold
    # if this folder contains a folder called 'constant', this is the case folder, so return this folder
    constfold = os.path.join(folder, 'constant')
    if os.path.exists(constfold):
        return folder
    else:
        return ''
    
def parseVTKSeries(folder):
    cf = casefolder(folder)
    bn = os.path.basename(cf)
    seriesfile = os.path.join(cf, 'VTK', bn+'.vtm.series')
    times = []
    if os.path.exists(seriesfile):
        with open(seriesfile, 'r') as f:
            for line in f:
                if 'name' in line:
                    times.append(float(re.split('time\" : | }', line)[1]))
        return times
    else:
        return []

def cleansession():
    Disconnect()
    Connect()

def csvfolder(folder, xpos):
    initialized = False
    if os.path.exists(folder):
        times = parseVTKSeries(folder)
        if len(times)>0:
            for time in tlist:
                xstr = '{:.1f}'.format(-2.412-xpos*1000)
                tstr = str(int(round(time*10)))
                ipfile = os.path.join(folder, "line_t_"+tstr+"_x_"+xstr+".csv")
                    # this is the file that all points for this time will be saved in
                if not os.path.exists(ipfile) or forceoverwrite: 
                    # only run this if the file hasn't been created already or we're being forced to
                    if not initialized: # if paraview hasn't already been initialized, initialize it
                        sv = initializeAll(folder, xpos)
                        sv.times = times
                        initialized = True
                    setTime(time, sv) 
                    
                    if sv.timeKeeper1.Time==time:
                        #raise Exception('Timekeeper broken. Start again.')
                        print('  time ', time, ', file: ', ipfile)
                        exportcsv(folder, ipfile, sv)
            ResetSession()
            cleansession()
            try:
                del sv
            except:
                return
    return

######################################################
####################### SCRIPT #######################
def mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)

for folder in folders:
    csvfolder(folder, xpos)
print('Done exporting csv files')
