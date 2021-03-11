import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from paraview_general import *
from typing import List, Dict, Tuple, Union, Any, TextIO
from vtkmodules.vtkCommonCore import vtkLogger
from folderparser import redoVTKSeriesNoLog, parseVTKSeries
vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_OFF)

FONT = 24

##############################################
######## SCREENSHOTS

####### file operations ######
    
# convert a time to a 3 digit string
def timestr(t:float) -> str:
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr
    return tstr

# get a full file name for an image
def fullFN(t:float, view:str, coloring:str, folder:str):
    fn = 't'+timestr(t)+'_'+view+'_'+coloring+'.png'
    filename = os.path.join(folder, 'images', fn)
    return filename

# get the full file name for the video
def vidFN(view:str, coloring:str, folder:str):
    fn = view+'_'+coloring+'.mp4'
    filename = os.path.join(folder, 'images', fn)
    return filename



######### show/hide/initializes ############
#### import the vtk series file

# initialize paraview, import series file, and create a time stamp
def initializeSV(sv:stateVars) -> None:
    # use this if you already have a stateVars object
    print('Screenshots: Initializing paraview for ', sv.folder)
    ResetSession()
    sv =  initializeP(sv) 
        # initialize Paraview
    sv =  initSeries(sv)
        # import the vtk series files
    sv = timeStamp(sv) 
        # add time stamp
    sv.initialized = True
    return 

        
#### initialize all of paraview
 

def initSeries(sv:stateVars) -> stateVars:
    casevtmseries = initSeries0(sv)  
    casevtmseriesDisplay = Show(casevtmseries, sv.renderView1, 'GeometryRepresentation')
    sv.renderView1.ResetCamera()
    sv.renderView1.Update()
    sv.casevtmseries = casevtmseries
    return sv    
 
# create a timestamp to add later
def timeStamp(sv:stateVars) -> stateVars:
    annotateTimeFilter1 = AnnotateTimeFilter(Input=sv.casevtmseries)
    annotateTimeFilter1.Format = '%2.1f s'
    annotateTimeFilter1Display = Show(annotateTimeFilter1, sv.renderView1, 'TextSourceRepresentation')
    annotateTimeFilter1Display.FontSize = FONT
    sv.timeStamp = annotateTimeFilter1
    return sv


###############
# camera operations

# go to a specific view area
def resetCam(sv:stateVars, time:float) -> None:
    setTime(time, sv)
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)

# go to a specific view point
def setView(st:str, sv:stateVars) -> None:
    if st=="z":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, 0, 10]
        sv.renderView1.CameraViewUp = [0,1,0]
    elif st=="y":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, -10, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="x":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [-10, 0.001, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="a":
        sv.renderView1.CameraFocalPoint = [0.0006, -0.0002, 0.0005]
        sv.renderView1.CameraPosition = [2, 1, 2]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="b":
        sv.renderView1.CameraFocalPoint = [0.001, 0, 0.001]
        sv.renderView1.CameraPosition = [-1,-1,1]
        sv.renderView1.CameraViewUp = [0,0,1]

# go to a specific viewpoint and view area
def setAndUpdate(st:str, sv:stateVars) -> None:
    setView(st, sv)
    sv.renderView1.Update()
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)


#### annotations

# put the color bar in the bottom
def positionCB(ColorBar) -> None:
    ColorBar.AutoOrient = 0
    ColorBar.Orientation = 'Horizontal'
    ColorBar.WindowLocation = 'LowerCenter'
    ColorBar.ScalarBarLength = 0.7
    ColorBar.UseCustomLabels = 1
    ColorBar.TitleFontSize = FONT
    ColorBar.LabelFontSize = FONT
    ColorBar.ComponentTitle = ''

# velocity magnitude color bar
def uColorBar(renderView1, display, umax:float) -> None:
    # set scalar coloring
    #ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction('U')
    uPWF = GetOpacityTransferFunction('U')

    uLUT.RescaleTransferFunction(0.0, umax)
    uPWF.RescaleTransferFunction(0.0, umax)
    if umax>0.015:
        du = 0.005
    else:
        du = 0.002

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = np.arange(0, umax+du, du)
    uLUTColorBar.Title = '|Velocity| (m/s)'
    uLUTColorBar.RangeLabelFormat = '%-#0.3f'
    
# velocity magnitude color bar
def nuColorBar(renderView1, display, color) -> None:
    # set scalar coloring
    

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction(color)
    uPWF = GetOpacityTransferFunction(color)

    # rescale to full range of newt, HB experiments
    uLUT.RescaleTransferFunction(10**-6, 10**2)
    uPWF.RescaleTransferFunction(10**-6, 10**2)
    
    # convert to log space. If we've already converted to log space during creation of another image, skip this
    if not uLUT.UseLogScale==1:
        uLUT.MapControlPointsToLogSpace()
        uLUT.UseLogScale = 1

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = [10**i for i in range(-6, 3)]
    uLUTColorBar.Title = 'Viscosity (MPa*s)'
 #     uLUTColorBar.RangeLabelFormat = '%-#2.0e'
    
# volume fraction color bar
def alphaColorBar(renderView1, display) -> None:
    # set scalar coloring
    ColorBy(display, ('CELLS', 'alpha.ink'))
    ColorBy(display, ('POINTS', 'alpha.ink'))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(renderView1, True)
    alphainkLUT = GetColorTransferFunction('alphaink')
    alphainkLUT.RescaleTransferFunction(0.0, 1)
    alphainkPWF = GetOpacityTransferFunction('alphaink')

    aLUTColorBar = GetScalarBar(alphainkLUT, renderView1)
    positionCB(aLUTColorBar)
    aLUTColorBar.CustomLabels = np.arange(0, 1.25, 0.25)
    aLUTColorBar.Title = 'Volume fraction ink'
    aLUTColorBar.RangeLabelFormat = '%-#0.2f'


################ types of surfaces #########################

# filament surface with velocity coloring
def vol(sv:stateVars) -> None:
    hideAll()
    alphasurface(sv, 'U')

# filament surface with no coloring
def plainalphasurface(sv:stateVars) -> None:
    hideAll()
    alphasurface(sv, 'None')
    
def setDisplayColor(display, color:str):
    try:
        display.SetScaleArray = ['POINTS', color]
        display.OpacityArray = ['POINTS', color]
        ColorBy(display, ('POINTS', color))
    except:
        pass
    
def viscColor(display, visc:str) -> None:
    colordict = {0.000001:[59, 76, 192],\
                 0.00001:[97, 128, 233],\
                 0.0001:[140, 175, 254],\
                 0.001:[183, 207, 249], \
                 0.01:[220, 221, 221], \
                 0.1:[244, 197, 173], \
                 1:[244, 153, 122],\
                 10:[222, 97, 77],\
                 100:[180, 4, 38]}
    try:
        c = colordict[float(visc)]
    except:
        raise Exception('Viscosity not in list')
    c = [i/255 for i in c]
    display.AmbientColor = c
    display.DiffuseColor = c
    
    
def inkclip(sv:stateVars, clipinput, colVar:str, invert:int, **kwargs) -> None:
    clip = Clip(Input=clipinput)
    clip.Scalars = ['POINTS', 'alpha.ink']
    if colVar=='nu1' or colVar=='inknu':
        clip.Value = 0.9
    elif colVar=='nu2' or colVar=='supnu':
        clip.Value = 0.1
    elif 'clipVal' in kwargs:
        clip.Value = kwargs['clipVal']
    else:
        clip.Value = 0.5
    clip.ClipType = 'Scalar'
    clip.Invert = invert

    # show data in view
    clipDisplay = Show(clip, sv.renderView1, 'UnstructuredGridRepresentation')
    clipDisplay.Representation = 'Surface'
    
    if colVar=="U":
        setDisplayColor(clipDisplay, 'U')
        
    elif colVar=='nu1' or colVar=='nu2':
        setDisplayColor(clipDisplay, colVar) 
        nuColorBar(sv.renderView1, clipDisplay, colVar)
    elif colVar=='inknu' or colVar=='supnu':
        if colVar=='inknu':
            viscColor(clipDisplay, sv.inknu)
        else:
            viscColor(clipDisplay, sv.supnu)
    else:
#         clipDisplay.SetScaleArray = 'None'
#         ColorBy(clipDisplay, None)
        c = [0.9, 0.9, 0.9]
        clipDisplay.AmbientColor = c
        clipDisplay.DiffuseColor = c
        
    Hide(clipinput, sv.renderView1)
    
    return clipDisplay

# filament surface
# colVar is the variable being colored ('U' or 'None')
def alphasurface(sv:stateVars, colVar:str) -> None:
    clipDisplay = inkclip(sv, sv.casevtmseries, colVar, 0)
    uColorBar(sv.renderView1, clipDisplay, 0.01)
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate('y', sv)
    sv.renderView1.Update()
    return 

def sliceMake(sv:stateVars, origin:List[float], normal:List[float]):
    slice1 = Slice(Input=sv.casevtmseries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Origin = origin
    slice1.SliceType.Normal = normal # look down the y axis
    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.SetScaleArray = ['POINTS', 'U']
    slice1Display.OpacityArray = ['POINTS', 'U']
    Hide(sv.casevtmseries, sv.renderView1)
    return slice1

def viscSlice(sv:stateVars, origin:List[float], normal:List[float], view:str) -> None:
    out = sv.readVisc()
    if out>0:
        return
    slice1 = sliceMake(sv, origin, normal)
    inkclip(sv, slice1, sv.inkfunc, 0)
    inkclip(sv, slice1, sv.supfunc, 1)

    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def viscy(sv:stateVars) -> None:
    hideAll()
    viscSlice(sv, [0,0,0], [0, -1, 0], 'y')
    
def viscx(sv:stateVars) -> None:
    hideAll()
    viscSlice(sv, [-0.001, 0, 0], [1,0,0], 'x')
    
def uSlice(sv:stateVars, origin:List[float], normal:List[float], view:str) -> None:
    slice1 = sliceMake(sv, origin, normal)
    d1 = inkclip(sv, slice1, 'U', 0, clipVal = 0.9)
    d2 = inkclip(sv, slice1, 'U', 1, clipVal = 0.1)
    for d in [d1, d2]:
        uColorBar(sv.renderView1, d, 0.02)
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def uslicey(sv:stateVars) -> None:
    hideAll()
    uSlice(sv, [0,0,0], [0, -1, 0], 'y')
    
def uslicex(sv:stateVars) -> None:
    hideAll()
    uSlice(sv, [-0.001, 0, 0], [1,0,0], 'x')  
 
    

# cross-section with mesh overlay and volume fraction coloring
def mesh(sv:stateVars) -> None:
    hideAll()
    slice1 = Slice(Input=sv.casevtmseries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Normal = [0, -1, 0] # look down the y axis

    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.SetScaleArray = ['POINTS', 'U']
    slice1Display.OpacityArray = ['POINTS', 'U']

    # show mesh edges
    slice1Display.SetRepresentationType('Surface With Edges')
    setAndUpdate('y', sv)
    alphaColorBar(sv.renderView1, slice1Display)    
    return 
 
# plain slice
def plainslice(sv:stateVars):
    # create a new 'Slice'
    slice2 = Slice(Input=sv.casevtmseries)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice2.SliceType)

    # Properties modified on slice1.SliceType
    slice2.SliceType.Normal = [0.0, -1.0, 0.0]

    # show data in view
    slice2Display = Show(slice2, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.ColorArrayName = [None, '']
    # turn off scalar coloring
    ColorBy(slice2Display, None)

    # change solid color
    slice2Display.AmbientColor = [1.0, 1.0, 1.0]
    slice2Display.DiffuseColor = [1.0, 1.0, 1.0]
    
    setAndUpdate('y', sv)
    
    return slice2
    
def yarrows(slice2, sv):
    # create a new 'Glyph'
    glyph1 = Glyph(Input=slice2, GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'U']
    glyph1.ScaleArray = ['POINTS', 'U']
    glyph1.ScaleFactor = 0.0006029999814927578
    glyph1.GlyphTransform = 'Transform2'
    glyph1.MaximumNumberOfSamplePoints = 2000
    glyph1.ScaleFactor = 0.03

    # show data in view
    glyph1Display = Show(glyph1, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    uColorBar(sv.renderView1, glyph1Display, 0.02)
    setAndUpdate('y', sv)
    
def sliceandyarrows(sv):
    hideAll()
    slice2 =  plainslice(sv)
    yarrows(slice2, sv)
    
def tube(sv):
    streamTracer1 = StreamTracer(Input=sv.casevtmseries, SeedType='High Resolution Line Source')
    streamTracer1.Vectors = ['POINTS', 'U']
    streamTracer1.MaximumStreamlineLength = 0.006
    streamTracer1.SeedType.Point1 = [-0.003014999907463789, -0.0021104998886585236, sv.tubeh]
    streamTracer1.SeedType.Point2 = [0.003014999907463789, 0.0021104998886585236, sv.tubeh]
    streamTracer1.SeedType.Resolution = 50
    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=streamTracer1.SeedType)
    # show data in view
    streamTracer1Display = Show(streamTracer1, sv.renderView1, 'GeometryRepresentation')
    # trace defaults for the display properties.
    streamTracer1Display.ColorArrayName = [None, '']
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    streamTracer1Display.ScaleTransferFunction.Points = [-3.2023908359124786, 0.0, 0.5, 0.0, 3.5858197468415374, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    streamTracer1Display.OpacityTransferFunction.Points = [-3.2023908359124786, 0.0, 0.5, 0.0, 3.5858197468415374, 1.0, 0.5, 0.0]
    # update the view to ensure updated data information
    sv.renderView1.Update()
    # create a new 'Tube'
    tube1 = Tube(Input=streamTracer1)
    tube1.Scalars = ['POINTS', 'AngularVelocity']
    tube1.Vectors = ['POINTS', 'Normals']
    tube1.Radius = 1e-05
    # show data in view
    tube1Display = Show(tube1, sv.renderView1, 'GeometryRepresentation')
    ColorBy(tube1Display, ('POINTS', 'U', 'Magnitude'))

    Hide(streamTracer1, sv.renderView1)
    sv.renderView1.Update()

    uColorBar(sv.renderView1, tube1Display, 0.01)

    # update the view to ensure updated data information
    sv.renderView1.Update()
    
def surfandtube(sv):
    hideAll()
    plainalphasurface(sv)
    tube(sv)
    
    
################ scripting #########################



class ssVars:
    def __init__(self, tag, tlist, **kwargs):
        if tag=='volumes':
            if 'volviewlist' not in kwargs:
                self.vollist = ['y']
            else:
                self.vollist = kwargs['volviewlist']
            self.coloring = 'umag'
            self.function = vol
        elif tag=='meshes':
            self.coloring = 'mesh'
            self.function = mesh
            self.vollist = ['y']
        elif tag=='vectors':
            self.coloring = 'vecs'
            self.function = sliceandyarrows
            self.vollist = ['y']
        elif tag=='tubes':
            if 'volviewlist' not in kwargs:
                self.vollist = ['a']
            else:
                self.vollist = kwargs['volviewlist']
            if 'tubeh' not in kwargs:
                self.tubeh = 0.001
            else:
                self.tubeh = kwargs['tubeh']
            self.coloring = 'stre'
            self.function = surfandtube
        elif tag=='viscy':
            self.coloring='viscy'
            self.function=viscy
            self.vollist = ['y']
        elif tag=='viscx':
            self.coloring = 'viscx'
            self.function = viscx
            self.vollist = ['x']
        elif tag=='uslicey':
            self.coloring = 'uslicey'
            self.function = uslicey
            self.vollist = ['y']
        elif tag=='uslicex':
            self.coloring = 'uslicex'
            self.function = uslicex
            self.vollist = ['x']
        self.tlist = tlist
        self.tag = tag

def runthrough(v:ssVars, sv:stateVars) -> None:
    if v.tag=='tubes':
        sv.tubeh = v.tubeh
    
    if len(v.tlist)>0:
        tlist = [int(round(10*i)) for i in v.tlist]
    else:
        tlist = [int(round(10*i)) for i in sv.times]
    
    tlist2 = []
    fnlist = []
    viewlist = []
    # get a list of files to create
    
    for view in v.vollist:
        for t in tlist:
            fn = fullFN(t, view, v.coloring, sv.folder)
                # find the full file name for this set of circumstances
            if (not os.path.exists(fn)) and t/10 in sv.times:
                # if the file does not exist, add this set to the list
                tlist2.append(t)
                fnlist.append(fn)
                viewlist.append(view)
                
    # if we already have all the files we need, don't run this
    if len(tlist2)==0:
        return sv

    if not sv.initialized:
        try:
            initializeSV(sv)
        except Exception as e:
            print('InitializeSV exception:', e)
            raise

    #print(tlist, v.vollist, v.function)
    try:
        f = v.function
        f(sv)                                  # create the image
        Show(sv.timeStamp, sv.renderView1)     # add time stamp
    except Exception as e:
        print('Create image exception:', e)
        return sv

    setTime(sv.times[-1]/10, sv) # for some reason we need to set this twice, or it will position the imagewrong
    
    # iterate through times and take snapshots of the surfaces
    
    for i,t in enumerate(tlist2):
        print(f'--t {t/10}, view: {viewlist[i]}, file: {fnlist[i]}' )
        setTime(t/10, sv)
        setView(viewlist[i], sv)
        SaveScreenshot(fnlist[i], sv.renderView1)
 
    return sv
         
def folderscript(folder:str, ssvList:List[ssVars]):
    
    try:
        if not os.path.exists(folder):
            return
        sv = stateVars(folder)
        sv.times = readTimes(folder)
        # sv.times = parseVTKSeries(folder)
        # if len(sv.times)==0:
        #     generateVTKseries(folder, False)
        #     sv.times = parseVTKSeries(folder)
        
        for ssv in ssvList:
            sv = runthrough(ssv, sv)
    except Exception as e:
        print(e)
        return
    cleansession()
    return    
