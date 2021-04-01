import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from paraview_general import *
from typing import List, Dict, Tuple, Union, Any, TextIO

FONT = 24

##############################################
######## SCREENSHOTS

# # this holds objects related to the paraview display
# class stateVars():
#     def __init__(self, folder):
#         mksubdirs(folder)
#         self.folder = folder        
#         self.casevtmseries = ""
#         self.renderView1 = ""
#         self.animationScene1 = ""
#         self.timeKeeper1 = ""
#         self.timeStamp = ""
#         self.times = ""
#         self.initialized = False


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
    ColorBy(display, ('POINTS', 'U', 'Magnitude'))

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
    alphasurface(sv, 'U')

# filament surface with no coloring
def plainalphasurface(sv:stateVars) -> None:
    alphasurface(sv, 'None')

# filament surface
# color is the variable being colored ('U' or 'None')
def alphasurface(sv:stateVars, color:str) -> Tuple:
    # sv = stateVars object
    # color = 'U' or 'None'

    clip1 = Clip(Input=sv.casevtmseries)
    clip1.Scalars = ['POINTS', 'alpha.ink']
    clip1.Value = 0.5
    clip1.ClipType = 'Scalar'
    clip1.Invert = 0

    # show data in view
    clip1Display = Show(clip1, sv.renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = [None, '']
    clip1Display.OSPRayScaleArray = 'U'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.ScaleFactor = 0.00018132181139662862
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'None'
    clip1Display.GaussianRadius = 9.066090569831431e-06
    
    if color=="U":
        clip1Display.SetScaleArray = ['POINTS', 'U']
        clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
        clip1Display.OpacityArray = ['POINTS', 'U']
        clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
        clip1Display.DataAxesGrid = 'GridAxesRepresentation'
        clip1Display.PolarAxes = 'PolarAxesRepresentation'
        clip1Display.ScalarOpacityUnitDistance = 5.84412949766512e-05
        clip1Display.ExtractedBlockIndex = 1
        clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
        clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
        uColorBar(sv.renderView1, clip1Display, 0.01)
    else:
        clip1Display.SetScaleArray = 'None'
        ColorBy(clip1Display, None)
        
    Hide(sv.casevtmseries, sv.renderView1)

    resetCam(sv, (sv.times[-1])) 
    setAndUpdate('y', sv)
    sv.renderView1.Update()
    return [clip1, clip1Display]


def viscSlice(sv:stateVars) -> Tuple:
    out = sv.readVisc()
    if out>0:
        
    
    

# cross-section with mesh overlay and volume fraction coloring
def mesh(sv:stateVars) -> Tuple:
    hideAll()
    slice1 = Slice(Input=sv.casevtmseries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Normal = [0, -1, 0] # look down the y axis

    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.OSPRayScaleArray = 'U'
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.0006029999814927578
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 3.0149999074637892e-05
    slice1Display.SetScaleArray = ['POINTS', 'U']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = ['POINTS', 'U']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'
    slice1Display.ScaleTransferFunction.Points = [-0.003019144991412759, 0.0, 0.5, 0.0, 0.010859496891498566, 1.0, 0.5, 0.0]
    slice1Display.OpacityTransferFunction.Points = [-0.003019144991412759, 0.0, 0.5, 0.0, 0.010859496891498566, 1.0, 0.5, 0.0]

    # show mesh edges
    slice1Display.SetRepresentationType('Surface With Edges')
    setAndUpdate('y', sv)
    alphaColorBar(sv.renderView1, slice1Display)    
    return [slice1, slice1Display]
 
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
    slice2Display.Representation = 'Surface'
    slice2Display.ColorArrayName = [None, '']
    slice2Display.OSPRayScaleArray = 'U'
    slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice2Display.SelectOrientationVectors = 'None'
    slice2Display.ScaleFactor = 0.0006029999814927578
    slice2Display.SelectScaleArray = 'None'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'None'
    slice2Display.GaussianRadius = 3.0149999074637892e-05
    slice2Display.SetScaleArray = ['POINTS', 'U']
    slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice2Display.OpacityArray = ['POINTS', 'U']
    slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice2Display.DataAxesGrid = 'GridAxesRepresentation'
    slice2Display.PolarAxes = 'PolarAxesRepresentation'
    slice2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    slice2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
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
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'U'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.03
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 3.018014947883785e-05
    glyph1Display.SetScaleArray = ['POINTS', 'U']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'U']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'
    glyph1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    glyph1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]   
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
    streamTracer1.SeedType.Point1 = [-0.003014999907463789, -0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Point2 = [0.003014999907463789, 0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Resolution = 50

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=streamTracer1.SeedType)

    # show data in view
    streamTracer1Display = Show(streamTracer1, sv.renderView1, 'GeometryRepresentation')


    # trace defaults for the display properties.
    streamTracer1Display.Representation = 'Surface'
    streamTracer1Display.ColorArrayName = [None, '']
    streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
    streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    streamTracer1Display.SelectOrientationVectors = 'Normals'
    streamTracer1Display.ScaleFactor = 0.001738997269421816
    streamTracer1Display.SelectScaleArray = 'AngularVelocity'
    streamTracer1Display.GlyphType = 'Arrow'
    streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
    streamTracer1Display.GaussianRadius = 8.69498634710908e-05
    streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
    streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
    streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
    streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
    streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
    streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

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

    # trace defaults for the display properties.
    tube1Display.Representation = 'Surface'
    tube1Display.ColorArrayName = [None, '']
    tube1Display.OSPRayScaleArray = 'AngularVelocity'
    tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tube1Display.SelectOrientationVectors = 'Normals'
    tube1Display.ScaleFactor = 0.001094136107712984
    tube1Display.SelectScaleArray = 'AngularVelocity'
    tube1Display.GlyphType = 'Arrow'
    tube1Display.GlyphTableIndexArray = 'AngularVelocity'
    tube1Display.GaussianRadius = 5.4706805385649205e-05
    tube1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
    tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tube1Display.OpacityArray = ['POINTS', 'AngularVelocity']
    tube1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tube1Display.DataAxesGrid = 'GridAxesRepresentation'
    tube1Display.PolarAxes = 'PolarAxesRepresentation'
    tube1Display.ScaleTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]
    tube1Display.OpacityTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]

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

        

def runthrough(tlist, viewlistin, coloring, sv, f, mode):
    # tlist is the list of times to go through
    # viewlist is the list of views to go through
    # coloring is the coloring label
    # subfolder is 'frames' or 'videos'
    # sv is a stateVars object
    # f is a function to create the rendering
    tlist2 = []
    fnlist = []
    viewlist = []

    for view in viewlistin: # for each view point
        for t in tlist: # tlist is *10
            # go through each time step
            fn = fullFN(t, view, coloring, sv.folder)
                # find the full file name for this set of circumstances
            if (not os.path.exists(fn)) and t/10 in sv.times:
                # if the file does not exist, add this set to the list
                tlist2.append(t)
                fnlist.append(fn)
                viewlist.append(view)
    
    # if we already have all the files we need, don't run this
    if len(tlist2)>0:
        if not sv.initialized:
            try:
                initializeSV(sv)
            except Exception as e:
                print(e)
                raise

        try:
            f(sv)
            Show(sv.timeStamp, sv.renderView1)     # add time stamp

    #         # iterate through times and take snapshots of the surfaces
            for [i,t] in enumerate(tlist2):
                # if t-round(t)==0:
                print(f't = {t}, view = {viewlist[i]}, filename = {fnlist[i]}' )
                setTime(t/10, sv)
                setView(viewlist[i], sv)
                SaveScreenshot(fnlist[i], sv.renderView1)
        except Exception as e:
            print(e)
            raise
    return sv
 


def folderscript(folder, modes, volumes, volviewlist, meshes, vectors, tubes, tlist):
    
    try:
        if not os.path.exists(folder):
            return
        sv = stateVars(folder)
        sv.times = parseVTKSeries(folder)
        for mode in modes:
            if mode==0:
                trange = [int(round(10*i)) for i in tlist]
            else:
                trange = [int(round(10*i)) for i in sv.times]

            coloringlist = []
            flist = []
            vollist = []

            # ink-support interfaces
            if volumes:
                coloringlist.append('umag')
                flist.append(vol)
                vollist.append(volviewlist)

            # mesh from the y view
            if meshes:
                coloringlist.append('mesh')
                flist.append(mesh)
                vollist.append(['y'])

            # flow vectors from the y view
            if vectors:
                coloringlist.append('vecs')
                flist.append(sliceandyarrows)
                vollist.append(['y'])

            # streamlines with interface
            if tubes:
                coloringlist.append('stre')
                flist.append(surfandtube)
                vollist.append(volviewlist)

            for i in range(len(coloringlist)):
                try:
                    sv = runthrough(trange, vollist[i], coloringlist[i], sv, flist[i], mode)
                except Exception as e:
                    print(e)
                    return
        cleansession()
        return
    except Exception as e:
        print(e)
        cleansession()
        return