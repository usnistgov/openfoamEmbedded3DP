#!/usr/bin/env pvpython
'''Functions for collecting line traces through the bath in vtk files'''

# external packages
import os
import logging
import csv
from paraview.simple import * # import the simple module from the paraview

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(currentdir)
sys.path.append(parentdir)
from paraview_general import *
import folderparser as fp
from pythonscripts.pvCleanup import addUnits

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"


#----------------------------------------------------------------

def initializeAll(folder:str, constDir:str, pos:float):
    '''initialize all of paraview. constDir is direction that has const value ('x' or 'z'). pos is the value that stays constant, in x or z (e.g. -0.001)'''
    print('Folder ', folder)
    print(' Initializing paraview ')
    sv = stateVars(folder)        # make subdirectories
    sv =  initializeP(sv)        # initialize Paraview
    le = fp.legendUnique(folder)

    if constDir=='x':
        btc = float(le['bath_top_coord'])/1000 * 0.97
        bbc = float(le['bath_bottom_coord'])/1000 * 0.97
        sv = initSeries(sv, x0=pos, xf=pos, z0=btc, zf=bbc)       # import the vtk series files  
    elif constDir=='z':
        blc = float(le['bath_left_coord'])/1000 * 0.97
        brc = float(le['bath_right_coord'])/1000 * 0.97
        sv = initSeries(sv, x0=blc, xf=brc, z0=pos, zf=pos)
    return sv


def initSeries(sv:stateVars, x0:float=-0.004, xf:float=0.004, z0:float=0.0021104998886585236, zf:float=-0.0021104998886585236) -> stateVars:
    '''Initialize the series file and get a line trace'''
    caseVTMSeries = initSeries0(sv)
    sv.caseVTMSeries = caseVTMSeries
    # sv.times = caseVTMSeries.TimestepValues
    shearRate = computeShearRate(sv)
    plotOverLine2 = PlotOverLine(Input=shearRate, Source='High Resolution Line Source')
    # plotOverLine2.Source.Point1 = [0.002588, 0, 0.0021104998886585236]
    # plotOverLine2.Source.Point2 = [0.002588, 0, -0.0021104998886585236]
    
    plotOverLine2.Source.Point1 = [x0, 0, z0]
    plotOverLine2.Source.Point2 = [xf, 0, zf]

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
    Hide(caseVTMSeries, sv.renderView1)
    sv.caseVTMSeries = caseVTMSeries
    return sv


def exportCsv(file:str, sv:stateVars) -> None:
    '''Export the trace to a csv'''
    SaveData(file, proxy=sv.plotOverLine2, PointDataArrays=['U', 'alpha.ink', 'arc_length', 'p', 'p_rgh'], AddTime=1)
    addUnits(file)   
    
def exportEmpty(file:str) -> None:
    '''Export an empty line trace to a csv. This is useful if the simulation didn't produce any results for the line scan, but we don't want to re-check it every time we run the script.'''
    
    with open(file, mode='w') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['time', 'vx', 'vy', 'vz', 'alpha', 'p', 'rAU', 'vtkValidPointMask', 'arc_length', 'x', 'y', 'z'])
        writer.writerow(['s', 'm/s', 'm/s', 'm/s', '', 'kg/(m*s^2)', 'm^3*s/kg', '', 'm', 'm', 'm', 'm'])
        
def convertToRelative(x:float, folder:str) -> float:
    '''convert an absolute x position in m to its position relative to the nozzle in mm'''
    le = fp.legendUnique(folder)
    return x*1000+float(le['nozzle_center_x_coord'])

def csvfolder(folder:str, constDir:str, pos:float, tlist:List[float], forceOverwrite:bool=False) -> None:
    '''Export line trace for the folder. 
    constDir is the direction that is constant, either 'x' or 'z'
    pos is the constant x or z position of the trace, in absolute coordinates.
    tlist is the list to times at which to take the trace
    forceOverwrite to overwrite existing files'''
    initialized = False
    if os.path.exists(folder):
        times = fp.parseVTKSeries(folder)
        if len(times)>0:
            for time in tlist:
                if time in times:
                    # if constDir=='x':
                    #     xstr = '{:.1f}'.format(convertToRelative(pos, folder))
                    # else:
                    #     xstr = '{:.1f}'.format(1000*pos)
                    xstr = '{:.1f}'.format(pos)
                    tstr = str(int(round(time*10)))
                    ipfile = os.path.join(folder, f"line_t_{tstr}_{constDir}_{xstr}_di.csv")
                        # this is the file that all points for this time will be saved in
                    if not os.path.exists(ipfile) or forceOverwrite: 
                        # only run this if the file hasn't been created already or we're being forced to
                        if not initialized: # if paraview hasn't already been initialized, initialize it
                            le = fp.legendUnique(folder)
                            posabs = (float(le['nozzle_center_x_coord']) + pos*float(le['nozzle_inner_width']))/1000
                            sv = initializeAll(folder, constDir, posabs)
                            sv.times = times
                            initialized = True
                        setTime(time, sv) 

                        if sv.timeKeeper1.Time==time:
                            #raise Exception('Timekeeper broken. Start again.')
                            print('  time ', time, ', file: ', ipfile)
                            exportCsv(ipfile, sv)
            ResetSession()
            cleanSession()
            try:
                del sv
            except:
                return
    return