#!/usr/bin/env pvpython
'''Functions for collecting line traces through the bath in vtk files'''

import os
import logging

from paraview.simple import * # import the simple module from the paraview

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from paraview_general import *
import folderparser as fp
from interfacemetrics import addUnits

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


#----------------------------------------------------------------

def initializeAll(folder, xpos):
    '''initialize all of paraview'''
    print('Folder ', folder)
    print(' Initializing paraview ')
    sv = stateVars(folder)        # make subdirectories
    sv =  initializeP(sv)        # initialize Paraview
    sv =  initSeries(sv, xpos)        # import the vtk series files  
    return sv


def initSeries(sv:stateVars, xpos:float) -> stateVars:
    '''Initialize the series file and get a line trace'''
    caseVTMSeries = initSeries0(sv)
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


def exportCsv(file:str, sv:stateVars) -> None:
    '''Export the trace to a csv'''
    SaveData(file, proxy=sv.plotOverLine2, PointDataArrays=['U', 'alpha.ink', 'arc_length', 'p', 'p_rgh'], AddTime=1)
    addUnits(file)   
    

def csvfolder(folder:str, xpos:float) -> None:
    '''Export line trace for the folder'''
    initialized = False
    if os.path.exists(folder):
        times = fp.parseVTKSeries(folder)
        if len(times)>0:
            for time in tlist:
                xstr = '{:.1f}'.format(xpos*1000+2.412)
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
                        exportCsv(ipfile, sv)
            ResetSession()
            cleansession()
            try:
                del sv
            except:
                return
    return