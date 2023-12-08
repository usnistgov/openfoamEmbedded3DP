#!/usr/bin/env python
'''Functions for generating legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import os
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket, sys
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import file.file_handling as fh
from file.backwards_read import fileReadBackwards
from file.file_export import *
from scrape_tools import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)



#-------------------------------------------------------------------------------------------------  



class scrape:
    '''the scrape class is used to store variables scraped from the log files
        scrape contains a function scrape.table(), which converts all the data in the class into a 2 column table'''
        
    def __init__(self, folder:str):
        '''folder is a folder name (string) which should point either to the case folder or one level above it'''
        
        self.fold = folder # full path
        self.fh = fh.folderHandler(folder)
        self.meshfold = self.fh.meshFolder()
        self.casefold = self.fh.caseFolder()
        self.folder = ['folder', os.path.basename(folder)] # just the folder name
        self.compareto = ['compare_to', '', '']

        self.initRates()
        self.initGeo()
                        
        # these lists store values ripped from the dictionaries
        self.initMesh()
        self.initDMD()                
        self.initTransport()
        self.initControl()                
        self.initFV()
    
    #-------------------------
    def initRates(self):
        '''initial (empty) values for simulation rates'''
        
        # these values are ripped from log files
        self.shmtimes = ['snappyHexMesh_time', '', 's'] 
            # how long it takes snappyHexMesh to run
        self.shmtimem = ['snappyHexMesh_time', '', 'min']
        self.iftimes = ['interFoam_time_s', '', 's'] 
            # how long it takes interfoam to run
        self.iftimehr = ['interFoam_time_hr', '', 'hr']
        self.simTime = ['simulation_time', '', 's'] 
            # how many seconds within the simulation we ran
        self.simrate = ['simulation_rate', '',  'hr/s'] 
            # how many real hours it takes to run each second within the simulation
        self.version = ['openfoam_version', '', '']

        
        
    #-------------------------
    def initGeo(self):
        '''these values measure the geometry of the system and are stored in geometry.csv, which is created by ncreate3d.py/noz3dscript.ipynb'''
        self.GEOniw = ['nozzle_inner_width', '0.603', 'mm']
        self.GEOnt = ['nozzle_thickness', '0.152', 'mm']
        self.GEObw = ['bath_width', '', 'mm']
        self.GEObd = ['bath_depth', '', 'mm']
        self.GEOnl = ['nozzle_length', '', 'mm']
        self.GEOblc = ['bath_left_coord', '', 'mm'] # x coord
        self.GEObrc = ['bath_right_coord', '', 'mm'] # x coord
        self.GEObfc = ['bath_front_coord', '', 'mm'] # y
        self.GEObbackc = ['bath_back_coord', '', 'mm'] # y
        self.GEObbotc = ['bath_bottom_coord', '', 'mm'] # z
        self.GEObtc = ['bath_top_coord', '', 'mm'] # z
        self.GEOnbc = ['nozzle_bottom_coord', '', 'mm'] # z
        self.GEOncxc = ['nozzle_center_x_coord', '', 'mm']
        self.GEOncyc = ['nozzle_center_y_coord', '', 'mm']
        self.GEOna = ['nozzle_angle', '0', 'degrees'] # RG
        self.GEOhoriz = ['horizontal', False, '']
        self.GEOadj = ['adjacent_filament_orientation', '', '']
        self.GEOdst = ['adjacent_filament_offset', '', 'mm']
        self.GEObathv = ['bath_velocity', '', 'm/s']
        self.GEOinkv = ['ink_velocity', '', 'm/s']
    
    #-------------------------
    def initMesh(self):
        '''initial values for mesh variables'''
        
        self.SHMlist = [['castellatedMesh', '', ''],\
                        ['snap', '', ''],\
                        ['addLayers', '', '']]
        self.CMClist = [['maxLocalCells', '', ''], \
                        ['maxGlobalCells', '', ''], \
                        ['minRefinementCells', '', ''], \
                        ['nCellsBetweenLevels', '', ''],\
                        ['resolveFeatureAngle', '', ''],\
                        ['locationInMesh', '', ''], \
                        ['allowFreeStandingZoneFaces', '', '']]
        self.CMCfixedWallsLevel = ['fixedWalls_level', '', '']
        self.SClist = [['nSmoothPatch', '', ''], \
                       ['tolerance', '', ''], \
                       ['nSolveIter', '', ''], \
                       ['nRelaxIter', '', ''],\
                       ['nFeatureSnapIter', '', ''], \
                       ['implicitFeatureSnap', '', ''],\
                       ['explicitFeatureSnap', '', ''], \
                       ['multiRegionFeatureSnap', '', '']]
        self.ALClist = [['relativeSizes', '', ''], \
                        ['expansionRatio', '', ''], \
                        ['finalLayerThickness', '', ''],\
                        ['minThickness', '', ''], \
                        ['nGrow', '', ''], \
                        ['featureAngle', '', ''], \
                        ['nRelaxIter', '', ''], \
                        ['nSmoothSurfaceNormals', '', ''], \
                        ['nSmoothNormals', '', ''], \
                        ['nSmoothThickness', '', ''],\
                        ['maxFaceThicknessRatio', '', ''],\
                        ['maxThicknessToMedialRatio', '', ''],\
                        ['minMedialAxisAngle', '', ''],\
                        ['nBufferCellsNoExtrude', '', ''], \
                        ['nLayerIter', '', ''], \
                        ['nRelaxedIter', '', '']]
        self.ALCfixedWallsLayers = ['fixedWalls_nSurfaceLayers', '', '']
        self.MQClist = [['nSmoothScale', '', ''], \
                        ['errorReduction', '', ''], \
                        ['maxNonOrtho', '']]
        self.SHMmergeTolerance = ['mergeTolerance', '', '']
                        
        self.blocksdims = ['blocks_dims', '', '']
    
    #-------------------------
    def initDMD(self):
        '''initial values for dynamic mesh dictionary'''
        self.DMDlist = [['refineInterval', '', 'steps'], \
                        ['lowerRefineLevel', '', ''], \
                        ['upperRefineLevel', '', ''], \
                        ['unrefineLevel', '', 'cells'],\
                        ['nBufferLayers', '', 'cells'], \
                        ['maxRefinement', '', 'cells'], \
                        ['maxCells', '', 'cells']]
    
    #-------------------------
    def initTransport(self):
        '''initial values for transport model'''
        self.TPinklist = [['ink', '', ''],\
                        ['ink_transportModel', '', ''],\
                          ['ink_nu', '', 'm^2/s'], \
                          ['ink_nu0', '', 'm^2/s'], \
                          ['ink_tau0', '', 'm^2/s^2'], \
                          ['ink_k', '', 'm^2*s^(n-2)'],\
                          ['ink_n', '', ''], \
                          ['ink_rho', '', 'kg/m^3']]
        self.TPsuplist = [['sup', '', ''],\
                        ['sup_transportModel', '', ''],\
                          ['sup_nu', '', 'm^2/s'], \
                          ['sup_nu0', '', 'm^2/s'], \
                          ['sup_tau0', '', 'm^2/s^2'], \
                          ['sup_k', '', 'm^2*s^(n-2)'],\
                          ['sup_n', '', ''], \
                          ['sup_rho', '', 'kg/m^3']]
        self.inkLabel = ''
        self.supLabel = ''
        self.TPsigma = ['sigma', '', 'J/m^2']
    
    #-------------------------
    def initControl(self):
        '''initial values for controlDict'''
        self.controlDictList =  [['application', '', ''],\
                                 ['endTime', '', 's'], \
                                 ['deltaT', '', 's'], \
                                 ['adjustTimeStep', '', ''], \
                                 ['maxCo', '', ''],\
                                 ['maxAlphaCo', '', ''],\
                                 ['maxDeltaT', '', 's'], \
                                 ['writeControl', '', ''], \
                                 ['writeInterval', '', 's'], \
                                 ['purgeWrite', '', ''],\
                                 ['runTimeModifiable', '', '']]
    
    #-------------------------
    def initFV(self):
        '''initial values for fvSolution, fvSchemes'''
        self.fvSailist = [['nAlphaSubCycles', '', ''],\
                          ['cAlpha', '', ''], \
                          ['nAlphaCorr', '', ''], \
                          ['MULESCorr', '', ''], \
                          ['nLimiterIter', '', ''], \
                          ['solver', '', ''], \
                          ['smoother', '', ''],\
                          ['tolerance', '', ''],\
                          ['relTol', '', ''], \
                          ['isofaceTol', '', ''], \
                          ['surfCellTol', '', ''], \
                          ['nAlphaBounds', '', ''], \
                          ['snapTol', '', ''],\
                          [ 'clip', '', '']]                
        self.fvSpcorrlist = [['solver', '', ''], \
                             ['preconditioner', '', ''], \
                             ['tolerance', '', ''], \
                             ['relTol', '', '']]
        self.fvSprghlist = [['solver', '', ''], \
                            ['preconditioner', '', ''],\
                            ['tolerance', '', ''],\
                            ['relTol', '', '']]
        self.fvSprghfinallist = [['relTol', '', '']]
        self.fvSUlist = [['solver', '', ''],\
                         ['smoother', '', ''], \
                         ['tolerance', '', ''], \
                         ['relTol', '', '']]
        self.fvSPIMPLElist = [['momentumPredictor', '', ''], \
                              ['nOuterCorrectors', '', ''], \
                              ['nCorrectors', '', ''], \
                              ['nNonOrthogonalCorrectors', '', '']]
        
        
    #-------------------------
    
    def scrapeAll(self):
        self.scrapeGeo()
        self.scrapeBlockMeshDict()
        self.scrapeSetFieldsDict()
        self.scrapeU()
        self.scrapeSHM()
        self.scrapeDMD()
        self.scrapeTP()
        self.scrapeLabels()
        self.scrapeCD()
        self.scrapeFV()

    def table(self) -> List[List[str]]:
        '''this function takes all the values we ripped and stored as variables in the scrape object and turns them into a table with 2 columns: the variable name and the variable value'''
        col = [] # start with an empty table, then add to it sequentially
        for i in [self.folder, self.version, self.compareto, self.shmtimes, \
                  self.shmtimem, self.iftimes, self.iftimehr, \
                  self.simTime, self.simrate]: # RG
            col.append(i)
        ca(col, ['', '', 'mesh','GEOMETRY'],\
           [self.GEOniw, self.GEOnt, self.GEObw, self.GEObd, \
            self.GEOnl, self.GEOblc, self.GEObrc, self.GEObfc,\
            self.GEObbackc, self.GEObbotc, self.GEObtc, self.GEOnbc,\
            self.GEOncxc, self.GEOncyc, self.GEOna, self.GEOhoriz,\
            self.GEOadj, self.GEOdst, self.GEObathv, self.GEOinkv]) # RG
        ca(col, ['', 'SYSTEM'], [])    
        ca(col, ['snappyHexMeshDict'], self.SHMlist)
        ca(col, ['castellatedMeshControls'], self.CMClist)
        col.append(self.CMCfixedWallsLevel)
        ca(col, ['snapControls'], self.SClist)
        ca(col, ['addLayersControls'], self.ALClist)
        col.append(self.ALCfixedWallsLayers)
        ca(col, ['meshQualityControls'], self.MQClist)
        col.append(self.SHMmergeTolerance)
        ca(col, ['', 'blockMeshDict'], [self.blocksdims])
        ca(col, ['', 'case', 'CONSTANT', 'dynamicMeshDict'], self.DMDlist)
        ca(col, ['', 'transportProperties'], self.TPinklist)
        ca(col, [], self.TPsuplist)
        col.append(self.TPsigma)
        ca(col, ['', 'SYSTEM', 'controlDict'], self.controlDictList)
        ca(col, ['', 'fvSolution', 'alpha.ink'], self.fvSailist)
        ca(col, ['', 'pcorr'], self.fvSpcorrlist)
        ca(col, ['', 'p_rgh'], self.fvSprghlist)
        ca(col, ['', 'p_rghFinal'], self.fvSprghfinallist)
        ca(col, ['', 'U'], self.fvSUlist)
        ca(col, ['', 'PIMPLE'], self.fvSPIMPLElist)
        return col

    
    #---------------------------------------
    
    def scrapeSHMLog(self) -> None:
        '''scrape the snappyHexMesh log file
        s is a scrape object'''
        shmlog = os.path.join(self.meshfold, 'log_snappyHexMesh')
        if os.path.exists(shmlog):
            # it's useful to go backwards here because we're just looking for the time reported at the end of the file
            for line in fileReadBackwards(shmlog):
                if line.startswith('Finished meshing'):
                    strs = re.split('Finished meshing in = | s', line)
                    shmtime = float(strs[1])
                    self.shmtimes[1] = "%.2f" % shmtime
                    self.shmtimem[1] = "%.2f" % (shmtime/60)
                    return 



    def scrapeIFLog(self) -> None:
        '''scrape the interFoam log file
            Because interFoam can take hours to days, sometimes runs get split into pieces. 
            Each interFoam run adds onto the existing log file. 
            s is a scrape object  '''
        ifLog = os.path.join(self.casefold, 'log_interFoam')
        if not os.path.exists(ifLog):
            ifLog = os.path.join(self.fold, 'log_interFoam')
        if not os.path.exists(ifLog):
            return
        iftime = 0 # this variable adds up all the times for separate interFoam runs
        waitfortop = False # this variable tells us what to look for, so we only collect the last reported time from each run
        simtime = 0
        # it's useful to go backwards here because we're just looking for the time reported at the end of the run
        
        # extract OpenFOAM version number
        with open(ifLog) as f:
            line = ''
            lii = 0
            while not 'Version' in line and lii<6:
                line = f.readline()
            if 'Version' in line:
                spl = re.split('Version:| ', line)
                vers = spl[-1]
                try:
                    vers = int(vers)
                except:
                    pass
                self.version[1] = vers
        
        # extract simulation time
        for line in fileReadBackwards(ifLog):
            if simtime==0 and line.startswith('Time = '):
                strs = re.split('Time = ', line)
                simtime = float(strs[1])
            if (waitfortop and line.startswith('fileModificationChecking')): 
                # when we hit fileModificationChecking, it's the end of the run, so now we should look for the next reported time
                waitfortop = False
            if (not waitfortop and line.startswith('ExecutionTime')):
                # read the last reported ExecutionTime
                strs = re.split('ExecutionTime = | s', line)
                iftime+=float(strs[1])
                waitfortop = True # now that we've read the time, look for the next end of run
        # store the time in the scrape object in seconds and hr
        self.iftimes[1] = "%.2f" % iftime 
        self.iftimehr[1] = "%.2f" % (iftime/60/60)
        try:
            if simtime==0:
                simtime = "%.2f" % float(self.simTime[1])
            if simtime>0:
                # if we've already measured a simulation time, calculate the simulation speed
                self.simrate[1] = "%.2f" % (iftime/60/60/simtime)
        except:
            pass
        return

    def scrapeRunTime(self) -> None:
        '''Get the simulation time from the folder'''
        ti = self.fh.times()
        if len(ti)>0:
            self.simTime[1] = str(max(ti))
        else:
            self.simTime[1] = '0'
        return


    def scrapeLogs(self) -> None:
        '''scrape all of the times (run time, simulation time, etc.)
        s is a scrape object'''
        self.scrapeRunTime()
        self.scrapeSHMLog()
        self.scrapeIFLog()


    def scrapeBlockMeshDict(self) -> None:
        '''scrape blockMeshDict
        s is a scrape object'''
        bm = os.path.join(self.meshfold, 'system', 'blockMeshDict')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                while not line.startswith('vertices'):
                    line = f.readline()
                # now we have reached the list of vertices
                # because we establish a basic mesh with blockmeshDict and refine with snappyHexMesh, this list only contains 8 vertices
                for i in range(2):
                    line = f.readline() # read coords from the first list of points
                strs = re.split('\(|\)| ', line)
                self.GEOblc[1] = strs[1] # bath left coord
                self.GEObbackc[1] = strs[2] # bath back coord
                self.GEObbotc[1] = strs[3] # bath bottom coord
                for i in range(7):
                    line = f.readline() # read coords from the last list of points
                strs = re.split('\(|\)| ', line)
                self.GEObrc[1] = strs[1] # bath right coord
                self.GEObfc[1] = strs[2] # bath front coord
                self.GEObtc[1] = strs[3] # bath top coord
                self.GEObw[1] = str(float(self.GEObrc[1]) - float(self.GEOblc[1])) # bath width
                self.GEObd[1] = str(float(self.GEObfc[1]) - float(self.GEObbackc[1])) # bath depth
                while not line.startswith('blocks'):
                    line = f.readline()    
                # now we have reached the list of blocks
                # again, there is only one block because we're using snappyHexMesh
                for i in range(2):
                    line = f.readline() # read coords from the first list of blocks
                strs = re.split('\) \(|\) s', line)
                self.blocksdims[1] = strs[1] # number of cells in the blocks: this should look like (# # #)
                
                return


    def scrapeSetFieldsDict(self) -> None:
        '''scrape setFieldsDict
        s is a scrape object'''
        bm = os.path.join(self.casefold, 'system', 'setFieldsDict')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                while not line.startswith('\t\tp1'):
                    line = f.readline()
                # now we have reached the bottom point of the nozzle
                strs = re.split('\(|\)| ', line) # RG
                while '' in strs:
                    strs.remove('')
                self.GEOncxc[1] = str(1000*float(strs[1])) # nozzle center x
                self.GEOncyc[1] = str(1000*float(strs[2])) # nozzle center y
                self.GEOnbc[1] = (1000*float(strs[3])) # nozzle bottom 
                try:
                    btc = float(self.GEObtc[1]) # bath top coord
                except:
                    self.GEOnl[1] = ""
                else:
                    self.GEOnl[1] = str(btc - self.GEOnbc[1]) # nozzle length
                self.GEOnbc[1] = str(self.GEOnbc[1])
                return


    def scrapeU(self) -> None:
        '''scrape 0/U
        s is a scrape object'''
        bm = os.path.join(self.casefold, '0', 'U')
        
        if os.path.exists(bm):
            with open (bm, "r") as f:
                # get geometry
                if len(self.GEOnl[1])==0 or len(self.GEOna[1])==0:
                    self.scrapeGeo()
                horiz = self.GEOhoriz[1]
                if type(horiz) is str:
                    if 'true' in horiz.lower():
                        horiz = True
                    else:
                        horiz = False
                
                # read bath speed
                line = f.readline()
                while not line.startswith('\tbathFlow'):
                    line = f.readline()
                # now we have reached the bath flow section
                for i in range(3):
                    line = f.readline()
                if horiz:
                    strs = re.split('\(|\)| |-', line)
                    self.GEObathv[1] = strs[5] # bath velocity: this will be reported in m/s
                else:
                    strs = re.split('\(|\)| ', line)
                    self.GEObathv[1] = strs[2] # bath velocity: this will be reported in m/s
                
                # read ink speed
                while not line.startswith('\tinkFlow'):
                    line = f.readline()
                # now we have reached the ink flow section
                for i in range(3):
                    line = f.readline()
                strs = re.split('\(|\)| |-', line)
                
                theta = np.radians(float(self.GEOna[1]))
                vt = strs[5]
                if theta==0:
                    self.GEOinkv[1] = vt
                else:
                    ri = float(self.GEOniw[1])/2
                    l = float(self.GEOnl[1])
                    self.GEOinkv[1] = float(vt)*(ri+l*np.tan(theta))**2/ri**2 # ink velocity: this will be reported in m/s
                return
        else:
            logging.info(f'path {bm} does not exist')


    def scrapeSHM(self) -> None:
        '''scrape snappyHexMeshDict
        here we use listLevel because snappyHexMeshDict has a lot of sections and we're scraping all of the data
        to change which fields we collect, go back to the scrape class definition
        we're collecting data for SHMlist, CMClist, CMCfixedWallsLevel, SClist, ALClist, ALCfixedWallsLayers, MQClist, and SHMmergeTolerance
        s is a scrape object'''
        bm = os.path.join(self.meshfold, 'system', 'snappyHexMeshDict')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                # first determine if castellatedMesh, snapping, and layers are active
                line = listLevel('castellatedMesh', 'geometry', line, f, self.SHMlist)
                # then read in castellated mesh controls
                if self.SHMlist[0][1]=='true':
                    # only collect castellatedMesh variables if we're using a castellated mesh
                    line = listLevel('castellatedMeshControls', '\tfeatures', line, f, self.CMClist)
                    while not line.startswith('\t\t\tfile\t\"fixed'):
                        line = f.readline()
                    line = f.readline()
                    strs = re.split(';|\t', line)
                    self.CMCfixedWallsLevel[1] = strs[4]
                if self.SHMlist[1][1]=='true':
                    # only collect snap variables if we're using snapping
                    line = listLevel('snapControls', 'addLayersControls', line, f, self.SClist)
                if self.SHMlist[2][1]=='true':
                    # only collect layers variables if we're using layers
                    listLevel('addLayersControls', '\tlayers', line, f, self.ALClist)
                    while not line.startswith('\t\tfixed'):
                        line = f.readline()
                    for i in range(2):
                        line = f.readline()
                    strs = re.split(';|\t', line)
                    self.ALCfixedWallsLayers[1] = strs[4]
                line = listLevel('meshQualityControls', 'writeFlags', line, f, self.MQClist)
                line = readLevel0('mergeTolerance', line, f, self.SHMmergeTolerance)
                return
        else:
            logging.warning(f'path {bm} does not exist')


    def scrapeDMD(self) -> None:
        '''scrape dynamicMeshDict
        s is a scrape object'''
        bm = os.path.join(self.casefold, 'constant', 'dynamicMeshDict')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                # just collect the variables in the dynamicRefineFvMeshCoeffs section
                line = listLevel('dynamicRefineFvMeshCoeffs', '\tcorrectFluxes', line, f, self.DMDlist)
                return
        else:
            logging.warning(f'path {bm} does not exist')


    def scrapeTP(self) -> None:
        '''scrape transportProperties
        s is a scrape object'''
        bm = os.path.join(self.casefold, 'constant', 'transportProperties')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                line = listLevel('ink', '}', line, f, self.TPinklist)
                line = listLevel('sup', '}', line, f, self.TPsuplist)
                line = readLevel0('sigma', line, f, self.TPsigma)
                return
        else:
            logging.warning(f'path {bm} does not exist')

    def scrapeLabels(self) -> None:
        '''scrape the labels.csv document'''
        bm = os.path.join(self.fold, 'labels.csv')
        if os.path.exists(bm):
            with open(bm, "r") as f:
                data = list(csv.reader(f))
                self.TPinklist[0][1]=data[0][1]
                self.TPsuplist[0][1]=data[1][1]


    def scrapeCD(self) -> None:
        '''scrape controlDict
        s is a scrape object'''
        bm = os.path.join(self.casefold, 'system', 'controlDict')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                line = listLevel('application', '//', line, f, self.controlDictList)
                return
        else:
            logging.warning(f'path {bm} does not exist')



    def scrapeFV(self) -> None:
        '''scrape fvSolution
            fvSolution files are broken into sections by the variable we're solving for, e.g. alpha, pcorr, p_rgh. 
            scrape each section into a different list stored in s
            s is a scrape object'''
        bm = os.path.join(self.casefold, 'system', 'fvSolution')
        if os.path.exists(bm):
            with open (bm, "r") as f:
                line = f.readline()
                line = listLevel('\t\"alpha', '\t}', line, f, self.fvSailist)
                line = listLevel('\t\"pcorr', '\t}', line, f, self.fvSpcorrlist)
                line = listLevel('\tp_rgh', '\t}', line, f, self.fvSprghlist)
                line = listLevel('\tp_rghFinal', '\t}', line, f, self.fvSprghfinallist)
                line = listLevel('\tU', '\t}', line, f, self.fvSUlist)
                line = listLevel('PIMPLE', '\t}', line, f, self.fvSPIMPLElist)
                return
        else:
            logging.warning(f'path {bm} does not exist')


    def scrapeGeo(self) -> None:
        '''scrape geometry.csv
            geometry.csv was created by ncreate3d.py/noz3dscript.ipynb when we generated the whole folder
            s is a scrape object'''
        bm = os.path.join(self.fold, 'geometry.csv')
        transfer = {'nozzle inner width':'niw', 'nozzle thickness':'nt', 'bath width':'bw', 'bath depth':'bd'
                    , 'nozzle length':'nl', 'bath left coord':'blc', 'bath right coord':'brc'
                    , 'bath front coord':'bfc', 'bath back coord':'bbackc', 'bath bottom coord':'bbotc', 'bath top coord':'btc'
                    , 'nozzle bottom coord':'nbc', 'nozzle center x coord':'ncxc', 'nozzle center y coord':'ncyc', 'nozzle angle':'na'
                    , 'horizontal':'horiz', 'adjacent filament orientation':'adj', 'adjacent filament offset':'dst'
                    , 'corresponding simulation':'cor', 'bath velocity':'bathv', 'ink velocity':'inkv'} # RG
        if os.path.exists(bm):
            with open(bm, "r") as f:
                data = list(csv.reader(f))
                for row in data:
                    # name of variable, remove units if there
                    s1 = re.split(' \(', row[0])[0]
                    # remove space from beginning of units
                    if row[2][0]==' ':
                        row[2]=row[2][1:]
                    # store variable
                    if s1 in transfer:
                        setattr(self, f'GEO{transfer[s1]}', row)
                    if s1 == 'corresponding simulation':
                        self.compareto[1] = row[1][1:] # RG
                return
        else:
            return
        
        
    