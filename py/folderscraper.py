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
from backwardsRead import fileReadBackwards
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from folderparser import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich", "Ross Gunther"]
__license__ = "NIST"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  


class scrape:
    '''the scrape class is used to store variables scraped from the log files
        scrape contains a function scrape.table(), which converts all the data in the class into a 2 column table'''
        
    def __init__(self, folder:str):
        '''folder is a folder name (string) which should point either to the case folder or one level above it'''
        
        self.fold = folder # full path
        self.meshfold = meshFolder(folder)
        self.casefold = caseFolder(folder)
        self.folder = ['folder', os.path.basename(folder)] # just the folder name
        self.compareto = ['compare to', '', '']

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
        self.shmtimes = ['snappyHexMesh time', '', 's'] 
            # how long it takes snappyHexMesh to run
        self.shmtimem = ['snappyHexMesh time', '', 'min']
        self.iftimes = ['interFoam time s', '', 's'] 
            # how long it takes interfoam to run
        self.iftimehr = ['interFoam time hr', '', 'hr']
        self.simTime = ['simulation time', '', 's'] 
            # how many seconds within the simulation we ran
        self.simrate = ['simulation rate', '',  'hr/s'] 
            # how many real hours it takes to run each second within the simulation
        self.version = ['openfoam_version', '', '']

        
        
    #-------------------------
    def initGeo(self):
        '''these values measure the geometry of the system and are stored in geometry.csv, which is created by ncreate3d.py/noz3dscript.ipynb'''
        self.GEOniw = ['nozzle inner width', '0.603', 'mm']
        self.GEOnt = ['nozzle thickness', '0.152', 'mm']
        self.GEObw = ['bath width', '', 'mm']
        self.GEObd = ['bath depth', '', 'mm']
        self.GEOnl = ['nozzle length', '', 'mm']
        self.GEOblc = ['bath left coord', '', 'mm'] # x coord
        self.GEObrc = ['bath right coord', '', 'mm'] # x coord
        self.GEObfc = ['bath front coord', '', 'mm'] # y
        self.GEObbackc = ['bath back coord', '', 'mm'] # y
        self.GEObbotc = ['bath bottom coord', '', 'mm'] # z
        self.GEObtc = ['bath top coord', '', 'mm'] # z
        self.GEOnbc = ['nozzle bottom coord', '', 'mm'] # z
        self.GEOncxc = ['nozzle center x coord', '', 'mm']
        self.GEOncyc = ['nozzle center y coord', '', 'mm']
        self.GEOna = ['nozzle angle', '0', 'degrees'] # RG
        self.GEOhoriz = ['horizontal', False, '']
        self.GEObathv = ['bath velocity', '', 'm/s']
        self.GEOinkv = ['ink velocity', '', 'm/s']
    
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
        self.CMCfixedWallsLevel = ['fixedWalls level', '', '']
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
        self.ALCfixedWallsLayers = ['fixedWalls nSurfaceLayers', '', '']
        self.MQClist = [['nSmoothScale', '', ''], \
                        ['errorReduction', '', ''], \
                        ['maxNonOrtho', '']]
        self.SHMmergeTolerance = ['mergeTolerance', '', '']
                        
        self.blocksdims = ['blocks dims', '', '']
    
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
                        ['transportModel', '', ''],\
                          ['nu', '', 'm^2/s'], \
                          ['nu0', '', 'm^2/s'], \
                          ['tau0', '', 'm^2/s^2'], \
                          ['k', '', 'm^2*s^(n-2)'],\
                          ['n', '', ''], \
                          ['rho', '', 'kg/m^3']]
        self.TPsuplist = [['sup', '', ''],\
                        ['transportModel', '', ''],\
                          ['nu', '', 'm^2/s'], \
                          ['nu0', '', 'm^2/s'], \
                          ['tau0', '', 'm^2/s^2'], \
                          ['k', '', 'm^2*s^(n-2)'],\
                          ['n', '', ''], \
                          ['rho', '', 'kg/m^3']]
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
                  self.simTime, self.simrate]:
            col.append(i) # RG
        ca(col, ['', '', 'mesh','GEOMETRY'],\
           [self.GEOniw, self.GEOnt, self.GEObw, self.GEObd, \
            self.GEOnl, self.GEOblc, self.GEObrc, self.GEObfc,\
            self.GEObbackc, self.GEObbotc, self.GEObtc, self.GEOnbc,\
            self.GEOncxc, self.GEOncyc, self.GEOna, self.GEOhoriz, self.GEObathv, self.GEOinkv])
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
        ti = times(self.fold)
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
                   , 'nozzle bottom coord':'nbc', 'nozzle center x coord':'ncxc', 'nozzle center y coord':'ncyc', 'nozzle angle':'na', 'horizontal':'horiz'
                   , 'bath velocity':'bathv', 'ink velocity':'inkv'}
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
                return
        else:
            return




def ca(col:List[List[str]], hlist:List[str], li:List[List[str]]) -> None:
    '''ca is used by scrape.table() to compile all of the data in the object into a table
    col is a table with 2 columns. 
    hlist is a list of headers that describe this chunk of data. they are added as rows with an empty value 
    li is a list of [variable name, value] to add to col'''
    for i in hlist:
        col.append([i, ''])
    for i in li:
        col.append(i)
        
 
 ###################### SCRAPING TOOLS ##################

        

def placeInList(l:List[List[str]], s:str, v:str) -> None: 
    '''placeInList puts a value v into a list with 2 columns, where s is the name of the variable
        this function will only place the value into the list if the variable name s is in the list and there isn't already a value there
        l is a list
        s is a string
        v is a value'''
    
    # find the index in the list where this value should go
    i = -1
    found=False
    while i<len(l)-1 and not found:
        i+=1
        if len(l[i])>0 and l[i][0]==s:
            found=True    
            
    # put the value in the list
    if not found:
        return
    else:
        l[i][1] = v # this puts the value in the list
    return


def cancelUnits(s:str) -> str:
    '''cancelUnits removes the units list (e.g. [ 0 2 -1 0 0 0 0]) from a value
    s is a string'''
    strs = re.split(' ', s)
    return strs[-1]


def listLevel(startString:str, endString:str, line:str, f:TextIO, l:List[List[str]]) -> str:
    '''listLevel is a tool that looks for sections of interest within files and scrapes values out of them
    startString is a string that tells us that we've reached the section of interest. It should be at the beginning of the line
    endString tells us that we've reached the end of the section of interest. It should be at the beginning of the line.
    line is the starting line
    f is a file stream created by open()
    l is a list of variable names and values that we're going to scrape values into
    returns the line we just read'''
    
    while not line.startswith(startString):
        line = f.readline()
    
    while not line.startswith(endString):
        strs = re.split(';|\t', line) # split the line at tabs and semicolons
        ii = 0
        si = 0
        while ii<len(strs) and len(strs[ii])==0:
            ii+=1 # find the first non-empty string in the list
        if ii+1<=len(strs)-1:
            s0 = strs[ii] # if there are enough entries in the split line to contain a name and value, place this entry into l
            s1 = cancelUnits(strs[ii+1])
            placeInList(l, s0, s1)
        line = f.readline()
    return line  


def readLevel0(s:str, line:str, f:TextIO, obj:List[List[str]]) -> str:
    '''readLevel0 is a simpler version of listLevel, where instead of placing many values in a list,
    we're looking for a single value. This function targets lines in files that have no tabs at the beginning, just "name\tvalue"
    s is a trigger string that tells us we've found the value we're looking for
    line is a starting line
    f is a file stream created by open()
    obj is a [1x2] list into which we'll store our value
    returns the line we just read'''
    while not line.startswith(s):
        line = f.readline()
    strs = re.split(';|\t', line) # split the line at ; and tabs
    obj[1] = strs[1] # the value will always be the second item
    return line


###################### SCRAPING FILES ##################





#----------------------------------------------------------------------------

def populate(folder:str, *varargin, readLogs:bool=True, overwrite:bool=False) -> List[List[str]]:
    '''populate scrapes all the data from the folder and exports it to a table called legend.csv
    folder is a full path name
    can also add strings that evaluate specific functions, e.g. scrapeTP, so that you can just scrape the transportproperties if there is already a legend'''
    if not isSimFolder(folder):
        raise Exception("Not a simulation folder")
    s = scrape(folder)   # create an object to store variables
    fn = os.path.join(folder, 'legend.csv')     # export file name
    s.compareto[1] = os.path.basename(os.path.dirname(folder))
    if readLogs:
        s.scrapeLogs()   # scrape the logs
        s.scrapeCD() # scrape the control dictionary
    if not overwrite and os.path.exists(fn):
        leOld,uOld = legendUnique(folder, units=True) # import old legend and units to dictionary
        t = s.table()
        leNew, uNew = legendTableToDict(t, units=True) # convert new legend and units to dictionary
        for key in leNew:
            if not key=='compare_to':
                # keep old values if not same, except for scraped timings
                if key in leOld and not leOld[key]==leNew[key]:
                    if not (('time' in key or 'rate' in key) and len(leNew[key])>0 and float(leNew[key])>0):
                        if leOld[key][0]==' ':
                            # remove space from beginning of unit
                            leOld[key] = leOld[key][1:]
                        leNew[key] = leOld[key]
                if key in uOld and not uOld[key]==uNew[key]:
                    if uOld[key][0]==' ':
                        # remove space from beginning of unit
                        uOld[key] = uOld[key][1:]
                    uNew[key] = uOld[key]
        t = [[key, leNew[key], uNew[key]] for key in leNew]                        
    else:
        # if there is no legend file, we have to go through all the files and scrape data
        s.scrapeAll()
        t = s.table() # generate a table from all the data we scraped
    exportCSV(fn, t) # export the updated table
    return t



def populateToTable(folder:str, repopulate:bool = False) -> List[List[str]]:
    '''take one folder and either import the existing legend or create a new one
    folder is a full path name
    repopulate is true if you want to overwrite existing legend.csv files'''
    if repopulate:
        # overwrite times in legend.csv, or if there is no file, create a new legend file
        t2 = populate(folder)
    else:
        # don't overwrite anything, just import the existing file
        t2 = importIf(folder)
    return t2


def populateList(liInit:List[str], exportFilename:str, repopulate:bool = False) -> None:
    '''populateList scrapes all data for all files in a list and creates a combined table
    liInit is a list of folders to scrape
    exportFilename is the destination to export the summary file
    repopulate is true if you want to overwrite existing legend.csv files'''
    li = []
    for folder in liInit:
        if os.path.exists(folder):
            li.append(folder)
    t1 = populateToTable(li[0], repopulate)   
        # first collect one table so you have the format
    tbig = [['' for i in range(len(li)+1)] for j in range(len(t1))] 
        # tbig combines all the legend files into one table
    for i,row in enumerate(t1):
        tbig[i][0] = row[0] 
            # import the variable names and the values 
            # for the first file into the big table
        tbig[i][1] = row[1]
    for j,l in enumerate(li[1:]): 
        # go through the rest of the files and put the values in the big table
        t2 = populateToTable(l, repopulate)
        for i,row in enumerate(t2):
            tbig[i][j+2] = row[1]
    exportCSV(exportFilename, tbig) 
        # export the combined table
    return


#------------------------------

# RETROFITTING FUNCTIONS


def updateGeoFile(folder:str) -> None:
    '''update the geometry file to put units in the third column'''
    fn = os.path.join(folder, 'geometry.csv')
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            geo = list(csv.reader(f))
    else:
        # no geometry file
        return
    if len(geo[0])==3:
        # geometry file already has units
        return
    new = []
    for row in geo:
        spl = re.split(' \(|\)', row[0])
        if len(spl)>1:
            newrow = [spl[0], row[1], spl[1]]
        else:
            newrow = [spl[0], row[1], '']
        new.append(newrow)
    exportCSV(fn, new)
    logging.info(f'Exported {fn}')
    
def addUnitsToLegend(folder:str) -> None:
    '''add units to the legend file'''
    updateGeoFile(folder)
    populate(folder)
        

