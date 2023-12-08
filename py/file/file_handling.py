#!/usr/bin/env python
'''Functions for handling files and folders'''

# external packages
import os
import sys
import numpy as np
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentdir)
from file_export import *


#--------------------


def simFolders(folder:str) -> list:
    '''list all case folders in the folder
    input folder should hold simulations, e.g. "C:\\HBHBsweep" '''
    if type(folder) is list:
        # iterate through topfolders
        flist = []
        for f in folder:
            flist = flist+simFolders(f)
        return flist
    else:
        # find all sim folders in this folder
        flist = []
        for f in os.listdir(folder):
            fold = os.path.join(folder,f)
            fh = folderHandler(fold)
            if fh.isSimFolder():
                flist.append(fold)
        flist.sort(key=lambda x: int(os.path.basename(x)[2:]))
        return flist
    

class folderHandler:
    '''for finding files within a simulation folder'''
    
    def __init__(self, folder:str):
        self.folder = folder
        
    def shortName(self) -> str:
        '''get a name for this folder that only uses the last two folders in the path name'''
        shortname = os.path.join(os.path.basename(os.path.dirname(self.folder)), os.path.basename(self.folder))
        return shortname
        
        
    def isSimFolder(self) -> bool:
        '''determines if the folder is a simulation folder
        for true, input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
        if os.path.basename(self.folder)=='mesh':
            return False
        cf = self.caseFolder()
        if not os.path.exists(cf):
            return False
        else:
            return True
        
    def setAndReturn(self, name:str, s:str) -> str:
        '''set the stored value of name to s, and return s'''
        setattr(self, name, s)
        return s

    def setCaseFold(self, s:str) -> str:
        '''set the case folder to this value and return'''
        return self.setAndReturn('caseFold', s)
    
    def caseFolder(self) -> str:
        '''find the path of the folder with the case files (e.g. constant, system, 0)
        input folder should be a simulation folder, e.g. "C:\\...\\nb30" '''
        # if there is a folder in this folder named 'case', return that folder
        if hasattr(self, 'caseFold'):
            return self.caseFold
        casefold = os.path.join(self.folder, 'case')
        if os.path.exists(casefold):
            return self.setCaseFold(casefold)
        # if this folder contains a folder called 'constant', this is the case folder, so return this folder
        constfold = os.path.join(self.folder, 'constant')
        if os.path.exists(constfold):
            return self.setCaseFold(self.folder)
        vtkfold = os.path.join(self.folder, 'VTK')
        if os.path.exists(vtkfold):
            return self.setCaseFold(self.folder)
        ipfold = os.path.join(self.folder, 'interfacePoints')
        if os.path.exists(ipfold):
            return self.setCaseFold(self.folder)
        legfold = os.path.join(self.folder, 'legend.csv')
        if os.path.exists(legfold):
            return self.setCaseFold(folder)
        else:
            return ''
        
    def meshFolder(self) -> str:
        '''find the path of the folder with the mesh files
        input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
        # if there is a folder in this folder called 'mesh', return that folder
        if hasattr(self, 'meshFold'):
            return self.meshFold
        mf = os.path.join(self.folder, 'mesh')
        if os.path.exists(mf):
            return self.setAndReturn('meshFold', mf)
        # if there is a folder in the parent folder called 'mesh', return that folder
        mf = os.path.join(os.path.dirname(self.folder), 'mesh')
        if os.path.exists(mf):
            return self.setAndReturn('caseFold', mf)
        else:
            return ''
        
    def VTKFolder(self) -> str:
        '''Find the path of the VTK folder .
        Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
        if hasattr(self, 'VTKFold'):
            return self.VTKFold
        
        if not os.path.exists(self.folder):
            raise Exception('Folder does not exist')
        vtkfold = os.path.join(self.folder, 'VTK')
        if os.path.exists(vtkfold):
            return self.setAndReturn('VTKFold', vtkfold)
        vtkfold = os.path.join(self.folder, 'case', 'VTK')
        if os.path.exists(vtkfold):
            return self.setAndReturn('VTKFold', vtkfold)
        else:
            return ''
        
        
    #-----------------------------------------------------------
    
 
        
    def series(self, loop:bool=True) -> str:
        '''Find the .vtm.series or .vtk.series file in the folder.
        Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
        if hasattr(self, 'seriesFile'):
            return self.seriesFile
        
        vtkfolder = self.VTKFolder()
        if not os.path.exists(vtkfolder):
            return ''

        for file in os.listdir(vtkfolder):
            if '.series' in file:
                return self.setAndReturn('seriesFile', os.path.join(vtkfolder, file))
                
        # if there is a vtk folder but no series file, generate one
        if loop:
            self.redoVTKSeriesNoLog()
            return self.series(loop=False)
        else:
            return ''
        
        
    def correctVTM(self, file:str) -> None:
        '''Sometimes the t=0 vtm file gets saved with the wrong time. This fixes the file. File should be a .vtm file'''
        folderlabel = ""
        time = ""
        st = ''
        if not os.path.exists(file):
            return
        with open(file, 'r') as f:
            line = f.readline()
            st = f'{st}{line}'
            while not ('DataSet name=' in line) and len(line)>0:
                line = f.readline()
                st = f'{st}{line}'
            if len(line)>0:
                strs = re.split('_|/internal', line)
                folderlabel = strs[1]
            while not ('TimeValue' in line) and len(line)>0:
                line = f.readline()
                st = f'{st}{line}'
            if len(line)>0:
                line = f.readline()
                time = line.replace('\n','')
                if folderlabel=='0' and not time=='0':
                    st = f'{st}0\n'
                else:
                    st = f'{st}{line}'
            while len(line)>0:
                line = f.readline()
                st = f'{st}{line}'
        exportFile(os.path.dirname(file), os.path.basename(file), st)


    def readVTM(self, file:str) -> Tuple[str, str]:
        '''Find the vtk folder name and corresponding time for a vtm file
        File should be a .vtm file'''
        folderlabel = ""
        time = ""
        if not os.path.exists(file):
            return folderlabel, time
        with open(file, 'r') as f:
            line = f.readline()
            while not ('DataSet name=' in line) and len(line)>0:
                line = f.readline()
            if len(line)>0:
                strs = re.split('_|/internal', line)
                folderlabel = strs[1]
            while not ('TimeValue' in line) and len(line)>0:
                line = f.readline()
            if len(line)>0:
                line = f.readline()
                time = line.replace('\n','')
        if folderlabel=='0' and not time=='0':
            self.correctVTM(file)
            time = '0'
        return folderlabel, time


    def generateVTKSeries(self, tlist:List[str], flist:List[str], ending:str, lastTime:float=0) -> None:
        '''This creates a new .vtk.series or .vtm.series file
        tlist is a list of times
        flist is a list of folders. The times and folders don't need to be sorted ahead of time
        cf is the casefolder, e.g. 'C:\\...\\nb64'
        ending is the file extension, either '.vtm' or '.vtk' 
        lastTime is the time at which the most recent vtm or vtk file was updated'''
        seriesfile = self.series(loop=False)
        if os.path.exists(seriesfile):
            seriesTime = os.path.getmtime(seriesfile)
            if seriesTime>lastTime:
                # the series file is up to date
                print(seriesTime, lastTime)
                return
            cfbasename = os.path.basename(seriesfile).replace(ending+'.series', '') # e.g. 'case'
        else:
            cf = self.caseFolder()
            cfbasename = os.path.basename(cf) # e.g. 'case' or 'nb64'
        if len(tlist)==0 or len(flist)==0 or not len(tlist)==len(flist):
            return
        l = (np.array([tlist, flist])).transpose() # this creates a combined time and folder name table
        l = l[np.argsort(l[:,0])]    # this sorts the time and folder name table by time

        # generate file
        st = '{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n' # opening line
        for time, folderlabel in l[0:-1]:
            st = st + '    { \"name\" : \"' + cfbasename + '_' + folderlabel
            st = st + ending + '\", \"time\" : ' + time + ' },\n' 
            # each vtk file gets a row in the file
        st = st + '    { \"name\" : \"' + cfbasename + '_' + l[-1,1]
        st = st + ending + '\", \"time\" : ' + l[-1,0] + ' }\n' 
            # last line contains no comma at the end
        st = st + '  ]\n}' 
            # closing line
        exportFile(os.path.join(cf, 'VTK'), f'{cfbasename}{ending}.series', st)
        return


    def redoVTKSeriesNoLog(self) -> None:
        '''rewrite the .vtm.series file to include all .vtm files in the vtm folder
        folder should be the folder for the simulation (e.g. 'C:\\...\\nb64')
        this uses the existing vtk and vtm files in the folder'''
        cf = self.caseFolder()
        if not os.path.exists(cf):
            return
        vtkfolder = self.VTKFolder()
        if not os.path.exists(vtkfolder):
            return
        flist = [] # folder numbers
        tlist = [] # times
        ending = '.vtm'
        lastTime = 0
        for file in os.listdir(vtkfolder):
            ffull = os.path.join(vtkfolder, file)
            if file.endswith('.vtm'):
                updatedTime = os.path.getmtime(ffull)
                if updatedTime>lastTime:
                    lastTime=updatedTime
                flabel, time = self.readVTM(ffull)
                if len(flabel)>0 and len(time)>0:
                    flist.append(flabel)
                    tlist.append(time)
            elif file.endswith('.vtk'):
                updatedTime = os.path.getmtime(ffull)
                if updatedTime>lastTime:
                    lastTime=updatedTime
                ending = '.vtk'
                flabel = int(re.split('\_|.v', file)[1])
                flist.append(flabel)
                if len(tlist)==0:
                    tlist.append(0)
                else:
                    tlist.append(tlist[-1]+0.1)   
        if ending=='.vtk':
            flist.sort()
            flist = [str(f) for f in flist]
            tlist = ['{:1.1f}'.format(t) for t in tlist]

        self.generateVTKSeries(tlist, flist, ending, lastTime=lastTime)
        return
    
    #----------------------------------------
    
    def vtkFiles(self) -> int:
        '''Determine how many .vtk or .vtm files there are. 
        Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
        vtkfolder = self.VTKFolder()
        if not os.path.exists(vtkfolder):
            return 0
        num = 0
        for file in os.listdir(vtkfolder):
            if file.endswith('.vtk') or file.endswith('.vtm'):
                num+=1
        return num


    def parseVTKSeries(self) -> List[float]:
        '''parseVTKSeries extracts a list of times from the .vtm.series files
        folder is a full file name that can point to the case folder or the folder above it
        input folder should be a simulation folder, e.g. "C:\\...\\nb30" or "C:\\...\\nb30\\case"'''
        seriesfile = self.series()
        times = []
        if os.path.exists(seriesfile):
            with open(seriesfile, 'r') as f:
                for line in f:
                    if 'name' in line:
                        times.append(float(re.split('time\" : | }', line)[1]))
        numvtkFiles = self.vtkFiles()
        if len(times)<numvtkFiles:
            self.redoVTKSeriesNoLog()
        return times
    
    #-----------------------------------------
    
    
    def timesFromFolder(self) -> List[float]:
        '''get a list of times simulated in the folder
            input folder should be a simulation folder, e.g. "C:\\...\\nb30"
            this will only tell you about OpenFOAM folders, e.g. "0.1".'''
        ff = os.listdir(self.caseFolder())
        slist = []
        for i in ff:
            try:
                s0 = float(i)
            except:
                pass
            else:
                slist.append(s0)
        return slist


    def times(self) -> List[float]:
        '''Get a list of simulated times. 
        Input folder should be a simulation folder, e.g. "C:\\...\\nb30"
        This gets the list of times from the vtk file and from the list of files in the folder. If the vtk file has fewer times, the folder list is returned. The vtk list is not updated because we might be in the middle of a run, in which case there will be more folders than vtm files. If the vtk file has more or the same amount of times, the vtk list is returned. If neither has any times, the vtk file gets updated and the new list from the vtk file is returned'''
        t1 = self.timesFromFolder()
        t2 = self.parseVTKSeries()
        if len(t1)>len(t2):
            return t1
        else:
            return t2