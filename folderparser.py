#!/usr/bin/env python
'''Functions for handling files and folders used in OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. folderparser identifies log files from interFoam, etc. and collects information into csv tables
'''

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
try:
    from config import cfg
except:
    pass

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  

def logFN(scriptFile:str) -> str:
    '''Get a log file name, given a script file name'''
    compname = socket.gethostname()
    base = os.path.splitext(os.path.basename(scriptFile))[0]
    dirpath = os.path.dirname(os.path.realpath(__file__))
    try:
        cfgbase = cfg.path.logs
    except:
        cfgbase = 'logs'
    logfolder = os.path.join(dirpath, cfgbase)
    if not os.path.exists(logfolder):
        os.mkdir(logfolder)
#         logfolder = os.path.join(os.path.dirname(dirpath), cfgbase)
#         if not os.path.exists(logfolder):
#             logfolder = dirpath
    return os.path.join(logfolder,f'{base}_{compname}.log')

def openLog(f:str, LOGGERDEFINED:bool, level:str="INFO", exportLog:bool=True) -> bool:
    '''this code lets you create log files, so you can track when you've moved files. f is the file name of the script calling the openLog function'''

    if not LOGGERDEFINED:
        loglevel = getattr(logging,level)
        root = logging.getLogger()
        if len(root.handlers)>0:
            # if handlers are already set up, don't set them up again
            return True
        root.setLevel(loglevel)

        # send messages to file
        if exportLog:
            logfile = logFN(f)
            filehandler = logging.FileHandler(logfile)
            filehandler.setLevel(loglevel)
            formatter = logging.Formatter("%(asctime)s/{}/%(levelname)s: %(message)s".format(socket.gethostname()), datefmt='%b%d/%H:%M:%S')
            filehandler.setFormatter(formatter)
            root.addHandler(filehandler)
            logging.info(f'Established log: {logfile}')

        # print messages
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(loglevel)
        formatter2 = logging.Formatter('%(levelname)s: %(message)s')
        handler.setFormatter(formatter2)
        root.addHandler(handler)
        LOGGERDEFINED = True
        
    return LOGGERDEFINED

def printCurrentTime() -> None:
    '''Print the current time'''
    now = datetime.now()
    current_time = "------ Current Time = " + now.strftime("%D, %H:%M:%S")
    logging.info(current_time)


    #-------------------------------------------------------------------------------------------------  
########### folder tools ##################

def shortName(folder:str) -> str:
    shortname = os.path.join(os.path.basename(os.path.dirname(folder)), os.path.basename(folder))
    return shortname


def caseFolder(folder:str) -> str:
    '''find the path of the folder with the case files (e.g. constant, system, 0)
    input folder should be a simulation folder, e.g. "C:\\...\\nb30" '''
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
    ipfold = os.path.join(folder, 'interfacePoints')
    if os.path.exists(ipfold):
        return folder
    legfold = os.path.join(folder, 'legend.csv')
    if os.path.exists(legfold):
        return folder
    else:
        return ''
        

def isSimFolder(folder:str) -> bool:
    '''determines if the folder is a simulation folder
        for true, input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
    if os.path.basename(folder)=='mesh':
        return False
    cf = caseFolder(folder)
    if not os.path.exists(cf):
        return False
    else:
        return True
    

def caseFolders(folder:str) -> List[str]:
    '''list all case folders in the folder
    input folder should hold simulations, e.g. "C:\\HBHBsweep" '''
    flist = []
    for f in os.listdir(folder):
        fold = os.path.join(folder,f)
        if isSimFolder(fold):
            flist.append(fold)
    return flist


def meshFolder(folder:str) -> str:
    '''find the path of the folder with the mesh files
    input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
    # if there is a folder in this folder called 'mesh', return that folder
    mf = os.path.join(folder, 'mesh')
    if os.path.exists(mf):
        return mf
    # if there is a folder in the parent folder called 'mesh', return that folder
    mf = os.path.join(os.path.dirname(folder), 'mesh')
    if os.path.exists(mf):
        return mf
    else:
        return ''

    #-------------------------------------------------------------------------------------------------  
########### VTK file tools ##################
    

def VTKFolder(folder:str) -> str:
    '''Find the path of the VTK folder .
    Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
    if not os.path.exists(folder):
        raise Exception('Folder does not exist')
    vtkfold = os.path.join(folder, 'VTK')
    if os.path.exists(vtkfold):
        return vtkfold
    vtkfold = os.path.join(folder, 'case', 'VTK')
    if os.path.exists(vtkfold):
        return vtkfold
    else:
        raise Exception('No VTK folder')


def series(folder:str, loop:bool=True) -> str:
    '''Find the .vtm.series or .vtk.series file in the folder.
    Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
    try:
        vtkfolder = VTKFolder(folder)
    except:
        return ''
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if '.series' in file:
                return os.path.join(vtkfolder, file)
        # if there is a vtk folder but no series file, generate one
        if loop:
            redoVTKSeriesNoLog(folder)
            return series(folder, loop=False)
        else:
            return ''
    return ""


def vtkFiles(folder:str) -> int:
    '''Determine how many .vtk or .vtm files there are. 
    Input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
    try:
        vtkfolder = VTKFolder(folder)
    except:
        return 0
    num = 0
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if file.endswith('.vtk') or file.endswith('.vtm'):
                num+=1
    return num


def parseVTKSeries(folder:str) -> List[float]:
    '''parseVTKSeries extracts a list of times from the .vtm.series files
    folder is a full file name that can point to the case folder or the folder above it
    input folder should be a simulation folder, e.g. "C:\\...\\nb30" or "C:\\...\\nb30\\case"'''
    seriesfile = series(folder)
    times = []
    if os.path.exists(seriesfile):
        with open(seriesfile, 'r') as f:
            for line in f:
                if 'name' in line:
                    times.append(float(re.split('time\" : | }', line)[1]))
    numvtkFiles = vtkFiles(folder)
    if len(times)<numvtkFiles:
        redoVTKSeriesNoLog(folder)
    return times


def timesFromFolder(folder:str) -> List[float]:
    '''get a list of times simulated in the folder
        input folder should be a simulation folder, e.g. "C:\\...\\nb30"
        this will only tell you about OpenFOAM folders, e.g. "0.1".'''
    ff = os.listdir(caseFolder(folder))
    slist = []
    for i in ff:
        try:
            s0 = float(i)
        except:
            pass
        else:
            slist.append(s0)
    return slist


def times(folder:str) -> List[float]:
    '''Get a list of simulated times. 
    Input folder should be a simulation folder, e.g. "C:\\...\\nb30"
    This gets the list of times from the vtk file and from the list of files in the folder. If the vtk file has fewer times, the folder list is returned. The vtk list is not updated because we might be in the middle of a run, in which case there will be more folders than vtm files. If the vtk file has more or the same amount of times, the vtk list is returned. If neither has any times, the vtk file gets updated and the new list from the vtk file is returned'''
    t1 = timesFromFolder(folder)
    t2 = parseVTKSeries(folder)
    if len(t1)>len(t2):
        return t1
    else:
        return t2
    

def currentTime(folder:str) -> float:
    '''Get the latest simulated time
    Input folder should be a simulation folder, e.g. "C:\\...\\nb30"
    If there is no vtk file or simulation folders, it takes the value from the legend'''
    rdict = currentRate(folder)
    endTime = rdict['end_time']
    t2 = times(folder)
    if len(t2)>0:
        simtime = max(t2)
    else:
        simtime = rdict['simulation_time']
        
    try:
        simtime = float(simtime)
    except:
        pass
    try:
        endTime = float(endTime)
    except:
        pass
    return simtime, endTime
            
def currentRate(folder:str) -> dict:
    '''get the current time, end time, and simulation rate from the folder'''
    t = legendUnique(folder)
    if len(t)==0:
        return {'simulation_time':'', 'end_time':'', 'rate':'', 'run_time':''}
    if 'endTime' in t:
        endTime = t['endTime']
    else:
        if 'endTime_(s)' in t:
            endTime = t['endTime_(s)']
        else:
            endTime = ''
    if 'simulation_time' in t:
        simTime = t['simulation_time']
    else:
        if 'simulation_time_(s)' in t:
            simTime = t['simulation_time_(s)']
        else:
            simTime = ''
    if '_interFoam_time' in t:
        runTime = t['_interFoam_time']
    else:
        if '_interFoam_time_(s)' in t:
            runTime = t['_interFoam_time_(s)']
        else:
            runTime = ''
    if 'simulation_rate' in t:
        rate = t['simulation_rate']
    else:
        if 'simulation_rate_(hr/s)' in t:
            rate = t['simulation_rate_(hr/s)']
        else:
            rate = ''
    return {'simulation_time':simTime, 'end_time':endTime, 'rate':rate, 'run_time':runTime}

#-------------------------------------------------------------------------------------------------  
    ##################################
##### overwriting vtk.series files



def correctVTM(file:str) -> None:
    '''Sometimes the t=0 vtm file gets saved with the wrong time. This fixes the file. File should be a .vtm file'''
    folderlabel = ""
    time = ""
    st = ''
    if not os.path.exists(file):
        return
    with open(file, 'r') as f:
        line = f.readline()
        st = st + line
        while not ('DataSet name=' in line) and len(line)>0:
            line = f.readline()
            st = st + line
        if len(line)>0:
            strs = re.split('_|/internal', line)
            folderlabel = strs[1]
        while not ('TimeValue' in line) and len(line)>0:
            line = f.readline()
            st = st + line
        if len(line)>0:
            line = f.readline()
            time = line.replace('\n','')
            if folderlabel=='0' and not time=='0':
                st = st + '0\n'
            else:
                st = st + line
        while len(line)>0:
            line = f.readline()
            st = st + line
    exportFile(os.path.dirname(file), os.path.basename(file), st)


def readVTM(file:str) -> Tuple[str, str]:
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
        correctVTM(file)
        time = '0'
    return folderlabel, time
    

def generateVTKSeries(tlist:List[str], flist:List[str], cf:str, ending:str, lastTime:float=0) -> None:
    '''This creates a new .vtk.series or .vtm.series file
    tlist is a list of times
    flist is a list of folders. The times and folders don't need to be sorted ahead of time
    cf is the casefolder, e.g. 'C:\\...\\nb64'
    ending is the file extension, either '.vtm' or '.vtk' 
    lastTime is the time at which the most recent vtm or vtk file was updated'''
    seriesfile = series(cf, loop=False)
    if os.path.exists(seriesfile):
        seriesTime = os.path.getmtime(seriesfile)
        if seriesTime>lastTime:
            # the series file is up to date
            return
        cfbasename = os.path.basename(seriesfile).replace(ending+'.series', '') # e.g. 'case'
    else:
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
    exportFile(os.path.join(cf, 'VTK'), cfbasename + ending + '.series', st)
    return


def redoVTKSeriesNoLog(folder:str) -> None:
    '''rewrite the .vtm.series file to include all .vtm files in the vtm folder
    folder should be the folder for the simulation (e.g. 'C:\\...\\nb64')
    this uses the existing vtk and vtm files in the folder'''
    cf = caseFolder(folder) # In older versions, sometimes the case file is below the folder and named 'case', 
                            # and sometimes it is the folder and named something like 'nb64'. 
                            # This allows flexibility for those two situations
    if cf=='':
        return
    cfbasename = os.path.basename(cf) # e.g. 'case' or 'nb64'
    vtkfolder = os.path.join(cf, 'VTK')
    flist = [] # folder numbers
    tlist = [] # times
    ending = '.vtm'
    lastTime = 0
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if file.endswith('.vtm'):
                updatedTime = os.path.getmtime(file)
                if updatedTime>lastTime:
                    lastTime=updatedTime
                flabel, time = readVTM(os.path.join(vtkfolder, file))
                if len(flabel)>0 and len(time)>0:
                    flist.append(flabel)
                    tlist.append(time)
            elif file.endswith('.vtk'):
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
        generateVTKSeries(tlist, flist, cf, ending, lastTime=lastTime)
    return

#-------------------------------------------------------------------------------------------------  
##### FILE HANDLING



def exportFile(folder:str, file:str, text:str) -> None:
    '''exportFile exports a text file
    folder is the folder to save the file to (e.g. 'C:\\...\\myfolder')
    file is a file base name (e.g. 'myfile.txt') within that folder
    text is the text to write to the file'''
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    logging.info("Exported file %s" % fn)


def exportCSV(fn:str, table:List[Any]) -> None:
    '''exportCSV exports a csv file
    fn is the full path of the file to export
    table is the table to export, as a list'''
    with open(fn, mode='w', newline='') as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in table:
            w.writerow(i)
    logging.info('Exported file %s' % fn)


def importIf(folder:str) -> List[List[str]]:
    '''importIf imports an existing legend file if it exists 
    folder is a full path, e.g. "C:\\...\\nb64"'''
    fn = os.path.join(folder, 'legend.csv')
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            return(list(csv.reader(f)))
    else:
        return []   

def legendUnique(folder:str, units:bool=False) -> Union[Tuple[dict,dict], dict]:
    '''legendUnique imports a legend, rewrites the variable names so they are all unique and ready to compile into a pandas dataframe, and outputs a dictionary'''
    t = importIf(folder)
    if len(t)==0:
        return {}
    headers = [i[0] for i in t]
    section = ''
    for i in range(len(t)):
        if t[i][0]=='':
            section = t[i+1][0]
        elif t[i][0] in ['sup', 'ink', 'controlDict', 'dynamicMeshDict']:
            section = t[i][0]
        else:
            if headers.count(t[i][0])>1:
                t[i][0] = section+'_'+t[i][0].replace(' ', '_')
            else:
                t[i][0] = t[i][0].replace(' ', '_')
        if len(t[i])==2:
            t[i] = t[i] +['']
    
    values = {a[0]:a[1] for a in t}
    if units:
        if len(t[2])==3:
            u = {a[0]:a[2] for a in t}
        else:
            u = legendUnique(cfg.path.legend_units, units=False)
        return values, u
    else:
        return values

def mkdirif(path:str) -> int:
    '''make a directory if it doesn't exist
    path is a full path name
    returns 1 if error, returns 0 for no error'''
    try:
        os.mkdir(path)
    except OSError as e:
        return 1
    else:
        logging.info ("Created directory %s" % path)
    return 0


def copy(src:str, dest:str) -> None:
    '''Copy directory src to directory dest. Both should be full folder paths'''
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            logging.error('Directory not copied. Error: %s' % e)
    return
 
    
######################################
#-------------------------------------------------------------------------------------------------  

    
def modifyControlDict(folder:str, tend:float) -> int:
    '''change the endtime in the controlDict
    returns 0 if the dictionary was changed, 1 if not'''
    cfi = os.path.join(folder, 'case')
    if os.path.exists(cfi):
        cf = cfi
    else:
        cf = folder
    cdfile = os.path.join(cf, 'system', 'controlDict')
    cdfilenew = os.path.join(cf, 'system', 'controlDict2')
    if not os.path.exists(cdfile):
        return 1
    retval = 0
    # if the endtime is already at the value, abort this loop and delete the new file
    with open(cdfile, 'r') as fold:
        with open(cdfilenew, 'w') as fnew:
            for line in fold:
                if line.startswith('endTime'):
                    linenew = 'endTime\t'+str(tend)+';\n'
                    if linenew==line:
                        retval = 1
                        break
                    else:
                        line = linenew
                fnew.write(line)
                
    # if we changed the endtime, overwrite the old file
    if retval==0:            
        os.remove(cdfile)
        os.rename(cdfilenew, cdfile)
        print('Set end time to '+str(tend)+' in '+cdfile)
    else:
        os.remove(cdfilenew)
    return retval

    
    
    
    
    ############## ARCHIVE
    
    
    

# def consolidateVTK(folder:str) -> None:
#     '''# sometimes you end up with a vtk folder in the main folder and a vtk folder inside the case folder, or just a vtk folder inside the case folder. If this happens, this consolidates all the files into the vtk folder inside the main folder
# # input folder should be a simulation folder, e.g. "C:\\...\\nb30"'''
#     case = os.path.join(folder, 'case', 'VTK')
#     if not os.path.exists(case):
#         return
#     vtk = os.path.join(folder, 'VTK')
#     if not os.path.exists(vtk):
#         mkdirif(vtk)
#     copyFilesInFolder(case, vtk)
#     shutil.rmtree(case)
#     try:
#         os.rmdir(os.path.join(folder, 'case'))
#     except:
#         pass
#     redoVTKSeriesNoLog(vtk)
# #     return

    

# def iPFileTime(file:str) -> int:
#     '''get the time of the interfacepoints file'''
#     f2num = int(re.split('_|.csv', file)[-2])
#     return f2num


# def rmIPFiles(ip:str, f2:str) -> None:
#     '''delete all files after this file in the interfacePoints folder
#     this is useful if there is an error generating interfacePoints files
#     ip is interfacepoints folder
#     f2 is the first bad file'''
#     f2num = iPFileTime(f2)
#     currFile = os.path.join(ip, 'interfacePoints_t_'+str(f2num)+'.csv')
#     while os.path.exists(currFile):
#         os.remove(currFile)
#         f2num+=1
#         currFile = os.path.join(ip, 'interfacePoints_t_'+str(f2num)+'.csv')
#     if os.path.exists(os.path.join(ip, 'temp.csv')):
#         os.remove(os.path.join(ip, 'temp.csv'))

        

# def sortIPFolder(f:str) -> List[str]:
#     '''sort files in interfacepoints folder by time and return a list of times
#     f is full path name '''
#     l1 = os.listdir(f)
#     l1.sort(key=iPFileTime)
#     return l1


# def secondLineInCSV(csvfile:str) -> str:
#     '''get the second line in a csv
#     csvfile needs to be a full path name'''
#     file = open(csvfile)
#     csv_reader = csv.reader(file)
#     next(csv_reader)
#     return next(csv_reader)

 
# def assessFolders(fromFolder:str) -> None:
#     '''for all folders inside fromFolder, determine if interface points were generated correctly   '''
#     logging.info('Parent folder: ', fromFolder)
#     for f in os.listdir(fromFolder):
#         try:
#             assessFolder(fromFolder, f)
#         except:
#             logging.info(f+' interface points error')
#     return



# def assessFolder(fromFolder:str, f:str) -> None:
#     '''determine if the interfacePoints files were correctly exported
#     fromFolder is parent folder, full path
#     f is simulation folder, just the basename'''
#     s1 = 0
#     s2 = 0
#     f2 = ''
#     excludelist = ['nb79', 'nb89']
#     if f.startswith('nb') and f not in excludelist:
#         ip = os.path.join(fromFolder, f, 'interfacePoints')
#         if os.path.exists(ip):
#             for f1 in sortIPFolder(ip):
#                 # if two consecutive files have the same size or 
#                 # the file size shrank by a lot, there was an error
#                 f1path = os.path.join(ip, f1)
#                 s1 = os.stat(f1path).st_size
#                 if (s1<0.5*s2 and s1<10000) or (s1==s2 and secondLineInCSV(f1path)==secondLineInCSV(f2)):
#                     rmIPFiles(ip, f2)
# #                     logging.info(f1path, s1)
#                     raise Exception(f+' interface points error')
#                 s2 = s1
#                 f2 = f1path
#     return


    
    

# def modifyControlDictRunTime(folder):
#     '''this modifies old controlDict files so runTimeModifiable is true
#     folder is full path, e.g. 'C:\...\nb64''''
#     cfi = os.path.join(folder, 'case')
#     if os.path.exists(cfi):
#         cf = cfi
#     else:
#         cf = folder
#     cdfile = os.path.join(cf, 'system', 'controlDict')
#     cdfilenew = os.path.join(cf, 'system', 'controlDict2')
#     if not os.path.exists(cdfile):
#         return
#     with open(cdfile, 'r') as fold:
#         with open(cdfilenew, 'w') as fnew:
#             for line in fold:
#                 if line.startswith('runTimeModifiable'):
#                     line = 'runTimeModifiable\tyes;\n'
#                 fnew.write(line)
#     os.remove(cdfile)
#     os.rename(cdfilenew, cdfile)
#     print(cdfile)