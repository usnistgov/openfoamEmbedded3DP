###############################
# folderparser identifies log files from interFoam, etc. and collects information into csv tables

import numpy as np
import os
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket, sys

#-------------------------------------------------------------------------------------------------  

# this code lets you create log files, so you can track when you've moved files
def openLog(f:str, LOGGERDEFINED):
    if not LOGGERDEFINED:
        root = logging.getLogger()
        root.setLevel(logging.INFO)

        # send messages to file
        filehandler = logging.FileHandler(f)
        filehandler.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s/{}/%(levelname)s: %(message)s".format(socket.gethostname()), datefmt='%b%d/%H:%M:%S')
        filehandler.setFormatter(formatter)
        root.addHandler(filehandler)

        # print messages
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        formatter2 = logging.Formatter('%(message)s')
        handler.setFormatter(formatter2)
        root.addHandler(handler)
        LOGGERDEFINED = True
    return LOGGERDEFINED

def printCurrentTime():
    now = datetime.now()
    current_time = "------ Current Time = " + now.strftime("%D, %H:%M:%S")
    logging.info(current_time)


    #-------------------------------------------------------------------------------------------------  
########### folder tools ##################



# find the path of the folder with the case files (e.g. constant, system, 0)
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def caseFolder(folder:str) -> str:
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
        
# determines if the folder is a simulation folder
# for true, input folder should be a simulation folder, e.g. 'C:\...\nb30'
def isSimFolder(folder:str) -> bool:
    b = os.path.basename(folder)
    if not b.startswith('nb'):
        return False
    cf = caseFolder(folder)
    if not os.path.exists(cf):
        return False
    else:
        return True
    
# list all case folders in the folder
# input folder should hold simulations, e.g. 'C:\HBHBsweep'
def caseFolders(folder:str) -> List[str]:
    flist = []
    for f in os.listdir(folder):
        fold = os.path.join(folder,f)
        if isSimFolder(fold):
            flist.append(fold)
    return flist

# find the path of the folder with the mesh files
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def meshFolder(folder:str) -> str:
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
    
# find the path of the VTK folder 
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def VTKFolder(folder:str) -> str:
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

# find the .vtm.series or .vtk.series file in the folder
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def series(folder:str) -> str:
    try:
        vtkfolder = VTKFolder(folder)
    except:
        return ''
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if '.series' in file:
                return os.path.join(vtkfolder, file)
    return ""

# determine how many .vtk or .vtm files there are
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def vtkfiles(folder:str) -> int:
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

# parseVTKSeries extracts a list of times from the .vtm.series files
    # folder is a full file name that can point to the case folder or the folder above it
    # input folder should be a simulation folder, e.g. 'C:\...\nb30' or 'C:\...\nb30\case'
def parseVTKSeries(folder:str) -> List[float]:
    seriesfile = series(folder)
    times = []
    if os.path.exists(seriesfile):
        with open(seriesfile, 'r') as f:
            for line in f:
                if 'name' in line:
                    times.append(float(re.split('time\" : | }', line)[1]))
    numvtkfiles = vtkfiles(folder)
    if len(times)<numvtkfiles:
        redoVTKSeriesNoLog(folder)
    return times

# get a list of times simulated in the folder
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
# this will only tell you about OpenFOAM folders, e.g. "0.1".
def timesFromFolder(folder:str) -> List[float]:
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

# get a list of simulated times
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
# this gets the list of times from the vtk file and from the list of files in the folder.
# if the vtk file has fewer times, the folder list is returned. 
    # The vtk list is not updated because we might be in the middle of a run, in which case there will be more folders than vtm files
# if the vtk file has more or the same amount of times, the vtk list is returned
# if neither has any times, the vtk file gets updated and the new list from the vtk file is returned
def times(folder:str) -> List[float]:
    t1 = timesFromFolder(folder)
    t2 = parseVTKSeries(folder)
    if len(t1)>len(t2):
        return t1
    else:
        return t2
    
# get the latest simulated time
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
# if there is no vtk file or simulation folders, it takes the value from the legend
def currentTime(folder:str) -> float:
    t = times(folder)
    if len(t)>0:
        return max(t)
    else:
        t = importIf(folder)
        if len(t)==0:
            return 0
        else:
            return t[6][1]   

#-------------------------------------------------------------------------------------------------  
    ##################################
##### overwriting vtk.series files

# sometimes you end up with a vtk folder in the main folder and a vtk folder inside the case folder, or just a vtk folder inside the case folder. If this happens, this consolidates all the files into the vtk folder inside the main folder
# input folder should be a simulation folder, e.g. 'C:\...\nb30'
def consolidateVTK(folder:str) -> None:
    case = os.path.join(folder, 'case', 'VTK')
    if not os.path.exists(case):
        return
    vtk = os.path.join(folder, 'VTK')
    if not os.path.exists(vtk):
        mkdirif(vtk)
    copyFilesInFolder(case, vtk)
    shutil.rmtree(case)
    try:
        os.rmdir(os.path.join(folder, 'case'))
    except:
        pass
    redoVTKSeriesNoLog(vtk)
    return

# sometimes the t=0 vtm file gets saved with the wrong time. This fixes the file.
# file should be a .vtm file
def correctVTM(file:str) -> None:
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

# find the vtk folder name and corresponding time for a vtm file
# file should be a .vtm file
# the folder label and time are both strings
def readVTM(file:str) -> Tuple[str, str]:
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
    
# this creates a new .vtk.series or .vtm.series file
# tlist is a list of times
# flist is a list of folders. The times and folders don't need to be sorted ahead of time
# cf is the casefolder, e.g. 'C:\\...\nb64'
# ending is the file extension, either '.vtm' or '.vtk'
def generateVTKSeries(tlist:List[str], flist:List[str], cf:str, ending:str) -> None:
    seriesfile = series(cf)
    if os.path.exists(seriesfile):
        cfbasename = os.path.basename(seriesfile).replace(ending+'.series', '') # e.g. 'case'
    else:
        cfbasename = os.path.basename(cf) # e.g. 'case' or 'nb64'
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

# rewrite the .vtm.series file to include all .vtm files in the vtm folder
# folder should be the folder for the simulation (e.g. 'C:\\...\nb64')
# this uses the existing vtk and vtm files in the folder
def redoVTKSeriesNoLog(folder:str) -> None:
    cf = caseFolder(folder) # sometimes the case file is below the folder and named 'case', 
                            # and sometimes it is the folder and named something like 'nb64'. 
                            # This allows flexibility for those two situations
    if cf=='':
        return
    cfbasename = os.path.basename(cf) # e.g. 'case' or 'nb64'
    vtkfolder = os.path.join(cf, 'VTK')
    flist = [] # folder numbers
    tlist = [] # times
    ending = '.vtm'
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if file.endswith('.vtm'):
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
        generateVTKSeries(tlist, flist, cf, ending)
    return

#-------------------------------------------------------------------------------------------------  
##### FILE HANDLING


# exportFile exports a text file
    # folder is a top folder
    # file is a file base name within that folder
    # text is the text to write to the file
def exportFile(folder:str, file:str, text:str) -> None:
    fn = os.path.join(folder, file)
    File_object = open(fn,"w")
    File_object.write(text)
    File_object.close()
    logging.info("Exported file %s" % fn)

# exportCSV exports a csv file
    # fn is the full path of the file to export
    # table is the table to export
    # this function needs import csv
def exportCSV(fn:str, table:List[Any]) -> None:
    with open(fn, mode='w', newline='') as f:
        w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in table:
            w.writerow(i)
    logging.info('Exported file %s' % fn)

# importIf imports an existing legend file if it exists 
    # folder is a full path, e.g. 'C:\\...\nb64'
def importIf(folder:str) -> List[List[str]]:
    fn = os.path.join(folder, 'legend.csv')
    if os.path.exists(fn):
        with open(fn, 'r') as f:
            return(list(csv.reader(f)))
    else:
#         return scrape(folder).table()
        return []
   
# legendUnique imports a legend, 
# rewrites the variable names so they are all unique and ready to compile into a pandas dataframe, 
# and outputs a dictionary
def legendUnique(folder:str) -> Dict:
    t = importIf(folder)
    if len(t)==0:
        return t
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
    return {a[0]:a[1] for a in t}
        

# make a directory if it doesn't exist
# path is a full path name
# returns 1 if error, returns 0 for no error
def mkdirif(path:str) -> int:
    try:
        os.mkdir(path)
    except OSError as e:
        return 1
    else:
        logging.info ("Created directory %s" % path)
    return 0

# copy directory src to directory dest
# both should be full folder paths
def copy(src:str, dest:str) -> None:
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

### interfacePoints
        
# get the time of the interfacepoints file
def iPFileTime(file:str) -> int:
    f2num = int(re.split('_|.csv', file)[-2])
    return f2num

# delete all files after this file in the interfacePoints folder
# this is useful if there is an error generating interfacePoints files
# ip is interfacepoints folder
# f2 is the first bad file
def rmIPFiles(ip:str, f2:str) -> None:
    f2num = iPFileTime(f2)
    currFile = os.path.join(ip, 'interfacePoints_t_'+str(f2num)+'.csv')
    while os.path.exists(currFile):
        os.remove(currFile)
        f2num+=1
        currFile = os.path.join(ip, 'interfacePoints_t_'+str(f2num)+'.csv')

        
# sort files in interfacepoints folder by time and return a list of times
# f is full path name 
def sortIPFolder(f:str) -> List[str]:
    l1 = os.listdir(f)
    l1.sort(key=iPFileTime)
    return l1

# get the second line in a csv
# csvfile needs to be a full path name
def secondLineInCSV(csvfile:str) -> str:
    file = open(csvfile)
    csv_reader = csv.reader(file)
    next(csv_reader)
    return next(csv_reader)


# determine if the interfacePoints files were correctly exported
# fromFolder is parent folder, full path
# f is simulation folder, just the basename
def assessFolder(fromFolder:str, f:str) -> None:
    s1 = 0
    s2 = 0
    f2 = ''
    excludelist = ['nb79', 'nb89']
    if f.startswith('nb') and f not in excludelist:
        ip = os.path.join(fromFolder, f, 'interfacePoints')
        if os.path.exists(ip):
            for f1 in sortIPFolder(ip):
                # if two consecutive files have the same size or 
                # the file size shrank by a lot, there was an error
                f1path = os.path.join(ip, f1)
                s1 = os.stat(f1path).st_size
                if (s1<0.5*s2 and s1<10000) or (s1==s2 and secondLineInCSV(f1path)==secondLineInCSV(f2)):
                    rmIPFiles(ip, f2)
#                     logging.info(f1path, s1)
                    raise Exception(f+' interface points error')
                s2 = s1
                f2 = f1path
    return
 
# for all folders inside fromFolder, determine if interface points were generated correctly    
def assessFolders(fromFolder:str) -> None:
    logging.info('Parent folder: ', fromFolder)
    for f in os.listdir(fromFolder):
        try:
            assessFolder(fromFolder, f)
        except:
            logging.info(f+' interface points error')
    return



   
# change the endtime in the controlDict
# returns 0 if the dictionary was changed, 1 if not
def modifycontroldict(folder:str, tend:float) -> int:
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
    
    
# this modifies old controlDict files so runTimeModifiable is true
# folder is full path, e.g. 'C:\...\nb64'
def modifycontroldictruntime(folder):
    cfi = os.path.join(folder, 'case')
    if os.path.exists(cfi):
        cf = cfi
    else:
        cf = folder
    cdfile = os.path.join(cf, 'system', 'controlDict')
    cdfilenew = os.path.join(cf, 'system', 'controlDict2')
    if not os.path.exists(cdfile):
        return
    with open(cdfile, 'r') as fold:
        with open(cdfilenew, 'w') as fnew:
            for line in fold:
                if line.startswith('runTimeModifiable'):
                    line = 'runTimeModifiable\tyes;\n'
                fnew.write(line)
    os.remove(cdfile)
    os.rename(cdfilenew, cdfile)
    print(cdfile)
