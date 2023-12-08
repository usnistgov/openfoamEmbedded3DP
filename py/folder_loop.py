#!/usr/bin/env python
'''Functions for handling files'''

# external packages
import os, sys
import re
import shutil
import time
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import pandas as pd
import subprocess
import time
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentdir)
sys.path.append(os.path.dirname(currentdir))
from tools.config import cfg
import file_handling as fh
from file.plainIm import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#----------------------------------------------

class folderLoop:
    '''loops a function over all printFolders in the topFolder. 
    the function needs to have only one arg, folder, and all other variables need to go in kwargs
    folders could be either the top folder to recurse into, or a list of folders'''
    
    def __init__(self, folders:Union[str, list], func, printTraceback:bool=False, printErrors:bool=True, folderDiag:int=0, **kwargs):
        if type(folders) is list:
            # list of specific folders
            self.folders = []
            for folder in folders:
                self.folders = self.folders + fh.simFolders(folder)
        elif not os.path.exists(folders):
            self.topFolder = ''
            self.folders = []
        else:
            # top folder, recurse
            self.topFolder = folders
            self.folders = fh.simFolders(folders)
        self.func = func
        self.kwargs = kwargs
        self.printTraceback = printTraceback
        self.printErrors = printErrors
        self.folderDiag = folderDiag
        
    def runFolder(self, folder:str) -> None:
        '''run the function on one folder'''
        if self.folderDiag>0:
            print(folder)
        try:
            self.func(folder, **self.kwargs)
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            self.folderErrorList.append({'folder':folder, 'error':e})
            if self.printErrors:
                print(e)
            if self.printTraceback:
                traceback.print_exc()

        
    def run(self) -> list:
        '''apply the function to all folders'''
        self.folderErrorList = []
        for folder in self.folders:
            self.runFolder(folder)
        return self.folderErrorList
    
    def testFolderError(self, i:int, **kwargs) -> None:
        '''test a single file that threw an error. i is the row in self.folderErrorList'''
        if i>len(self.folderErrorList):
            print(f'{i} is greater than number of error files ({len(self.folderErrorList)})')
            return
        row = self.folderErrorList[i]
        print(row)
        self.func(row['folder'], **kwargs)
        
    def exportErrors(self, fn:str) -> None:
        '''export the error list to file'''
        plainExp(fn, pd.DataFrame(self.folderErrorList), {'folder':'', 'error':''}, index=False)
    