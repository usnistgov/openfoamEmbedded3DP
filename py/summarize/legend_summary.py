#!/usr/bin/env python
'''Functions for generating legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import os
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging, platform, socket, sys

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import file.file_handling as fh
import folder_scraper as fs

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)



#-------------------------------------------------------------------------------------------------  

def legendSummary(topFolders:str, exportFN:str) -> None:
    '''scrape all legends into one table'''
    o = []
    for topfolder in topFolders:
        for f in fh.simFolders(topfolder):
            l = fs.legendUnique(f)
            if len(l)>0:
                o.append(l)
    p = pd.DataFrame(o)
    p.to_csv(exportFN)
    print(f'Exported {exportFN}')