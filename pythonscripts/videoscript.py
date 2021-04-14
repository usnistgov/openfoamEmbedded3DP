#!/usr/bin/env python
'''Script for combining images into videos. Loops continuously every 6 hours.'''

import os, sys
import time

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import folderparser as fp
from videofuncs import saveVid
from config import cfg

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

            
#-----------------------------------------------            
            
            
if __name__ == "__main__":
    serverfolder = os.path.join(cfg.path.server, 'yieldingsweep', 'HBHByielded')
    topfolders = [os.path.join(serverfolder,s) for s in ['k', 'n', 'tau0']]
    while True:
        try:
            for topfolder in topfolders:
                for folder in fp.caseFolders(topfolder):
                    for p in ['umag']:
                        for s in ['y']:
                            try:
                                saveVid(folder, s, p)      
                            except:
                                pass
        except:
            pass
        print('waiting 6 hrs for more files...')
        fp.printCurrentTime()
        time.sleep(60*60*6)