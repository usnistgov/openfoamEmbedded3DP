#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging

# local packages


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  


class DictList:
    '''DictList is an object that ensures uniform formatting of variables into files'''
        
    def __init__(self, title:str, form:int, proplist:List[str]):
        '''Inputs: title, form, proplist
        title is the title of the variable group
        form is 0 for bracketed entries e.g. group definitions, 1 for parentheses e.g. point lists
        proplist is a list of properties that we're storing within this group'''
        self.title = title
        self.proplist = proplist
        if form==0: 
            self.pl = "{" # parens (group) left
            self.pr = "}" # parens (group) right
            self.sp = "\t" # spacing
            self.rl = "" # row left
            self.rr = ";" # row right
        else:
            self.pl = "("
            self.pr = ");"
            self.sp = " "
            self.rl = "("
            self.rr = ")"
            
    def tabs(self, n:int) -> str:
        '''return n tabs in a row'''
        s = ""
        for i in range(n):
            s = s + "\t"
        return s
    
    
    def prnt(self, level:int) -> str:
        '''Format the title and variables into the format that OpenFOAM likes. 
        Input: the number of tabs before the title of the group. If you want no title, use level -1
        Output: formatted string'''
        if level<0:
            s = ""
        else:
            s = f'{self.tabs(level)}{self.title}\n{self.tabs(level)}{self.pl}\n'
        for p in self.proplist:
            if isinstance(p, DictList):  # if this property is a DictList, recurse
                s = s + p.prnt(level + 1)
            else: # row level
                s = s + self.tabs(level + 1)
                if isinstance(p, str): # single item
                    s = f'{s}{p}\n'
                else: # list of items
                    if len(p)>0:
                        # one row
                        s = s + self.rl
                        for pi in p[0:-1]:
                            s = f'{s}{pi}{self.sp}'
                        s = f'{s}{p[-1]}{self.rr}\n'
                    else:
                        s = f'{s}\n'
                if level<0:
                    s = f'{s}\n'
        if level>=0:
            s = f'{s}{self.tabs(level)}{self.pr}\n\n'
        return s