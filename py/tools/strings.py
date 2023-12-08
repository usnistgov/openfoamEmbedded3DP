#!/usr/bin/env python
'''Handling strings'''

# external packages
import sys
import os
import csv
import numpy as np
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#################################################################

class varNicknames:
    '''a class used to hold aliases for variable names'''
    
    def __init__(self):
        self.vars = {'nozzle_inner_width':'niw'
                     , 'nozzle_outer_width':'now'
                  , 'nozzle_thickness':'nt'
                  , 'bath_width':'bw'
                  , 'bath_depth':'bd'
                    , 'nozzle_length':'nl'
                  , 'bath_left_coord':'blc'
                  , 'bath_right_coord':'brx'
                    , 'bath_front_coord':'bfc'
                  , 'bath_back_coord':'bbackc'
                  , 'bath_bottom_coord':'bbotc'
                  , 'bath_top_coord':'btc'
                   , 'nozzle_bottom_coord':'nbz'
                  , 'nozzle_center_x_coord':'ncx'
                  , 'nozzle_center_y_coord':'ncy'
                  , 'nozzle_angle':'na'
                   , 'bath_velocity':'bv'
                  , 'ink_velocity':'iv'
                    , 'adjacent_filament_offset':'Ddis'
                    , 'adjacent_filament_orientation':'Ddir'}
        self.names = {'imsize':'s'
                 , 'fontsize':'fs'
                 , 'HerschelBulkley':'HB'
                , 'Newtonian':'N'
                , '_list':''
                     , 'transportModel':'tm'
                     , 'yvar':'y'
                     , 'xvar':'x'
                        
                     }
        self.headers = {'vx':'$x$ velocity (mm/s)', 
                        'vz':'$z$ velocity (mm/s)',
                        'nu':'Viscosity (Pa$\cdot$s)',
                        'p':'Pressure (Pa)'}
        self.plainEng = {'sigma':'Surface tension'
                         ,'nuink':'Ink viscosity'
                         ,'nusup':'Bath viscosity'
                         ,'nozzle_inner_width':'Nozzle inner diameter'
                         ,'nozzle_angle':'Nozzle angle'
                         , 'vink':'Flow speed'
                         , 'vsup':'Move speed'
                         , 'iv':'Flow speed'
                         , 'bv':'Move speed'
                         , 'ink_velocity':'Flow speed'
                         , 'sup_velocity':'Move speed'
                         , 'spacing':'Spacing'
                         , 'viscRatio':'Ink viscosity/bath viscosity'
                         , 'ReRatio':'Ink Re/bath Re'
                         , 'aspectration':'h/w/intended'
                         , 'aspectratio':'Height/width'
                         , 'speeddecay':'Speed/intended'
                         , 'vertdispn':'Bottom position/intended'
                         , 'arean':'Area/intended'
                         , 'adjacent_filament_orientation':'Shift direction'
                         , 'ink_transportModel':'Ink model'
                         , 'sup_transportModel':'Bath model'
                         , 'asymmetryh':'Left-heaviness'
                         , 'asymmetryv':'Bottom-heaviness'
                         , 'centeryn':'Right shift'
                         , 'centerzn':'Up shift'
                         , 'maxwidthn':'w/intended'
                         , 'maxheightn':'h/intended'
                         , 'xbehind':'Position'
                         , 'ralpha':'Volume Fraction Residual'
                         , 'simTime':'Time (s)'
                         , 'Shift direction=y':'In-plane'
                        , 'Shift direction=z':'Out-of-plane'
                      , 'Dir=y':'In-plane'
                        , 'Dir=z':'Out-of-plane'
                        , 'Shift direction: y':'In-plane'
                        , 'Shift direction: z':'Out-of-plane'
                        }
        self.symbols = {'sigma':'\u03C3'
                        ,'Surface tension':'\u03C3'
                        ,'Ink viscosity':'$\u03B7_{ink}$'
                        ,'Bath viscosity':'$\u03B7_{sup}$'
                        , 'Ink viscosity/bath viscosity':'$\u03B7_{ink}/\u03B7_{sup}$'
                        , 'viscosity':'\u03B7'
                        , 'nuink':'$\u03B7_{ink}$'
                        , 'nusup':'$\u03B7_{sup}$'
                        , 'nu':'\u03B7'
                        , 'Ink Re/bath Re': '$Re_{ink}/Re_{sup}$'
                        , 'Nozzle inner diameter':'$d_i$'
                        , 'niw':'$d_i$'
                        , 'nid':'$d_i$'
                        , 'now':'$d_o$'
                        , '^2':'$^2$'
                        , 'Nozzle angle':'$\\theta$'
                        , 'degrees':'$\degree$'
                        , '*':'$\cdot$'
                        , 'Flow speed':'$v_{flow}$'
                        , 'Move speed':'$v_{move}$'
                        , 'Height/width':'h/w'
                        , 'Bottom position/intended':'Pos/intended'
                        , 'Shift direction':'Dir'
                        , ' model':''
                        , 'Volume Fraction':'$f_{ink}$'
                        , 'Residual':'res'
                        , 'HerschelBulkley':'HB'
                        , 'Newtonian':'Newt'
                       }
        
    def replace(self, s:str, d:dict) -> str:
        if not type(s) is str:
            return s
        k = s
        for key,val in d.items():
            k = k.replace(key,val)
        return k
        
    def shorten(self, s:str) -> str:
        '''for shortening labels without using special characters'''
        return self.replace(s, {**self.vars, **self.names})
    
    def labDict(self, yvar:str) -> str:
        '''dictionary for variable headers in line traces'''
        if yvar in self.headers:
            return self.headers[yvar]
        else:
            return yvar
                         
    def toEnglish(self, s:str) ->str:
        return self.replace(s, self.plainEng)
        
    def symbolic(self, s:str) -> str:
        '''for shortening labels with special characters'''
        return self.replace(s, self.symbols)
    
    def shortSymbol(self, s:str) -> str:
        return self.symbolic(self.shorten(self.toEnglish(s)))
    
    def hasSpacing(self, s:str) -> bool:
        '''determine if a spacing variable is present in this variable name s'''
        if 'adjacent_filament_offset' in s or 'spacing' in s:
            return True
        else:
            return False
        
    def sigmaVelocity(self, s:str) -> bool:
        '''get a plain english value for a sigma,ink_velocity value'''
        spl = re.split(', ', s)
        sigma = spl[0]
        flow = {'0':'no flow', '10':'flow'}[spl[1]]
        return f'$\sigma=${sigma}\n{flow}'
    
#=====================================================================================
    
def expFormat(x:float) -> str:
    '''display a number in exponential format'''
    try:
        p0 = np.log10(x)
        p = round(p0)
        if abs(p-p0)<10**-8:
            return r'$10^{{{}}}$'.format(int(p))
        else:
            return x
    except:
        return x


def decideFormat(x:float) -> Tuple[Union[float, str], bool]:
    '''determine what kind of format to put the number in
    might return float or string. second value is whether it got turned to exponent'''
    if type(x) is tuple or type(x) is np.array or type(x) is list:
        return x, False
    if type(x) is str and ',' in x:
        return x, False
    if x==0:
        return x,False
    try:
        l = np.log10(x)
    except TypeError:
        return x,False
    if not l==round(l):
        return x,False
    else:
        return expFormat(x),True


def expFormatList(xlist:List[Any], returnPrecision:bool=False, forceFormat:bool=False, useExp:bool=True, prec:int=-1000) -> List[Any]:
    '''put the whole list in exponential or appropriate format
    returnPrecision=True to return an int that describes the precision of the list. if -1000, then not formatted. if 1000, then exponent format.
    forceFormat forces the list to be formatted even if there is only one element in the list
    useExp=True to use exponential form, otherwise do not allow exponential
    '''
    xout = []
    xlist = list(xlist)
    
    if prec>-1000 and prec<1000:
        # already determined precision
        xout = [round(float(x),prec) for x in xlist]
        if returnPrecision:
            return xout, prec
        else:
            return xout  

    if len(xlist)<2 and not forceFormat:
        # only one element in list, don't format it
        if returnPrecision:
            return xlist, -1000   # not formatted
        else:
            return xlist
    for x in xlist:
        if useExp:
            x1, useExp = decideFormat(x)
        else:
            x1 = x
        xout.append(x1)
        prec = 1000    # exponent format
    if not useExp:
        # not using exponents. this is for strings, floats, ints
        ints = True
        i = 0
        
        # determine if these are all ints
        while i<len(xlist) and ints:
            try:
                if not round(xlist[i])==xlist[i]:
                    ints=False
            except:
                if returnPrecision:
                    return xlist, -1000 # not a number, return list
                else:
                    return xlist
            i+=1
            
        # set precision
        if ints:
            # set all to ints
            prec = 0
            xout = [int(xi) for xi in xlist]
        else:
            # these are floats. determine precision of float
            xsort = xlist.copy()
            xsort.sort()
            diffs = [t - s for s, t in zip(xsort, xsort[1:])] # differences between steps
            diffs = [i for i in diffs if i>0]
            if len(diffs)==0:
                prec=1
                found = False
                while prec<20 and not found:
                    if abs(round(xlist[0], prec)-xlist[0])<10**-10:
                        found=True
                    else:
                        prec+=1
            else:
                prec = int(np.ceil(-np.log10(min(diffs))))+1
            xout = [round(x,prec) for x in xlist]
    if returnPrecision:
        return xout, prec
    else:
        return xout  


    
    

    