#!/usr/bin/env python
'''Functions for importing csv'''

# external packages
import os
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import numpy as np

# local packages


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#----------------------------------------------

def plainIm(file:str, ic:Union[int, bool]=0, checkUnits:bool=True) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
    '''import a csv to a pandas dataframe. ic is the index column. Int if there is an index column, False if there is none. checkUnits=False to assume that there is no units row. Otherwise, look for a units row'''
    if os.path.exists(file):
        try:
            toprows = pd.read_csv(file, index_col=ic, nrows=2)
            
            toprows = toprows.fillna('')
            toprows.columns = map(str.lower, toprows.columns) # set headers to lowercase
            row1 = list(toprows.iloc[0])
            if checkUnits and all([(type(s) is str or pd.isnull(s)) for s in row1]):
                # row 2 is all str: this file has units
                unitdict = dict(toprows.iloc[0])
                skiprows=[1]
            else:
                unitdict = dict([[s,'undefined'] for s in toprows])
                skiprows = []
            try:
                d = pd.read_csv(file, index_col=ic, dtype=float, skiprows=skiprows)
            except:
                d = pd.read_csv(file, index_col=ic, skiprows=skiprows)
        except Exception as e:
#             logging.error(str(e))
            return [],{}
        d.columns = map(str.lower, d.columns) # set headers to lowercase
        return d, unitdict
    else:
        return [], {}
    
    
def plainExp(fn:str, data:pd.DataFrame, units:dict) -> None:
    '''export the file'''
    if len(data)==0 or len(units)==0:
        return
    col = pd.MultiIndex.from_tuples([(k,units[k]) for k in data]) # index with units
    data = np.array(data)
    df = pd.DataFrame(data, columns=col)       
    df.to_csv(fn)
    logging.info(f'Exported {fn}')