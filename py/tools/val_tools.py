#!/usr/bin/env python
'''Functions for importing csv'''

# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#----------------------------------------------

def tryfloat(val:Any) -> Any:
    try:
        val = float(val)
    except:
        pass
    if type(val) is str:
        return val.strip()
    else:
        return val