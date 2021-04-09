#!/usr/bin/env python
'''Loads configuration settings'''

import yaml
from box import Box
import sys
import os
import logging.config

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#----------------------------------------------

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
configdir = os.path.join(currentdir, 'configs')
if not os.path.exists(configdir):
    configdir = os.path.join(parentdir, 'configs')
    if not os.path.exists(configdir):
        raise FileNotFoundError(f"No configs directory found")
    

with open(os.path.join(configdir, "config.yml"), "r") as ymlfile:
    cfg = Box(yaml.safe_load(ymlfile))
    

# # os.makedirs(cfg.path.logs, exist_ok=True)

# logConfigFile = os.path.join(currentdir, cfg.path.log_config)
# if not os.path.exists(logConfigFile):
#     logConfigFile = os.path.join(parentdir, cfg.path.log_config)
#     if not os.path.exists(logConfigFile):
#         logConfigFile = os.path.join(configdir, cfg.path.log_config)
#         if not os.path.exists(logConfigFile):
#             raise FileNotFoundError(f"Log yaml configuration file not found in {cfg.path.log_config}")

# with open(logConfigFile, "r") as ymlfile:
#     log_config = yaml.safe_load(ymlfile)

# # Set up the logger configuration
# logging.config.dictConfig(log_config)