import yaml
from box import Box
import sys
import os
import logging.config
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