import yaml
from box import Box

with open("config.yaml", "r") as ymlfile:
    cfg = Box(yaml.safe_load(ymlfile))