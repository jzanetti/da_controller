import os
import copy
from collections import Mapping
from collections import OrderedDict
import yaml
from datetime import datetime, timedelta
from copy import deepcopy
#from config_processor import namelist
import re
import ConfigParser
import ast

CONFIG_DIR = {'wrf': '/home/jzanetti/github/da_controller/config_processor/etc'}
CLUSTER_CONFIG = '{}_cluster_config.yml'

def nwp_config(model, process):
    config_file = os.path.join(CONFIG_DIR[process], model + '.py')
    if os.path.exists(config_file):
        return config_file

def get_model_config(model, process='wrf'):
    '''Recursively update config dictionary from eldest parent downwards
    (i.e. children overwrite their ancestors).'''
    if model is not None:
        config_namespace = {}
        with open(nwp_config(model, process=process)) as cfg:
            exec(cfg.read(), config_namespace)
        return update(get_model_config(config_namespace['parent_config']),
                          config_namespace['config'])
    else:
        return {}

def update(dict1, dict2):
    '''Recursively update dict1 with dict2.'''
    new_dict = copy.deepcopy(dict1)
    for k, v in dict2.items():
        if isinstance(v, Mapping):
            r = update(new_dict.get(k, {}), v)
            new_dict[k] = r
        else:
            new_dict[k] = dict2[k]
    return new_dict

def get_cluster_config(model, workflow='wrf'):
    with open(os.path.join(CONFIG_DIR['wrf'],
        CLUSTER_CONFIG.format(workflow.lower())), 'r') as config_file:
        cfg = yaml.load(config_file)
    try:
        model_cfg = cfg[model.split('-')[0]]
    except KeyError:
        model_cfg = cfg['default']
    if model in model_cfg['config']:
        model_cfg.pop('config')
    else:
        model_cfg = cfg['default']
        model_cfg.pop('config')
    return model_cfg
