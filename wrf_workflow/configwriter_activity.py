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
from config_processor.model_config import *
from config_processor.namelist import *

def WriteNamelist(config_file, namelist_path, 
                  process_stage, add_config=None):
    with open(config_file, 'r') as fin:
        cfg = yaml.load(fin)

    #nml = namelist.Namelist(cfg[process_stage], filetype = process_stage, 
    #                        **({'gsi_additional_config': gsi_obsinput} if gsi_obsinput else {}))
    nml = Namelist(cfg[process_stage])
    nml.write_to_file(namelist_path,
                      **({'filetype': process_stage, 
                          'add_config': add_config} 
                         if add_config else {}))


def write_common(model, base_time, nudging_cutoff_time, go_back_hours, out_file):
    cfg = get_model_config(model)
    for i, clength in enumerate(cfg['general_config']['domain_periods']):
        cfg['general_config']['domain_periods'][i] = clength + go_back_hours
    cluster_cfg = get_cluster_config(model)
    cluster_cfg.update(cfg['general_config'])
    nml = CommonNamelist(config=model,
                            base_time=base_time,
                            **cluster_cfg)

    with open(out_file, 'w') as outfile:
        yaml.dump({'common': nml.params}, outfile, default_flow_style=False)

    attributes = ['lat_0', 'lat_1', 'lat_2', 'llcrnrlat', 'llcrnrlon',
                      'lon_0', 'urcrnrlat', 'urcrnrlon', 'EARTH_RADIUS_M']
    basemap = {}
    for domain in cfg['basemap_parameters'].iterkeys():
        basemap[domain] = {
            x: getattr(cfg['basemap_parameters'][domain], x, None)
            for x in attributes}
    nml = CommonNamelist(**basemap)
    
    with open(out_file, 'a') as outfile:
        yaml.dump({'basemap': nml.params}, outfile, default_flow_style=False)
        
    fdda_setups = {'run_fdda': cfg['fdda_control']['run_fdda'],
                       'go_back_hours': cfg['fdda_control']['go_back_hours'],
                       'nudging_cutoff_time': str(nudging_cutoff_time)}
    
    nml = CommonNamelist(**fdda_setups)

    with open(out_file, 'a') as outfile:
        yaml.dump({'fdda_control': nml.params}, outfile, default_flow_style=False)


    var_setups = {'da_core': cfg['var_control']['da_core'],
                      'max_goback_hours': cfg['var_control']['max_goback_hours'],
                      'da_domain_id': cfg['var_control']['da_domain_id'],
                      'da_win_half': cfg['var_control']['da_win_half'],
                      'cycle_freq': cfg['var_control']['cycle_freq'],
                      'radar_da': cfg['var_control']['radar_da'],
                      'dfi_radar': cfg['var_control']['dfi_radar'],}
    nml = CommonNamelist(**var_setups)

    with open(out_file, 'a') as outfile:
        yaml.dump({'var_control': nml.params}, outfile, default_flow_style=False)

    wrf_archive_setups = {'archive_control': cfg['archive_control']['archive_prog']}
    
    nml = CommonNamelist(**wrf_archive_setups)

    with open(out_file, 'a') as outfile:
        yaml.dump({'archive_control': nml.params}, outfile, default_flow_style=False)


def write_wrf(model, base_time, work_dir, go_back_hours):
    nml = WRFNamelist(config=model,
                                base_time=base_time,
                                hourly_output_dir=work_dir,
                                go_back_hours=go_back_hours,
                                metgrid_path='met_em.d',
                                make_daily_subdir=False)
    
    return nml

def write_gsi(model):
    nml = GSINamelist(config=model)
    
    return nml
