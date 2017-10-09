import os
import copy
from collections import Mapping
from collections import OrderedDict
import yaml
from datetime import datetime, timedelta
from copy import deepcopy
from wrf_workflow import configwriter_activity
import re
import ConfigParser
import ast
import argparse
import gsi_process
import ConfigParser
from datetime import datetime
import gsi_control

gsiobs_config = ConfigParser.ConfigParser()
gsiobs = '/home/jzanetti/github/da_controller/etc/gsiobs.ini'
gsiobs_config.read(gsiobs)

def valid_datetime(timestamp):
    '''turn a timestamp into a datetime object'''
    try:
        return datetime.strptime(timestamp, "%Y%m%d%H%M")
    except ValueError:
        msg = "Not a valid date: '{}'.".format(timestamp)
    raise argparse.ArgumentTypeError(msg)


def setup_parser():
    """run data assimilation control"""
    PARSER = argparse.ArgumentParser(
        description='data assimilation control')

    PARSER.add_argument('model', type=str, help="model name")
    PARSER.add_argument('analysis_time', type=valid_datetime, help="analysis time")
    PARSER.add_argument('config_file', type=str, help='config file path')
    PARSER.add_argument('work_dir', type=str, help='work dir')
    PARSER.add_argument('gsi_base', type=str, help='gsi dir')
    PARSER.add_argument('crtm_base', type=str, help='crtm dir')
    PARSER.add_argument('background_path', type=str, help='background path')
    PARSER.add_argument('--prepbufr', type=str, dest='prepbufr', help="prepbufr path", default=None)
    PARSER.add_argument('--refbufr', type=str, dest='refbufr', help="ref. bufr path", default=None)
    PARSER.add_argument('--velbufr', type=str, dest='velbufr', help="vel. path", default=None)
    
    return PARSER.parse_args(['nz8kmN-ECMWF-var', '2017092100', 
                              '/home/jzanetti/github/da_controller/test_conf/nz8kmN-ECMWF-var_test.yaml',
                              '/home/jzanetti/github/da_controller/test',
                              '/home/jzanetti/programs/comGSIv3.5_EnKFv1.1/comGSIv3.5_EnKFv1.1_cloud',
                              '/home/jzanetti/programs/CRTM_2.2.3',
                              '/home/jzanetti/gsi_directory/practice_13/run/wrf_inout',
                              '--prepbufr', '/home/jzanetti/gsi_directory/practice_13/run/prepbufr'])
    return PARSER.parse_args()


if __name__ == '__main__':
    args = setup_parser()
    
    if os.path.exists(args.work_dir) == False:
        os.makedirs(args.work_dir)
    
    gsi_obstst = 'gsi_obstst'
    gsi_process.write_gsi_obsinput(gsiobs_config, gsi_obstst)
    configwriter_activity.WriteNamelist(args.config_file, os.path.join(args.work_dir, 'namelist.gsi'), 'gsi', add_config=gsi_obstst)


    gsi_control.gsi_control(args.work_dir, args.gsi_base, args.crtm_base,
               args.background_path, 
               prepbufr = args.prepbufr, refbufr = args.refbufr, velbufr = args.velbufr)

    print 'done'
    