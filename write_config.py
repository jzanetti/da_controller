import os
import copy
from collections import Mapping
from collections import OrderedDict
import yaml
from datetime import datetime, timedelta
from copy import deepcopy
from config_processor import namelist
import re
import ConfigParser
import ast

CONFIG_DIR = {'wrf': '/home/szhang/GitHub_branches/var_dir/config_processor_test/etc/wrf'}
CLUSTER_CONFIG = '{}_cluster_config.yml'
gsiobs = '/home/szhang/workspace/var_controller/etc/gsiobs.ini'


gsiobs_config = ConfigParser.ConfigParser()
gsiobs_config.read(gsiobs)
filename = 'gsi_obstst'


def write_gsi_obsinput(gsiobs_config, filename):

    class obsinput_type:
        dtype = dplat = dsis = dval = dthin = dsfcalc = None

    def return_obs_header_end(f_gsiobs, gsi_obsinput, obsout_format, 
                              header=True):
        if header:
            f_gsiobs.write('{}\n!  '.format('OBS_INPUT::'))
            for i in range(len(gsi_obsinput)):
                f_gsiobs.write(obsout_format[i].format(gsi_obsinput[i]))
            f_gsiobs.write('\n')
        else:
            f_gsiobs.write('::\n')
            

    def return_obs_key(req_list, dfile, dtype):
        if dtype in req_list[dfile].keys():
            return req_list[dfile][dtype]
        else:
            if req_list[dfile]['default'] == 'var':
                return dtype
            else:
                return req_list[dfile]['default']

    gsi_obsinput = ['dfile', 'dtype', 'dplat', 
                    'dsis', 'dval', 'dthin', 'dsfcalc']
    obsout_format = ['{0:15}', '{0:12}', '{0:10}', 
                     '{0:21}', '{0:8}', '{0:6}', '{0:7}']

    gsi_io = obsinput_type()
    for obs in gsi_obsinput[1:]:
        setattr(gsi_io, obs, ast.literal_eval(gsiobs_config.get("GSIOBS",obs)))
    
    with open(filename, 'w') as f_gsiobs:  
        return_obs_header_end(f_gsiobs, gsi_obsinput, obsout_format)
        
        for dfile in gsi_io.dtype.keys():
            for dtype in gsi_io.dtype[dfile]:
                obsout = [dfile, dtype]
                for obsin in gsi_obsinput[2:]:
                    obsout.append(return_obs_key(getattr(gsi_io, obsin), dfile, dtype))
                f_gsiobs.write('   ')
                for i in range(len(gsi_obsinput)):
                    f_gsiobs.write(obsout_format[i].format(str(obsout[i])))
                f_gsiobs.write('\n')
        
        return_obs_header_end(f_gsiobs, gsi_obsinput, obsout_format, header=False)
            
    f_gsiobs.close()

def WriteNamelist(config_file, namelist_path, 
                  process_stage, add_config=None):
    with open(config_file, 'r') as fin:
        cfg = yaml.load(fin)

    #nml = namelist.Namelist(cfg[process_stage], filetype = process_stage, 
    #                        **({'gsi_additional_config': gsi_obsinput} if gsi_obsinput else {}))
    nml = namelist.Namelist(cfg[process_stage])
    nml.write_to_file(namelist_path,
                      **({'filetype': process_stage, 
                          'add_config': add_config} 
                         if add_config else {}))

def dict_merge(a, b):
    # Merge two dict that are supposed to have the same architecture
    # It isn't a general deep merge but designed for efficiency considerations
    if not isinstance(b, dict):
        return
    for k, v in b.items():
        if isinstance(v, dict) and k in a:
            dict_merge(a[k], v)
        else:
            a[k] = v

def str_to_type(var):
    """From a string return an object of the appropriate type"""
    try:
        if int(var) == float(var):
            return int(var)
    except ValueError:
        try:
            float(var)
            return float(var)
        except ValueError:
            if ',' in var:
                return [str_to_type(x.strip()) for x in
                        var.rstrip(',').split(',')]
            else:
                if '.true.' in var:
                    return True
                elif '.false.' in var:
                    return False
                else:
                    return var.replace("'", "")

def type_to_str(var, precision):
    """Return a string representation suitable for including in namelist"""
    if isinstance(var, list):
        return ', '.join(type_to_str(x, precision) for x in var)
    if isinstance(var, str):
        return "'{}'".format(var)
    elif isinstance(var, float):
        return str('{:.{}f}'.format(var, precision))
    elif isinstance(var, bool):
        if var:
            return '.true.'
        else:
            return '.false.'
    else:
        return str(var)

class Namelist(object):

    """Object to hold information about a Fortran namelist"""
    precision = 6

    def __init__(self, namelist=None):
        self.namelist = namelist

    def read_from_file(self, filename):
        """Read data from a namelist file"""
        self.namelist = OrderedDict()
        with open(filename) as nml:
            current_section = None
            for line in nml:
                line = line.strip().split('=')
                if line[0].startswith('&'):
                    current_section = line[0].lstrip('&')
                    self.namelist[current_section] = OrderedDict()
                if len(line) > 1 and current_section != None:
                    entry = ''
                    values = line[1].strip().rstrip(',').split(',')
                    entry = []
                    for value in values:
                        if value != '':
                            entry.append(str_to_type(value.strip()))
                    self.namelist[current_section][line[0].strip()] = entry
        return self.namelist

    def write_to_file(self, filename):
        """Write data to a namelist file"""
        try:
            with open(filename, 'w') as nml:
                for section, values in self.namelist.items():
                    nml.write('&{}\n'.format(section))
                    for key, value in values.items():
                        nml.write(' {} = {}\n'.format(key,
                                            type_to_str(value, self.precision)))
                    nml.write('/\n\n')
        except (IOError, OSError):
            print 'failed'

    def update_var(self, section, key, var):
        """Check if environment variable has a value
        and if so update namelist"""
        if os.getenv(var):
            self.namelist[section][key] = str_to_type(os.getenv(var))

    def output(self):
        """Return a 'nice' version of the namelist suitable for printing"""
        output = ''
        for section in self.namelist.keys():
            output += '{}\n'.format(section)
            for key, value in self.namelist[section].items():
                output += '\t{}: {}\n'.format(key, value)
        return output

class GSINamelist(Namelist):
    default_settings = OrderedDict([
        ('setup', OrderedDict([
            ('miter', 3),
            ('niter(1)', 50), ('niter(2)', 50),
            ('write_diag(1)', True), ('write_diag(2)', False), ('write_diag(3)', True),
            ('gencode', 78), ('qoption', 2),
            ('factqmin', 0.0), ('factqmax', 0.0),
            ('iguess', -1),
            ('oneobtest', False), ('retrieval', False),
            ('nhr_assimilation', 3), ('l_foto', False),
            ('use_pbl', False),('lwrite_peakwt', True),
            ('lread_obs_save', False), ('lread_obs_skip', False),
            ('newpc4pred', True), ('adp_anglebc', True), ('angord', 4),
            ('passive_bc', True), ('use_edges', False), ('emiss_bc', True),
            ('diag_precon', True), ('step_start', 1.e-3), ('use_prepb_satwnd', True),
            ])),
        ('gridopts', OrderedDict([
           ('JCAP', 62),('JCAP_B', 62),('NLAT', 60),('NLON',60), ('nsig', 60), ('regional', True),
           ('wrf_nmm_regional', False), ('wrf_mass_regional', False),
           ('nems_nmmb_regional', False), ('nmmb_reference_grid', 'H'), ('diagnostic_reg', False),                                                            
           ('filled_grid', False), ('half_grid', True), ('netcdf', True)
            ])),
        ('bkgerr', OrderedDict([
           ('vs', 0),('hzscl', 0),('bw', 0),('fstat',True)])),
        ('obsqc', OrderedDict([
            ('dfact', 0.75),
            ('dfact1', 3.0), 
            ('noiqc', False), 
            ('c_varqc', 0.02), 
            ('vadfile', 'prepbufr')])),
        ('obs_input', OrderedDict([
            ('dmesh(1)', 120.0),
            ('dmesh(2)', 60.0), 
            ('dmesh(3)', 30), 
            ('time_window_max', 1.5), 
            ('ext_sonde', True)])),
        ('superob_radar', OrderedDict([
            ('del_azimuth', 5),
            ('del_elev', 0.25), 
            ('del_range', 5000), 
            ('del_time', 0.5), 
            ('elev_angle_max', 5),
            ('minnum', 50),
            ('range_max', 100000),
            ('l2superob_only', False),
            ])),
        ('hybrid_ensemble', OrderedDict([
           ('l_hyb_ens', False)])),
        ('rapidrefresh_cldsurf', OrderedDict([
           ('i_gsdcldanal_type', 1),
           ('i_lightpcp', 1),
           ])),   
        ])
    
    def __init__(self, config):
        ''' Example:
                config: the domain specific config name e.g nz8kmN-ECMWF
                domain: nz8kmN
        '''
        self.model_config_etc = get_model_config(config)
        gsi_conf = self.model_config_etc['gsi_namelist']
        self.namelist = deepcopy(self.default_settings)
        
        for section in self.namelist.keys():
            if section in gsi_conf.keys():
                for gsiopt in self.namelist[section]:
                    if gsiopt in gsi_conf[section].keys():
                        self.namelist[section][gsiopt] = gsi_conf[section][gsiopt]

class WRFNamelist(Namelist):

    """WRF namelist object"""
    precision = 6

    default_settings = OrderedDict([
        ('time_control', OrderedDict([
            ('run_days', 0),
            ('run_hours', 0),
            ('run_minutes', 0),
            ('run_seconds', 0),
            ('start_year', None),
            ('start_month', None),
            ('start_day', None),
            ('start_hour', None),
            ('start_minute', None),
            ('start_second', None),
            ('end_year', None),
            ('end_month', None),
            ('end_day', None),
            ('end_hour', None),
            ('end_minute', None),
            ('end_second', None),
            ('interval_seconds', None),
            ('input_from_file', True),
            ('history_interval', 60),
            ('frames_per_outfile', 1000),
            ('restart', False),
            ('restart_interval', 9000),
            ('io_form_history', 2),
            ('io_form_restart', 2),
            ('io_form_input', 2),
            ('io_form_boundary', 2),
            ('debug_level', 0),
            ('auxinput1_inname', None),
            ('auxinput11_interval_s', 120),
            ('auxinput11_end_h', 10),
            ('output_diagnostics', 1),
            ('auxhist3_outname', 'wrf_diag_d<domain>_dummy'),
            ('io_form_auxhist3', 2),
            ('auxhist3_interval', 60), # need to be hourly for wrf2grib
            ('frames_per_auxhist3', 1000),
            ('auxhist2_outname', None),
            ('io_form_auxhist2', 2),
            ('auxhist2_interval', 0),
            ('frames_per_auxhist2', 1),
            ('ignore_iofields_warning', False)])),
        ('domains', OrderedDict([
            ('time_step', None),
            ('time_step_fract_num', None),
            ('time_step_fract_den', None),
            ('max_dom', 1),
            ('s_we', 1),
            ('e_we', None),
            ('s_sn', 1),
            ('e_sn', None),
            ('s_vert', 1),
            ('e_vert', None),
            ('num_metgrid_levels', None),
            ('num_metgrid_soil_levels', None),
            ('dx', None),
            ('dy', None),
            ('grid_id', 1),
            ('parent_id', 1),
            ('i_parent_start', None),
            ('j_parent_start', None),
            ('parent_grid_ratio', 1),
            ('parent_time_step_ratio', 1),
            ('feedback', 1),
            ('smooth_option', 0),
            ('p_top_requested', 10000),
            ('numtiles', 1),
            ('nproc_x', -1),
            ('nproc_y', -1),
            ('rh2qv_wrt_liquid', True)])),
        ('physics', OrderedDict([
            ('mp_physics', 6),
            ('ra_lw_physics', 1),
            ('ra_sw_physics', 1),
            ('radt', 30),
            ('sf_sfclay_physics', 1),
            ('sf_surface_physics', 1),
            ('pxlsm_smois_init', 0),
            ('bl_pbl_physics', 1),
            ('bldt', 0),
            ('cu_physics', 1),
            ('cudt', 5),
            ('isfflx', 1),
            ('ifsnow', 0),
            ('icloud', 1),
            ('surface_input_source', 1),
            ('num_soil_layers', 5),
            ('sf_urban_physics', 0),
            ('mp_zero_out', 0),
            ('maxiens', 1),
            ('maxens', 3),
            ('maxens2', 3),
            ('maxens3', 16),
            ('ensdim', 144)])),
        ('fdda', OrderedDict([
            ('grid_fdda', 0),
            ('gfdda_inname', 'wrffdda_d<domain>'),
            ('gfdda_interval_m', 360),
            ('gfdda_end_h', 0),
            ('grid_sfdda', 0),
            ('sgfdda_inname', 'wrfsfdda_d<domain>'),
            ('sgfdda_interval_m', 60),
            ('sgfdda_interval_s', 3600),
            ('sgfdda_end_h', 0),
            ('io_form_gfdda', 2),
            ('fgdt', 0),
            ('if_no_pbl_nudging_uv', 0),
            ('if_no_pbl_nudging_t', 1),
            ('if_no_pbl_nudging_q', 1),
            ('if_zfac_uv', 0),
            ('k_zfac_uv', 10),
            ('if_zfac_t', 1),
            ('k_zfac_t', 10),
            ('if_zfac_q', 1),
            ('k_zfac_q', 10),
            ('guv', 0.0003),
            ('gt', 0.0003),
            ('gq', 0.0003),
            ('if_ramping', 0),
            ('dtramp_min', 0.0),
            ('obs_nudge_opt', 0),
            ('max_obs', 150000),
            ('fdda_start', 0.),
            ('fdda_end', 600),
            ('obs_nudge_wind', 1),
            ('obs_coef_wind', 6.E-4),
            ('obs_nudge_temp', 1),
            ('obs_coef_temp', 6.E-4),
            ('obs_nudge_mois', 1),
            ('obs_coef_mois', 6.E-4),
            ('obs_rinxy', 120.),
            ('obs_rinsig', 0.1),
            ('obs_twindo', 0.666666),
            ('obs_npfi', 10),
            ('obs_ionf', 2),
            ('obs_idynin', 0),
            ('obs_dtramp', 40.),
            ('obs_prt_freq', 10),
            ('obs_ipf_errob', True),
            ('obs_ipf_nudob', True),
            ('obs_ipf_in4dob', True)])),
        ('dynamics', OrderedDict([
            ('w_damping', 1),
            ('diff_opt', 1),
            ('km_opt', 4),
            ('diff_6th_opt', 0),
            ('diff_6th_factor', 0.12),
            ('base_temp', 290.),
            ('damp_opt', 0),
            ('zdamp', 5000.),
            ('dampcoef', 0.01),
            ('khdif', 0),
            ('kvdif', 0),
            ('non_hydrostatic', True),
            ('moist_adv_opt', 1),
            ('scalar_adv_opt', 1)])),
        ('bdy_control', OrderedDict([
            ('spec_bdy_width', 5),
            ('spec_zone', 1),
            ('relax_zone', 4),
            ('specified', True),
            ('nested', False)])),
        ('dfi_control', OrderedDict([
            ('dfi_opt', 0),
            ('dfi_nfilter', 7),
            ('dfi_write_filtered_input', True),
            ('dfi_write_dfi_history', False),
            ('dfi_cutoff_seconds', 3600),
            ('dfi_time_dim', 1000),
            ('dfi_bckstop_year', 2001),
            ('dfi_bckstop_month', 06),
            ('dfi_bckstop_day', 11),
            ('dfi_bckstop_hour', 11),
            ('dfi_bckstop_minute', 0),
            ('dfi_bckstop_second', 0),
            ('dfi_fwdstop_year', 2001),
            ('dfi_fwdstop_month', 06),
            ('dfi_fwdstop_day', 11),
            ('dfi_fwdstop_hour', 11),
            ('dfi_fwdstop_minute', 0),
            ('dfi_fwdstop_second', 0),
            #('dfi_radar',0) # not supported
            ])),
        ('grib2', OrderedDict()),
        ('namelist_quilt', OrderedDict([
            ('nio_tasks_per_group', 0),
            ('nio_groups', 1)]))
    ])

    max_dom_fields = {'time_control': {'input_from_file': None,
                                       'history_interval': None,
                                       'frames_per_outfile': None,
                                       'auxinput11_interval_s': None,
                                       'auxinput11_end_h': None,
                                       'auxhist3_interval': None,
                                       'frames_per_auxhist3': None,
                                       'auxhist2_interval': 0,
                                       'frames_per_auxhist2': 1
                                       },
                      'domains': {'s_we': None,
                                  's_sn': None,
                                  's_vert': None,
                                  },
                      'physics': {'mp_physics': None,
                                  'ra_lw_physics': None,
                                  'ra_sw_physics': None,
                                  'radt': None,
                                  'sf_sfclay_physics': None,
                                  'sf_surface_physics': None,
                                  'pxlsm_smois_init': None,
                                  'bl_pbl_physics': None,
                                  'bldt': None,
                                  'cu_physics': None,
                                  'cudt': None,
                                  },
                      'fdda': {'grid_fdda': 0,
                               'gfdda_interval_m': None,
                               'gfdda_end_h': 0,
                               'grid_sfdda': 0,
                               'sgfdda_interval_m': None,
                               'sgfdda_interval_s': None,
                               'sgfdda_end_h': None,
                               'fgdt': None,
                               'if_no_pbl_nudging_uv': None,
                               'if_no_pbl_nudging_t': None,
                               'if_no_pbl_nudging_q': None,
                               'if_zfac_uv': None,
                               'k_zfac_uv': None,
                               'if_zfac_t': 0,
                               'k_zfac_t': None,
                               'if_zfac_q': 0,
                               'k_zfac_q': None,
                               'guv': None,
                               'gt': None,
                               'gq': None,
                               'obs_nudge_opt': None,
                               'fdda_start': None,
                               'fdda_end': None,
                               'obs_nudge_wind': None,
                               'obs_coef_wind': None,
                               'obs_nudge_temp': None,
                               'obs_coef_temp': None,
                               'obs_nudge_mois': None,
                               'obs_coef_mois': None,
                               'obs_rinxy': None,
                               'obs_twindo': None,
                               'obs_ionf': None,
                               'obs_prt_freq': None,
                               },
                      'dynamics': {'zdamp': None,
                                   'dampcoef': None,
                                   'khdif': None,
                                   'kvdif': None,
                                   'non_hydrostatic': None,
                                   'moist_adv_opt': None,
                                   'scalar_adv_opt': None,
                                   },
                      'bdy_control': {'specified': False,
                                      'nested': True,
                                      }
                      }

    def __init__(self, config, base_time, hourly_output_dir,
                 go_back_hours,
                 fdda_offset=None, metgrid_levels=None, metgrid_soil_levels=None,
                 metgrid_interval=None, rh2qv_wrt_liquid=None, metgrid_path=None,
                 make_daily_subdir=True, restart_files=None):

        ''' Example:
                config: the domain specific config name e.g nz8kmN-ECMWF
                domain: nz8kmN
        '''
        self.model_config_etc = get_model_config(config)
        domain_conf = self.model_config_etc['general_config']
        #global_conf = global_data_config[domain_conf['global_source']]
        global_conf = {'levels': 10, 'soil_levels':10, 'rh2qv_wrt_liquid':1, 'data_interval':1}
        
        # Define the start date (for restart run)
        if restart_files is None:
            start_dates = ','.join(domain_conf['number_of_domains']*
                             [datetime.strftime(base_time, '%Y-%m-%d_%H:%M:%S')])
        else:
            for restart_file in restart_files:
                if 'wrfrst' in restart_file:
                    restart_date = '_'.join(restart_files[0].split('_')[-2:])
                    break
            start_dates = ','.join(domain_conf['number_of_domains']*[restart_date])

        end_dates = ','.join([datetime.strftime(base_time+timedelta(hours=(p + go_back_hours)),
                        '%Y-%m-%d_%H:%M:%S') for p in domain_conf['domain_periods']])

        super(WRFNamelist, self).__init__()
        domain = config
        for suffix in ['-ECMWF', '-NCEP', '-UKMO']:
            domain = domain.split(suffix)[0]

        self.domain_config_etc = get_model_config(domain)
        self.update_domain_param(config, domain)
        self.namelist = deepcopy(self.default_settings)
        self.get_geogrid_params(domain)
        self.set_dates(start_dates, end_dates)
        if self.model_config_etc['var_control']['dfi_radar']:
            self.set_dfi(start_dates)
        if fdda_offset is None:
            fdda_offset = 0
        self.update_hourly_section(
            config, start_dates.split(',')[0], fdda_offset, hourly_output_dir,
            make_daily_subdir)
        self.update_domains_section()
        self.update_physics_section()
        self.update_bdy_control_section()
        # self.update_from_env() ########################### Fix this
        if metgrid_levels is not None:
            self.namelist['domains']['num_metgrid_levels'] = metgrid_levels
        else:
            self.namelist['domains']['num_metgrid_levels'] = \
                                                global_conf['levels']
        if metgrid_soil_levels is not None:
            self.namelist['domains']['num_metgrid_soil_levels'] = \
                                                metgrid_soil_levels
        else:
            self.namelist['domains']['num_metgrid_soil_levels'] = \
                                            global_conf['soil_levels']
        if rh2qv_wrt_liquid is not None:
            self.namelist['domains']['rh2qv_wrt_liquid'] = rh2qv_wrt_liquid
        else:
            self.namelist['domains']['rh2qv_wrt_liquid'] = \
                                     global_conf['rh2qv_wrt_liquid']
        if metgrid_interval is not None:
            self.namelist['time_control']['interval_seconds'] = metgrid_interval
        else:
            self.namelist['time_control']['interval_seconds'] = \
                                    3600*global_conf['data_interval']
        if metgrid_path is not None:
            self.namelist['time_control']['auxinput1_inname'] = \
                '{}<domain>.<date>'.format(metgrid_path)
        if self.namelist['domains']['max_dom'] > 1:
            self.extend_params()

    def update_domain_param(self, config, domain):
        """ Change or Add domain parameters with config specific parameters.
        i.e different options for different driving moedl and same domain. """
        dict_merge(self.domain_config_etc, self.model_config_etc.get(
                                                        config, None))

    def get_geogrid_params(self, config):
        """Read required parameters from geogrid"""
        for key in self.namelist.keys():
            if key in self.model_config_etc['wrf_namelist']:
                self.namelist[key].update(self.model_config_etc['wrf_namelist'][key])
        for key in ('e_we', 'e_sn', 'dx', 'dy', 'parent_id', 'i_parent_start',
                    'j_parent_start', 'parent_grid_ratio'):
            self.namelist['domains'][key] = self.model_config_etc['wps_namelist'][key]
        self.namelist['domains']['max_dom'] = len(
            self.model_config_etc['wps_namelist']['e_we'])

    def update_domains_section(self):
        """Update domains section of namelist"""
        self.namelist['domains']['parent_time_step_ratio'] = \
            self.namelist['domains']['parent_grid_ratio']
        self.namelist['domains']['grid_id'] = 1
        if self.namelist['domains']['max_dom'] > 1:
            self.namelist['domains']['grid_id'] = [x + 1 for x in range(
                self.namelist['domains']['max_dom'])]
            for i in range(1, self.namelist['domains']['max_dom']):
                parent_id = self.namelist['domains']['parent_id'][i]
                grid_ratio = self.namelist['domains']['parent_grid_ratio'][i]
                self.namelist['domains']['dx'].append(
                    self.namelist['domains']['dx'][parent_id - 1] / grid_ratio)
                self.namelist['domains']['dy'].append(
                    self.namelist['domains']['dy'][parent_id - 1] / grid_ratio)

    def update_hourly_section(self, config, start_date, fdda_offset,
                              hourly_output_dir, make_daily_subdir):
        date = datetime.strptime(start_date, "%Y-%m-%d_%H:%M:%S")
        date = date + timedelta(hours=fdda_offset)
        dir_day = os.path.join(hourly_output_dir, date.strftime('%y%m%d%H'))
        if not os.path.isdir(dir_day) and make_daily_subdir:
            try:
                os.makedirs(dir_day)
            except OSError:
                # This can append if the nfs share isn't mounted on the box
                # If this append we still want CHAMP to run without hourly output
                # This can be removed once we have introduced hourly output on
                # every box.
                self.hourly = False
                print 'failed2'
        # Here we define the filename with a .temp at the end of the name.
        # All files contain a .temp suffix will be renamed by WRF (wrf-share.patch)
        # after the connection is closed. Making the writing atomic.
        # i.e: the file to process will match this naming convention with
        # the suffix removed
        self.namelist['time_control']['auxhist2_outname'] = \
            str(os.path.join(
                dir_day, 'wrf_hourly_{}_d<domain>_<date>.temp'.format(config)))

    def update_physics_section(self):
        """Update physics section of namelist"""
        self.namelist['physics']['radt'] = max(5,
                                               self.namelist['domains']['dx'][-1] / 1000)

    def update_bdy_control_section(self):
        """Update bdy_control section of namelist"""
        spec_zone = self.namelist['bdy_control']['spec_zone']
        relax_zone = self.namelist['bdy_control']['relax_zone']
        spec_bdy_width = self.namelist['bdy_control']['spec_bdy_width']
        if spec_zone + relax_zone != spec_bdy_width:
            print 'failed3'
            self.namelist['bdy_control']['spec_bdy_width'] = spec_zone + \
                relax_zone

    def set_dates(self, start_dates, end_dates):
        start_param = list(zip(*[re.split('[-_:]', x)
                            for x in start_dates.split(',') if x]))
        end_param = list(zip(*[re.split('[-_:]', x)
                          for x in end_dates.split(',') if x]))
        time_param = start_param + end_param
        time_vars = ['start_year', 'start_month', 'start_day', 'start_hour',
                     'start_minute', 'start_second', 'end_year', 'end_month',
                     'end_day', 'end_hour', 'end_minute', 'end_second']
        for i, var in enumerate(time_vars):
            self.namelist['time_control'][
                var] = str_to_type(','.join(time_param[i]))

    def set_dfi(self, start_dates):
        
        def dfi_length_update(start_dates, fwd_integr=30, bck_integr=60):
            new_start_dates = ''
            end_start_dates = ''
            start_dates = start_dates.split(',')[0]
            
            for start_vars in [start_dates]:
                start_var = datetime.strptime(start_vars, '%Y-%m-%d_%H:%M:%S') + timedelta(seconds=60*fwd_integr)
                end_var = datetime.strptime(start_vars, '%Y-%m-%d_%H:%M:%S') - timedelta(seconds=60*bck_integr)
                new_start_dates = new_start_dates + start_var.strftime('%Y-%m-%d_%H:%M:%S') + ','
                end_start_dates = end_start_dates + end_var.strftime('%Y-%m-%d_%H:%M:%S') + ','
            new_start_dates = new_start_dates[:-1]
            end_start_dates = end_start_dates[:-1]
            
            return new_start_dates, end_start_dates
        
        start_dates, end_dates = dfi_length_update(start_dates)
        start_param = list(zip(*[re.split('[-_:]', x)
                            for x in start_dates.split(',') if x]))
        end_param = list(zip(*[re.split('[-_:]', x)
                          for x in end_dates.split(',') if x]))
        time_param = start_param + end_param
        time_vars = ['dfi_fwdstop_year', 'dfi_fwdstop_month', 'dfi_fwdstop_day', 'dfi_fwdstop_hour', 'dfi_fwdstop_minute', 'dfi_fwdstop_second',
                     'dfi_bckstop_year', 'dfi_bckstop_month', 'dfi_bckstop_day', 'dfi_bckstop_hour', 'dfi_bckstop_minute', 'dfi_bckstop_second']
                
        for i, var in enumerate(time_vars):
            self.namelist['dfi_control'][var] = str_to_type(','.join(time_param[i]))

    # def update_from_env(self):
    #     """Update parameters from environment variables"""
    #     self.update_var('time_control', 'interval_seconds', 'pregrid_interval')
    #     if os.getenv('wps_output'):
    #         self.namelist['time_control']['auxinput1_inname'] = \
    #             '{}<domain>.<date>'.format(os.getenv('wps_output'))

    #     domains_env_vars = ['num_metgrid_levels', 'num_metgrid_soil_levels',
    #                         'rh2qv_wrt_liquid']
    #     for env_var in domains_env_vars:
    #         self.update_var('domains', env_var, env_var)

    #     fdda_env_vars = ['grid_sfdda_1', 'grid_sfdda_2']
    #     for env_var in fdda_env_vars:
    #         self.update_var('fdda', env_var, env_var)
    #     self.update_var('fdda', 'grid_fdda', 'fdda_grid')
    #     self.update_var('fdda', 'gfdda_end_h', 'fdda_grid_end')
    #     self.update_var('fdda', 'obs_nudge_opt', 'fdda_obs')

    def update_from_file(self, filename):
        """Update parameters using an external file"""
        with open(filename) as fln:
            for line in fln:
                section, key, value = line.split(None, 2)
                if section not in self.namelist:
                    self.namelist[section] = OrderedDict()
                if ',' in value:
                    self.namelist[section][key] = []
                    for entry in value.split(','):
                        self.namelist[section][key].append(str_to_type(entry))
                else:
                    self.namelist[section][key] = value

    def extend_params(self):
        """For those parameters that need an entry per domain, create a list
        max_dom long and fill with the given value or repeat the first entry
        if none is given"""
        fields = deepcopy(self.max_dom_fields)
        for section, params in fields.items():
            for key, value in params.items():
                try:
                    if len(self.namelist[section][key]) == self.namelist['domains']['max_dom']:
                        continue
                except TypeError:
                    pass
                self.namelist[section][key] = [self.namelist[section][key]]
                for _ in range(1, self.namelist['domains']['max_dom']):
                    if value != None:
                        self.namelist[section][key].append(value)
                    else:
                        self.namelist[section][key].append(
                            self.namelist[section][key][-1])


class CommonNamelist(object):

    """Common namelist object"""

    def __init__(self, **kwargs):
        self.params = OrderedDict()
        for key in kwargs:
            self.params[key] = kwargs[key]
        self.update_common_params()

    def update_common_params(self):
        return self.params

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

def write_common(nudging_cutoff_time, go_back_hours, out_file):
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
    
        
if __name__ == '__main__':
    model = 'nz8kmN-ECMWF-var'
    go_back_hours = 0
    base_time = datetime.strptime('2017092100','%Y%m%d%H')
    nudging_cutoff_time = '2017092100'
    out_file = 'test.yaml' 
    work_dir = ''
    
    # write common
    write_common(nudging_cutoff_time, go_back_hours, out_file)

    # write gsi
    nml_gsi = write_gsi(model)
    with open(out_file, 'a') as outfile:
        yaml.dump({'gsi': nml_gsi.namelist}, outfile, default_flow_style=False)
    write_gsi_obsinput(gsiobs_config, filename)

    # write wrf
    nml_wrf = write_wrf(model, base_time,work_dir,go_back_hours)
    with open(out_file, 'a') as outfile:
        yaml.dump({'wrf': nml_wrf.namelist}, outfile, default_flow_style=False)

    WriteNamelist(out_file, 'namelist.input', 'wrf')
    WriteNamelist(out_file, 'namelist.gsi', 'gsi', add_config=filename)


    print 'done at {}'.format(out_file)
    