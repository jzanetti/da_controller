parent_config = 'nz8kmN'

config = {
    'general_config': {
        'number_of_domains': 2,
        'domain_periods': [6, 6],
        'global_source': 'ECMWF_ANZ0125',
        'constant_global_sources': [],
        'constant_intermediate_data': ['lsm-ECMWF_ANZ0125']
        },
    'wps_namelist': {
        'geog_data_res': ['gtopo_10m+usgs_10m+nesdis_greenfrac+10m', 'gtopo_2m+usgs_2m+nesdis_greenfrac+2m'],
        },
    'wrf_namelist': {
        'domains': {
            'time_step': 60,
            },
        'physics':{
            'num_land_cat': 24,
            },
        'time_control': {
            'auxhist2_interval': [60, 60]
            },
        },
     'gsi_namelist': {
        'setup': {
            'miter': 3,
            'niter(1)': 50,
            'niter(2)': 50,
            'write_diag(1)': True,
            'write_diag(2)': False,
            'write_diag(3)': True,
            'gencode': 78,
            'qoption': 2,
            'factqmin': 0.0,
            'factqmax': 0.0,
            'iguess': -1,
            'oneobtest': False,
            'retrieval': False,
            'nhr_assimilation': 3,
            'l_foto': False,
            'use_pbl': False,
            'lwrite_peakwt': True,
            'lread_obs_save': False,
            'lread_obs_skip': False,
            'newpc4pred': True,
            'adp_anglebc': True,
            'angord': 4,
            'passive_bc': True,
            'use_edges': False,
            'emiss_bc': True,
            'diag_precon': True,
            'step_start': 1.e-3},
           },
     'archive_control':{
           'archive_type': ['wrfout'],
           'archive_domain': [1,2],
           'archive_prog': '1,3',
        },
     'var_control': {
           'da_core': 'gsi',
           'cycle_freq': 3,
           'max_goback_hours': 24,
           'da_domain_id': [2],
           'da_win_half': [0.5],
           'use_nudge': True,
           'dfi_radar': True,
           }
    }

