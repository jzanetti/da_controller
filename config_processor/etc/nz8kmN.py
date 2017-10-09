from config_processor.basemap import LambertConformalBasemapParameters as LCBP

parent_config = 'default'

config = {
    'basemap_parameters': {
        'd1': LCBP(lat_0=-41.00000763,
                   lat_1=-30.,
                   lat_2=-60.,
                   llcrnrlat=-55.4167823792,
                   llcrnrlon=135.009307861,
                   lon_0=167.5,
                   urcrnrlat=-21.3899841309,
                   urcrnrlon=-174.317016602),
        'd2': LCBP(lat_0=-40.74848557,
                   lat_1=-30.,
                   lat_2=-60.,
                   llcrnrlat=-48.7435379028,
                   llcrnrlon=164.509033203,
                   lon_0=167.5,
                   urcrnrlat=-32.3397026062,
                   urcrnrlon=179.364990234)
        },
    'wps_namelist': {
        'parent_id': [1, 1],
        'parent_grid_ratio': [1, 3],
        'i_parent_start': [0, 74],
        'j_parent_start': [0, 48],
        'e_we': [165, 166],
        'e_sn': [165, 214],
        'dx': [24000],
        'dy': [24000],
        'geog_data_res': ['10m', '2m'],
        },
    'wrf_namelist': {
        'time_control': {
            'auxhist2_interval': [0, 60]
            },
        'domains': {
            'e_vert': [51, 51],
            'time_step': 144,
            'time_step_fract_num': 0,
            'time_step_fract_den': 1
            },
        'fdda': {
            'gfdda_interval_m': 180,
            },
        'namelist_quilt': {
            'nio_tasks_per_group': 2
            }
        }
    }
