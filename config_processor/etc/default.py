from collections import OrderedDict

parent_config = None

config = OrderedDict([
    ('wrf_namelist', OrderedDict([
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
            ('history_interval', 0),
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
            ('auxhist3_interval', 60),
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
        ('grib2', OrderedDict()),
        ('namelist_quilt', OrderedDict([
            ('nio_tasks_per_group', 2),
            ('nio_groups', 1)]))
        ])),
    ('var_control', OrderedDict([
           ('da_core', None),
           ('max_goback_hours', None),
           ('da_domain_id', None),
           ('da_win_half', None),
           ('cycle_freq', None),
           ('radar_da', None),
           ('use_nudge', None),
           ('dfi_radar', None)
           ])),
    ('fdda_control', OrderedDict([
        ('run_fdda', False),
        ('go_back_hours', 6),
        ])),
    ])




