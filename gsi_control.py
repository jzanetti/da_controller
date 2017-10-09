import argparse
from datetime import datetime
from glob import glob
import inspect
import os
import shutil
import ntpath

GSIPROC = 3
run_command = 'mpirun -np {} '.format(GSIPROC)

BK_CORE = 'ARW'
BKCV_OPT = 'global'
BYTE_ORDER = 'Big_Endian'

EMIS_BIN_LIST = [
    'Nalli.IRwater.EmisCoeff.bin',
    'NPOESS.IRice.EmisCoeff.bin',
    'NPOESS.IRland.EmisCoeff.bin',
    'NPOESS.IRsnow.EmisCoeff.bin',
    'NPOESS.VISice.EmisCoeff.bin',
    'NPOESS.VISland.EmisCoeff.bin',
    'NPOESS.VISsnow.EmisCoeff.bin',
    'NPOESS.VISwater.EmisCoeff.bin',
    'FASTEM6.MWwater.EmisCoeff.bin',
    'AerosolCoeff.bin',
    'CloudCoeff.bin']

SATBIAS_LIST = {
    'coef':'gdas1.t00z.abias.20150617',
    'passive_coef': 'gdas1.t00z.abias_pc.20150617'}


def gsi_control(work_dir, gsi_base, crtm_base,
               background_path, 
               prepbufr = None, refbufr = None, velbufr = None):
    #args = setup_parser()
    
    da_startend(work_dir, gsi_base, start=True, end=False)
    
    link_gsi_bk_and_obs(background_path, work_dir,
                         prepbufr = prepbufr, 
                         refbufr = refbufr, 
                         velbufr = velbufr)
    
    link_gsi_err(os.path.join(gsi_base, 'fix'), BYTE_ORDER, work_dir)

    link_gsi_info(os.path.join(gsi_base, 'fix'), 
                  work_dir, bkcv_opt=BKCV_OPT,
                  byte_order=BYTE_ORDER)

    link_ctrm_emis_coeff(crtm_base, work_dir, BYTE_ORDER,
                          EMIS_BIN_LIST)
    sattype_list = retrieve_sattype(os.path.join(work_dir, 'satinfo'))
    link_sat_coef(crtm_base, work_dir, BYTE_ORDER, sattype_list)
    link_single_obs_table(os.path.join(gsi_base, 'fix'), work_dir)
    link_satbias(os.path.join(gsi_base, 'fix'), work_dir, SATBIAS_LIST)


def link_satbias(fix_root, work_dir, satbias_list):
    for coeff_type in satbias_list.keys():
        f = 'satbias_in' if coeff_type == 'coef' else 'satbias_pc'
        shutil.copyfile(os.path.join(fix_root, satbias_list[coeff_type]),
                    os.path.join(work_dir, f))


def link_single_obs_table(fix_root, work_dir):
    prepobs_table = os.path.join(fix_root, 'prepobs_prep.bufrtable')

    shutil.copyfile(prepobs_table, 
                    os.path.join(work_dir, 'prepobs_prep.bufrtable'))

def retrieve_sattype(satinfo_path):
    sattype_list = []
    with open(satinfo_path) as f:
        for line in f:
            sat = filter(None, line.split(' '))[0]
            if (not sat.startswith('!')) and ((sat in sattype_list) == False):
                sattype_list.append(sat)
    
    return sattype_list
                

def link_gsi_bk_and_obs(bkpath, work_root, prepbufr=None, refbufr=None, velbufr=None):

    obsname_map = {'prepbufr': 'prepbufr', 
                   'refbufr': 'refInGSI',
                   'velbufr': 'radarbufr',
                   'bkpath': 'wrf_inout'}
    
    if os.path.exists(bkpath) == False:
        return False

    for cbufr in inspect.getargspec(link_gsi_bk_and_obs).args:
        if eval(cbufr):
            if  os.path.exists(eval(cbufr)) == False:
                return False
            else:
                if cbufr in obsname_map.keys():
                    shutil.copyfile(eval(cbufr), 
                               os.path.join(work_root, obsname_map[cbufr]))
    return True

def link_gsi_info(fix_root, work_root, bkcv_opt = 'global',
                   msanavinfo = None, byte_order = 'Big_Endian'):
    anavinfo_table = {'global': 'anavinfo_ndas_netcdf_glbe', 'nam': 'nam_nmmstat_na.gcv', 
                      'ms': msanavinfo, 
                      'byte_order': byte_order, 'map': 'anavinfo'}
    info_table = {
        'satinfo': 'global_satinfo.txt',
        'convinfo': 'global_convinfo_user_defined.txt',
        'ozinfo': 'global_ozinfo.txt',
        'pcpinfo': 'global_pcpinfo.txt'}
                                   
    os.symlink(os.path.join(fix_root, anavinfo_table['byte_order'], anavinfo_table[bkcv_opt]),
                    os.path.join(work_root, anavinfo_table['map']))

    for info_type in info_table.keys():                                    
        os.symlink(os.path.join(fix_root, info_table[info_type]),
                    os.path.join(work_root, info_type))

def link_ctrm_emis_coeff(crtm_root, work_root, byte_order,
                          emis_bin_list):
    ctrm_root_order = os.path.join(crtm_root, byte_order)
    
    for coef_bin in emis_bin_list:                                    
        os.symlink(os.path.join(ctrm_root_order, coef_bin),
                    os.path.join(work_root, coef_bin))

def link_sat_coef(crtm_root, work_dir, byte_order, sattype_list):
    spc_tau_list = ['SpcCoeff.bin', 'TauCoeff.bin']
    binfile = []
    for spc_tau in spc_tau_list:
        for sat in sattype_list:
            binfile.extend(glob(os.path.join(crtm_root, byte_order, '{}*.{}'.format(sat, spc_tau))))
        
    for bfline in binfile:
        os.symlink(bfline, os.path.join(work_dir, ntpath.basename(bfline)))

def link_sgleobs_table(fix_root, work_root):
    os.symlink(os.path.join(fix_root, 'prepobs_prep.bufrtable'),
                os.path.join(work_root, 'prepobs_prep.bufrtable'))

def link_gsi_err(fix_root, byte_order, work_root, 
                  bk_core = 'global',
                  msberror = None, msoberr = None):

    err_table = {
        'berror':   {'global': 'nam_glb_berror.f77.gcv', 'nam': 'nam_nmmstat_na.gcv', 'ms': msberror, 
                     'byte_order': byte_order, 'map': 'berror_stats'},
        'oberror':  {'global': 'prepobs_errtable.global', 'nam': 'nam_errtable.r3dv', 'ms': msoberr,
                     'byte_order': '', 'map': 'errtable'}}

    for err_type in err_table.keys():
        os.symlink(os.path.join(fix_root, err_table[err_type]['byte_order'], err_table[err_type][bk_core]),
                        os.path.join(work_root, err_table[err_type]['map']))

def diag_process(work_root, analysis_obstime, 
                 diag_type_list = ['conv']):
    for diag_type in diag_type_list:
        if len(glob(os.path.join(work_root, '*', diag_type))) == 0:
            continue
        cmd01 = 'cat pe*.{}_01  > diag_{}_{}.{}'.format(diag_type, diag_type, 'ges', analysis_obstime) 
        cmd03 = 'cat pe*.{}_03  > diag_{}_{}.{}'.format(diag_type, diag_type, 'anl', analysis_obstime)
        os.system(cmd01)
        os.system(cmd03)

def da_startend(work_dir, gsi_base, start=True, end=False):
    if start:
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        os.symlink(os.path.join(gsi_base, 'run', 'gsi.exe'), os.path.join(work_dir, 'gsi.exe'))
