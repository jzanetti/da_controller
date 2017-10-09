import ast
from config_processor import namelist
import yaml

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
        setattr(gsi_io, obs, ast.literal_eval(gsiobs_config.get("GSIOBS", obs)))
    
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