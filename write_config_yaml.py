from wrf_workflow import configwriter_activity
from datetime import datetime
import yaml
import os

model = 'nz8kmN-ECMWF-var'
base_time = '2017092100'
work_dir = '/home/jzanetti/github/da_controller/test_conf'

#####################################
# config_file
#####################################
nudging_cutoff_time = base_time
config_filepath = model + '.py'
base_time = datetime.strptime(base_time,'%Y%m%d%H')
out_file = os.path.join(work_dir, model + '_test.yaml')
go_back_hours = 0

if os.path.exists(work_dir) == False:
    os.makedirs(work_dir)

#####################################
# write yaml file
#####################################
# write common
configwriter_activity.write_common(model, base_time, nudging_cutoff_time, go_back_hours, out_file)

# write gsi
nml_gsi = configwriter_activity.write_gsi(model)
with open(out_file, 'a') as outfile:
    yaml.dump({'gsi': nml_gsi.namelist}, outfile, default_flow_style=False)
#write_gsi_obsinput(gsiobs_config, filename)

# write wrf
nml_wrf = configwriter_activity.write_wrf(model, base_time,work_dir,go_back_hours)
with open(out_file, 'a') as outfile:
    yaml.dump({'wrf': nml_wrf.namelist}, outfile, default_flow_style=False)

#####################################
# write namelist.input
#####################################
configwriter_activity.WriteNamelist(out_file, os.path.join(work_dir, 'namelist.input'), 'wrf')

print 'write {}/namelist.input'.format(work_dir)
print 'write {}'.format(out_file)
