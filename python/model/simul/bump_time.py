import os
from write import *

path = os.getcwd()

replace_global('IF_INI_COND', 0) 
replace_global('INI_COND_ID', 1) 

replace_global('IF_TRIALS', 0) 
replace_global('TRIAL_ID', 4) 

# replace_global('PHI_CUE', 0.25)
# replace_global('PHI_CUE', 'TRIAL_ID/N_TRIALS') 

replace_global('folder', "'bump_on'") 
exec(open( path + '/m1_phi_time.py').read()) 

replace_global('folder', "'bump_off'") 
exec(open( path + '/m1_phi_time.py').read()) 

replace_global('IF_INI_COND', 0) 
replace_global('INI_COND_ID', 1) 

replace_global('IF_TRIALS', 0) 
replace_global('TRIAL_ID', 1) 

replace_global('PHI_CUE', 0.25) 
