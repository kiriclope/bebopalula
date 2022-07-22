import os
from write import *

path = os.getcwd()

replace_global('IF_INI_COND', 0)
replace_global('INI_COND_ID', 1)

replace_global('IF_TRIALS', 0) 
replace_global('TRIAL_ID', 1) 

replace_global('folder', "'noise_500_on'") 
exec(open( path + '/spatial.py').read()) 

replace_global('folder', "'noise_500_off'") 
exec(open( path + '/spatial.py').read()) 
    
replace_global('IF_TRIALS', 0) 
replace_global('TRIAL_ID', 1) 
