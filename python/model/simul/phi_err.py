import sys, importlib

import numpy as np
import matplotlib.pyplot as plt

import get_m1
importlib.reload(sys.modules['get_m1']) 

from get_m1 import * 
from utils import * 
from write import *

importlib.reload(sys.modules['params'])


gv.folder = 'albert_off'
gv.IF_INI_COND = 0
gv.IF_TRIALS = 0
gv.init_param()

path = gv.path

def get_diffusion(path):
    phi_trial = []
    for i_trial in range(1,26):
    
        phi_ini = []
        for i_ini in range(1,11):
            gv.path = path
            gv.path += '/trial_%d' % i_trial ; 
            gv.path += '/ini_cond_%d' % i_ini ; 
            print(gv.path)
        
            time, rates = get_time_rates(path=gv.path) 
            phi = get_phi(rates) 
            # phi_ini.append( phi[0] - (1.0 - i_trial/10) * np.pi ) 
            phi_trial.append( phi[0] - (1.0 - i_trial/25) * np.pi ) 
        
        # phi_trial.append(phi_ini)
        
    phi_trial = np.asarray(phi_trial) 
    # print('phi', phi_trial.shape) 
    
    phi_trial = np.mean(phi_trial[..., -8:], axis=-1) # average over time 
    # print('phi', phi_trial.shape) 
    
    return phi_trial * 180 / np.pi 

# phi_stim = np.arange(1,11)/10 * 180 

phi_off = get_diffusion(path)
# diff_off = np.mean(phi_off, axis=-1) # average over initial conditions 
diff_off = phi_off # average over initial conditions 
diff_off[diff_off>90] -= 180 
diff_off[diff_off<-90] += 180 

print('diff_off', diff_off.shape) 
print(diff_off)

mean_off = np.mean( diff_off, axis=-1) # average over trials 
std_off = np.std( diff_off, axis=-1)

print('mean_off', mean_off, 'std_off', std_off) 

plt.hist(diff_off, histtype='step', color='b') 
plt.xlabel('Error (deg)') 

# path = path.replace('off', 'off_2') # change dirname 
path = path.replace('off', 'on') # change dirname 

phi_on = get_diffusion(path) 
# diff_on = np.mean(phi_on, axis=-1) # average over initial conditions 
diff_on = phi_on # average over initial conditions 
diff_on[diff_on>90] -= 180 
diff_on[diff_on<-90] += 180 
print('diff_on', diff_on.shape) 
print(diff_on) 

mean_on = np.mean( diff_on, axis=-1) # average over trials 
std_on = np.std( diff_on, axis=-1) 

print('mean_on', mean_on, 'std_on', std_on) 

plt.hist(diff_on, histtype='step', color='r') 
plt.xlabel('Error (deg)') 

# plt.boxplot([diff_off, diff_on], patch_artist=True, labels=['off', 'on'], showfliers=False, notch=True)
# plt.ylabel('Error (deg)')
