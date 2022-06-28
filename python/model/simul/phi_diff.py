import sys, importlib

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat

from get_m1 import * 
from utils import * 
from write import *

importlib.reload(sys.modules['params']) 
importlib.reload(sys.modules['get_m1']) 

gv.IF_INI_COND = 0
gv.IF_TRIALS = 0
gv.N_TRIALS = 10 

gv.init_param()

path = gv.path

def get_diffusion(path):
    phi_trial = []
    for i_trial in range(1, gv.N_TRIALS+1): 
        phi_ini = []
        for i_ini in range(1, 10+1):
            gv.path = path
            gv.path += '/trial_%d' % i_trial ;        
            gv.path += '/ini_cond_%d' % i_ini ; 
            print(gv.path)

            try:
                time, rates = get_time_rates(path=gv.path)
                _, phi = decode_bump(rates[:,0]) 
                phi_ini.append(phi)
                phi_ext = ( 1-i_trial / gv.N_TRIALS) * 180
                Dphi = phi[10] * 180 / np.pi - phi_ext 
                
                print('phi', phi[10] * 180 / np.pi,
                    'phi_ext', phi_ext,
                    'Dphi', Dphi)
            
            except:
                print('error')
                phi_ini.append(np.nan*np.zeros(40))
                pass
            
        phi_trial.append(phi_ini) 
    
    phi_trial = np.asarray(phi_trial) 
    print('phi', phi_trial.shape) 
    
    return phi_trial * 180 / np.pi 

phi_off = get_diffusion(path)
phi_off_avg = stat.circmean(phi_off[..., 12:16], high=180, nan_policy='omit', axis=-1) 
diff_off = np.hstack( phi_off_avg - stat.circmean(phi_off_avg, high=180, nan_policy='omit', axis=-1)[:,np.newaxis] )

diff_off[diff_off>90] -= 180 
diff_off[diff_off<-90] += 180 

print('diff_off', diff_off.shape) 

# figname = gv.folder + '_on_' + 'diff_hist'
figname = 'off_on_' + 'diff_hist'
plt.figure(figname)
plt.hist(2*diff_off, histtype='step', color='b') 
plt.xlabel('Diffusion (deg)') 
plt.ylabel('Count') 

path  = path.replace('off', 'on') 

phi_on = get_diffusion(path)
phi_on_avg = stat.circmean(phi_on[..., 12:16], high=180, nan_policy='omit', axis=-1) 
diff_on = np.hstack( phi_on_avg - stat.circmean(phi_on_avg, high=180, nan_policy='omit', axis=-1)[:,np.newaxis] ) 
print('diff_on', diff_on.shape)

diff_on[diff_on>90] -= 180 
diff_on[diff_on<-90] += 180 

plt.hist(2*diff_on, histtype='step', color='r') 

plt.savefig(figname + '.svg', dpi=300)
# plt.boxplot([diff_off, diff_on], patch_artist=True, labels=['off', 'on'])
# plt.ylabel('Diffusion (deg)')

