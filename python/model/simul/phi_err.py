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

gv.N_TRIALS = 20
gv.init_param()

path = gv.path

def get_diffusion(path):
    Dphi_trial = []
    for i_trial in range(1, gv.N_TRIALS + 1):
        gv.path = path
        gv.path += '/trial_%d' % i_trial ; 
        print(gv.path)
        try:
            time, rates = get_time_rates(path=gv.path) 
            _, phi = decode_bump(rates[:,0]) 
            # phi_ini.append( phi[0] - (1.0 - i_trial/gv.N_TRIALS) * np.pi )                
            
            print('phi', phi.shape)
            Dphi = ( phi - (1-gv.PHI_DIST) * np.pi )
            
            phi_trial.append(phi) 
            print('phi', phi[28] * 180 / np.pi,
                  'phi_dist', 180-gv.PHI_DIST * 180,
                  'Dphi', Dphi[28] * 180 / np.pi) 
        except:
            phi_trial.append(np.nan*np.zeros(40)) 
            print('error') 
            pass
                
    phi_trial = np.asarray(phi_trial) 
    print('phi_trial', phi_trial.shape) 
    
    return phi_trial * 180 / np.pi 

Dphi_off = get_diffusion(path)

Dphi_off[Dphi_off>90] -= 180 
Dphi_off[Dphi_off<-90] += 180 

# Dphi_off[np.abs(Dphi_off)>10] = np.nan

drift_off_avg = stat.circmean(Dphi_off[..., 28:32], high=90, low=-90, axis=-1, nan_policy='omit') 

# figname = gv.folder + 'off_on_' + 'drift_hist'
figname = 'off_on_' + 'drift_hist'

plt.figure(figname)
plt.hist(2*drift_off_avg, histtype='step', color='b') 
plt.xlabel('Angular Deviation (deg)') 
plt.ylabel('Count')

path = path.replace('off', 'on') # change dirname 

Dphi_on = get_diffusion(path) 

Dphi_on[Dphi_on>90] -= 180 
Dphi_on[Dphi_on<-90] += 180 

# Dphi_on[np.abs(Dphi_on)>10] = np.nan 

drift_on_avg = stat.circmean(Dphi_on[..., 28:32], high=90, low=-90, axis=-1, nan_policy='omit') 
plt.hist(2*drift_on_avg, histtype='step', color='r') 

plt.savefig(figname + '.svg', dpi=300)

# plt.boxplot([drift_off, drift_on], patch_artist=True, labels=['off', 'on'], showfliers=False, notch=True)
# plt.ylabel('Error (deg)')
