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
    for i_trial in range(1, gv.N_TRIALS + 1):
    
        phi_ini = []
        for i_ini in range(1, 1 + 1):
            gv.path = path
            gv.path += '/trial_%d' % i_trial ; 
            gv.path += '/ini_cond_%d' % i_ini ; 
            print(gv.path)
            try:
                time, rates = get_time_rates(path=gv.path) 
                _, phi = decode_bump(rates[:,0]) 
                # phi_ini.append( phi[0] - (1.0 - i_trial/gv.N_TRIALS) * np.pi )                
                # Dphi = ( phi - (1-i_trial/gv.N_TRIALS) * np.pi ) # remember 1st
                Dphi = ( phi - (i_trial/gv.N_TRIALS) * np.pi ) # remember 2nd
                
                phi_ini.append(Dphi) 
                print('phi', phi[10] * 180 / np.pi,
                      'phi_ext', (1-i_trial/gv.N_TRIALS)*180,
                      'Dphi', Dphi[10] * 180 / np.pi) 
            except:
                phi_ini.append(np.nan*np.zeros(40)) 
                print('error') 
                pass
            
        phi_trial.append(phi_ini)
    
    phi_trial = np.asarray(phi_trial) 
    # print('phi', phi_trial.shape) 
    
    return phi_trial * 180 / np.pi 

Dphi_off = get_diffusion(path)

Dphi_off[Dphi_off>90] -= 180 
Dphi_off[Dphi_off<-90] += 180 

drift_off = stat.circmean(Dphi_off[..., 24:28], high=90, low=-90, axis=-1, nan_policy='omit') # over time
drift_off_avg = np.nanmean(drift_off, axis=-1) # over realisations

# figname = gv.folder + 'off_on_' + 'drift_hist'
figname = 'off_on_' + 'drift_hist'
plt.figure(figname)

trials = np.arange(1, 11)
distance = np.sqrt( np.abs( ( (1.0-trials/10 )*180 )**2 - (trials/10 * 180)**2 ) ) 

print(distance.shape, drift_off_avg.shape) 
plt.plot(distance, 2*drift_off_avg, color='b')

plt.ylabel('Angular Deviation (deg)') 
plt.xlabel('Distance')

path = path.replace('off', 'on') # change dirname 

Dphi_on = get_diffusion(path) 

Dphi_on[Dphi_on>90] -= 180 
Dphi_on[Dphi_on<-90] += 180 

drift_on = stat.circmean(Dphi_on[..., 24:28], high=90, low=-90, axis=-1, nan_policy='omit') # over time
drift_on_avg = np.nanmean(drift_on, axis=-1) # over realisations

plt.plot(distance, 2*drift_on_avg, color='r')

plt.savefig(figname + '.svg', dpi=300)

# plt.boxplot([drift_off, drift_on], patch_artist=True, labels=['off', 'on'], showfliers=False, notch=True)
# plt.ylabel('Error (deg)')
