import sys, importlib

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat 
from get_m1 import * 
from utils import * 
from write import *

importlib.reload(sys.modules['params'])
importlib.reload(sys.modules['get_m1']) 

gv.folder = 'christos_off_2'
gv.IF_INI_COND = 0
gv.IF_TRIALS = 0

gv.N_TRIALS = 10 
gv.init_param()

path = gv.path

def get_diffusion(path):
    phi_trial = []
    for i_trial in range(1, 1 + 1):
    
        phi_ini = []
        for i_ini in range(1, 1 + 1):
            gv.path = path
            # gv.path += '/trial_%d' % i_trial ; 
            # gv.path += '/ini_cond_%d' % i_ini ; 
            print(gv.path)
            try:
                time, rates = get_time_rates(path=gv.path) 
                phi = get_phi(rates) 
                # phi_ini.append( phi[0] - (1.0 - i_trial/gv.N_TRIALS) * np.pi ) 
                Dphi = ( phi - phi[60] ) 
                # Dphi = ( phi - (i_trial/gv.N_TRIALS) * np.pi ) 
                phi_ini.append(Dphi) 
                print('phi', phi[10] * 180 / np.pi,
                      'phi_ext', (i_trial/gv.N_TRIALS)*180,
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

print('Dphi_off', Dphi_off.shape)

Sig_off = np.sqrt(np.nanmean(np.vstack(Dphi_off**2), axis=0)) 

# path = path.replace('christos_off', 'on_2') # change dirname 
path = path.replace('off', 'on') # change dirname 

Dphi_on = get_diffusion(path) 

Dphi_on[Dphi_on>90] -= 180 
Dphi_on[Dphi_on<-90] += 180 

print('Dphi_on', Dphi_on.shape)

Sig_on = np.sqrt(np.nanmean(np.vstack(Dphi_on**2), axis=0)) 

# figname = gv.folder + 'off_on_' + 'drift_hist'
figname = 'off_on_' + 'drift'

plt.figure(figname)
time = np.linspace(0, 14, Sig_off.shape[-1])
plt.plot(time, Sig_off, 'b')
plt.plot(time, Sig_on, 'g')

plt.xlabel('Time (s)') 
plt.ylabel('$\Sigma(t)$') 

plt.savefig(figname + '.svg', dpi=300)

# plt.boxplot([drift_off, drift_on], patch_artist=True, labels=['off', 'on'], showfliers=False, notch=True)
# plt.ylabel('Error (deg)')
