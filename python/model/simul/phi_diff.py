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
gv.IF_CHRISTOS=0 

gv.init_param()
              
path = gv.path

def get_diffusion(path):
    phi_trial = []
    for i_trial in range(1, 10+1): 
        phi_ini = []
        for i_ini in range(1, 20+1): 
            gv.path = path
            
            gv.path += '/christos'
            phi_cue = i_trial / gv.N_TRIALS
            
            gv.path += '/cue_A_%.2f_eps_%.2f_phi_%.3f' % (gv.A_CUE, gv.EPS_CUE, phi_cue)
            gv.path += '/dist_A_%.2f_eps_%.2f_phi_%.3f' % (gv.A_DIST, gv.EPS_DIST, 1-phi_cue)
            
            gv.path += '/trial_%d' % i_trial ; 
            gv.path += '/ini_cond_%d' % i_ini ; 
            print(gv.path)

            try:
                time, rates = get_time_rates(path=gv.path)
                _, phi = decode_bump(rates[:,0]) 

                phi_ini.append(phi)
                
                # print('phi', phi.shape)
                if(phi.shape[0]!=160):
                    print('simul not done')
                    phi = np.nan*np.zeros(160) 
                
                phi_cue = ( 1.0 - i_trial / gv.N_TRIALS) * 180 
                Dphi = phi[int(3/gv.T_WINDOW)] * 180 / np.pi - phi_cue 
                
                print('phi', phi[int(3/gv.T_WINDOW)] * 180 / np.pi,
                    'phi_cue', phi_cue,
                    'Dphi', Dphi ) 
                
            except:
                print('error')
                phi_ini.append(np.nan*np.zeros(160))
                pass
            
        phi_trial.append(phi_ini) 
    
    phi_trial = np.asarray(phi_trial) 
    print('phi', phi_trial.shape) 
    
    return phi_trial * 180 / np.pi 

bins =  [int(4/gv.T_WINDOW), int(5/gv.T_WINDOW)] 

phi_off = get_diffusion(path) # trials x ini x time
print('phi_off', phi_off.shape)

mean_phi_off = stat.circmean(phi_off, high=180, low=0, nan_policy='omit', axis=1) # average over initial conditions
diff_off = phi_off - mean_phi_off[:,np.newaxis]
print('diff_off', diff_off.shape) 

diff_off[diff_off>90] -= 180 
diff_off[diff_off<-90] += 180 

diff_off_avg = np.nanmean(diff_off[..., bins[0]:bins[1]], axis=-1) # averag
# diff_off_avg = stat.circmean(diff_off[..., bins[0]:bins[1]], high=180, low=0, nan_policy='omit', axis=-1) # average over time 
diff_off_stack = np.hstack( diff_off_avg ) # stack trials and ini together 

path  = path.replace('off', 'on') 

phi_on = get_diffusion(path)  # trials x ini x time 
print('phi_on', phi_off.shape)

mean_phi_on = stat.circmean(phi_on, high=180, low=0, nan_policy='omit', axis=1) # average over initial conditions
diff_on = phi_on - mean_phi_on[:,np.newaxis]
print('diff_on', diff_on.shape) 

diff_on[diff_on>90] -= 180 
diff_on[diff_on<-90] += 180 

diff_on_avg = np.nanmean(diff_on[..., bins[0]:bins[1]], axis=-1) # average
# diff_on_avg = stat.circmean(diff_on[..., bins[0]:bins[1]], high=180, low=0, nan_policy='omit', axis=-1) # average over time 
diff_on_stack = np.hstack( diff_on_avg ) # stack trials and ini together

figname = gv.folder + '_on_' + 'diff_hist'
plt.figure(figname)
plt.hist(2*diff_off_stack, histtype='step', color='b') 
plt.hist(2*diff_on_stack, histtype='step', color='r') 
plt.xlabel('Diffusion (°)') 
plt.ylabel('Count') 

plt.savefig(figname + '.svg', dpi=300)

figname = gv.folder + '_on_' + 'circular'
plt.figure(figname)
phi_off_avg = np.nanmean(phi_off[..., bins[0]:bins[1]], axis=-1)
x_off = np.cos(2*phi_off_avg / 180 * np.pi)
y_off= np.sin(2*phi_off_avg/ 180 * np.pi)

plt.plot(x_off, y_off, 'bx', ms=4)

phi_on_avg = np.nanmean(phi_on[..., bins[0]:bins[1]], axis=-1)
x_on = np.cos(2*phi_on_avg / 180 * np.pi)
y_on= np.sin(2*phi_on_avg/ 180 * np.pi)

plt.plot(x_on, y_on, 'rx', ms=4)

plt.savefig(figname + '.svg', dpi=300)

