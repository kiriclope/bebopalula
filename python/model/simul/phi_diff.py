import sys, importlib

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
from joblib import Parallel, delayed

import progressbar as pgb 
from get_m1 import * 
from utils import * 
from write import *

importlib.reload(sys.modules['params']) 
importlib.reload(sys.modules['get_m1']) 

gv.IF_INI_COND = 0 
gv.IF_TRIALS = 0 
gv.N_TRIALS = 10 
gv.N_INI = 50
gv.IF_CHRISTOS = 0 
gv.folder = 'noise_200_off'
gv.init_param()
              
path = gv.path

def parloop(path, i_trial, i_ini, N_TRIALS=gv.N_TRIALS, A_CUE=gv.A_CUE, EPS_CUE=gv.EPS_CUE, A_DIST=gv.A_DIST, EPS_DIST=gv.EPS_DIST):

    phi_cue = i_trial / N_TRIALS 
    
    path += '/christos'
    path += '/cue_A_%.2f_eps_%.2f_phi_%.3f' % (A_CUE, EPS_CUE, phi_cue)
    path += '/dist_A_%.2f_eps_%.2f_phi_%.3f' % (A_DIST, EPS_DIST, 1-phi_cue)
    
    path += '/trial_%d' % i_trial ; 
    path += '/ini_cond_%d' % i_ini ; 
    
    try:
        _, rates = get_time_rates(path=path) 
        m1, phi = decode_bump(rates[:,0])
        
        if(phi.shape[0]!=160):
            m1 = np.nan*np.zeros(160) 
            phi = np.nan*np.zeros(160) 
    except:
        m1 = np.nan*np.zeros(160)
        phi = np.nan*np.zeros(160)
        print(path)
        pass

    return m1, phi 

with pgb.tqdm_joblib( pgb.tqdm(desc='phi off', total= gv.N_INI*gv.N_TRIALS) ) as progress_bar: 
    
    m1_off, phi_off = zip( *Parallel(n_jobs=-1)(delayed(parloop)(path, i_trial, i_ini) 
        for i_trial in range(1, gv.N_TRIALS+1) 
        for i_ini in range(1, gv.N_INI+1) ) )

m1_off = np.asarray(m1_off).reshape(gv.N_INI, gv.N_TRIALS, 160) 
phi_off = np.asarray(phi_off).reshape(gv.N_INI, gv.N_TRIALS, 160) * 180 / np.pi 

print(phi_off.shape) 
# print(phi_off[0, : , int(3/gv.T_WINDOW)])

path  = path.replace(gv.folder, 'bump_off') 

with pgb.tqdm_joblib( pgb.tqdm(desc='phi on', total= gv.N_INI *gv.N_TRIALS) ) as progress_bar: 
    
    m1_on, phi_on = zip( *Parallel(n_jobs=-1)(delayed(parloop)(path, i_trial, i_ini) 
        for i_trial in range(1, gv.N_TRIALS+1) 
        for i_ini in range(1, gv.N_INI+1) ) )

m1_on = np.asarray(m1_on).reshape(gv.N_INI, gv.N_TRIALS, 160) 
phi_on = np.asarray(phi_on).reshape(gv.N_INI, gv.N_TRIALS, 160) * 180 / np.pi 

print(phi_on.shape) 

bins =  [int(4/gv.T_WINDOW), int(5/gv.T_WINDOW)] 

mean_phi_off = stat.circmean(phi_off, high=180, low=0, nan_policy='omit', axis=1) # average over initial conditions
diff_off = phi_off - mean_phi_off[:,np.newaxis]
print('diff_off', diff_off.shape) 

diff_off[diff_off>90] -= 180 
diff_off[diff_off<-90] += 180 

diff_off_avg = np.nanmean(diff_off[..., bins[0]:bins[1]], axis=-1) # averag
diff_off_stack = np.hstack( diff_off_avg ) # stack trials and ini together 

mean_phi_on = stat.circmean(phi_on, high=180, low=0, nan_policy='omit', axis=1) # average over initial conditions
diff_on = phi_on - mean_phi_on[:,np.newaxis]
print('diff_on', diff_on.shape) 

diff_on[diff_on>90] -= 180 
diff_on[diff_on<-90] += 180 

diff_on_avg = np.nanmean(diff_on[..., bins[0]:bins[1]], axis=-1) # average
diff_on_stack = np.hstack( diff_on_avg ) # stack trials and ini together

figname = gv.folder + '_on_' + 'diff_hist'
plt.figure(figname)
plt.hist(2*diff_off_stack, histtype='step', color='b') 
plt.hist(2*diff_on_stack, histtype='step', color='r') 
plt.xlabel('Diffusion (°)') 
plt.ylabel('Count') 

plt.savefig(figname + '.svg', dpi=300)

figname = gv.folder + '_on_' + 'circular'
plt.figure(figname, figsize=(3.5, 3.5))

m1_off_avg = np.nanmean(m1_off[..., bins[0]:bins[1]], axis=-1)
phi_off_avg = np.nanmean(phi_off[..., bins[0]:bins[1]], axis=-1)

x_off = np.cos(2*phi_off_avg / 180 * np.pi)
y_off= np.sin(2*phi_off_avg/ 180 * np.pi)

plt.plot(m1_off_avg * x_off, m1_off_avg * y_off, 'bx', ms=5)

m1_on_avg = np.nanmean(m1_on[..., bins[0]:bins[1]], axis=-1)
phi_on_avg = np.nanmean(phi_on[..., bins[0]:bins[1]], axis=-1)

x_on = np.cos(2*phi_on_avg / 180 * np.pi)
y_on= np.sin(2*phi_on_avg/ 180 * np.pi)

plt.plot(m1_on_avg * x_on, m1_on_avg * y_on, 'rx', ms=5)
plt.axis('off')

phi_cues = np.linspace(0.1, 1, gv.N_TRIALS) * 180

x_cues = np.cos(2*phi_cues / 180 * np.pi)
y_cues= np.sin(2*phi_cues / 180 * np.pi)

m1_cues = np.nanmean( m1_off_avg + m1_on_avg ) / 2

plt.plot(m1_cues * x_cues, m1_cues * y_cues, 'k+', ms=20)

plt.savefig(figname + '.svg', dpi=300)

