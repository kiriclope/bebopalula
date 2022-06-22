import sys, os, importlib
from importlib import reload
from scipy.signal import savgol_filter 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import params as gv 
importlib.reload(sys.modules['params'])
import get_m1
importlib.reload(sys.modules['get_m1']) 

from get_m1 import * 
from utils import * 
from rank_utils import * 

alpha = [1,.05]
# if gv.IF_SPEC == 1:
#     if len(sys.argv)>0:
#         gv.KAPPA = float(sys.argv[1])
#         gv.KAPPA_1 = float(sys.argv[1])
    
# if gv.IF_INI_COND==1: 
#     if len(sys.argv)>0:
#         if gv.IF_SPEC:
#             gv.INI_COND_ID = int(sys.argv[2])
#         else:
#             gv.INI_COND_ID = int(sys.argv[1])

#     gv.MAP_SEED=gv.INI_COND_ID

# if gv.IF_TRIALS==1:    
#     if len(sys.argv)>0:
#         if gv.IF_SPEC:
#             gv.TRIAL_ID = int(sys.argv[2])
#         else:
#             gv.TRIAL_ID = int(sys.argv[1])

gv.init_param() 
        
time, rates = get_time_rates(path=gv.path)
m0 = np.nanmean(rates, axis=-1).T
avg_m0 = np.nanmean(m0, axis=-1) 

m1, phi, smooth_rates = get_m1_phi_smooth_rates(rates, osc=0) 

print('m0', avg_m0, 'm1', np.mean(m1, axis=-1), 'smooth rates', smooth_rates.shape ) 

if(gv.RANK==2):
    time, rates_perm = get_time_rates(MAP=1, path=gv.path, con_path=gv.con_path) 
    print(rates_perm.shape) 
    m1_perm, phi_perm, smooth_rates_perm = get_m1_phi_smooth_rates(rates_perm) 
    print('smooth_rates_perm', smooth_rates_perm.shape) 
    
# figtitle = 'm1_phi_time' + gv.folder 
# fig = plt.figure(figtitle, figsize=(1.618*2*1.5*2, 1.618*2*gv.RANK)) 

figtitle = 'm1_phi_time' 
fig = plt.figure(figtitle, figsize=(1.618*2*1.5*2, 1.618*2*gv.RANK)) 

# ax = plt.subplot(gv.RANK,2,1) 
    
# for i_pop in range(gv.n_pop): 
#     theta = np.linspace(0, np.pi, gv.n_size[i_pop]) 
#     ax.plot(theta* 180/np.pi, np.nanmean(smooth_rates[i_pop, -4:, :gv.n_size[i_pop]], axis=0), '-', color=gv.pal[i_pop], alpha=alpha[i_pop] ) 

# theta = np.linspace(0, np.pi, gv.n_size[0]) 
# plt.plot(theta * 180/np.pi, np.mean(m0[0][-4:],axis=-1) + np.mean(m1[0][-4:], axis=-1) * np.cos(2*(theta + np.mean(phi[0][-4:], axis=-1) ) ), 'r--') 

# if(gv.IF_TRIALS): 
#     plt.vlines( (1-gv.TRIAL_ID) * 180 / 25.0 % 180, 0, 20, ls='--', lw=2) 
# else: 
#     plt.vlines( (1-gv.PHI_EXT) * 180, 0, 20, ls='--', lw=2) 

# print(np.mean(m0[0][-4:], axis=-1), np.mean(m1[0][-4:], axis=-1), np.mean(phi[0][-4:], axis=-1)) 

# # ax.set_xlabel('$\\theta_0$ (rad)')
# ax.set_xlabel('Prefered Orientation')
# ax.set_ylabel('Rates (Hz)') 
# plt.xlim([0, 180])
# plt.xticks([0, 45, 90, 135, 180])
# # plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
# #            ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])

ax = plt.subplot(gv.RANK,2,1)

for i_pop in range(gv.n_pop) : 
    plt.plot(time, m1[i_pop], '-', color=gv.pal[i_pop], alpha=alpha[i_pop]) 

plt.xlabel('Time (s)') 
plt.ylabel('Bump Amplitude (Hz)')

# plt.ylabel('$\\nu^{(1)}_0 / \\nu^{(0)}_0}$ (Hz)') 
plt.xlim([0, 10]) 

ax = plt.subplot(gv.RANK,2,2) 

for i_pop in range(gv.n_pop) :
    # plt.plot(time, phi[i_pop], color=gv.pal[i_pop], alpha=alpha[i_pop]) 
    # time2, phi2 = phitoPi(time, phi[i_pop]) 
    plt.plot(time, phi[i_pop]*180/np.pi, color=gv.pal[i_pop], alpha=alpha[i_pop])

    if(gv.IF_TRIALS):
        plt.hlines( gv.TRIAL_ID * 180 / gv.N_TRIALS % 180, 0, 10, ls='--') 
        print('phi', np.mean(phi[i_pop, -4:]) * 180 / np.pi, 'phi_trial', gv.TRIAL_ID / gv.N_TRIALS * 180) 
    else:
        plt.hlines( (gv.PHI_EXT) * 180, 0, 10, ls='--') 
        print('phi', np.mean(phi[i_pop, -4:]) * 180 / np.pi , 'phi_ext', (gv.PHI_EXT) * 180 % 180) 
       
plt.xlabel('Time (s)')
plt.ylabel('Bump Phase (°)')
plt.ylim([0, 180])
plt.yticks([0, 45, 90, 135, 180])
# plt.yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
#            ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
# plt.xlim([0, 10])

if(gv.RANK==2):

    ax = plt.subplot(2,3,4) 
    
    for i_pop in range(gv.n_pop): 
        theta = np.linspace(0, np.pi, gv.n_size[i_pop]) 
        ax.plot(theta, np.nanmean(smooth_rates_perm[i_pop, :, :gv.n_size[i_pop]], axis=0), '-', color=gv.pal[i_pop], alpha=alpha[i_pop]) 
    
    ax.set_xlabel('$\\theta_1$ (rad)') 
    ax.set_ylabel('$\\nu(\\theta_1)$ (Hz)') 
    plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
               ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$']) 
    
    ax = plt.subplot(2,3,5) 
    
    for i_pop in range(gv.n_pop) : 
        plt.plot(time, m1_perm[i_pop]/m0[i_pop], '-', color=gv.pal[i_pop]) 
    
    plt.xlabel('Time (s)') 
    plt.ylabel('$\\nu^{(1)}_1 / \\nu^{(0)}_1$ (Hz)') 
    # plt.xlim([0, 10])
    ax = plt.subplot(2,3,6) 

    for i_pop in range(gv.n_pop) : 
        plt.plot(time, phi_perm[i_pop], color=gv.pal[i_pop], alpha=alpha[i_pop]) 
    
    plt.xlabel('Time (s)') 
    plt.ylabel('$\phi_1$ (rad)') 
    plt.yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
               ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
    # plt.xlim([0, 10])
        
# plt.show()
