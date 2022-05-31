import sys, os, importlib
from importlib import reload
from scipy.signal import savgol_filter 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import params as gv 
importlib.reload(sys.modules['params']) 

from get_m1 import *
from utils import *

gv.init_param()
    
filter_rates = pd.read_csv(gv.path + '/filter_rates.dat', sep='\s+').to_numpy()
# filter_rates = filter_rates[20:]
time = filter_rates[:,0] / 1000
rates = np.delete(filter_rates, [0], axis=1) 

if gv.n_pop!=1:
    rates = pad_along_axis(rates, gv.n_size[0]-gv.n_size[1]) 
    n_neurons = int(rates.shape[1]/2) 
    rates = np.reshape(rates, (rates.shape[0], 2, n_neurons))  
else: 
    n_neurons = rates.shape[0] 
    
mean_rates = np.nanmean(rates, axis=-1) 
avg_mean_rates = np.nanmean(mean_rates[20:,:], axis=0) 

figtitle = 'spatial_profile_' + gv.folder 
figtitle = 'spatial_profile_1' 
fig = plt.figure(figtitle, figsize=(1.25*1.618*1.5*2, 1.618*1.25*2))

ax = fig.add_subplot(2,2,1) 
for i_pop in range(gv.n_pop): 
    plt.plot(time, mean_rates[:,i_pop], lw=1, color=gv.pal[i_pop]) 
    
plt.title('$m_E^0=$%.2f, $m_I^0=$%.2f' % (avg_mean_rates[0], avg_mean_rates[1]) ) 
plt.xlabel('Time (s)') 
plt.ylabel('Rates (Hz)') 
add_vlines()

ax = fig.add_subplot(2,2,2)

avg_m1=[] 
for i_pop in range(gv.n_pop) : 
    
    theta = np.linspace(0, np.pi, gv.n_size[i_pop]) 
    pop_rates = rates[:, i_pop, : gv.n_size[i_pop]]
    
    smooth_rates = circular_convolution(pop_rates, int(pop_rates.shape[-1]*.1) ) # over neurons 
    
    m1 = compute_m1(smooth_rates) 
    print('m1', m1.shape) 
    # m1 = circular_convolution(m1, int( m1.shape[0]*.01 ), axis=0 ) # over time 
    
    m0 = np.nanmean(pop_rates, axis=1) 
    plt.plot(time, m1, '-', color=gv.pal[i_pop]) 

    avg_m1.append(np.nanmean(m1)) 
    
plt.title('$m^1_E=$%.2f, $m^1_I=$%.2f' % (avg_m1[0], avg_m1[1]) ) 
plt.xlabel('Time (s)') 
plt.ylabel('$m_1$ (Hz)') 
add_vlines()

ax = fig.add_subplot(2,2,3) 

avg_phi = [] 
for i_pop in range(gv.n_pop) : 
    
    theta = np.linspace(0, np.pi, gv.n_size[i_pop]) 
    pop_rates = rates[:, i_pop, : gv.n_size[i_pop]]
        
    smooth_rates = circular_convolution(pop_rates, int(pop_rates.shape[-1]*.01) ) # over neurons
    
    phi = compute_phi(smooth_rates) 
    # phi = circular_convolution(phi, int( phi.shape[0]*.1 ), axis=0 ) 
    
    if(i_pop==0): 
        plt.plot(time, phi, color=gv.pal[i_pop]) 
    
    avg_phi.append(np.nanmean(phi))
    
plt.title('$\phi_E=$%.2f, $\phi_I=$%.2f' % (avg_phi[0] *180/np.pi, avg_phi[1]*180/np.pi) ) 

plt.xlabel('Time (s)') 
plt.ylabel('$\phi$ (rad)') 
plt.yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
           ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{2}$', r'$\pi$'])

add_vlines()    
    
ax = fig.add_subplot(2,2,4) 

avg_rates = np.nanmean(rates, axis=0) # over time 
for i_pop in range(gv.n_pop):

    pop_rates = avg_rates[i_pop]
    pop_rates = pop_rates[~np.isnan(pop_rates)]
    
    smooth_avg_rates = circular_convolution(pop_rates, int(pop_rates.shape[0]*.1 ) ) # over neurons 
    smooth_avg_rates = np.flip(smooth_avg_rates, axis=-1) 
    
    print(smooth_avg_rates.shape) 
    avg_m1 = compute_m1(smooth_avg_rates) 
    avg_phi = compute_phi(smooth_avg_rates) 
    
    print('[<m1>]', avg_m1, '[<phi>]', avg_phi) 
    
    theta = np.linspace(0, np.pi, gv.n_size[i_pop]) 
    cos_func =  avg_mean_rates[i_pop] + avg_m1 * np.cos(2*theta - avg_phi ) 
    
    # print(theta.shape, cos_func.shape) 
    
    plt.plot(theta, smooth_avg_rates, color=gv.pal[i_pop] ) 
    # plt.plot(theta, cos_func, color=gv.pal[i_pop] ) 
    

plt.xlabel('$\\theta$ (rad)') 
plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
           ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{2}$', r'$\pi$'])

plt.ylabel('Rates (Hz)') 

plt.show()
