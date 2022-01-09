import sys, os, importlib
from importlib import reload
from scipy.signal import savgol_filter

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import constants as gv 
importlib.reload(sys.modules['constants']) 

gv.init_param()
    
filter_rates = pd.read_csv(gv.path + 'filter_rates.dat', sep='\s+').to_numpy() 
# print(filter_rates.shape)

time = filter_rates[:,0] 
# print('time', time.shape)

rates = np.delete(filter_rates, [0], axis=1)
# print('rates', rates.shape) 

if gv.n_pop!=1:
    n_neurons = int(rates.shape[1]/2)
    rates = np.reshape(rates, (rates.shape[0], 2, n_neurons))
else: 
    n_neurons = rates.shape[0] 
    
print('rates', rates.shape)

mean_rates = np.mean(rates, axis=-1) 

# print(mean_rates.shape) 

avg_mean_rates = np.mean(mean_rates, axis=0) 
print('avg_mean_rates', avg_mean_rates) 

figtitle = 'rates_' + gv.folder 
fig = plt.figure(figtitle, figsize=(1.25*1.618*1.5*2, 1.618*1.25)) 
ax = fig.add_subplot(int('121'))

# for _ in range(5) :
#     i_neuron = np.random.randint(0, n_neurons) 
#     plt.plot(time, rates[..., i_neuron], alpha=0.25) 

plt.plot(time, mean_rates[:,0], lw=2, color='r') 
plt.plot(time, mean_rates[:,1], lw=2, color='b')  
   
plt.xlabel('Time (ms)') 
plt.ylabel('Rates (Hz)') 
# plt.ylim([0, 100]) 

ax = fig.add_subplot(int('122')) 
theta = np.linspace(0, np.pi, gv.n_size)

filter_rates = savgol_filter(rates, int(np.ceil(1000./2.0) * 2 + 1), polyorder = 1, deriv=0, axis=-1, mode='mirror') 

fft_rates = np.fft.rfft(filter_rates, axis=-1) / gv.n_size 

print(fft_rates.shape)
m1 = fft_rates[...,1].real 

# OSI = m1 / mean_rates ;

for i_pop in range(gv.n_pop):
    plt.plot(time, m1[:,i_pop], color=gv.pal[i_pop]) 

plt.xlabel('Time (ms)') 
plt.ylabel('$m_1$') 

# avg_rates = np.mean(rates, axis=0) 

# filter_rates = savgol_filter(avg_rates, int(np.ceil(500./2.0) * 2 + 1), polyorder = 0, deriv=0, axis=-1, mode='mirror') 
# # print(avg_rates.shape) 

# for i_pop in range(gv.n_pop):
#     plt.plot(theta, filter_rates[i_pop], color=gv.pal[i_pop]) 

# plt.xlabel('Theta (rad)')
# plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
#            ['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{2}$', r'$\pi$'])
# plt.ylabel('Rates (Hz)') 

plt.show()
