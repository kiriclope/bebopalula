import sys, os, importlib
from importlib import reload

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import constants as gv 
importlib.reload(sys.modules['constants']) 

gv.init_param()

from balance_inputs_dist import inputs_dist, vec_Phi

# mean_field = inputs_dist() 
# mean_var = mean_field.x 

# if gv.n_pop==2 :
#     mf_inputs_E = np.random.normal(loc=mean_var[0], scale=np.sqrt(mean_var[2]), size=gv.n_size) 
#     mf_inputs_I = np.random.normal(loc=mean_var[1], scale=np.sqrt(mean_var[3]), size=gv.n_size) 
    
#     mf_rates_E = vec_Phi(mf_inputs_E)*1000
#     mf_rates_I = vec_Phi(mf_inputs_I)*1000 
# else :
#     mf_inputs = np.random.normal(loc=mean_var[0], scale=np.sqrt(mean_var[1]), size=gv.n_size) 
#     mf_rates = vec_Phi(mf_inputs) 
    
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
    
# print('rates', rates.shape)

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

avg_rates = np.mean(rates, axis=0) 
# print(avg_rates.shape)

for i_pop in range(gv.n_pop):
    plt.hist(avg_rates[i_pop], histtype='step') 

# if gv.n_pop==2:
#     plt.hist(mf_rates_I, histtype='step', ls='--', color='b') 
#     plt.hist(mf_rates_E, histtype='step', ls='--', color='r') 
# else:
#     plt.hist(mf_rates, histtype='step', ls='--', color='b') 

plt.xlabel('Rates (Hz)')
plt.ylabel('Count')

plt.show()
