import sys, os, importlib 
from importlib import reload 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from joblib import Parallel, delayed 
import progressbar as pgb 

import constants as gv 
importlib.reload(sys.modules['constants']) 

gv.init_param()

# raw_spike_times = pd.read_csv('/homecentral/alexandre.mahrach/bistable_mongilo/src/data/spikes.txt', sep='\s+').to_numpy() 
raw_spike_times = pd.read_csv(gv.path + 'spike_times.dat', sep='\s+').to_numpy() 
print('data', raw_spike_times.shape) 

neurons_id = raw_spike_times[:,0] 
print('neurons_id', neurons_id.shape, neurons_id[500:505]) 

spike_times = raw_spike_times[:,1] 
print('spike_times', spike_times.shape, spike_times[500:505]) 

figtitle = 'raster' 
fig = plt.figure(figtitle, figsize=(1.25*1.618*1.5*10, 1.618*1.25*2)) 

excitatory_idx = np.where(neurons_id<gv.n_size) 
inhibitory_idx = np.where(neurons_id>=gv.n_size) 

plt.scatter(spike_times[excitatory_idx], neurons_id[excitatory_idx], marker='|', alpha=0.25, color='r') 
plt.scatter(spike_times[inhibitory_idx], neurons_id[inhibitory_idx], marker='|', alpha=0.25, color='b') 

plt.xlabel('Time (ms)') 
plt.ylabel('Neuron #') 

# plt.xlim([0, 2000]) 
