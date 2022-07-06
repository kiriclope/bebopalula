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

gv.N_INI = 10 
gv.init_param()

path = gv.path

def parloop_spkC(neurons_id, i_neuron): 

    idx_neuron = neurons_id==i_neuron
    spike_count = np.sum(idx_neuron) 
    
    return spike_count 

def get_spike_count(path):
    
    raw_spike_times = pd.read_csv(path + '/spike_times.dat', sep='\s+').to_numpy() 
    # print('data', raw_spike_times.shape)     
    neurons_id = raw_spike_times[:,0] 
    # print('neurons_id', neurons_id.shape, neurons_id[500:505]) 
    
    spike_times = raw_spike_times[:,1] 
    print('spike_times', spike_times.shape, spike_times[500:505]) 
    
    neurons_id = neurons_id[spike_times<2000] # time in ms 
    
    n_neurons = gv.n_neurons * 10000 
    
    with pgb.tqdm_joblib( pgb.tqdm(desc='spike_count', total=n_neurons) ) as progress_bar: 
        spike_counts = Parallel(n_jobs=-1, backend='multiprocessing')(delayed(parloop_spkC)(neurons_id, i_neuron) for i_neuron in range(n_neurons) ) 
    
    return spike_counts

def get_fano_factor(path):
    spike_count_inis = []
    
    for i_ini in range(1, gv.N_INI + 1):
        gv.path = path
        gv.path += '/ini_cond_%d' % i_ini ; 
        print(gv.path)
        
        spike_count_inis.append(get_spike_count(gv.path)) 
    
    spike_count_inis = np.asarray(spike_count_inis) 
    print('spike_count', spike_count_inis.shape, spike_count_inis[0, :10]) 
    
    mean_spike_count = np.nanmean(spike_count_inis, axis=0) # over inis 
    var_spike_count = np.nanvar(spike_count_inis, axis=0) 
    
    fano_factor = var_spike_count / mean_spike_count 
    
    print('fano_factor', fano_factor.shape, fano_factor[:10]) 
    
    return fano_factor, spike_count_inis 

# fano_off, spike_off = get_fano_factor(path)
path = path.replace('off', 'on') # change dirname 
fano_on, spike_on = get_fano_factor(path)
