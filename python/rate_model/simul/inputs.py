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
# else : 
#     mf_inputs = np.random.normal(loc=mean_var[0], scale=np.sqrt(mean_var[1]), size=gv.n_size) 
    
try :
    inputs = pd.read_csv(gv.path + 'inputs.dat', sep='\s+').to_numpy()
except:
    print('file not found  ...  \n')
    raise

time = inputs[:,0] 
print('time', time.shape)

inputs = np.delete(inputs, [0], axis=1)
print('inputs', inputs.shape) 

if gv.n_pop!=1:
    n_neurons = int(inputs.shape[1]/4)
    inputs = np.reshape(inputs, (inputs.shape[0], 2, 2, n_neurons))
    
    print('inputs per pop', inputs.shape)

mean_inputs = np.mean(inputs, axis=-1) 
    
print('mean_inputs', mean_inputs.shape) 

figtitle = 'inputs' 
fig = plt.figure(figtitle, figsize=(1.25*1.618*1.5*2, 1.618*1.25*2)) 

if gv.n_pop==1:
    
    ext_inputs = np.sqrt(gv.K) * gv.ext_inputs * np.ones((len(time), inputs.shape[-1]) ) 
    net_inputs = inputs + ext_inputs
        
    ax = fig.add_subplot(int('121'))
    
    plt.plot(time, mean_inputs, 'b', lw=2)
    plt.plot(time, np.mean(ext_inputs, axis=-1), 'r', lw=2) 
    plt.plot(time, np.mean(net_inputs, axis=-1), 'k', lw=2) 
    
    plt.xlabel('time (ms)') 
    plt.ylabel('inputs (mA)') 

    ax = fig.add_subplot(int('122')) 
        
    plt.hist( np.mean( ext_inputs, axis=0), color='r', histtype='step') 
    plt.hist( np.mean( inputs, axis=0), color='b', histtype='step') 
    plt.hist( np.mean( net_inputs, axis=0), color='k', histtype='step') 
    
    plt.hist( mf_inputs, color='k', ls='--', histtype='step') 
    
else:
    for i_pop in range(gv.n_pop):
        ax = fig.add_subplot( int('22%d'% (i_pop+1) ) ) 
        net_inputs = np.sqrt(gv.K) * gv.ext_inputs[i_pop] + mean_inputs[..., 0, i_pop] + mean_inputs[..., 1, i_pop]
        
        plt.plot(time, np.sqrt(gv.K) * gv.ext_inputs[i_pop] + mean_inputs[..., 0, i_pop], 'r', lw=2) 
        plt.plot(time, net_inputs, 'k', lw=2) 
        plt.plot(time, mean_inputs[..., 1, i_pop], 'b', lw=2) 
        
        plt.xlabel('Time (ms)') 
        plt.ylabel('Inputs (mA)')
        if i_pop==0:
            plt.title('Excitatory Pop.')
        else:
            plt.title('Inhibitory Pop.')
            
        ax = fig.add_subplot(int('22%d' %(i_pop+3))) 
        
        net_inputs = np.sqrt(gv.K) * gv.ext_inputs[i_pop] + np.mean( inputs[:, 0, i_pop], axis=0) + np.mean( inputs[:, 1, i_pop], axis=0) 
        plt.hist( np.sqrt(gv.K) * gv.ext_inputs[i_pop] + np.mean( inputs[:, 0, i_pop], axis=0) , color='r', histtype='step') 
        plt.hist( np.mean( inputs[:, 1, i_pop], axis=0), color='b', histtype='step') 
        plt.hist( net_inputs, color='k', histtype='step') 

        # if i_pop==0: 
        #     plt.hist( mf_inputs_E, color='k', ls='--', histtype='step') 
        # else: 
        #     plt.hist( mf_inputs_I, color='k', ls='--', histtype='step') 
        
        plt.xlabel('Inputs (mA)') 
        plt.ylabel('Count') 

plt.show()
