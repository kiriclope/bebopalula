import sys, os, importlib 
from importlib import reload 
from scipy.signal import savgol_filter 

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from joblib import Parallel, delayed

import params as gv 
importlib.reload(sys.modules['params']) 
import progressbar as pgb 

from utils import *

gv.IF_SPEC=0 
gv.IF_LOW_RANK=0 
gv.IF_INI_COND=0 
gv.IF_TRIALS=0 
gv.IF_DPA = 0
gv.IF_DUAL = 0
gv.IF_DRT = 0
gv.IF_ODR = 0 

gv.IF_STP=1

gv.RANK=2
gv.init_param() 

kappa_list = np.arange(3.5, 4.05, .05) 
trial_list = np.arange(0, 8, 1) 

MEAN_KSI = gv.MEAN_KSI 
VAR_KSI = gv.VAR_KSI 
print(gv.path) 

def parloop(kappa, trial, n_pop=gv.n_pop, path=gv.path, MAP=gv.MAP, ksi_path=gv.ksi_path, RANK=gv.RANK, verbose=0, n_size=gv.n_size, con_path=gv.con_path) : 
    
    if(RANK==1): 
        ksi_path = con_path + '/low_rank/rank_1/seed_ksi_3'
        path += '/low_rank/kappa_%.2f' % kappa
    if(RANK==2): 
        ksi_path = con_path + '/low_rank/rank_2/seed_ksi_3'
        path += '/low_rank/kappa_%.2f_kappa_1_%.2f' % (kappa, kappa)
        
    path += '/trial_%d' % trial ; 
    if verbose: 
        print(path)
        print(ksi_path)
    
    # if 0==0:
    try: 
        time, rates = get_time_rates(0, path, con_path)

        rates = rates[-10:]
        if(verbose):
            print('time', time.shape, 'rates', rates.shape)
        
        avg_rates = np.nanmean(rates.copy(), axis=0) # avg over time 
        mean_rates = np.nanmean(avg_rates.copy(), axis=-1) # mean over neurons 
        
        if verbose: 
            print('rates', mean_rates, '<rates>', avg_rates.shape) 
        
        overlap = get_overlap(rates.copy(), ksi_path=ksi_path, MAP=0, n_size=n_size) 
        overlap = np.nanmean(overlap, axis=-1) 
        
        if RANK==2:
            overlap_1 = get_overlap(rates.copy(), ksi_path=ksi_path, MAP=1, n_size=n_size) 
            overlap_1 = np.nanmean(overlap_1, axis=-1) 
        else:
            overlap_1 = overlap.copy() 
        
        if verbose:
            print('overlap', overlap) 
            
    except :
        mean_rates = np.zeros(n_pop) * np.nan 
        overlap = np.zeros(n_pop) * np.nan 
        overlap_1 = np.zeros(n_pop) * np.nan 

        if verbose:
            print('error') 
    
    return mean_rates, overlap, overlap_1

parloop(10, 0, verbose=1) 

with pgb.tqdm_joblib( pgb.tqdm(desc='computing rates and overlap', total=int(kappa_list.shape[0] * trial_list.shape[0]) ) ) as progress_bar: 
        
    rates_list, overlap_list, overlap_1_list = zip( *Parallel(n_jobs=-1, backend="multiprocessing")(delayed(parloop)( kappa, trial) 
                                                                                                    for kappa in kappa_list 
                                                                                                    for trial in trial_list ) ) 

rates_list = np.array(rates_list).T 
overlap_list = np.array(overlap_list).T 
overlap_1_list = np.array(overlap_1_list).T 

rates_list = rates_list.reshape(gv.n_pop, kappa_list.shape[0], trial_list.shape[0]) 
overlap_list = overlap_list.reshape(gv.n_pop, kappa_list.shape[0], trial_list.shape[0]) 
overlap_1_list = overlap_1_list.reshape(gv.n_pop, kappa_list.shape[0], trial_list.shape[0]) 

print('rates', rates_list.shape, 'overlap', overlap_list.shape, 'overlap_1', overlap_1_list.shape)
# print(rates_list)

mean_rates = np.nanmean(rates_list, axis=-1) 
mean_overlap = np.nanmean(overlap_list, axis=-1) 
mean_overlap_1 = np.nanmean(overlap_1_list, axis=-1) 

std_rates = np.nanstd(rates_list, axis=-1) 
std_overlap = np.nanstd(overlap_list, axis=-1) 
std_overlap_1 = np.nanstd(overlap_list, axis=-1) 

print('rates', mean_rates) 
print('overlap', mean_overlap) 

figname = 'rates_overlap_kappa_' + gv.folder 
fig, axis = plt.subplots(1, gv.RANK+1, figsize=(1.25*1.618*1.5*(gv.RANK+1), 1.618*1.5), num=figname) 
    
axis[0].set(ylim=(0 , np.ceil( np.nanmax(rates_list)*12 ) /10 ) ) 
axis[0].set_xlabel('$var(\\xi)$') 
axis[0].set_ylabel('Rates (Hz)') 

# axis[1].set(ylim=(np.ceil( -np.nanmax(np.abs(overlap_list))*12 ) /10 , np.ceil( np.nanmax(np.abs(overlap_list))*12 ) /10 ) ) 
axis[1].set_xlabel('$var(\\xi^{Sample})$') 
axis[1].set_ylabel('Sample Overlap (a.u.)') 

if gv.RANK==2:
    # axis[2].set(ylim=(np.ceil( -np.nanmax(np.abs(overlap_1_list))*12 ) /10 , np.ceil( np.nanmax(np.abs(overlap_1_list))*12 ) /10 ) ) 
    axis[2].set_xlabel('$var(\\xi^{Dist})$') 
    axis[2].set_ylabel('Dist. Overlap (a.u.)') 

for i_pop in range(gv.n_pop):
    axis[0].plot(kappa_list, mean_rates[i_pop], '-o', color = gv.pal[i_pop]) 
    axis[0].plot(kappa_list, rates_list[i_pop,:], 'x', color = gv.pal[i_pop], alpha=0.2) 
    axis[0].fill_between(kappa_list, mean_rates[i_pop]-std_rates[i_pop], mean_rates[i_pop]+std_rates[i_pop], alpha=.1, color=gv.pal[i_pop])

def overlap_signed(x, sign=1):
    x_signed = x.copy()
    # print(x.shape)
    
    for kappa in range(x.shape[0]):
        for trial in range(x.shape[1]):        
            if x_signed[kappa,trial]>0 and sign==-1: 
                x_signed[kappa,trial] = np.nan
            if x_signed[kappa,trial]<0 and sign==1: 
                x_signed[kappa,trial] = np.nan
                
    return x_signed

overlap_plus = overlap_signed(overlap_list[0])
# overlap_minus = np.nanmean( -np.abs(overlap_list[0]), axis=-1) 
# print(overlap_plus)

# for i_pop in range(gv.n_pop):
# axis[1].plot(kappa_list, mean_overlap[0], '-o', color = gv.pal[0], label=gv.label[0]) 
# axis[1].fill_between(kappa_list, mean_overlap[i_pop]-std_overlap[i_pop], mean_overlap[i_pop]+std_overlap[i_pop], alpha=.1, color=gv.pal[i_pop])
# axis[1].plot(kappa_list, overlap_plus, '-o', color = gv.pal[0], alpha=1) 
axis[1].plot(kappa_list, overlap_list[0], 'o', color = gv.pal[0], alpha=1) 

if gv.RANK==2:
    axis[2].plot(kappa_list, overlap_1_list[0], 'o', color = gv.pal[0], alpha=1) 

plt.show() 
