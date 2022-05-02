import os, sys, importlib 

import numpy as np
import matplotlib 
matplotlib.use('GTK3cairo') 

import matplotlib.pyplot as plt 
from scipy.ndimage import uniform_filter1d
from scipy.signal import savgol_filter

from sklearn.preprocessing import StandardScaler

from joblib import Parallel, delayed 
from oasis.functions import deconvolve
from sklearn.feature_selection import VarianceThreshold 

from . import fct_facilities as fac 

fac.SetPlotParams() 

from . import constants as gv 
from . import preprocessing as pp 
from . import plotting as pl 
from . import utils as data 
from . import featureSel as fs
from . import progressbar as pg

from models.glms import *

def getAllDays(day=None): 

    data.get_days() 
    
    if day is None:
        day_list = gv.days 
        print('all days') 
        gv.all_days = 1 
    elif day=='first': 
        print('first days') 
        gv.first_days = 1 
        day_list = gv.days[:3] 
    elif day=='last': 
        print('last days') 
        gv.last_days = 1 
        day_list = gv.days[3:] 
    elif day=='all': 
        gv.all_days = 1 
        print('all days') 
        day_list = gv.days 
    elif isinstance(day, int): 
        print('single day') 
        day_list = [day] 
        
    print('days', day_list)
        
    for i_day, gv.day in enumerate(day_list):
        
        X_all, y_all = data.get_fluo_data() 
        if i_day==0:
            X_days = np.zeros( ( len(day_list), 3, 2 , int(gv.n_trials/2), gv.n_neurons, len(gv.time) ) ) 

        X_trials, _ = data.get_X_S1_S2(X_all, y_all) 
        
        X_days[i_day, :,:,:,0:X_trials.shape[3]] = X_trials 
        
    X_days = np.swapaxes(X_days, 2, 3) 
    X_days = np.hstack(X_days) 
    X_days = np.swapaxes(X_days, 1, 2) 
    
    print('X_days',X_days.shape) 
    
    return X_days

def create_figdir(**kwargs):
    globals().update(kwargs) 
    pl.figDir()
    
    try:
        gv.figdir = gv.figdir.replace("/z_score_bl", "")
    except:
        pass
    
    try:
        gv.figdir = gv.figdir.replace("/EDvsLD", "") 
    except:
        pass 
    
    if gv.all_days == 1 :        
        gv.figdir = gv.figdir + '/fluo_traces/%s/%s/all_days' % (gv.tasks[kwargs['i_task']], gv.mouse) 
    elif gv.first_days == 1 :
        gv.figdir = gv.figdir + '/fluo_traces/%s/%s/first_days' % (gv.tasks[kwargs['i_task']], gv.mouse)
    elif gv.last_days == 1 :
        gv.figdir = gv.figdir + '/fluo_traces/%s/%s/last_days' % (gv.tasks[kwargs['i_task']], gv.mouse)
    else:
        gv.figdir = gv.figdir + '/fluo_traces/%s/%s/day_%d' % (gv.tasks[kwargs['i_task']], gv.mouse, gv.day)
        
    if not os.path.isdir(gv.figdir): 
        os.makedirs(gv.figdir) 
        print(gv.figdir) 
def FluoTrace(day=None, mean=0, i_neuron=None, seed=None, **kwargs): 
    # set options and get globals 
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    data.get_days() 
    data.get_stimuli_times() 
    data.get_delays_times() 
    
    X_tasks = getAllDays(day=day)
    create_figdir(**options) # must come after getAllDays() 
    
    X_tasks = np.swapaxes(X_tasks, 0, 1) 
    X_tasks = np.hstack(X_tasks) 

    if mean:
        obj_str = 'mean'         
        title = 'population average'
        figtitle = 'population_average_trace'
    else:
        obj_str = 'neuron'
        title = 'neuron %d' % i_neuron
        figtitle = 'single_trace_neuron_%d' % i_neuron
        
    plt.figure(figtitle, figsize=pl.set_size(160)) 
    plt.title(title) 
    plt.ylabel('Fluo. (a.u.)') 
    plt.xlabel('Time (s)') 
    plt.xticks([0,2,4,6,8,10,12,14])
    pl.add_vlines()
    
    preprocess_X(X):
    X = X_tasks[i_task]
    X = preprocess_X(X) 
    
    if mean:
        X = np.mean(X, axis=1) # avg over neurons
    
    print('mouse', gv.mouse, 'day', gv.day, 'task', gv.task_str[options['i_task']], 'X', X.shape) 
    
    # set trial, neuron and seed 
    if seed is None:
        seed = np.random.randint(0, 1e6) 
        np.random.seed(seed)
    
    if i_neuron is None:
        i_neuron = np.random.randint(0, X.shape[1])
                
    print('X', X.shape) 
                    
                
    # for i_trial in range(5): 
    #     # i_trial = np.random.randint(0, X.shape[0]) 
    #     if mean==0: 
    #         plt.plot(gv.time, X[i_trial, i_neuron], 'k', alpha=0.25, color=gv.pal[i_task]) 
    #     else: 
    #         plt.plot(gv.time, X[i_trial], 'k', alpha=0.25, color=gv.pal[i_task]) 
            
    # plot median over trials and percentiles 
    X_avg = np.mean(X, axis=0)
    X_std = np.std(X, axis=0)

    lower = X_avg - X_std
    upper = X_avg + X_std
        
    # alpha = 0.95 
    # p = ((1.0-alpha)/2.0)*100 
    # lower = np.percentile(X, p, axis=0) 
            
    # p = ( alpha + (1.0-alpha)/2.0)*100 
    # upper = np.percentile(X, p, axis=0)
        
    if mean==0: 
        plt.plot(gv.time,  X_avg[i_neuron], color=gv.pal[i_task], label='raw') 
        # plt.fill_between(gv.time, lower[i_neuron], upper[i_neuron], color=gv.pal[i_task], alpha=.25)
    else: 
        plt.plot(gv.time,  X_avg, color=gv.pal[i_task], label='raw') 
        plt.fill_between(gv.time, lower, upper, alpha=.1, color=gv.pal[i_task]) 
        
    plt.ylim([-.05,.15])
    plt.yticks([-.05,0,.05,.1, .15])
    
    pl.save_fig(figtitle) 

def popFluoTrace(day=None, mean=0, i_neuron=None, seed=None, **kwargs): 
    # set options and get globals 
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    data.get_days() 
    data.get_stimuli_times() 
    data.get_delays_times() 
    
    X_tasks = getAllDays(day=day)
    create_figdir(**options) # must come after getAllDays() 
    
    X_tasks = np.swapaxes(X_tasks, 0, 1) 
    X_tasks = np.hstack(X_tasks) 

    if mean:
        obj_str = 'mean'         
        title = 'population average'
        figtitle = 'population_average_trace'
    else:
        obj_str = 'neuron'
        title = 'neuron %d' % i_neuron
        figtitle = 'single_trace_neuron_%d' % i_neuron
        
    plt.figure(figtitle, figsize=pl.set_size(160)) 
    plt.title(title) 
    plt.ylabel('Fluo. (a.u.)') 
    plt.xlabel('Time (s)') 
    plt.xticks([0,2,4,6,8,10,12,14])
    pl.add_vlines()
        
    for i_task in range(3):
        X = X_tasks[i_task]
        
        if mean:
            X = np.mean(X, axis=1) # avg over neurons
    
        print('mouse', gv.mouse, 'day', gv.day, 'task', gv.task_str[options['i_task']], 'X', X.shape) 
    
        # set trial, neuron and seed 
        if seed is None:
            seed = np.random.randint(0, 1e6) 
            np.random.seed(seed)
    
        if i_neuron is None:
            i_neuron = np.random.randint(0, X.shape[1])

        if gv.SAVGOL:
            X = savgol_filter(X, int(np.ceil(gv.frame_rate) * 2 + 3), polyorder = gv.SAVGOL_ORDER, deriv=0, axis=-1, mode='mirror') 
                
        print('X', X.shape) 
                    
                
        # for i_trial in range(5): 
        #     # i_trial = np.random.randint(0, X.shape[0]) 
        #     if mean==0: 
        #         plt.plot(gv.time, X[i_trial, i_neuron], 'k', alpha=0.25, color=gv.pal[i_task]) 
        #     else: 
        #         plt.plot(gv.time, X[i_trial], 'k', alpha=0.25, color=gv.pal[i_task]) 
            
        # plot median over trials and percentiles 
        X_avg = np.mean(X, axis=0)
        X_std = np.std(X, axis=0)

        lower = X_avg - X_std
        upper = X_avg + X_std
        
        # alpha = 0.95 
        # p = ((1.0-alpha)/2.0)*100 
        # lower = np.percentile(X, p, axis=0) 
            
        # p = ( alpha + (1.0-alpha)/2.0)*100 
        # upper = np.percentile(X, p, axis=0)
        
        if mean==0: 
            plt.plot(gv.time,  X_avg[i_neuron], color=gv.pal[i_task], label='raw') 
            # plt.fill_between(gv.time, lower[i_neuron], upper[i_neuron], color=gv.pal[i_task], alpha=.25)
        else: 
            plt.plot(gv.time,  X_avg, color=gv.pal[i_task], label='raw') 
            plt.fill_between(gv.time, lower, upper, alpha=.1, color=gv.pal[i_task]) 
        
    plt.ylim([-.05,.15])
    plt.yticks([-.05,0,.05,.1, .15])
    
    pl.save_fig(figtitle) 

def FluoDist(day=None, mean=0, i_neuron=None, seed=None, **kwargs): 
    # set options and get globals 
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    data.get_days() 
    data.get_stimuli_times() 
    data.get_delays_times() 
    
    X = getAllDays(day=day)
    create_figdir(**options) # must come after getAllDays() 
    
    X = np.swapaxes(X, 0, 1) 
    X = np.hstack(X)
    
    X = X[i_task] 
    
    print('mouse', gv.mouse, 'day', gv.day, 'task', gv.task_str[options['i_task']], 'X', X.shape) 

    if gv.SAVGOL:
        X = savgol_filter(X, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = gv.SAVGOL_ORDER, deriv=0, axis=-1)
    
    BL = X[..., gv.bins_BL] 
    F0 = np.mean(BL, axis=-1) 
    
    print(gv.bins_BL) 
    
    trial = np.random.randint(0, F0.shape[0])
    
    print(trial, X.shape, X[trial, 0:5, 0])
    print(trial, F0.shape, F0[trial, 0:5])
    
    plt.figure('single_trial_F0')
    plt.hist(F0[trial], int(F0.shape[1]/10)) 
    plt.axvline(np.mean(F0[trial]), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('F0') 

    avgF0 = np.mean(F0, axis=0)
    
    plt.figure('trial_averaged_F0')
    plt.hist(avgF0, int(avgF0.shape[0]/10)) 
    plt.axvline(np.mean(avgF0), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('$<F0>_{trials}$') 

    gv.DCV_THRESHOLD=0 
    rates = pp.deconvolveFluo(X) 
    ratesBL = np.mean(rates[..., gv.bins_BL], axis=-1) 

    print(trial, rates.shape, rates[trial, 0:5, 0])
    print(trial, ratesBL.shape, ratesBL[trial, 0:5])
    
    plt.figure('single_trial_dcv')
    plt.hist(ratesBL[trial], int(ratesBL.shape[1]/10)) 
    plt.axvline(np.mean(ratesBL[trial]), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('$r_0$') 

    avgRatesBL= np.mean(ratesBL, axis=0)
    
    plt.figure('trial_averaged_dcv')
    plt.hist(avgRatesBL, int(avgRatesBL.shape[0]/10)) 
    plt.axvline(np.mean(avgRatesBL), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('$<r_0>_{trials}$') 
    
