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

        X_trials = data.get_X_S1_S2(X_all, y_all) 
        
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

    X = getAllDays(day=day)
    create_figdir(**options) # must come after getAllDays() 
    
    X = np.swapaxes(X, 0, 1) 
    X = np.hstack(X)
    
    X = X[i_task] 
    print('X', X.shape)
    
    # # get raw data 
    # X_all, y_all = data.get_fluo_data() 
    
    # # get given task trials 
    # task_labels = [0, 13, 14] 
    # bool_task = (y_all[4]==task_labels[options['i_task']]) & (y_all[8]==0) 
    # task_idx = np.where( bool_task )[0] 
    
    # X = X_all[task_idx] 
    print('mouse', gv.mouse, 'day', gv.day, 'task', gv.task_str[options['i_task']], 'X', X.shape) 

    # set trial, neuron and seed 
    if seed is None:
        seed = np.random.randint(0, 1e6) 
    np.random.seed(seed)
    
    if i_neuron is None:
        i_neuron = np.random.randint(0, X.shape[1])
        
    if gv.Z_SCORE or gv.Z_SCORE_BL: 
        X = pp.z_score(X)
        
        if mean==0:
            figtitle = 'z_score_neuron_%d' % i_neuron 
            plt.figure(figtitle, figsize=(2.1*1.25*1.5, 1.85*1.25)) 
            plt.title('neuron %d' % i_neuron ) 
        else:
            X = np.mean(X, axis=1) # avg over neurons
            
            figtitle = 'population_z_score'
            plt.figure(figtitle, figsize=(2.1*1.25*1.5, 1.85*1.25))         
            plt.title('population average') 
            
        plt.ylabel('Z-score')
        
    else:
        
        if mean ==0:
            figtitle = 'single_trace_neuron_%d' % i_neuron
            plt.figure(figtitle, figsize=(2.1*1.25*1.5, 1.85*1.25))         
            plt.title('neuron %d' % i_neuron ) 
        else:
            X = np.mean(X, axis=1) # avg over neurons
            
            figtitle = 'population_trace' 
            plt.figure(figtitle, figsize=(2.1*1.25*1.5, 1.85*1.25)) 
            plt.title('population average') 
            
        plt.ylabel('Fluorescence (a.u.)') 
        
    plt.xlabel('Time (s)') 
    plt.xticks([0,2,4,6,8,10,12,14])    
    pl.add_vlines()
        
    for _ in range(5):
        i_trial = np.random.randint(0, X.shape[0]) 
        if mean==0: 
            plt.plot(gv.time, X[i_trial, i_neuron], 'k', alpha=0.25) 
        else: 
            plt.plot(gv.time, X[i_trial], 'k', alpha=0.25) 
            
    # plot median over trials and percentiles 
    X_avg = np.mean(X, axis=0)
    
    alpha = 0.95 
    p = ((1.0-alpha)/2.0)*100 
    lower = np.percentile(X, p, axis=0) 

    p = ( alpha + (1.0-alpha)/2.0)*100 
    upper = np.percentile(X, p, axis=0)
    
    # X_avg = np.mean(X, axis=0) 
    # ci = 1.96 * np.std(X, axis=0)/np.mean(X, axis=0) 
    # lower = X_avg - ci
    # upper = X_avg + ci 
    
    if mean==0: 
        plt.plot(gv.time,  X_avg[i_neuron], 'k', label='raw') 
        plt.fill_between(gv.time, lower[i_neuron], upper[i_neuron], color='c', alpha=.25)
        
        if gv.Z_SCORE or gv.Z_SCORE_BL: 
            plt.ylim([-10,30]) 
        else:
            plt.ylim([-1,1])
    else: 
        plt.plot(gv.time,  X_avg, 'k', label='raw')
        plt.fill_between(gv.time, lower, upper, color='c', alpha=.25)
        
        # if gv.Z_SCORE or gv.Z_SCORE_BL: 
        #     # plt.ylim([-.2,1]) 
        #     plt.ylim([-.5,5]) 
        # else:
        #     plt.ylim([-.15,.25]) 
        
    pl.save_fig(figtitle) 
    
def populationTrace():
    scaler = StandardScaler()
    gv.mouse=gv.mice[2] 

    gv.n_days = 6 
    gv.day=1
    
    data.get_days() 
    data.get_stimuli_times() 
    data.get_delays_times() 

    gv.data_type = 'raw' 
    gv.NORMALIZE = 0 
    
    gv.Z_SCORE = 0 
    gv.Z_SCORE_BL = 1
    
    gv.F0_THRESHOLD = None 
    
    gv.DECONVOLVE = 0 
    gv.DCV_THRESHOLD = 0 
    
    gv.DETREND = 0
    gv.DETREND_ORDER = 7
    
    X_all, y_labels = data.get_fluo_data()

    i_task = 0 
    task_labels = [0, 13, 14] 
    bool_task = (y_labels[4]==task_labels[i_task]) & (y_labels[8]==0) 
    task_idx = np.where( bool_task )[0] 
    
    X = X_all[task_idx]
    
    print('mouse', gv.mouse, 'session', gv.session, 'data X', X.shape) 
    
    data.get_delays_times() 
    data.get_bins() 
    
    if gv.mouse in [gv.mice[0]]: 
        gv.n_trials = 40 
    else: 
        gv.n_trials = 32 
        
    time = np.linspace(0, gv.duration, X.shape[-1]) 
    
    figtitle = 'population_trace' 
    plt.figure(figtitle, figsize=(2.1*1.25*2, 1.85*1.25))
    plt.title('%s'% gv.task_str[i_task])
    seed = np.random.randint(0, 1e6) 
    np.random.seed(seed)

    X = np.mean(X, axis=1) # avg over neurons
    for _ in range(5):
        i_trial = np.random.randint(0, X.shape[0])
        plt.plot(time, X[i_trial], 'k', label='raw', alpha=0.25) 
        
    plt.xlabel('Time (s)') 
    plt.xticks([0,2,4,6,8,10,12,14]) 
    plt.ylabel('Fluorescence (a.u.)') 
        
    pl.add_vlines() 
    X_avg = np.median(X, axis=0) # average over trials 
    lower = np.percentile(X, 25, axis=0) 
    upper = np.percentile(X, 75, axis=0) 
    
    plt.plot(time,  X_avg, 'k', label='raw') 
    plt.fill_between(time, lower, upper, color='b', alpha=.1) 
    
    pl.save_fig(figtitle) 
    
    if gv.Z_SCORE or gv.Z_SCORE_BL: 
        print('Z_SCORE')         
        X = pp.z_score(X) 
        
        figtitle = 'population_zscore' 
        plt.figure(figtitle, figsize=(2.1*1.25*2, 1.85*1.25)) 
        plt.title('%s'% gv.task_str[i_task])
        
        np.random.seed(seed) 
        
        for _ in range(5): 
            i_trial = np.random.randint(0, X.shape[0]) 
            plt.plot(time, X[i_trial], 'k', label='raw', alpha=0.25) 
            
        X_avg = np.median(X, axis=0) 
        lower = np.percentile(X, 25, axis=0) 
        upper = np.percentile(X, 75, axis=0) 
        
        plt.plot(time,  X_avg, 'k', label='raw') 
        plt.fill_between(time, lower, upper, color='b', alpha=.1) 
        
        plt.xlabel('Time (s)') 
        plt.xticks([0,2,4,6,8,10,12,14]) 
        plt.ylabel('z-score (a.u.)')                 
        pl.add_vlines() 
        
        pl.save_fig(figtitle) 

def X_trace():
    scaler = StandardScaler()
    gv.mouse=gv.mice[2] 
    
    gv.n_days = 6 
    gv.day=1
    
    data.get_days() 
    data.get_stimuli_times() 
    data.get_delays_times() 
    
    gv.data_type = 'raw' 
    gv.NORMALIZE = 1 
    
    gv.Z_SCORE = 0
    gv.Z_SCORE_BL = 1
    
    gv.F0_THRESHOLD = 0
    
    gv.DECONVOLVE = 0
    gv.DCV_THRESHOLD = 0
    
    gv.DETREND = 0
    gv.DETREND_ORDER = 7
    
    X, y = data.get_fluo_data() 
    
    print('mouse', gv.mouse, 'session', gv.session, 'data X', X.shape,'y', y.shape) 
    
    data.get_delays_times() 
    data.get_bins() 
    
    if gv.mouse in [gv.mice[0]]: 
        gv.n_trials = 40 
    else: 
        gv.n_trials = 32 
    
    trial = np.random.randint(0, X.shape[0]) 
    neuron = np.random.randint(0, X.shape[1]) 
    # trial = 5 
    # neuron = 10 
    
    time = np.linspace(0, gv.duration, X.shape[-1]) 
    
    plt.figure('X') 

    # gv.SAVGOL = 0 
    # gv.SAVGOL_ORDER = 1 
    
    # if gv.SAVGOL: 
    #     X = savgol_filter(X, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = gv.SAVGOL_ORDER, deriv=0, axis=-1) 
    
    # idx = fs.featSel.var_fit_transform(np.mean(X[..., gv.bins_STIM], axis=-1), threshold=0.1) 
    # print(idx.shape) 
    # X = X[:,idx,:]
    
    X_avg = np.mean(X, axis=1) 
    X_avg = scaler.fit_transform(X_avg.T).T
    plt.plot(time, X_avg[trial], 'k')
    
    if gv.DETREND :
        X_det = pp.detrend_X(X, order=gv.DETREND_ORDER) 
        X_avg = np.mean(X_det, axis=1) 
        plt.plot(time, X_avg[trial], 'k') 

    if gv.F0_THRESHOLD is not None:
        print('F0')
        dFF = pp.dFF0_remove_silent(X) 
        X_avg = np.mean(dFF, axis=1) 
        X_avg = scaler.fit_transform(X_avg.T).T

        plt.plot(time, X_avg[trial], 'g') 
        
        if gv.DETREND :
            dFF_det = pp.detrend_X(dFF, order=gv.DETREND_ORDER) 
            X_avg = np.mean(dFF_det, axis=1)                        
            plt.plot(time, X_avg[trial], 'g') 

    if gv.DECONVOLVE: 
        print('DCV') 
        X_dcv = pp.deconvolveFluo(X) 
        X_avg = np.mean(X_dcv, axis=1) 
        X_avg = scaler.fit_transform(X_avg.T).T
        
        plt.plot(time, X_avg[trial], 'k') 
        
        if gv.DETREND:
            X_det = pp.detrend_X(X_dcv, order=gv.DETREND_ORDER) 
            X_avg = np.mean(X_det, axis=1) 
            plt.plot(time, X_avg[trial], 'k') 
            
    if gv.NORMALIZE: 
        print('NORM')
        X_norm = pp.normalize(X) 
        X_avg = np.mean(X_norm, axis=1) 
        # X_avg = scaler.fit_transform(X_avg.T).T
        
        plt.plot(time, X_avg[trial], 'r') 
        
        if gv.DETREND :
            X_norm = pp.detrend_X(X_norm, order=gv.DETREND_ORDER) 
            X_avg = np.mean(X_norm, axis=1) 
            plt.plot(time, X_avg[trial], 'r') 
        
    if gv.Z_SCORE or gv.Z_SCORE_BL:
        print('Z_SCORE')
        X_score = pp.z_score(X)
        
        X_avg = np.mean(X_score, axis=1)        
        X_avg = scaler.fit_transform(X_avg.T).T
        plt.plot(time, X_avg[trial], 'b')

        if gv.DETREND :
            X_score = pp.detrend_X(X_score, order=gv.DETREND_ORDER) 
            X_avg = np.mean(X_score, axis=1) 
            plt.plot(time, X_avg[trial], 'b') 

    plt.ylabel('X (a.u.)') 
    plt.xlabel('time (s)') 
    
    pl.add_vlines() 
        
def F0_dist():
    gv.mouse=gv.mice[2] 

    gv.n_days = 6 
    gv.day=1
    
    data.get_sessions_mouse() 
    data.get_stimuli_times() 
    data.get_delays_times() 

    gv.type = 'dF'
    X, y = data.get_fluo_data()
    
    X_score = np.mean( pp.z_score(X), axis=1)
    X = pp.normalize(X) 
    
    print('mouse', gv.mouse, 'session', gv.session, 'data X', X.shape,'y', y.shape) 
    
    data.get_delays_times() 
    data.get_bins()
    
    if gv.mouse in [gv.mice[0]]: 
        gv.n_trials = 40 
    else: 
        gv.n_trials = 32
        
    F0 = np.mean(X[...,gv.bins_BL],axis=-1)
    print('F0', F0.shape)
    
    # trial = np.random.randint(0,F0.shape[0])    
    # neuron = np.random.randint(0,F0.shape[1])
    trial = 0 
    neuron = 10
    
    plt.figure('F0 neuron')
    plt.hist(F0[trial], int(F0.shape[1]/2)) 
    plt.axvline(np.mean(F0[trial]), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('F0') 

    plt.figure('F0 trials') 
    x = F0.T[neuron] 
    print('x', x.shape) 
    plt.hist(x, int(F0.shape[0]/2)) 
    plt.axvline(np.mean(x), c='k', ls='-') 
    plt.ylabel('# trials') 
    plt.xlabel('F0') 
    
    plt.figure('<F0>') 
    plt.hist(np.mean(F0,axis=0), int(F0.shape[1]/2)) 
    plt.axvline(np.mean(np.mean(F0)), c='k', ls='-') 
    plt.ylabel('# neurons') 
    plt.xlabel('<F0>_{trials}') 

    print('trial', trial, 'F0', np.mean(F0[trial]), '<F0>', np.mean(F0) ) 

    
    F0 = F0[..., np.newaxis]         
    
    # removing silent neurons
    X_avg = np.mean(X, axis=1)
    
    idx = np.where(F0<=0.2)
    F0 = np.delete(F0, idx, axis=-2) 
    X = np.delete(X, idx, axis=-2)
    
    print('X', X.shape, 'F0', F0.shape)
    
    dFF = (X-F0) / (F0 + gv.eps) 

    bins = np.arange(0, dFF.shape[-1])
    dFF_avg = np.mean(dFF, axis=1)
    dFF_std = np.std(dFF, axis=1) 
    
    plt.figure('dFF') 
    # plt.plot(bins, X_avg[trial], 'k')
    # plt.plot(bins, X_score[trial], 'r')
    gv.duration = X.shape[-1]/gv.frame_rate 
    time = np.linspace(0,gv.duration,X.shape[-1]) 
    
    # plt.errorbar(time, dFF_avg[trial], yerr=dFF_std[trial],c='b')
    plt.plot(time, dFF_avg[trial])
    pl.add_vlines()
    plt.ylabel('<dFF> (a.u.)') 
    plt.xlabel('time (s)') 
    print(gv.duration, len(gv.bins_BL), gv.frame_rate)
    print(gv.bins_BL)
    print(gv.bins_STIM)
    
def test_deconvolve(threshold=0.5):
    
    gv.mouse=gv.mice[2] 
    
    data.get_sessions_mouse() 
    data.get_stimuli_times() 
    data.get_delays_times() 
        
    gv.session=gv.sessions[-1]
    gv.day=1 
    
    X, y = data.get_fluo_data() 
    print('mouse', gv.mouse, 'session', gv.session, 'data X', X.shape,'y', y.shape) 
    
    data.get_delays_times() 
    data.get_bins()
    
    if gv.mouse in [gv.mice[0]]: 
        gv.n_trials = 40 
    else: 
        gv.n_trials = 32
        
    F0 = np.mean(X[...,gv.bins_BL],axis=-1)
    # F0 = np.mean( np.mean(X[...,gv.bins_BL],axis=-1), axis=0 ) 
    # F0 = F0[np.newaxis] 
    # F0 = np.percentile(X, 15, axis=-1) 
    X_noise = X - uniform_filter1d( X, int(gv.frame_rate/4) ) # gv.frame = 1s , gv.frame/2=0.5s , gv.frame/4=0.25s 

    # def F0_loop(X, n_trial, n_neuron, bins): 
    #     X_ij = X[n_trial, n_neuron] 
    #     c, s, b, g, lam = deconvolve(X_ij, penalty=1) 
    #     return b 
    
    # # loop over trials and neurons 
    # with pg.tqdm_joblib(pg.tqdm(desc='F0', total=X.shape[0]*X.shape[1])) as progress_bar: 
    #     F0 = Parallel(n_jobs=gv.num_cores)(delayed(F0_loop)(X, n_trial, n_neuron, gv.bins_BL) 
    #                                         for n_trial in range(X.shape[0]) 
    #                                         for n_neuron in range(X.shape[1]) ) 
    # F0 = np.array(F0).reshape( (X.shape[0], X.shape[1]) )
    
    def X_loop(X, F0, n_trial, n_neuron):
        X_ij = X[n_trial, n_neuron]
        F0_ij = F0[n_trial, n_neuron]
        # F0_ij = F0[0, n_neuron] 
        
        # c, s, b, g, lam = deconvolve(X_ij, penalty=1, b_nonneg=False) 
        c, s, b, g, lam = deconvolve(X_ij, penalty=1, b=F0_ij, b_nonneg=False) 
        # c, s, b, g, lam = deconvolve(X_ij,  g=(None, None), penalty=1, b=F0_ij) 
        return c
        
    def S_loop(X, F0, n_trial, n_neuron):
        X_ij = X[n_trial, n_neuron] 
        F0_ij = F0[n_trial, n_neuron] 
        # F0_ij = F0[0, n_neuron] 
        
        # c, s, b, g, lam = deconvolve(X_ij, penalty=1, b_nonneg=False) 
        c, s, b, g, lam = deconvolve(X_ij, penalty=1, b=F0_ij) 
        # c, s, b, g, lam = deconvolve(X_ij, g=(None, None), penalty=1, b=F0_ij) 
        return s 
    
    def corr_loop(X, X0, n_trial, n_neuron):
        X_ij = X[n_trial, n_neuron] 
        X0_ij = X0[n_trial, n_neuron] 
        corr_ij = np.corrcoef(X_ij, X0_ij)[0,1] 
        return corr_ij 
    
    # loop over trials and neurons 
    with pg.tqdm_joblib(pg.tqdm(desc='denoise', total=X.shape[0]*X.shape[1])) as progress_bar: 
        X_dcv = Parallel(n_jobs=gv.num_cores)(delayed(X_loop)(X, F0, n_trial, n_neuron) 
                                            for n_trial in range(X.shape[0]) 
                                            for n_neuron in range(X.shape[1]) ) 
    X_dcv = np.array(X_dcv).reshape(X.shape) 
    
    # print(np.corrcoef(X[0,0], X_dcv[0,0])) 
    
    # with pg.tqdm_joblib(pg.tqdm(desc='corr', total=X.shape[0]*X.shape[1])) as progress_bar: 
    #     corr = Parallel(n_jobs=gv.num_cores)(delayed(corr_loop)(X_dcv, X, n_trial, n_neuron) 
    #                                          for n_trial in range(X.shape[0]) 
    #                                          for n_neuron in range(X.shape[1]) )        
    # corr =  np.array(corr).reshape(X.shape)
    # mean_corr = np.mean(corr)
    
    # print("Mean corr of deconvolved activity  with ground truth ('Fluo') : %.4f" %  mean_corr) 
    
    with pg.tqdm_joblib(pg.tqdm(desc='deconvolve', total=X.shape[0]*X.shape[1])) as progress_bar: 
        S_dcv = Parallel(n_jobs=gv.num_cores)(delayed(S_loop)(X, F0, n_trial, n_neuron) 
                                              for n_trial in range(X.shape[0]) 
                                              for n_neuron in range(X.shape[1]) ) 
        
    S_dcv = np.array(S_dcv).reshape(X.shape)

    # def R_loop(S, n_trial, n_neuron):
    #     spiketrain = S[n_trial, n_neuron] 
    #     R = instantaneous_rate(spiketrain, sampling_period=50*ms)
    #     return R
    
    # with pg.tqdm_joblib(pg.tqdm(desc='firing rate', total=X.shape[0]*X.shape[1])) as progress_bar: 
    #     R_dcv = Parallel(n_jobs=gv.num_cores)(delayed(R_loop)(X, n_trial, n_neuron) 
    #                                           for n_trial in range(X.shape[0]) 
    #                                           for n_neuron in range(X.shape[1]) ) 

    trial = np.random.randint(0,X.shape[0])
    neuron = np.random.randint(0,X.shape[1])
    # trial = 0 
    # neuron = 10
    plt.figure('F0')
    plt.hist(F0[trial], int(F0.shape[1]/2)) 
    plt.ylabel('bins') 
    plt.xlabel('F0') 

    plt.figure('Traces')        
    plt.plot(gv.time, X[trial,neuron], '-b') 
    plt.plot(gv.time, X_dcv[trial,neuron], '-r') 
    plt.plot(gv.time, X_noise[trial,neuron], '-c') 
    plt.plot(gv.time, S_dcv[trial,neuron], '-g') 
    plt.xlabel('bins')
    plt.ylabel('fluo') 
    pl.add_vlines()
    
    # plt.axhline(np.percentile(X[trial,neuron], 15), c='k', ls='--') 
    plt.axhline(F0[0,neuron], c='k', ls='--') 
    
    plt.figure('spikes') 
    plt.plot(X_noise[trial,neuron], '-c') 
    plt.plot(S_dcv[trial,neuron], '-g') 
    plt.xlabel('bins')
    plt.ylabel('P(s)') 
    
    def threshold_spikes(S_dcv, threshold):
        S_th = S_dcv
        S_th[S_th<threshold] = 0 
        S_th[S_th>threshold] = 1 
        S_th = uniform_filter1d( S_dcv, int(gv.frame_rate/4) ) # gv.frame = 1s , gv.frame/2=0.5s , gv.frame/4=0.25s
        return S_th 
    
    S_th = threshold_spikes(S_dcv, threshold) 
    plt.plot(S_th[trial,neuron], 'k') 
    
    S_bl = np.mean(S_th[...,gv.bins_BL], axis=-1) 

    print('trial', trial, 'neuron', neuron, 'corrFluo', np.corrcoef(X[trial,neuron], X_dcv[trial,neuron])[0,1],
          'idv bl mean', S_bl[trial, neuron], 'population average bl', np.mean(S_bl[trial], axis=0) )
    
    S_bl = np.mean(S_bl, axis=0) 

    S_avg = np.mean(S_th[...,gv.bins_BL],axis=-1) 
    S_avg = S_avg[..., np.newaxis]
    
    print('X_avg', np.mean(S_avg))
    # removing silent neurons 
    idx = np.where(S_avg<=0) 
    S_th = np.delete(S_th, idx, axis=1)
    
    print('X_dcv', S_th.shape) 
    
    plt.figure('FR dist') 
    plt.hist(S_bl, int(S_bl.shape[0]/2)) 
    plt.xscale('log')
    plt.ylabel('bins')
    plt.xlabel('neuron')
    plt.xlim([0.1, 10])
    
    # S_flt = uniform_filter1d( S_dcv, int(gv.frame_rate/2) ) 
    # S_flt = savgol_filter(S_dcv, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = 5, deriv=0) 
    # plt.plot(S_flt[trial,neuron], 'k') 
    
    # def scaler_loop(S, n_trial, bins_BL): 
    #     S_i = S[n_trial] 
    #     scaler = StandardScaler() 
    #     scaler.fit(S_i[:,bins_BL].T) 
    #     return scaler.transform(S_i.T).T 
            
    # with pg.tqdm_joblib(pg.tqdm(desc='standardize', total=X.shape[0])) as progress_bar: 
    #     S_scaled = Parallel(n_jobs=gv.num_cores)(delayed(scaler_loop)(S_flt, n_trial, gv.bins_BL) 
    #                                              for n_trial in range(X.shape[0]) ) 
        
    # S_scaled = np.array(S_scaled) 
    # S_scaled = savgol_filter(S_scaled, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = 3, deriv=0) 
    # plt.plot(S_scaled[trial,neuron], 'y') 
