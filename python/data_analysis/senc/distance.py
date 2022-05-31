#!/usr/bin/env python3
from importlib import reload
import inspect, sys
import numpy as np 
import matplotlib.pyplot as plt

sys.path.insert(0, '../')

# plt.off()

import utils.constants as gv 
reload(gv)
from utils.options import *

import utils.get_data as data
reload(data)

from utils.get_days import * 

import utils.preprocessing as pp
reload(pp)
import utils.plot_utils as pl 
reload(pl)

import senc.utils
reload(senc.utils)
from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 
from senc.sel_time import *

def get_X_S1_X_S2(**options): 

    data.get_days() # do not delete that !!        
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'],
                                         stimulus=options['stimulus'],
                                         task=options['task'],
                                         trials=options['trials'])
    
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=options['scaler_BL'],
                                         center=options['center_BL'],
                                         scale=options['scale_BL'],
                                         avg_mean=options['avg_mean_BL'],
                                         avg_noise=options['avg_noise_BL'],
                                         unit_var=options['unit_var']) 
    
    return X_S1, X_S2

def get_coding_dir(**kwargs):

    options = kwargs.copy()
    options['trials'] = 'correct'
    
    if options['stimulus'] == 'distractor':
        options['task'] = 'Dual'
        options['tasks'] = ['Dual']
    else:
        options['task'] = 'all'
        options['tasks'] = ['all'] 
    
    # get data
    X1, X2 = get_X_S1_X_S2(**options) 
    
    # average over given epoch 
    X1_avg = pp.avg_epochs(X1, options['bins']) 
    X2_avg = pp.avg_epochs(X2, options['bins']) 
    
    # get coding direction of the sample memory 
    coefs = get_coefs(X1_avg, X2_avg, return_clf=0, **options) 

    return coefs
    
def get_proj(X, coefs, model):
    
    # scaler = model.named_steps['scaler'] 
    # X_rescaled = (X - scaler.mean_[np.newaxis, : , np.newaxis]) / scaler.var_[np.newaxis, : , np.newaxis]    
    # coefs *= scaler.var_
    
    X_proj = np.zeros( X.shape[-1] ) 
    
    for i_epoch in range(X.shape[-1]): 
        X_proj[i_epoch] = np.mean( np.dot(coefs, X[..., i_epoch].T), axis=0) 
        
        # for i_trials in range(X.shape[0]):
        #     X_proj[i_epoch] += np.dot(coefs, X[i_trials, :, i_epoch]) 
    
    # X_avg = np.mean(X, axis=0)     
    # X_proj = np.dot(coefs,  X_avg )
    
    return -X_proj 

def get_avg_proj(X_S1, X_S2, coefs, model):
    
    X_proj = np.zeros( X_S1.shape[-1] ) 
    
    for i_epoch in range(X_S1.shape[-1]): 
        # X_proj[i_epoch] = np.sum( np.dot(coefs, X_S1[..., i_epoch].T), axis=0 ) 
        # X_proj[i_epoch] -= np.sum( np.dot(coefs, X_S2[..., i_epoch].T), axis=0 ) 

        a = np.mean( np.dot(coefs, X_S2[..., i_epoch].T), axis=0 )
        b = np.mean( np.dot(coefs, X_S1[..., i_epoch].T), axis=0 ) 
        X_proj[i_epoch] = np.sqrt((a-b)**2)
    
    return X_proj 

# return X_proj / (X_S1.shape[0]+X_S2.shape[0]) 

def get_projections(**options):
    
    coefs = get_coding_dir(**options) 
    print('coefs', coefs.shape)

    eps = 1
    if options['stimulus'] == 'distractor':
        eps = -1
        options['task'] = 'Dual'
    
    # options['stimulus'] = 'sample' 
    # get projection of sample trials onto sample axis 
    print(options['stimulus'], options['task'])
    X1, X2 = get_X_S1_X_S2(**options) 
    models=[1,2]
    
    P1 = get_proj(X1, coefs[0], models[0]) 
    P2 = get_proj(X2, coefs[0], models[0]) 
    P12 = get_avg_proj(X1, X2, coefs[0], models[0]) 

    P1_ci = []
    _, P1_ci = my_bootstraped_ci(X1, X2,
                                 statfunction=lambda x, y: get_proj(x, coefs[0], models[0]) 
                                 , n_samples=options['n_samples'] ) 

    P2_ci = []
    _, P2_ci = my_bootstraped_ci(X1, X2,
                                 statfunction=lambda x, y: get_proj(y, coefs[0], models[0]) 
                                 , n_samples=options['n_samples'] ) 

    P12_ci = []
    _, P12_ci = my_bootstraped_ci(X1, X2,
                                  statfunction=lambda x, y: get_avg_proj(x, y, coefs[0],
                                                                         models[0])
                                  , n_samples=options['n_samples'] ) 

    P12_shuffle = []
    P12_shuffle = shuffle_stat(X1, X2,
                               lambda x, y: get_avg_proj(x, y, coefs[0], models[0]),
                               n_samples=options['n_shuffles'])
    
    return P1, P2, P12, P1_ci, P2_ci, P12_ci, P12_shuffle

if __name__ == '__main__':
    
    kwargs = dict() 
    
    kwargs['ci'] = 0 
    kwargs['n_samples'] = 1000 
    kwargs['shuffle'] = 0 
    kwargs['n_shuffles'] = 1000 
    
    kwargs['T_WINDOW'] = 0.5 
    kwargs['sample'] = 'S1' 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'robust' 
    kwargs['avg_mean_BL'] = 1 
    kwargs['avg_noise_BL'] = 1 
    kwargs['unit_var'] = 1
    
    kwargs['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo']) 
    # kwargs['tasks'] = ['DPA', 'Dual'] 
    # kwargs['clf']='logitnetAlphaCV' 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
        # kwargs['stimulus'] = sys.argv[5]
        # kwargs['sample'] = sys.argv[7]

    kwargs['n_days'] = 9 
    
    kwargs['clf_name'] = 'LogisticRegressionCV'
    # kwargs['clf_name'] = 'LDA'
    # kwargs['clf_name'] = 'bolasso'
    # kwargs['penalty'] = 'l1'
    kwargs['n_lambda'] = 100 
    
    kwargs['prescreen'] = True  
    kwargs['pval']= .05 # .05, .01, .001  
    
    kwargs['standardize'] = False 
    kwargs['fit_intercept'] = True 
    
    # kwargs['inner_score']= 'roc_auc' 
    # kwargs['inner_score']= 'neg_log_loss' 
    kwargs['inner_score']= 'accuracy' 
    kwargs['in_fold'] = 'loo'
    kwargs['n_in'] = 5 
    kwargs['alpha'] = 1.0 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    options['bins'] = ['ED']
    options['stimulus'] = 'sample' 
    S1, S2, S12, S1_ci, S2_ci, S12_ci, S12_shuffle = get_projections(**options) 

    S12_mean_shuffle  = np.nanmean(S12_shuffle, axis=0) 
    S12_perc_shuffle = np.nanpercentile(S12_shuffle, [2.5, 97.5], axis=0) 
    
    options['bins'] = ['MD']
    options['stimulus'] = 'distractor'
    D1, D2, D12, D1_ci, D2_ci, D12_ci, D12_shuffle = get_projections(**options) 

    D12_mean_shuffle  = np.nanmean(D12_shuffle, axis=0) 
    D12_perc_shuffle = np.nanpercentile(D12_shuffle, [2.5, 97.5], axis=0) 
    
    # plot figure
    create_figdir(**options) 
    
    options['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo', 'Dual']) 
    options['task'] = sys.argv[2] 
    options['i_task'] = np.argwhere(options['tasks']==options['task'])[0][0] 
    
    # figname = 'overlap_sample_dist_' + options['task'] + '_day_' + options['day'] + '_' + options['trials'] 
    figname = 'distance_' + options['task'] + '_day_' + options['day'] + '_' + options['trials'] 
    
    # if options['task']=='DPA': 
    fig, axis = plt.subplots(1, 2, figsize=(1.25*1.618*1.5*2, 1.618*1.5), num=figname) 
    pl.add_vlines(axis[0]) 
    pl.add_vlines(axis[1]) 
    
    axis[0].set_xlabel('Time (s)') 
    axis[0].set_ylabel('Sample Readout') 
    axis[0].set_xticks([0,2,4,6,8,10,12,14]) 
    axis[0].set_ylim([-.25,5]) 
    
    axis[1].set_xlabel('Time (s)') 
    axis[1].set_ylabel('Distractor Readout') 
    axis[1].set_xticks([0,2,4,6,8,10,12,14]) 
    axis[1].set_ylim([-.25,8]) 
    # else:
    #     axis = plt.gcf().get_axes() 
    
    # axis[0].plot(gv.time, S1, color=gv.pal[options['i_task']])
    # axis[0].fill_between(gv.time, S1-S1_ci[:,0], S1+S1_ci[:,1],
    #                      alpha=0.1, color=gv.pal[options['i_task']]) 
    
    # axis[0].plot(gv.time, S2, color=gv.pal[options['i_task']], ls='--')
    # axis[0].fill_between(gv.time, S2-S2_ci[:,0], S2+S2_ci[:,1],
    #                      alpha=0.1, color=gv.pal[options['i_task']], ls='--') 
    
    axis[0].plot(gv.time, S12, color=gv.pal[options['i_task']])
    axis[0].fill_between(gv.time, S12-S12_ci[:,0], S12+S12_ci[:,1],
                         alpha=0.1, color=gv.pal[options['i_task']]) 
    
    axis[0].plot(gv.time, S12_mean_shuffle, '--' , color=gv.pal[options['i_task']]) 
    axis[0].fill_between(gv.time, S12_perc_shuffle[0], S12_perc_shuffle[1], color=gv.pal[options['i_task']], alpha=.1) 
    
    # axis[1].plot(gv.time, D1, color=gv.pal[options['i_task']])
    # axis[1].fill_between(gv.time, D1-D1_ci[:,0], D1+D1_ci[:,1],
    #                      alpha=0.1, color=gv.pal[options['i_task']]) 

    # axis[1].plot(gv.time, D2, color=gv.pal[options['i_task']], ls='--')
    # axis[1].fill_between(gv.time, D2-D2_ci[:,0], D2+D2_ci[:,1],
    #                      alpha=0.1, color=gv.pal[options['i_task']], ls='--') 
    
    axis[1].plot(gv.time, D12, color=gv.pal[options['i_task']])
    axis[1].fill_between(gv.time, D12-D12_ci[:,0], D12+D12_ci[:,1],
                         alpha=0.1, color=gv.pal[options['i_task']]) 
    
    axis[1].plot(gv.time, D12_mean_shuffle, '--' , color=gv.pal[options['i_task']]) 
    axis[1].fill_between(gv.time, D12_perc_shuffle[0], D12_perc_shuffle[1], color=gv.pal[options['i_task']], alpha=.1) 
    
    pl.save_fig(figname,1) 
    
    plt.show()
    plt.close('all')
    
