from importlib import reload

import inspect, sys
import numpy as np 
import matplotlib.pyplot as plt

sys.path.insert(0, '../')

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

from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 

def get_X_S1_X_S2(**options): 

    data.get_days() # do not delete that !!        
    
    X1, X2 = get_X_S1_X_S2_days_task(day=options['day'],
                                         stimulus=options['stimulus'],
                                         task=options['task'],
                                         trials=options['trials'])

    print('X_S1', X1.shape, 'X_S2', X2.shape) 
    
    X1, X2 = pp.preprocess_X_S1_X_S2(X1, X2,
                                     scaler=options['scaler_BL'],
                                     center=options['center_BL'],
                                     scale=options['scale_BL'],
                                     avg_mean=options['avg_mean_BL'],
                                     avg_noise=options['avg_noise_BL'],
                                     unit_var=options['unit_var']) 
    
    X1 = pp.avg_epochs(X1, ['STIM', 'ED', 'DIST', 'MD', 'CUE',  'LD'])
    X2 = pp.avg_epochs(X2, ['STIM', 'ED', 'DIST', 'MD', 'CUE', 'LD'])
    
    return X1, X2

def get_overlap(X_D1, X_D2, Delta_sample, **kwargs): 
    
    # X_S1, X_S2 = scale_data(X_S1, X_S2, scaler=standardize) 
    # Delta_sample = get_coding_direction(X_S1, X_S2, pval) 
    
    # X_D1, X_D2 = scale_data(X_D1, X_D2, scaler=kwargs['scaler'])
    Delta_distractor = -get_coding_direction(X_D1, X_D2, **kwargs) 
    
    overlap = np.zeros( Delta_sample.shape[-1] )
    
    if kwargs['bins'] is not None: 
        # Delta_bins = np.nanmean( Delta_sample[:, kwargs['bins']], axis=-1) 
        Delta_bins = Delta_sample[:,1]
        for i_epoch in range(Delta_sample.shape[-1]): 
            overlap[i_epoch] = cos_between(Delta_bins, Delta_distractor[:, i_epoch]) 
    else: 
        for i_epoch in range(Delta_sample.shape[-1]): 
            overlap[i_epoch] = cos_between(Delta_sample[:,i_epoch], Delta_distractor[:, i_epoch]) 
    
    return overlap 

def overlap_time(**options):

    print('get sample axis')
    options['stimulus'] = 'sample'
    
    dum = options['task'] 
    options['task'] = 'all'
    X_S1, X_S2 = get_X_S1_X_S2(**options)
    
    dX_sample = get_coding_direction(X_S1, X_S2, **options) 
    
    print('get distractor axis')
    options['stimulus'] = 'distractor'
    options['task'] = dum 
    X_D1, X_D2 = get_X_S1_X_S2(**options)
    
    print('get overlap sample/dist')

    if(options['bins']=='ED'):
        options['bins'] = gv.bins_ED        
    elif(options['bins']=='MD'):
        options['bins'] = gv.bins_MD 
    elif(options['bins']=='STIM'):
        options['bins'] = gv.bins_STIM 
    elif(options['bins']=='DIST'):
        options['bins'] = gv.bins_DIST 
    else:
        options['bins'] = None 
    print('bins', options['bins'])
    
    overlap = get_overlap(X_D1, X_D2, dX_sample, **options) 
    
    ci = []
    if options['ci']:
        print('get bootstrap ci')
        mean, ci = my_bootstraped_ci(X_D1, X_D2,
                                     statfunction=lambda x, y: get_overlap(x, y, dX_sample, **options),
                                     n_samples=options['n_samples'] )     
        print('ci', ci.shape) 

    shuffle = [] 
    if options['shuffle']:
        print('get shuffle statistics') 
        shuffle = shuffle_stat(X_D1, X_D2,
                               lambda x,y: -get_overlap(x, y, dX_sample, **options), 
                               n_samples=options['n_shuffles'] )     
        print('shuffle', shuffle.shape) 

    return overlap, ci, shuffle

def plot_overlap(overlap, overlap_ci, overlap_shuffle, **options):
    
    if(options['bins'] ==None):
        figtitle = 'overlap_equal_time_%s_day_%s_%s_trials' % ( options['task'], str(options['day']), options['trials'] )
    else:
        figtitle = 'overlap_ED_%s_day_%s_%s_trials' % ( options['task'], str(options['day']), options['trials'] )
    
    fig = plt.figure(figtitle) 

    gv.time = [0,1,2,3,4,5]
    plt.plot(gv.time, overlap,'ko') 
    plt.xlabel('Time (s)')
    
    if options['bins'] is not None : 
        plt.ylabel('Cosine\n' r'Early Sample vs. Dist') 
    else: 
        plt.ylabel('Cosine\n' r'Sample vs. Dist') 
    
    plt.ylim([-.3,.3])
    plt.yticks([-.3,-.2,-.1,0,.1,.2,.3]) 
    # plt.xlim([0, 14]) 
    # plt.xticks([0,2,4,6,8,10,12,14]) 
    # plt.xlim([0,3]) 
    
    if options['ci']:
        plt.errorbar(gv.time, overlap, yerr=overlap_ci.T, ls=None) 
    
    if options['shuffle']:
        mean = np.nanmean(overlap_shuffle, axis=0) 
        std = np.nanpercentile(overlap_shuffle, [2.5, 97.5], axis=0) 
        plt.plot(gv.time, mean, 'k--') 
        plt.fill_between(gv.time, std[0], std[1], alpha=.1, color='k') 
    
    plt.ylim([-.25, 0.5]) 
    plt.yticks([-0.25, 0, .25, .5]) 
    # plt.xticks([0,2,4,6,8,10,12,14]) 
    plt.xticks( [0,1,2,3,4,5], ['Sample', 'Early delay', 'Distractor',
                    'Middle delay', 'Cue', 'Late delay'])
    
    # pl.add_vlines() 
    # plt.text(2.5, 0.55, 'Sample', horizontalalignment='center', fontsize=10) 
    # plt.text(5, 0.55, 'Dist.', horizontalalignment='center', fontsize=10) 
    # plt.text(7, 0.55, 'Cue', horizontalalignment='center', fontsize=10) 
    # plt.text(9.5, 0.55, 'Test', horizontalalignment='center', fontsize=10) 
    
    pl.save_fig(figtitle) 
    plt.close('all') 
    
if __name__ == '__main__':
    
    kwargs = dict() 
    kwargs['T_WINDOW']= 0.5 
    kwargs['ci'] = 1
    kwargs['n_samples']= 1000 
    kwargs['shuffle'] = 0
    kwargs['n_shuffles']= 1000 
    
    kwargs['n_days'] = 9 
    kwargs['prescreen'] = True 
    kwargs['pval']= .05  # .05, .01, .001 
    
    kwargs['standardize'] = False 
    kwargs['fit_intercept'] = True 
    
    # kwargs['clf_name'] = 'LDA' 
    # kwargs['shrinkage'] = 'auto' 
    kwargs['clf_name'] = 'LogisticRegressionCV' 
    kwargs['penalty'] = 'l2' 
    kwargs['n_lambda'] = 100 
    kwargs['inner_score']= 'accuracy' 
    kwargs['in_fold'] = 'loo' 
    kwargs['n_in'] = 5 
    kwargs['alpha'] = 1.0 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3]
        kwargs['trials'] = sys.argv[4] 
        kwargs['bins'] = sys.argv[5] 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    create_figdir(**options)
    
    overlap, overlap_ci, overlap_shuffle = overlap_time(**options) 
    plot_overlap(overlap, overlap_ci, overlap_shuffle, **options)
    
