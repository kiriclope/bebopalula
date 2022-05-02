#!/usr/bin/env python3
import inspect, sys 
from importlib import reload 
sys.path.insert(0, '../')

import time 
import numpy as np 
import matplotlib.pyplot as plt

from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=False)

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 
from sklearn.pipeline import Pipeline 
from sklearn.preprocessing import StandardScaler, RobustScaler 

plt.ion() 

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

import matplotlib
matplotlib.use('gtk3Cairo')

def dum():

    count = 0
    data.get_days() 
    day_list = gv.days

    perf = []
    
    for i_day, day in enumerate(day_list):
        count = 0
    
        gv.day=day
        X_all, y_all = data.get_fluo_data()
        
        gv.task = "DualGo_S1_unpair_correct" 
        y = data.which_trials(y_all) 
        count += y.shape[0]
        print(y.shape)
    
        gv.task = "DualGo_S2_unpair_correct" 
        y = data.which_trials(y_all)     
        count += y.shape[0]
        print(y.shape)

        perf.append(count/16)
        
    print('CR', perf) 
    
def performance_days(**kwargs): 
    options = set_options(**kwargs) 
    set_globals(**options)
    
    create_figdir(**options) 
    
    data.get_days() 
    
    day_list = gv.days

    figtitle = 'perf_days' 
    fig = plt.figure(figtitle) 
    
    for i_task in range(3):
        perf_S1_S2 = np.zeros(len(day_list)) 
        for i_day, day in enumerate(day_list):
            X_trials, y_trials = get_X_y_day(day=day) 
            
            # print('X_trials', X_trials.shape, 'y_trials', y_trials.shape) 
            perf_S1_S2[i_day] = np.count_nonzero(~np.isnan(y_trials[i_task][1])) / y_trials[i_task][0].shape[0] / 16 
            
            print('day', day_list[i_day], 'performance: S1_S2', perf_S1_S2[i_day]) 

        plt.plot(day_list, perf_S1_S2, '-o', color=gv.pal[i_task]) 
    
    plt.xlabel('Day')
    plt.ylabel('Performance') 
    plt.ylim([0.25, 1.25]) 
    
    # pl.save_fig(figtitle) 
    # plt.close('all') 

if __name__ == '__main__':

    kwargs = dict() 
    
    kwargs['ci']=1
    kwargs['shuffle']=1
    kwargs['perm_test']=1
    
    kwargs['n_samples'] = 1000
    kwargs['n_shuffles'] = 1000 
    kwargs['T_WINDOW'] = .5 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'robust' 
    kwargs['avg_mean_BL'] = 1 
    kwargs['avg_noise_BL'] = 1 
    kwargs['unit_var'] = 1 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
        kwargs['obj'] = sys.argv[5] 
        kwargs['stimulus'] = sys.argv[6] 
        kwargs['t_train'] = sys.argv[7] 
    
    kwargs['tasks'] = ['DPA', 'DualGo', 'DualNoGo'] 
    kwargs['IF_SAVE']=0 
        
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    # options['clf'] = get_clf(**options) 
    # print('clf', options['clf']) 

    # dum()
    performance_days(**options)
    
    
