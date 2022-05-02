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

def my_stat(sel_perm_first_tasks, sel_perm_last_tasks, sel):
    
    stat =(sel_perm_first_tasks[0] - sel_perm_last_tasks[0])**2 + (sel_perm_first_tasks[1] - sel_perm_last_tasks[1])**2 + (sel_perm_first_tasks[2] - sel_perm_last_tasks[2])**2
    
    obs = (sel[0][0] - sel[1][0])**2 + (sel[0][1] - sel[1][1])**2 + (sel[0][2] - sel[1][2])**2 
    
    pval = np.mean( np.abs(stat) > np.abs(obs), axis=0 ) 

    print('stat', stat.shape, 'obs', obs.shape)
    
    return pval

def cols_high(obj='cos'):
    cols = [-.1,0,.1]
    if obj=='frac': 
        high = [0.375, 0.3] 
        low = [-.0,-.0,-.0] 
        corr = [-0.025, -0.025, -0.025] 
    if obj=='norm': 
        high = [1.1, 0.95] 
        low = [-.15,-.15,-.15] 
        corr = [-0.025, -0.025, -0.025] 
        # high = [2.1, 1.9] 
        # # high = [9-.25, 8-.25] 
        # low = [-.95,-.95,-.95] 
        # # corr = [-0.25, -0.25, -0.25] 
        # corr = [-0.025, -0.025, -0.025] 
    if obj=='cos': 
        high = [1.1, 0.95] 
        low = [-.15,-.15,-.15] 
        corr = [-0.025, -0.025, -0.025] 
    if obj=='proj': 
        high = [5, 4.5] 
        low = [-.25,-.25,-.25] 
        corr = [-0.1, -0.1, -0.1] 
    if obj=='score': 
        high = [1.05, 0.95] 
        low = [.45,.45,.45] 
        corr = [-0.025, -0.025, -0.025] 
    if obj=='non_zero': 
        high = [0.375, 0.3] 
        low = [-.0,-.0,-.0] 
        corr = [-0.025, -0.025, -0.025] 
    
    return cols, high, low, corr

def add_pval_shuffle(pval_shuffle, obj='cos'):
    cols, high, low, corr = cols_high(obj)
    
    for i_task in range(pval_shuffle.shape[0]): 
        for i_epoch in range(pval_shuffle.shape[1]): 
            
            if pval_shuffle[i_task, i_epoch]<=0.001: 
                plt.text( i_epoch * .4 + cols[i_task], low[i_task], "***",
                          ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_shuffle[i_task, i_epoch]<=.01: 
                plt.text( i_epoch*.4 + cols[i_task], low[i_task], "**",
                          ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_shuffle[i_task, i_epoch]<=.05: 
                plt.text( i_epoch*.4 + cols[i_task], low[i_task], "*",
                          ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_shuffle[i_task, i_epoch]>.05: 
                plt.text( i_epoch*.4 + cols[i_task], low[i_task], "n.s.",
                          ha='center', va='bottom', color='k', fontsize=10) 
                
def add_pval_perm(pval_perm, obj='cos'):
    cols, high, low, corr = cols_high(obj)

    for i_task in range(pval_perm.shape[0]): 
        for i_epoch in range(pval_perm.shape[1]): 
            plt.plot( [(i_epoch)*.4 + cols[0], (i_epoch)*.4  + cols[i_task+1]] , [high[i_task], high[i_task]] , lw=1, c='k') 
            
            if pval_perm[i_task, i_epoch]<=.001: 
                plt.text(( 2*(i_epoch)*.4 +cols[0] + cols[i_task+1])*.5, high[i_task] + corr[i_task], "***", 
                         ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_perm[i_task, i_epoch]<=.01: 
                plt.text(( 2*(i_epoch)*.4 +cols[0] + cols[i_task+1])*.5, high[i_task] + corr[i_task], "**", 
                         ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_perm[i_task, i_epoch]<=.05: 
                plt.text(( 2*(i_epoch)*.4 +cols[0] + cols[i_task+1])*.5, high[i_task] + corr[i_task], "*",
                         ha='center', va='bottom', color='k', fontsize=10) 
            elif pval_perm[i_task, i_epoch]>.05: 
                plt.text(( 2*(i_epoch)*.4 +cols[0] + cols[i_task+1])*.5, high[i_task], "n.s.",
                         ha='center', va='bottom', color='k', fontsize=10) 
    
def get_delta_all_epochs(**options): 
    
    data.get_days() # do not delete that !!        
    X_S1_all, X_S2_all = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task='all',
                                                 trials=options['trials']) 
    
    X_S1_all, X_S2_all = pp.preprocess_X_S1_X_S2(X_S1_all, X_S2_all,
                                                 scaler=options['scaler_BL'],
                                                 center=options['center_BL'], scale=options['scale_BL'],
                                                 avg_mean=options['avg_mean_BL'], avg_noise=options['avg_noise_BL'],
                                                 unit_var=options['unit_var'])
        
    X_S1_all = pp.avg_epochs(X_S1_all, ['ED']) 
    X_S2_all = pp.avg_epochs(X_S2_all, ['ED']) 
    
    Delta_all = get_coding_direction(X_S1_all, X_S2_all, **options) 
    print('Delta0', Delta_all.shape) 
    
    return Delta_all 

def get_sel_epochs(**options): 

    X_S1_days = []
    X_S2_days = []
    
    data.get_days() # do not delete that !!
    
    for day in ['first', 'last']:
        print(day)

        X_S1_tasks = [] 
        X_S2_tasks = [] 
        
        for i_task in range(len(options['tasks'])): 
            options['task'] = options['tasks'][i_task] 
            X_S1, X_S2 = get_X_S1_X_S2_days_task(day=day, stimulus=options['stimulus'],
                                                 task=options['task'], trials=options['trials'])
            
            print('X_S1', X_S1.shape) 
            
            X_S1_tasks.append(X_S1) 
            X_S2_tasks.append(X_S2) 
            
            X_S1_tasks[i_task], X_S2_tasks[i_task] = pp.preprocess_X_S1_X_S2(X_S1_tasks[i_task], X_S2_tasks[i_task],
                                                                             scaler=options['scaler_BL'],
                                                                             center=options['center_BL'],
                                                                             scale=options['scale_BL'],
                                                                             avg_mean=options['avg_mean_BL'],
                                                                             avg_noise=options['avg_noise_BL'],
                                                                             unit_var=options['unit_var']) 
            
            X_S1_tasks[i_task] = pp.avg_epochs(X_S1_tasks[i_task], options['epochs']) 
            X_S2_tasks[i_task] = pp.avg_epochs(X_S2_tasks[i_task], options['epochs']) 
        
        # X_S1_tasks = np.array(X_S1_tasks) 
        # X_S2_tasks = np.array(X_S2_tasks) 

        # print('X_S1_tasks', X_S1_tasks.shape)
        
        X_S1_days.append(X_S1_tasks) 
        X_S2_days.append(X_S2_tasks) 
    
    # X_S1_days = np.array(X_S1_days) 
    # X_S2_days = np.array(X_S2_days) 
    
    print('X_S1_days', X_S1_days[0][0].shape) 
        
    sel = [] 
    sel_ci = [] 
    sel_shuffle = [] 
    pval_shuffle = [] 
    
    options['bins'] = None 
    options['bins_ED'] = gv.bins_ED
    options['bins_MD'] = gv.bins_MD
    options['bins_LD'] = gv.bins_LD

    print('bins', options['bins_ED'], options['bins_MD'] )
    
    Delta_tasks = [] 
    alpha_tasks = [] 
    lbd_tasks = [] 
    tuned_clfs = [] 
    
    for i_day in range(2): 
        sel_task = [] 
        for i_task in range(len(options['tasks'])): 
            options['task'] = options['tasks'][i_task] 
            print('task', options['task'], 'X_S1', X_S1_days[i_day][i_task].shape, 'X_S2', X_S2_days[i_day][i_task].shape) 
            
            # sample selectivity 
            options['n_jobs'] = -10
            startbuild = time.time() 
            sel_task.append( get_sel( X_S1_days[i_day][i_task], X_S2_days[i_day][i_task], **options) ) 
            endbuild = time.time() 
            build_time = endbuild - startbuild 
            print("runtime: %.2f s" % build_time) 
            options['n_jobs'] = None 
        
        sel.append(sel_task) 
    
    sel = np.array(sel) 
    print('sel', sel) 
    
    # permutation test 
    pval_perm = []     
    options2 = options 

    sel_perm_first_tasks = [] 
    sel_perm_last_tasks = [] 

    if options['perm_test']:
        for i_task in range(len(kwargs['tasks'])): 
            
            sel_perm_first, sel_perm_last = get_sel_perm( X_S1_days[0][i_task], X_S2_days[0][i_task],
                                                          X_S1_days[1][i_task], X_S2_days[1][i_task],
                                                          lambda x,y: get_sel(x, y, **options),
                                                          lambda x,y: get_sel(x, y, **options2),
                                                          n_samples=options['n_samples'] ) 

            sel_perm_first_tasks.append(sel_perm_first)
            sel_perm_last_tasks.append(sel_perm_last)
            
            print('sel_perm', sel_perm_first.shape, sel_perm_last.shape) 
            
            pval = np.mean( np.abs( sel_perm_first - sel_perm_last) > np.abs(sel[0][i_task] - sel[1][i_task]), axis=0 ) 
            print(options['tasks'][i_task], 'perm_test pval', pval ) 
            
            pval_perm.append(pval)
            
        pval_stat = my_stat(sel_perm_first_tasks, sel_perm_last_tasks, sel) 
        print('new stat', 'perm_test pval', pval_stat ) 
          
        pval_perm = np.array(pval_perm) # get rid of Early delay 
        print('pval_perm', pval_perm.shape, pval_perm)
    
    return sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm 

def plot_sel_epochs(sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm, **options):
    create_figdir(**options) 
    
    figtitle = '%s_%s_epochs_%s_day_%s' % (options['mouse_name'][options['i_mice']], options['obj'], options['trials'], str(options['day']))
    if options['off_diag']:
        figtitle += '_'+ options['t_train'] +'_train' 
    
    fig = plt.figure(figtitle)
    
    cols = [-.1,0,.1] 
    
    for i_task, task in enumerate(options['tasks']):         
        
        # if options['obj']=='cos': 
        #     plt.plot([cols[i_task], .4 + cols[i_task]] , sel[i_task,1:3], 'o', color=gv.pal[i_task], ms=2) 
        # else: 
        plt.plot([cols[i_task], .4 + cols[i_task], .8 + cols[i_task]] , sel[i_task], 'o', color=gv.pal[i_task], ms=2) 
                    
        if options['ci']: 
            # if options['obj']=='cos': 
            #     plt.errorbar([cols[i_task], .4 + cols[i_task]] , sel[i_task,1:3], yerr=sel_ci[i_task, :, 1:3],
            #                  ecolor=gv.pal[i_task], color=gv.pal[i_task], ls='none') 
            # else: 
            plt.errorbar([cols[i_task], .4 + cols[i_task], .8 + cols[i_task]] , sel[i_task], yerr=sel_ci[i_task],
                         ecolor=gv.pal[i_task], color=gv.pal[i_task], ls='none') 
            
    if options['shuffle']: 
        add_pval_shuffle(pval_shuffle, obj=options['obj']) 
    
    if options['perm_test']: 
        add_pval_perm(pval_perm, obj=options['obj']) 
    
    if options['obj']!='cos': 
        if options['off_diag'] : 
            if options['t_train'] == 'ED': 
                plt.xticks([0,.4,.8], ['Early vs.\n' r"Early delay", 'Early vs.\n' r"Middle delay", 'Early vs.\n' r"Late delay"]) 
            else: 
                plt.xticks([0,.4,.8], ['Late vs.\n' r"Early delay", 'Late vs.\n' r"Middle delay", 'Late vs.\n' r"Late delay"]) 
        else: 
            plt.xticks([0,.4,.8], ['Early vs.\n' r"Early delay", 'Middle vs.\n' r"Middle delay", 'Late vs.\n' r"Late delay"]) 
        plt.xlim([-0.25, 1.05]) 
    else: 
        # plt.xticks([0,.4], ['Early vs.\n' r"Middle delay", 'Early vs.\n' r"Late delay"]) 
        # plt.xlim([-0.25, .65]) 

        plt.xticks([0,.4,.8], ['Early vs.\n' r"Middle delay", 'Early vs.\n' r"Late delay", 'Middle vs.\n' r"Late delay"]) 
        plt.xlim([-0.25, 1.05]) 
    
    if options['obj']=='norm': 
        plt.ylabel('Sample Sel.') 
        plt.ylim([-.25, 1.25]) 
        plt.yticks([0,.25,.5,.75,1]) 
        # plt.ylim([0, 2.5]) 
        # plt.yticks([0,2,4,6,8,10]) 
        # plt.yticks([0,.25,.5,.75,1])
        # plt.yticks([-1,-.5,0,.5,1,1.5,2, 2.5]) 
    
    if options['obj']=='cos':  
        # plt.ylabel('Overlap\n' r'Early Sample vs. Sample') 
        plt.ylabel('Cosine') 
        plt.ylim([-.25, 1.25]) 
        plt.yticks([0,.25,.5,.75,1]) 
        
    if options['obj']=='frac': 
        plt.ylabel('Frac. Selective') 
        plt.ylim([-.1, 0.4]) 
    
    if options['obj']=='proj': 
        plt.ylabel('S1/S2 memory axis') 
        plt.ylim([-.5, 6]) 
    
    if options['obj']=='score':
        plt.ylabel('Score')
        plt.ylim([0.4, 1.15]) 
        plt.yticks([0.5,.75,1.0]) 

    if options['obj']=='non_zero':  
        plt.ylabel('Non zero weights') 
        plt.ylim([-.1, 0.4]) 
        plt.yticks([0,.1,.2,.3,.4]) 

    # plt.savefig(figtitle) 
    pl.save_fig(figtitle, 1) 
    
if __name__ == '__main__':
    
    kwargs = dict() 
    
    kwargs['ci']=0 
    kwargs['shuffle']=0
    kwargs['perm_test']=1
    
    kwargs['n_samples'] = 1000
    kwargs['n_shuffles'] = 1000 
    kwargs['T_WINDOW'] = .5 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'robust' 
    kwargs['avg_mean_BL'] = 0 
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
    if kwargs['obj']=='norm':
        kwargs['epochs'] = ['ED','MD','LD','BL'] 
    else: 
        kwargs['epochs'] = ['ED','MD','LD'] 
    
    kwargs['epochs'] = ['ED','MD','LD'] # needed if comparing X_Si to X_Si BL 
    kwargs['off_diag'] = 1
    # kwargs['t_train'] = 'LD' 
    
    kwargs['add_vlines']=0 
    kwargs['IF_SAVE']=0 
    
    kwargs['n_days'] = 6 
    
    kwargs['clf_name'] = 'LogisticRegressionCV' 
    kwargs['clf_name'] = 'LDA' 
    kwargs['penalty'] = 'l2' 
    kwargs['n_lambda'] = 10 
    
    kwargs['out_fold'] = 'stratified' 
    kwargs['n_out'] = 10 
    kwargs['outer_score']= 'accuracy' 
    
    kwargs['in_fold'] = 'stratified' 
    kwargs['n_in'] = 10
    kwargs['inner_score']= 'accuracy' 
    
    kwargs['prescreen'] = True 
    kwargs['pval']= 0.05 
    kwargs['standardize'] = True  
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    # options['Delta0'] = get_delta_all_epochs(**options) 
    sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm = get_sel_epochs(**options) 
    # plot_sel_epochs(sel[1]-sel[0], sel_shuffle, pval_shuffle, sel_ci, pval_perm, **options) 
    # plot_sel_epochs(sel[1], sel_shuffle, pval_shuffle, sel_ci, pval_perm, **options) 
    
    
