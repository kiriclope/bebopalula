#!/usr/bin/env python3
import inspect, sys 
from importlib import reload 
sys.path.insert(0, '../')

import time 
import numpy as np 
import matplotlib.pyplot as plt

from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=0)

plt.ioff() 

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
    
    return cols, high, low, corr

def add_pval_shuffle(pval_shuffle, obj='cos'):
    cols, high, low, corr = cols_high(obj)
    
    for i_task in range(3): 
        for i_epoch in range(2): 
            
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

    for i_task in range(2): 
        for i_epoch in range(2): 
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
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task='all', trials=options['trials']) 
    
    # compute Delta ED for pulled tasks
    X_S1_all = np.vstack(X_S1) 
    X_S2_all = np.vstack(X_S2) 
    
    X_S1_S2_all = np.vstack( (X_S1_all, X_S2_all) ) 
    
    X_BL = X_S1_S2_all[..., gv.bins_BL] 
    
    X_S1_S2_all = np.vstack( (X_S1_all, X_S2_all) ) 
    X_all_scale, center, scale = pp.standard_scaler_BL(X_S1_S2_all, avg_mean=1, avg_noise=1) 
    
    X_S1_all = X_S1_S2_all[:X_S1_all.shape[0]] 
    X_S2_all = X_S1_S2_all[X_S1_all.shape[0]:] 
    
    X_S1_all = pp.avg_epochs(X_S1_all, gv.epochs) 
    X_S2_all = pp.avg_epochs(X_S2_all, gv.epochs) 
    
    print('X_S1_S2_all', X_S1_S2_all.shape, 'X_S1_all', X_S1_all.shape, 'X_S2_all', X_S2_all.shape) 
    
    Delta_all = get_coding_direction(X_S1_all, X_S2_all, **options) 
    print('Delta0', Delta_all.shape) 

    return Delta_all

def get_sel_epochs(**options): 

    X_S1_tasks = []
    X_S2_tasks = []
    
    data.get_days() # do not delete that !!
    
    for i_task in range(len(options['tasks'])):
        options['task'] = options['tasks'][i_task] 
        X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task=options['task'], trials=options['trials']) 
        X_S1_tasks.append(X_S1) 
        X_S2_tasks.append(X_S2) 
        
        X_S1_tasks[i_task], X_S2_tasks[i_task] = pp.preprocess_X_S1_X_S2(X_S1_tasks[i_task], X_S2_tasks[i_task],
                                                                         scaler=options['scaler_BL'],
                                                                         center=options['center_BL'], scale=options['scale_BL'],
                                                                         avg_mean=options['avg_mean_BL'], avg_noise=options['avg_noise_BL']
                                                                         , unit_var=options['unit_var']) 
        
        X_S1_tasks[i_task] = pp.avg_epochs(X_S1_tasks[i_task], options['epochs']) 
        X_S2_tasks[i_task] = pp.avg_epochs(X_S2_tasks[i_task], options['epochs']) 
    
    sel = [] 
    sel_ci = [] 
    sel_shuffle = [] 
    pval_shuffle = [] 
    
    options['bins'] = None 
    
    Delta_tasks = []  
    alpha_tasks = [] 
    lbd_tasks = []
    tuned_clfs = []
    
    for i_task in range(len(options['tasks'])): 
        options['task'] = options['tasks'][i_task] 
        print('task', options['task'], 'X_S1', X_S1_tasks[i_task].shape, 'X_S2', X_S2_tasks[i_task].shape)
        
        # sample selectivity 
        if options['obj']=='cos': 
            sel_task, Delta = get_sel( X_S1_tasks[i_task], X_S2_tasks[i_task], return_Delta=1, **options) 
            options['Delta0'] = Delta # fixes Delta for stats. 
            Delta_tasks.append(Delta) 
        # elif options['obj'] == 'score':
        #     sel_task, alpha, lbd = get_sel( X_S1_tasks[i_task], X_S2_tasks[i_task], return_hyper=1, **options) 
        #     alpha_tasks.append(alpha)
        #     lbd_tasks.append(lbd)
        #     print('alpha', alpha, 'lbd', lbd)
        #     options['alpha']= alpha 
        #     options['lbd']= lbd[0] 
        #     # options['l1_ratio'] = alpha 
        #     # options['C'] = 1.0/lbd[0]
        elif options['obj']=='score':
            options['IF_TUNE']=1 
            options['IF_CI']=1 
            # tuned_clfs.append(tuned_clf) 
            startbuild = time.time() 
            sel_task = get_sel( X_S1_tasks[i_task], X_S2_tasks[i_task], **options) 
            endbuild = time.time() 
            build_time = endbuild - startbuild 
            print("Building time: %.2f s" % build_time) 
            # sel_task, lower, upper = get_sel( X_S1_tasks[i_task], X_S2_tasks[i_task], **options) 
            # sel_ci.append(np.vstack((lower, upper)).T) 
        else: 
            sel_task = get_sel( X_S1_tasks[i_task], X_S2_tasks[i_task], **options) 

        print('sel', sel_task) 
        sel.append(sel_task) 
        
        # bootstrapped confidence interval 
        if options['ci'] :
            # options['IF_TUNE']=0
            # options['clf'] = tuned_clf 
            sel_ci.append( my_bootstraped_ci(X_S1_tasks[i_task], X_S2_tasks[i_task], statfunction=lambda x,y: get_sel(x, y, **options),
                                             n_samples=options['n_samples']) ) 
        
        if options['shuffle']: 
            # options['IF_TUNE']=0 
            # options['clf'] = tuned_clf 
            # shuffle statistics
            sel_shuffle.append(shuffle_stat(X_S1_tasks[i_task], X_S2_tasks[i_task], lambda x,y: get_sel(x, y, **options),
                                            n_samples=options['n_shuffles'] ).T )
            
            # # p value with respect to shuffle 
            pval_shuffle.append( 2.0*np.amin( np.stack( [np.mean( sel_shuffle[i_task] >= sel[i_task][..., np.newaxis], axis=-1 ), 
                                                         np.mean( sel_shuffle[i_task] <= sel[i_task][..., np.newaxis], axis=-1 ) ] ) 
                                              , axis=0 ) 
            )
    
    sel = np.array(sel) 
    print('sel', sel) 
    
    if options['ci']:
        sel_ci = np.array(sel_ci)
        # sel_ci[..., 0] = np.abs(sel - sel_ci[..., 0]) 
        # sel_ci[..., 1] = np.abs(-sel + sel_ci[..., 1]) 
        
        print('sel_ci', sel_ci.shape) 
        print(sel_ci)
        
    if options['shuffle']:
        sel_shuffle = np.array(sel_shuffle) 
        pval_shuffle = np.array(pval_shuffle)[:,1:] 
        print('shuffle', sel_shuffle.shape) 
        print('pval', pval_shuffle.shape, pval_shuffle) 
    
    # permutation test 
    pval_perm = [] 
    if options['obj']=='cos': 
        options['Delta0'] = Delta_tasks[0]    
    options2 = options 
    
    if options['perm_test']:
        for i_task in range(1, len(kwargs['tasks'])): # Dual Go and Dual No Go 
            
            if options['obj']=='cos': 
                options2['Delta0'] = Delta_tasks[i_task]
            
            sel_perm_DPA, sel_perm_Other = get_sel_perm( X_S1_tasks[0], X_S2_tasks[0], X_S1_tasks[i_task], X_S2_tasks[i_task],
                                                         lambda x,y: get_sel(x, y, **options), lambda x,y: get_sel(x, y, **options2),
                                                         n_samples=options['n_samples'] ) 
                        
            print('sel_perm', sel_perm_DPA.shape, sel_perm_Other.shape) 
            
            # compare Delta_perm = DPA_perm - Other_perm vs Delta = DPA - Other 
            pval_perm.append( 2.0*np.amin( np.stack( [np.mean( sel_perm_DPA - sel_perm_Other >= sel[0] - sel[i_task], axis=0 ),
                                                      np.mean( sel_perm_DPA - sel_perm_Other <= sel[0] - sel[i_task], axis=0 )] ) 
                                           , axis=0 ) 
            ) 
        
        pval_perm = np.array(pval_perm) # get rid of Early delay 
        print('pval_perm', pval_perm.shape, pval_perm) 
        pval_perm = pval_perm[:,1:] # get rid of Early delay 
    
    return sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm 

def plot_sel_epochs(sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm, **options):
    create_figdir(**options) 
    
    figtitle = '%s_%s_epochs_%s_day_%s' % (options['mouse_name'][options['i_mice']], options['obj'], options['trials'], str(options['day'])) 
    fig = plt.figure(figtitle) 
    cols = [-.1,0,.1]
    
    for i_task, task in enumerate(options['tasks']):         
        if options['obj']=='norm':
            
            BL = sel[i_task, -1] 
            ED = sel[i_task, 0] 
            
            sel[i_task] /= ED
            
            if options['ci']:
                sel_ci[i_task] /= ED 

        # convert percentile to errorbars 
        # sel_ci[i_task, :, 0] = sel[i_task] - sel_ci[i_task, :, 0] 
        # sel_ci[i_task, :, 1] = - sel[i_task] + sel_ci[i_task, :, 1] 
        
        plt.plot([cols[i_task], .4 + cols[i_task]] , sel[i_task,1:3], 'o', color=gv.pal[i_task], ms=2) 
        
        if options['ci']:
            plt.errorbar([cols[i_task], .4 + cols[i_task]], sel[i_task, 1:3], yerr=sel_ci[i_task, :, 1:3], 
                         ecolor=gv.pal[i_task], color=gv.pal[i_task], ls='none') 
    
    if options['shuffle']:
        add_pval_shuffle(pval_shuffle, obj=options['obj'])
    
    if options['perm_test']:
        add_pval_perm(pval_perm, obj=options['obj'])
    
    plt.xticks([0,.4], ['Early vs.\n' r"Middle delay", 'Early vs.\n' r"Late delay"]) 
    plt.xlim([-0.25, .65]) 
    
    if options['obj']=='norm': 
        plt.ylabel('Sample Sel.') 
        plt.ylim([-.25, 1.25]) 
        plt.yticks([0,.25,.5,.75,1]) 
        # plt.ylim([0, 2.5]) 
        # plt.yticks([0,2,4,6,8,10]) 
        # plt.yticks([0,.25,.5,.75,1])
        # plt.yticks([-1,-.5,0,.5,1,1.5,2, 2.5]) 
    
    if options['obj']=='cos':  
        plt.ylabel('Overlap\n' r'Early Sample vs. Sample') 
        # plt.ylabel('Cosine') 
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
        
    # plt.savefig(figtitle) 
    pl.save_fig(figtitle) 
    
if __name__ == '__main__':

    kwargs = dict()
    
    kwargs['ci']=1
    kwargs['shuffle']=1
    kwargs['perm_test']=1
    
    kwargs['pval']= 0.05 
    
    kwargs['n_samples'] = 1000
    kwargs['n_shuffles'] = 1000
    kwargs['T_WINDOW'] = 0.5 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    
    # kwargs['clf_name']='LogisticRegression' 
    # kwargs['penalty'] = 'l1' 
    # kwargs['solver'] = 'liblinear' 
    
    kwargs['clf_name'] = 'logitnetCV' 
    kwargs['alpha'] = 1.0 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
        kwargs['obj'] = sys.argv[5] 
        kwargs['stimulus'] = sys.argv[6] 
    
    kwargs['tasks'] = ['DPA', 'DualGo', 'DualNoGo'] 
    if kwargs['obj']=='norm':
        kwargs['epochs'] = ['ED','MD','LD','BL'] 
    else:
        kwargs['epochs'] = ['ED','MD','LD']
    
    kwargs['epochs'] = ['ED','MD','LD'] # needed if comparing X_Si to X_Si BL 
    
    kwargs['add_vlines']=0 
    kwargs['IF_SAVE']=0 
    
    if kwargs['obj']=='cos':
        kwargs['n_iter'] = 1 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    if options['obj'] == 'score' :
        options['clf'] = get_clf(**options) 
        print('clf', options['clf']) 
    
    # options['IF_TUNE']=1 
    # # options['Delta0'] = get_delta_all_epochs(**options) 
    sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm = get_sel_epochs(**options) 
    plot_sel_epochs(sel, sel_shuffle, pval_shuffle, sel_ci, pval_perm, **options) 

