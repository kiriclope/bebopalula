from importlib import reload
import sys
import gc 
import numpy as np 
import matplotlib.pyplot as plt

from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=False) 

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 
from sklearn.pipeline import Pipeline 
from sklearn.preprocessing import StandardScaler, RobustScaler 

# plt.ioff() 
plt.ion() 

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

import senc.utils
reload(senc.utils)
from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 

def get_delta_time(**options): 
    
    X_S1_tasks=[]
    X_S2_tasks=[]
    
    print(options['tasks'])
    
    X_S1_all, X_S2_all = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task='all', trials='correct') 
    
    # X_S1_all = np.vstack(X_S1_tasks) 
    # X_S2_all = np.vstack(X_S2_tasks) 
    
    X_S1_all, X_S2_all = pp.preprocess_X_S1_X_S2(X_S1_all, X_S2_all,
                                                 scaler=options['scaler_BL'],
                                                 center=options['center_BL'], scale=options['scale_BL'],
                                                 avg_mean=options['avg_mean_BL'], avg_noise=options['avg_noise_BL'],
                                                 unit_var=options['unit_var']
    ) 
    
    X_S1_all = pp.avg_epochs(X_S1_all, ['ED']) 
    X_S2_all = pp.avg_epochs(X_S2_all, ['ED']) 
    
    print('X_S1_all', X_S1_all.shape, 'X_S2_all', X_S2_all.shape)
    
    # X_S1_all, X_S2_all, center, scale = scale_data(X_S1, X_S2, scaler=options['scaler'], return_center_scale=1)
    
    # Delta_all = get_coding_direction(X_S1_all, X_S2_all, **options)[0]
    Delta_all = -get_coefs(X_S1_all, X_S2_all, **options)[0] 
    print('Delta_all', Delta_all.shape)  
    
    if options['return_center_scale']: 
        return Delta_all, center, scale
    else:
        return Delta_all
    
def sel_time(**options): 
    
    data.get_days() # do not delete that !!
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task=options['task'], trials=options['trials']) 
    # that must be before bins    
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=options['scaler_BL'],
                                         center=options['center_BL'], scale=options['scale_BL'],
                                         avg_mean=options['avg_mean_BL'], avg_noise=options['avg_noise_BL'], unit_var=options['unit_var']) 
    
    # epochs=['STIM','ED','DIST','MD','CUE','LD'] 
    # X_S1_all = pp.avg_epochs(X_S1_all, gv.epochs) 
    # X_S2_all = pp.avg_epochs(X_S2_all, gv.epochs) 
    
    if options['stimulus']=='sample': 
        if options['bins']=='ED': 
            options['bins'] = gv.bins_ED
            options['t_train'] = 'ED'
        if options['bins']=='MD': 
            options['bins'] = gv.bins_MD 
            options['t_train'] = 'MD'
        if options['bins']=='LD': 
            options['bins'] = gv.bins_LD 
            options['t_train'] = 'LD'
        if options['bins']=='Dist': 
            options['bins'] = gv.bins_DIST 
        if options['bins']=='Sample': 
            options['bins'] = gv.bins_STIM 
    else: 
        options['bins'] = gv.bins_DIST
    
    print('bins', options['bins']) 
    # options['bins'] = None 
    
    print('task', options['task'], 'X_S1', X_S1.shape, 'X_S2', X_S2.shape) 
    
    options['n_jobs'] = -10
    
    if options['obj']=='cos' or options['obj']=='proj' : 
        # options['Delta0'] = None # this fixes delta for the stats. 
        sel, Delta = get_sel(X_S1, X_S2, return_Delta=1, **options) 
        # options['Delta0'] = Delta # this fixes delta for the stats. 
    else:
        sel = get_sel(X_S1, X_S2, **options)
    
    options['n_jobs'] = None 
    
    print('sel', sel.shape, sel[:10]) 
    
    sel_ci = None 
    if options['ci']==1: 
        print('bootstraped ci') 
        _, sel_ci = my_bootstraped_ci(X_S1, X_S2, statfunction=lambda x, y: get_sel(x, y, **options) 
                                   , n_samples=options['n_samples'] ) 
        print('ci', sel_ci.shape) 
    
    sel_shuffle = None 
    if options['shuffle']==1:
        print('shuffle statistics') 
        sel_shuffle = shuffle_stat(X_S1, X_S2, lambda x,y: get_sel(x, y,**options), n_samples=options['n_shuffles']) 
        mean_shuffle = np.nanmean(sel_shuffle, axis=0) 
        perc_shuffle = np.nanpercentile(sel_shuffle, [2.5, 97.5], axis=0) 
        print('sel_shuffle', sel_shuffle.shape, 'mean', mean_shuffle.shape, 'perc', perc_shuffle.shape) 

    if options['obj']=='norm':
        # bins = np.arange(0, 12) 
        # BL = np.nanmean( sel[bins] ) 
        # sel -= BL 
        max_ = np.amax(sel)
        
        sel /= max_ 
        # sel /= BL
        
        if options['ci']: 
            # sel_ci -= BL 
            sel_ci /= max_ 
            # sel_ci /= BL 
        if options['shuffle']: 
            # sel_shuffle -= BL 
            sel_shuffle /= max_ 
            # sel_shuffle /= BL
    
    return sel, sel_ci, sel_shuffle
 
def plot_sel_time(sel, sel_ci=None, sel_shuffle=None, **options): 
    create_figdir(**options)
    
    if(len(options['task'])=='Dual'): 
        figtitle = '%s_time_%s_Dual_%s_day_%s_trials' % ( options['obj'], options['stimulus'], str(options['day']), options['trials'] ) 
    else: 
        figtitle = '%s_time_%s_day_%s_%s_trials' % ( options['obj'], options['stimulus'], str(options['day']), options['trials'] ) 

    if options['obj']=='score':
        if options['off_diag']:
            figtitle += '_'+ options['t_train'] +'_train'             
    
    fig = plt.figure(figtitle) 
    
    plt.plot(gv.time, sel, color=gv.pal[options['i_task']]) 
    if options['ci'] == 1: 
        plt.fill_between(gv.time, sel-sel_ci[:,0], sel+sel_ci[:,1], alpha=0.1) 
    if options['shuffle']==1: 
        mean_shuffle = np.nanmean(sel_shuffle, axis=0) 
        perc_shuffle = np.nanpercentile(sel_shuffle, [2.5, 97.5], axis=0) 
        
        plt.plot(gv.time, mean_shuffle, '--' , color=gv.pal[options['i_task']]) 
        plt.fill_between(gv.time, perc_shuffle[0], perc_shuffle[1], color=gv.pal[options['i_task']], alpha=.1) 
    
    if options['add_vlines']==1:    
        pl.add_vlines()
        if options['obj']=='proj':
            plt.text(2.5, 6.1, 'Sample', horizontalalignment='center', fontsize=10) 
            plt.text(5, 6.1, 'Dist.', horizontalalignment='center', fontsize=10) 
            plt.text(7, 6.1, 'Cue', horizontalalignment='center', fontsize=10) 
            plt.text(9.5, 6.1, 'Test', horizontalalignment='center', fontsize=10) 
        else:
            plt.text(2.5, 1.1, 'Sample', horizontalalignment='center', fontsize=10) 
            plt.text(5, 1.1, 'Dist.', horizontalalignment='center', fontsize=10) 
            plt.text(7, 1.1, 'Cue', horizontalalignment='center', fontsize=10) 
            plt.text(9.5, 1.1, 'Test', horizontalalignment='center', fontsize=10) 
            
    plt.xlabel('Time (s)') 
    plt.xticks([0, 2, 4, 6, 8, 10, 12, 14]) 
    plt.xlim([0, 14]) 
    
    if options['obj']=='norm':
        if options['stimulus']=='distractor': 
            plt.ylabel('Distractor Sel.') 
        else: 
            plt.ylabel('Sample Sel.')         
        plt.ylim([-.25, 1.1]) 
        plt.yticks([-0.25, 0, .25, .5, .75, 1.0]) 
        
    if options['obj']=='cos':
        plt.ylabel('Cosine') 
        # plt.ylabel('Overlap\n' r'Early Sample vs. Sample') 
        plt.ylim([-.25, 1.05]) 
        plt.yticks([-0.25, 0, .25, .5, .75, 1.0]) 
    
    if options['obj']=='frac':
        plt.ylabel('Frac. Selective') 
        plt.ylim([-.05, 0.3]) 
        plt.yticks([0, .1, .2, .3])
        
    if options['obj']=='proj': 
        
        plt.ylabel('Sample Memory Axis') 
        plt.ylim([-0.5, 6.0]) 
        plt.yticks([-0.5, 0.0, 2.0, 4.0, 6.0]) 
        
    if options['obj']=='score': 
        plt.ylabel('Score') 
        plt.ylim([0.25, 1.1]) 
        plt.yticks([.25, .5, 0.75, 1.0]) 
    
    if(options['IF_SAVE']==1):
        pl.save_fig(figtitle, 1)

    plt.close('all')
    plt.show()
    
if __name__ == '__main__':
    
    kwargs = dict() 
    
    kwargs['T_WINDOW'] = 0.5 
    
    kwargs['ci'] = 0 
    kwargs['n_samples'] = 100 
    kwargs['shuffle'] = 0 
    kwargs['n_shuffles'] = 100 
    
    kwargs['scaler'] = 'standard' #'standard' # if not standardized gives strange results for norm 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    kwargs['unit_var'] = 1 
    kwargs['n_days'] = 9
    kwargs['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo', 'Dual']) 
    # kwargs['tasks'] = ['DPA', 'Dual'] 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
        kwargs['obj'] = sys.argv[5] 
        kwargs['stimulus'] = sys.argv[6] 
        kwargs['bins'] = sys.argv[7] 

    kwargs['off_diag'] = 1
    if kwargs['bins'] == 'DIAG':
        kwargs['off_diag'] = 0
        
    # kwargs['bins'] = 'LD' 
    # kwargs['t_train'] = 'LD' 
    
    # kwargs['clf_name'] = 'logitnetCV' 
    # kwargs['lbd'] = 'lambda_min' 
    
    kwargs['clf_name'] = 'LogisticRegressionCV'
    # kwargs['clf_name'] = 'LDA'
    # kwargs['shrinkage'] = 1
    
    kwargs['penalty'] = 'l2'
    # kwargs['solver'] = 'saga'
    kwargs['n_lambda'] = 20
    
    kwargs['in_fold'] = 'stratified'
    kwargs['n_in'] = 10
    kwargs['inner_score'] = 'neg_log_loss'
    
    kwargs['out_fold'] = 'stratified'
    kwargs['n_out'] = 10
    kwargs['outer_score'] = 'roc_auc' 
    
    kwargs['prescreen'] = True 
    kwargs['pval'] = 0.05 
    kwargs['standardize'] = True 
        
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    print(options['task']) 
    
    if options['task']=='all': 
        options['add_vlines']=0 
        options['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo']) 
        
        if options['obj']=='proj':
            options['Delta0'] = get_delta_time(**options) 
        
        for options['i_task'] in range(len(options['tasks'])): 
            options['task'] = options['tasks'][options['i_task']] 
            
            if options['i_task']==len(options['tasks'])-1: 
                options['add_vlines']=1 
                options['IF_SAVE']=1 
            
            sel, sel_ci, sel_shuffle = sel_time(**options) 
            plot_sel_time(sel, sel_ci, sel_shuffle, **options) 

        plt.close('all')
    else: 
        # options['tasks'] = np.array(['DualGo', 'DualNoGo']) 
        # options['tasks'] = np.array(['all']) 
        # options['stimulus'] = 'sample' 
        # options['return_center_scale'] = 0 
        # options['Delta0'] = get_delta_time(**options) 
        
        options['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo', 'Dual']) 
        # options['tasks'] = np.array(['DPA','Dual']) 
        # options['bins'] = 'ED' 
        options['i_task'] = np.argwhere(options['tasks']==options['task'])[0][0] 
        print(options['tasks'], options['task'], options['i_task']) 
        
        options['add_vlines'] = 1 
        sel, sel_ci, sel_shuffle = sel_time(**options) 
        plot_sel_time(sel, sel_ci, sel_shuffle, **options) 
        
