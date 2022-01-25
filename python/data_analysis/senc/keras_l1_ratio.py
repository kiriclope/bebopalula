from importlib import reload
import inspect, sys, os 
import numpy as np 
import matplotlib.pyplot as plt

from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=0)

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses
    
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

# from utils.dum import *

from utils.plot_settings import SetPlotParams
SetPlotParams()

from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 

def cv_score_l1_ratio(**kwargs):
    
    data.get_days() # do not delete that !! 
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=kwargs['day'], stimulus='sample', task=kwargs['task'], trials=kwargs['trials']) 
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=kwargs['scaler_BL'],
                                         center=kwargs['center_BL'], scale=kwargs['scale_BL'],
                                         avg_mean=kwargs['avg_mean_BL'], avg_noise=kwargs['avg_noise_BL']) 
    
    X_S1 = pp.avg_epochs(X_S1.copy(), epochs=['ED', 'MD', 'LD']) 
    X_S2 = pp.avg_epochs(X_S2.copy(), epochs=['ED', 'MD', 'LD']) 
    
    print('X_S1', X_S1.shape,'X_S2', X_S2.shape) 
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    alphas = np.linspace(0, 1, kwargs['n_alpha']) 
    
    # def func(alpha=0): 
    #     return get_model(alpha, input_dim=X_S1_S2.shape[-2]) 
    
    # clf = KerasClassifier(build_fn=func, epochs=10, batch_size=10, verbose=0) 
    
    clf = kwargs['clf']  

    score_mat = []     
    for i_epochs in range(X_S1.shape[-1]): 
        
        scores = [] 
        X_t_train = X_S1_S2[...,i_epochs] 
        
        for j_epochs in range(X_S1.shape[-1]):
            
            X_t_test = X_S1_S2[...,j_epochs] 
                
            score, grid = temp_nested_cv(clf, X_t_train, X_t_test, y, kwargs['param_grid'], scaling='standard', 
                                         in_fold=kwargs['in_fold'], n_in=kwargs['n_in'], in_score=kwargs['inner_score'], 
                                         out_fold=kwargs['out_fold'], n_out=kwargs['n_out'], out_score=kwargs['outer_score'], 
                                         fix_hyper=0, n_jobs=-1) 
            
            print(grid.best_params_, kwargs['inner_score'], grid.best_score_, kwargs['outer_score'], score) 

            scores.append(score) 
        score_mat.append(scores) 
    
    print(score_mat) 
    fig = plt.figure()
    ax = plt.gca()
    
    img = ax.imshow(score_mat, origin='lower', cmap='jet', vmin=0.5, vmax=1) 
    ax.grid(False) 
    cbar = fig.colorbar(img, ticks=[.5, .6, .7 , .8, .9 , 1]) 
    # plt.clim([0.5,1]) 
    cbar.set_label('AUC') 
    
    plt.xlabel('Test Delay') 
    plt.ylabel('Train Delay') 
    ax.set_yticks([0,1,2]) 
    ax.set_xticks([0,1,2]) 
    ax.set_yticklabels(['Early', 'Middle', 'Late']) 
    ax.set_xticklabels(['Early', 'Middle', 'Late']) 
    # plt.plot(time, scores) 
    # pl.add_vlines() 
        
if __name__ == '__main__': 
    
    kwargs = dict() 
    
    kwargs['ci'] = 0 
    kwargs['n_samples'] = 1000 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    
    kwargs['clf_name'] = 'logitnet' 
    kwargs['T_WINDOW'] = 0.5 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
    
    kwargs['n_alpha'] = 10 
    kwargs['n_out'] = 10 
    kwargs['n_in'] = 10 
    kwargs['out_fold'] = 'stratified' 
    kwargs['in_fold'] = 'stratified' 
    kwargs['n_lambda'] = 10 
    
    kwargs['lbd'] = 'lambda_1se' 
    kwargs['inner_score']= 'neg_mean_squared_error' 
    kwargs['outer_score']= 'roc_auc' 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    cv_score_l1_ratio(**options)     
    
