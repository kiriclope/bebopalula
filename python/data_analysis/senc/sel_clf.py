import warnings
warnings.filterwarnings("ignore")

from importlib import reload
import inspect, sys
sys.path.insert(0, '../')

from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=0)

import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from joblib import Parallel, delayed 

# from keras.wrappers.scikit_learn import KerasClassifier
# from utils.keras_wrapper import build_keras_model

import utils.progressbar as pg 
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

from utils.plot_settings import SetPlotParams
SetPlotParams()

from utils.glmnet_wrapper import logitnet

def get_logitnet(alpha=1, n_lambda=20, lbd=1, scoring='accuracy',
                 prescreen=False, standardize=True, fit_intercept=True, thresh=1e-4 , maxit=1e6):
    
    return logitnet(alpha, n_lambda, lbd, scoring,
                    prescreen, standardize, fit_intercept, thresh , maxit)


def cv_score(**kwargs):
        
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
    
    # func = lambda alpha: build_keras_model(alpha, input_dim=X_S1_S2.shape[1]) 
    
    # if kwargs['clf'] =='keras_logreg':
    #     clf = KerasClassifier(build_fn=func, epochs=1, batch_size=1, verbose=0) 
    # else:
    clf = get_clf(**kwargs)
        
    print('clf', clf)
    
    # alphas = np.linspace(0, 1, 10)    
    # param_grid = dict(clf__alpha=alphas) # this is important if using pipeline
    
    Cs = np.logspace(-2, 2, 100) 
    param_grid = dict(clf__C=Cs) # this is important if using pipeline 
    
    # param_grid = dict(clf__C=Cs, clf__l1_ratio=alphas) # this is important if using pipeline 
    # param_grid = dict(clf__lbd=Cs) # this is important if using pipeline 
    # param_grid = dict(clf__lbd=Cs, clf__alpha=alphas) # this is important if using pipeline 
    
    # X_S1_S2 = scipy.signal.resample(X_S1_S2, axis=-1, num=int(X_S1_S2.shape[-1]/2) ) 
    # time = np.linspace(0, 14, X_S1_S2.shape[-1]) 
    
    score_mat = [] 
    for i_epochs in range(X_S1_S2.shape[-1]): 
        X_t_train = X_S1_S2[..., i_epochs] 
        scores = []
        for j_epochs in range(X_S1_S2.shape[-1]): 
            
            X_t_test = X_S1_S2[..., j_epochs] 
            
            # # nested cv with hyperparam tuning 
            # score = nested_cv(clf, X, y, param_grid, n_jobs=-1)
            
            score = temp_nested_cv(clf, X_t_train, X_t_test, y, param_grid) 
            print('i_epoch', i_epochs, 'j_epoch', j_epochs, '<nested_cv_score>', score) 
            
            scores.append(score) 
        score_mat.append(scores) 
    
    score_mat = np.array(score_mat) 
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
    kwargs['pval']= 0.05  
    kwargs['n_samples']= 1000

    kwargs['clf'] = 'LogisticRegression' 
    kwargs['penalty'] = 'l1' 
    kwargs['solver'] = 'liblinear' 

    kwargs['T_WINDOW'] = 0.5 
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    
    # kwargs['clf'] = 'keras_logreg'
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3]
        kwargs['trials'] = sys.argv[4] 
    
    options = set_options(**kwargs) 
    set_globals(**options)    
    
    cv_score(**options) 
        
    
