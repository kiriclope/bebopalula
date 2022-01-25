import warnings
warnings.filterwarnings('ignore')

from importlib import reload
import inspect, sys
sys.path.insert(0, '../')

import time 
from sklearnex import patch_sklearn
patch_sklearn(global_patch=True, verbose=0)

import numpy as np 
import matplotlib.pyplot as plt

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

from utils.dum import *

def modulation_dist(**kwargs):
        
    data.get_days() # do not delete that !!
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=kwargs['day'], stimulus='sample', task=kwargs['task'], trials=kwargs['trials'])
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=kwargs['scaler_BL'],
                                         center=kwargs['center'], scale=kwargs['scale'],
                                         avg_mean=kwargs['avg_mean'], avg_noise=kwargs['avg_noise']) 
    
    # X_S1 = pp.avg_epochs(X_S1.copy(), epochs=['ED', 'LD'])[...,-1] 
    # X_S2 = pp.avg_epochs(X_S2.copy(), epochs=['ED', 'LD'])[...,-1] 

    X_S1 = pp.avg_epochs(X_S1.copy(), epochs=['ED', 'LD'])
    X_S2 = pp.avg_epochs(X_S2.copy(), epochs=['ED', 'LD']) 
    
    print('X_S1', X_S1.shape,'X_S2', X_S2.shape) 
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    # clf = get_clf(**kwargs) 
    # Cs = np.logspace(-2, 2, 10) 
    # alphas = np.logspace(-2, 2, 10) 
    # param_grid = dict(clf__C=Cs, clf__l1_ratio=alphas) # this is important if using pipeline 
    
    def func(alpha=0): 
        return get_model(alpha, input_dim=X_S1_S2.shape[-2]) 
    
    clf = KerasClassifier(build_fn=func, epochs=1, batch_size=1, verbose=0) 
    # clf.fit(X, y) 
    
    # single cv loop 
    # score = outer_cv(clf, X, y) 
    # print('<cv_score>', score) 
    
    # nested cv with hyperparam tuning 
    alphas = np.linspace(0, 1, 10) 
    param_grid = dict(clf__alpha=alphas) # this is important if using pipeline 
    # score = nested_cv(clf, X[...,0], y, param_grid) 
    # # score = temp_nested_cv(clf, X[...,0], X[...,-1], y, param_grid) 
    # # score = inner_cv(clf, X, y, param_grid) 
    
    # print('<nested_cv_score>', score)
    
    score_mat = [] 
    for i_epochs in range(X_S1_S2.shape[-1]): 
        X_t_train = X_S1_S2[..., i_epochs] 
        
        # score = nested_cv(clf, X_t_test, y, param_grid, n_jobs=-1)
        
        scores = []
        for j_epochs in range(X_S1_S2.shape[-1]): 
            
            X_t_test = X_S1_S2[..., j_epochs] 
            
            startbuild = time.time() 
            score = temp_nested_cv(clf, X_t_train, X_t_test, y, param_grid, n_jobs=-1) 
            endbuild = time.time() 
            build_time = endbuild - startbuild 
            print("Building time: %.2f s" % build_time) 
            print('i_epoch', i_epochs, 'j_epoch', j_epochs, '<nested_cv_score>', score) 
            
            scores.append(score) 
        score_mat.append(scores) 
    
    score_mat = np.array(score_mat) 
    print(score_mat) 
    plt.figure() 
    plt.imshow(score_mat, origin='lower') 
    plt.grid(False) 
    plt.colorbar() 
    plt.clim([0.5,1]) 
    
if __name__ == '__main__':
    
    kwargs = dict() 
    kwargs['pval']= 0.05  
    kwargs['n_samples']= 1000
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3]
        kwargs['trials'] = sys.argv[4] 

    options = set_options(**kwargs) 
    set_globals(**options)    

    modulation_dist(**options) 
    
