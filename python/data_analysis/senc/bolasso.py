from importlib import reload
import inspect, sys, os 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats

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

from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 

from utils.plot_settings import SetPlotParams
SetPlotParams()

def cv_score_l1_ratio(**kwargs):

    clf = kwargs['clf']
    
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

    scores = [] 
    conf_int = [] 
        
    random_state = int(np.random.rand()*1000) 
    clf.random_state = random_state # same random_state for each values of alpha

    coefs = []
    
    for i_epochs in range(X_S1.shape[-1]): 
            
        X = X_S1_S2[:,:,i_epochs] 
        
        clf.alpha = 0.5 
        boots_coefs = bootstrap_coefs(clf, X, y, n_boots=1000, n_jobs=-1) 
        mean_coefs = np.mean(boots_coefs, axis=0) 
        print(boots_coefs.shape, mean_coefs[:10]) 
    
        _, p_val = stats.ttest_1samp(boots_coefs, 0, axis=0) 
    
        print(p_val.shape, p_val[:10]) 
        print(np.sum(p_val>=0.05)) 
        
        X_fs = X[:, p_val>=0.05] 
        print(X_fs.shape) 
        
        clf.alpha = 0 
        clf.fit(X_fs, y) 

        coefs.append(clf.coef_)
        
        print('lbd min', clf.lbd_min_, 'non zeros', clf.non_zero_min_) 

        
if __name__ == '__main__':
    
    kwargs = dict()
    
    kwargs['ci'] = 0
    kwargs['n_samples'] = 1000 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    
    kwargs['T_WINDOW'] = 0.5 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
    
    kwargs['n_alpha'] = 20 
    kwargs['n_lambda'] = 40 
    kwargs['lbd'] = 'lambda_min' 
    
    kwargs['out_fold'] = 'loo' 
    kwargs['n_out'] = 10 
    kwargs['outer_score']= 'accuracy' 
    
    kwargs['in_fold'] = 'stratified' 
    kwargs['n_in'] = 10 
    kwargs['inner_score']= 'neg_log_loss' 
    
    kwargs['clf_name'] = 'logitnetCV' 
    
    # kwargs['lbds'] = np.exp(np.linspace(-4, 4, kwargs['n_lambda']) ) 
    kwargs['lbds'] = np.logspace(-4, 2, kwargs['n_lambda'])    
    kwargs['param_grid'] = dict(lbd=kwargs['lbds']) 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    cv_score_l1_ratio(**options) 
    
