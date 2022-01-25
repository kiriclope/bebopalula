from importlib import reload
import inspect, sys, os 
import numpy as np 
import matplotlib.pyplot as plt

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

def cv_score(**kwargs):
        
    data.get_days() # do not delete that !!
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'], stimulus='sample', task=options['task'], trials=options['trials'])
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=options['scaler_BL'],
                                         center=options['center'], scale=options['scale'],
                                         avg_mean=options['avg_mean'], avg_noise=options['avg_noise']) 
    
    X_S1 = pp.avg_epochs(X_S1.copy(), epochs=['ED', 'LD']) 
    X_S2 = pp.avg_epochs(X_S2.copy(), epochs=['ED', 'LD']) 
    
    print('X_S1', X_S1.shape,'X_S2', X_S2.shape) 
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 

    clf = get_clf(**kwargs) 
    Cs = np.logspace(-2, 2, 10)

    kwargs['penalty'] = 'l1'
    
    if kwargs['penalty']== 'elasticnet':
        kwargs['solver'] = 'saga'
        alphas = np.logspace(-2, 2, 10) 
        param_grid = dict(clf__C=Cs, clf__l1_ratio=alphas) # this is important if using pipeline 
    else: 
        kwargs['solver'] = 'liblinear' 
        param_grid = dict(clf__C=Cs)
    
    scores = [] 
    for i_epochs in range(X_S1.shape[-1]): 
        X = X_S1_S2[:,:,i_epochs] 
        
        # single cv loop
        # score = outer_cv(clf, X, y) 
        # print('<cv_score>', score) 
    
        # nested cv with hyperparam tuning 
        scores.append( nested_cv(clf, X, y, param_grid) ) 
        # score = inner_cv(clf, X, y, param_grid) 
        
    print('<nested_cv_scores>', scores) 
    
if __name__ == '__main__':
    
    kwargs = dict() 
    kwargs['pval']= 0.05  
    kwargs['n_samples']= 1000
    
    kwargs['clf'] = 'LogisticRegression'
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3]
        kwargs['trials'] = sys.argv[4] 
    
    options = set_options(**kwargs) 
    set_globals(**options)    
    
    cv_score(**options) 
    
