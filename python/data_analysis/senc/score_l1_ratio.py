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
        
    for i_epochs in range(X_S1.shape[-1]): 
        
        scores = [] 
        conf_int = [] 
        
        X = X_S1_S2[:,:,i_epochs] 
        random_state = int(np.random.rand()*1000) 
        
        clf.random_state = random_state 
        
        for i_alpha in range(len(alphas)): 
            
            # plt.figure() 
            clf.alpha = alphas[i_alpha] 
            # clf.alpha = 1 
            clf.n_jobs = -1 
            
            # clf.fit(X, y) 
            # clf.lasso_path() 
            # model_ = clf.model_ 
            # idx_lbd_min = model_['lambdau'] == model_['lambda_min'] 
            
            # lbd_min = model_['lambdau'][idx_lbd_min] 
            # score = model_['cvm'][idx_lbd_min] 
            
            # print('epoch', i_epochs, 'alpha', alphas[i_alpha], 'score', score, 'lbd', lbd_min) 
            
            # scores.append(score) 

            # grid = inner_cv(clf, X, y, kwargs['param_grid'], scaling='standard', n_in=kwargs['n_in'], n_jobs=-1) 
            
            score = outer_cv(clf, X, y,
                             n_out = kwargs['n_out'], folds=kwargs['fold_type'],
                             inner_score=kwargs['inner_score'], outer_score = kwargs['outer_score'],
                             random_state = random_state) 
            
            scores.append( score ) 
            print('epoch', i_epochs, 'alpha', alphas[i_alpha], 'score', score) 
            
            #     if kwargs['ci']:
            #         clf.n_jobs = None
            #         ci = my_bootstraped_ci(X_S1[..., i_epochs], X_S2[..., i_epochs],
            #                                statfunction=lambda x1,x2: outer_cv(clf, np.vstack((x1,x2)), y,
            #                                                                    n_out = kwargs['n_out'],
            #                                                                    folds=kwargs['fold_type'],
            #                                                                    inner_score=kwargs['inner_score'],
            #                                                                    outer_score=kwargs['outer_score'],
            #                                                                    random_state = random_state),
            #                                n_samples=options['n_samples']) 
            
            #         # shuffle.append(shuffle_stat(X_S1, X_S2, lambda x1,x2: outer_cv(clf, np.vstack((x1,x2)), y,
            #         #                                                                folds=kwargs['fold_type']),
            #         #                             n_samples=options['n_shuffles'] ).T )
            
            
            #         ci_list = [item for sublist in ci for item in sublist] 
            
            #         conf_int.append( ci_list ) 
            
            #         print('epoch', i_epochs, 'alpha', alphas[i_alpha], 'score', score, 'ci', ci_list) 
        
        plt.figure('score') 
        plt.plot(alphas, scores) 
        
        max_idx = np.argmax(scores) 
        max_alpha = alphas[max_idx] 
        plt.axvline(max_alpha, ls='dashed')
        
        if kwargs['ci']:
            conf_int = np.array(conf_int) 
            print(conf_int.shape) 
            plt.fill_between(alphas, scores-conf_int[:,0], scores+conf_int[:,1], alpha=0.25)         
        
        plt.xlabel('l1 ratio')
        plt.ylabel('Score') 
        
if __name__ == '__main__':
    
    kwargs = dict()
    
    kwargs['ci'] = 0
    kwargs['n_samples'] = 1000 
    
    kwargs['scaler'] = 'standard' 
    kwargs['scaler_BL'] = 'standard' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1 
    
    kwargs['clf_name'] = 'logitnetCV' 
    kwargs['T_WINDOW'] = 0.5 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3]
        kwargs['trials'] = sys.argv[4] 
    
    kwargs['n_alpha'] = 20 
    kwargs['n_out'] = 40 
    kwargs['n_in'] = 10 
    kwargs['n_lambda'] = 20 
    
    kwargs['lbd'] = 'lambda_min' 
    kwargs['fold_type'] = 'stratified' 
    kwargs['inner_score']= 'accuracy' 
    kwargs['outer_score']= 'accuracy' 
    
    options = set_options(**kwargs) 
    set_globals(**options)    
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    cv_score_l1_ratio(**options)     
    
