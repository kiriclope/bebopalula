import os, sys, importlib 
from importlib import reload

import warnings
warnings.filterwarnings("ignore")

import numpy as np

import scipy.stats as stats 
from sklearn.preprocessing import StandardScaler 

from joblib import Parallel, delayed 
import multiprocessing

import utils.constants as gv 
reload(gv) 
from utils.options import *

from . import bootstrap 
reload(bootstrap) 
from .bootstrap import bootstrap, crossVal 

import utils.glms 
reload(utils.glms) 
from utils.glms import * 

''' 
We want to estimate the statistical significance of a model,
to do so, we bootstrap with replacement, fit the model 
on each bootstrap sample, and compute the quantity of interest (mean, score, corr, ...) for each boot.
'''

def getAlphaLambda(X, y, **options):
    
    clf_ = logitnetAlphaIterCV(lbd=options['lbd'], n_alpha=options['n_alpha'], n_lambda=options['n_lambda'], n_splits=options['inner_splits'],
                               standardize=True, fit_intercept=options['fit_intercept'], prescreen=options['prescreen'],
                               fold_type=options['fold_type'],
                               scoring=options['inner_scoring'], thresh=options['tol'] , maxit=options['max_iter'],
                               n_jobs=-1, verbose=0) 
    
    # clf = logitnetIterCV(lbd=options['lbd'], alpha=options['alpha'], n_lambda=options['n_lambda'], n_splits=options['n_splits'],
    #                           standardize=options['standardize'], fit_intercept=options['fit_intercept'], prescreen=options['prescreen'],
    #                           fold_type=options['fold_type'],
    #                           scoring=options['inner_scoring'], thresh=options['tol'] , maxit=options['max_iter'],
    #                           n_jobs=-1, verbose=options['verbose']) 
    
    # scaler = StandardScaler() 
    # X = scaler.fit_transform(X) 
    clf_.fit(X,y) 
    
    return clf_.alpha_, clf_.lbd_min_, clf_.non_zero_min 

def get_coefs(X_S1, X_S2, return_clf=0, **kwargs):
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    model = deepcopy(kwargs['clf'])
    
    try :
        model.named_steps['clf'].n_jobs = kwargs['n_jobs'] 
    except:
        model.n_jobs = kwargs['n_jobs']
    
    # if kwargs['t_train'] == 'LD':
    #     X_t_train = X_S1_S2[..., -1].copy()
    #     X_S1_S2[..., -1] = X_S1_S2[..., 0].copy() 
    #     X_S1_S2[..., 0] = X_t_train
    
    coefs = np.zeros( (X_S1.shape[-1], X_S1.shape[-2]) ) 
    
    models = []
    for i_epochs in range(X_S1.shape[-1]): 
        
        # model = deepcopy(kwargs['clf'])
        X = X_S1_S2[:,:,i_epochs].copy()
        
        try:
            model.fit(X,y) 
            # model = model.best_estimator_
            
            models.append(model) 
            # coefs.append(coef) # N_samples x N_neurons 
            
            try:
                pval = model.named_steps['filter'].pvalues_ 
                coef = model.named_steps['clf'].coef_ 
                idx = pval<=kwargs['pval'] 
                # print('non zero', np.sum(idx), np.sum(coef!=0) ) 
                coefs[i_epochs, idx] = coef 
            except:
                try:
                    coef = model.named_steps['clf'].coef_ 
                    coefs[i_epochs] = coef 
                except: 
                    coefs[i_epochs] = model.coef_ 
        
        except:
            pass
        
    if return_clf:
        return coefs, models 
    else:
        return coefs  

def get_non_zero(X_S1, X_S2, return_clf=0, **kwargs):
        
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    model = deepcopy(kwargs['clf'])
    try :
        model.named_steps['clf'].n_jobs = kwargs['n_jobs'] 
    except:
        model.n_jobs = kwargs['n_jobs']
    
    coefs = np.zeros( (X_S1.shape[-1], X_S1.shape[-2]) ) 
    
    models = []
    for i_epochs in range(X_S1.shape[-1]): 
        
        X = X_S1_S2[:,:,i_epochs]

        try:
            model.fit(X,y) 
            models.append(model)
            # coefs.append(coef) # N_samples x N_neurons 
        
            try:
                pval = model.named_steps['filter'].pvalues_ 
                coef = model.named_steps['clf'].coef_ 
                idx = pval<=kwargs['pval'] 
                # print('non zero', np.sum(idx), np.sum(coef!=0) ) 
                coefs[i_epochs, idx] = coef 
            except:
                try:
                    coef = model.named_steps['clf'].coef_ 
                    coefs[i_epochs] = coef 
                except: 
                    coefs[i_epochs] = model.coef_ 
        
        except:
            pass
        
    if return_clf:
        return coefs, models 
    else:
        return coefs  

def get_hyperparam(X_S1, X_S2, **kwargs): 
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    # get estimation of best (alpha, lambda) on the concatenated trials for the train epoch 
    if i_epochs==0 : 
        print('fix_alpha_lbd') 
        alpha, lbd, non_zero = getAlphaLambda(X, y, **kwargs) 
        print('alpha', alpha, 'lambda', lbd, 'non_zero', non_zero) 
        gv.clf.alpha = alpha 
        gv.clf.lbd = lbd 
    
def get_score(X_S1, X_S2, return_hyper=0, **kwargs):
    
    clf = get_clf(**kwargs) 
    # if kwargs['verbose'] :
    #     print('alpha', clf.alpha, 'lbd', clf.lbd) 
    # print('l1_ratio', clf.l1_ratio, 'C', clf.C) 
    
    # print('clf', clf) 
    model = crossVal(clf, method=kwargs['fold_type'], n_iter=kwargs['n_iter'], scoring=kwargs['scoring']
                     , scaler=kwargs['scaler'], n_jobs=None, verbose=kwargs['verbose']) 
    
    # model=clf 
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    score = [] 
    
    # X_train = np.mean( X_S1_S2[..., [24,25,26]], axis=-1)
    X_train = X_S1_S2[..., -1].copy() 
    
    for i_epochs in range(X_S1.shape[-1]): 
        # X_train = X_S1_S2[..., i_epochs].copy() 
        X_test = X_S1_S2[..., i_epochs].copy() 
        
        dum = model.get_scores(X_train, X_test, y) 
        score.append(dum) # N_samples x N_neurons 
    
    if return_hyper:
        return np.array(score), alpha, lbd 
    else:
        return np.array(score)        
    
def get_cv_score(X_S1, X_S2, return_hyper=0, **kwargs):
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    cv_scores = [] 
    lower = [] 
    upper = [] 
    
    random_state = int(np.random.rand()*1000) 
    
    model = deepcopy(kwargs['clf'])

    try :
        model.named_steps['clf'].n_jobs = kwargs['n_jobs'] 
        model.named_steps['clf'].random_state = random_state 
    except:
        model.n_jobs = kwargs['n_jobs']
        model.random_state = random_state
    
    if kwargs['off_diag'] :
        if X_S1.shape[-1]!=84 :
            if kwargs['t_train'] == 'ED':
                X_t_train = X_S1_S2[..., 0] 
            else:
                X_t_train = X_S1_S2[..., -1]
                X_S1_S2[..., -1] = X_S1_S2[..., 0].copy() 
                X_S1_S2[..., 0] = X_t_train.copy() 
        else:
            X_t_train = np.mean(X_S1_S2[..., kwargs['bins']], axis=-1) 
    
    # print('X_t_train', X_t_train.shape) 
    
    for i_epochs in range(X_S1_S2.shape[-1]): 
        # X = X_S1_S2[..., i_epochs] 
        if kwargs['off_diag']==False :
            X_t_train = X_S1_S2[..., i_epochs] 
        X_t_test = X_S1_S2[..., i_epochs] 
        
        # cv_score = nested_cv(kwargs['clf'], X, y, kwargs['param_grid'],
        #                      in_fold=kwargs['in_fold'], n_in=kwargs['n_in'], in_score=kwargs['inner_score'],
        #                      out_fold=kwargs['out_fold'], n_out=kwargs['n_out'],  out_score=kwargs['outer_score'],
        #                      fix_hyper=1, n_jobs=kwargs['n_jobs']) 
        
        # cv_score = outer_cv(model, X, y,
        #                     n_out = kwargs['n_out'], folds=kwargs['out_fold'],
        #                     inner_score=kwargs['inner_score'], outer_score = kwargs['outer_score'],
        #                     random_state = random_state, n_jobs=kwargs['n_jobs']) 
        
        cv_score = outer_temp_cv(model, X_t_train, X_t_test, y,
                                 n_out = kwargs['n_out'], folds=kwargs['out_fold'],
                                 inner_score=kwargs['inner_score'], outer_score = kwargs['outer_score'],
                                 random_state = random_state, n_jobs=kwargs['n_jobs']) 
        
        # print(grid.best_params_, kwargs['inner_score'], grid.best_score_, kwargs['outer_score'], cv_score) 
        
        cv_scores.append(cv_score) 
    
    return np.array(cv_scores) 
    
def get_boots_coefs(X_S1, X_S2, **kwargs):
            
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 
    
    boots_coefs = [] 
    
    random_state = int(np.random.rand()*1000) 
    kwargs['clf'].random_state = random_state
    
    for i_epochs in range(X_S1_S2.shape[-1]): 
        X = X_S1_S2[..., i_epochs] 

        kwargs['clf'].fit(X,y)
        boots_coefs.append(kwargs['clf'].coef_)
        # coefs = bootstrap_coefs(kwargs['clf'], X, y, n_jobs=-1) 
        # boots_coefs.append(coefs) 
    
    return np.array(boots_coefs) 

def get_proj_coefs(X_S1, X_S2, return_Delta=0, **kwargs): 
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) ) 
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) )) 

    model = deepcopy(kwargs['clf'])
    
    if return_Delta:
        coefs = np.zeros( X_S1.shape[-2] ) 
    
        X_S1_S2_bins = np.nanmean( X_S1_S2[..., kwargs['bins']], axis=-1) 
        X_S1_S2_avg = np.nanmean(X_S1_S2, axis=-1) 
        
        model.fit(X_S1_S2_avg, y) 
        pval = model.named_steps['filter'].pvalues_
        idx = pval<=kwargs['pval']
        coef = model.named_steps['clf'].coef_ 
        coefs[idx] = coef 
        scaler = model.named_steps['scaler']
    
    else:
        coefs = kwargs['Delta0'] 
        scaler = kwargs['fit_scaler'] 
    
    X_S1_rescaled = (X_S1 - scaler.mean_[np.newaxis, : , np.newaxis]) / scaler.var_[np.newaxis, : , np.newaxis]
    X_S2_rescaled = (X_S2 - scaler.mean_[np.newaxis, : , np.newaxis]) / scaler.var_[np.newaxis, : , np.newaxis]
    
    # X_proj = ( np.dot(coefs,  X_S1_rescaled ) +  np.dot(coefs,  X_S2_rescaled ) ) / 2.0 
    X_S1_proj = np.zeros( X_S1.shape[-1] ) 
    X_S2_proj = np.zeros( X_S2.shape[-1] ) 
    
    for i_epoch in range(X_S1.shape[-1]): 
        X_S1_proj[i_epoch] = np.mean( np.dot(coefs, X_S1_rescaled[..., i_epoch].T), axis=0 )  
        X_S2_proj[i_epoch] = np.mean( np.dot(coefs, X_S2_rescaled[..., i_epoch].T), axis=0 )  
    
    X_proj = (X_S1_proj - X_S2_proj) / 2.0 

    if return_Delta:
        return X_proj, coefs
    else:
        return X_proj
