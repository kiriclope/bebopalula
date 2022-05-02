import numpy as np
from copy import deepcopy

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler

from sklearn.svm import LinearSVC 
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import GridSearchCV, cross_val_score, cross_validate, cross_val_predict
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut, RepeatedStratifiedKFold
from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 

# from sklearn.ensemble import BaggingClassifier

from joblib import Parallel, delayed, parallel_backend

from . import constants as gv 
from .options import * 
# from .glmnet_wrapper import logitnet, logitnetCV, logitnetIterCV, logitnetAlphaCV, logitnetAlphaIterCV 
from .bolasso_sklearn import bolasso 

import utils.progressbar as pg 

# import keras
# from .keras_wrapper import keras_logreg

# from .dot_clf import dot_clf

from . cross_temp_utils import temp_cross_val_score, temp_cross_validate

def get_clf(**kwargs):
    
    # options = set_options(**kwargs)
    # globals().update(options) 

    # if 'dot' in kwargs['clf_name']:
    #     clf = dot_clf(pval=kwargs['pval']) 
    
    # sklearn    
    if 'LDA' in kwargs['clf_name']: 
        clf = LinearDiscriminantAnalysis(tol=kwargs['tol'], solver='lsqr', shrinkage='auto')
        # clf = LinearDiscriminantAnalysis(tol=kwargs['tol'], solver='svd', shrinkage=None)
        
        # clf = LinearDiscriminantAnalysis() 
        
        if kwargs['standardize'] and kwargs['prescreen']: 
            clf = Pipeline([ ('scaler', StandardScaler()), 
                             ('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                             ('clf', clf) ] )             
        elif kwargs['standardize']: 
            clf = Pipeline([('scaler', StandardScaler()), ('clf', clf)]) 
        elif kwargs['prescreen']: 
            clf = Pipeline([('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                            ('clf', clf)])         
            # clf = Pipeline([('filter', SelectPercentile(f_classif, percentile=int(100*kwargs['pval'])) ), 
            #                 ('clf', clf)])         

        # shrinks = np.linspace(0, 1, 10)         
        # param_grid = dict(clf__shrinkage=shrinks) 
        # clf = GridSearchCV(clf, param_grid=param_grid, cv=kwargs['n_in'], n_jobs=1) 
    
    # if 'LinearSVC' in kwargs['clf_name']:
    #     clf = LinearSVC(C=C, penalty=penalty, loss=loss, dual=False,
    #                     tol=tol, max_iter=int(max_iter), multi_class='ovr',
    #                     fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
    #                     class_weight=None, verbose=0, random_state=None) 
    
    if 'LogisticRegressionCV' in kwargs['clf_name']:
    
        if kwargs['in_fold'] == 'stratified': 
            cv = StratifiedKFold(n_splits=kwargs['n_in'], shuffle=True, random_state=kwargs['random_state']) # outer cv loop for scoring 
        if kwargs['in_fold'] == 'loo': 
            cv = LeaveOneOut() 
        if kwargs['in_fold'] == 'repeated':
            cv = RepeatedStratifiedKFold(n_splits=kwargs['n_in'], n_repeats=100, random_state=kwargs['random_state'])
        
        if kwargs['penalty'] != 'elasticnet': 
            kwargs['l1_ratios']=None
        else :
            kwargs['l1_ratios']= np.linspace(0,1, kwargs['n_alpha'])            
            
        clf = LogisticRegressionCV(Cs=np.logspace(-4, 4, kwargs['n_lambda']), solver=kwargs['solver'], 
                                   penalty=kwargs['penalty'], l1_ratios=kwargs['l1_ratios'], 
                                   tol=kwargs['tol'], max_iter=int(kwargs['max_iter']), scoring=kwargs['inner_score'], 
                                   fit_intercept=kwargs['fit_intercept'], intercept_scaling=kwargs['intercept_scaling'], 
                                   cv=cv, n_jobs=kwargs['n_jobs'], verbose=0) 
        
        if kwargs['standardize'] and kwargs['prescreen']:
            
            clf = Pipeline([('scaler', StandardScaler() ), 
                            ('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                            ('clf', clf) ] ) 
        
        elif kwargs['standardize']: 
            clf = Pipeline([('scaler', StandardScaler()), ('clf', clf)]) 
        elif kwargs['prescreen']: 
            clf = Pipeline([('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                            ('clf', clf)])         
        
    elif 'LogisticRegression' in kwargs['clf_name']:
        if kwargs['penalty'] != 'elasticnet': 
            kwargs['l1_ratio']=None 
            
        clf = LogisticRegression(C=kwargs['C'], solver=kwargs['solver'], penalty=kwargs['penalty'],
                                 tol=kwargs['tol'], max_iter=int(kwargs['max_iter']),
                                 fit_intercept=kwargs['fit_intercept'],  intercept_scaling=kwargs['intercept_scaling'],
                                 l1_ratio=kwargs['l1_ratio'], n_jobs=kwargs['n_jobs']) 
    
        if kwargs['standardize'] and kwargs['prescreen']:
            
            clf = Pipeline([('scaler', StandardScaler() ), 
                            ('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                            ('clf', clf) ] ) 
        
        elif kwargs['standardize']: 
            clf = Pipeline([('scaler', StandardScaler()), ('clf', clf)]) 
        elif kwargs['prescreen']: 
            clf = Pipeline([('filter', SelectFpr(f_classif, alpha=kwargs['pval']) ), 
                            ('clf', clf)])         
    
    # glmnet_python
    # if 'logitnetAlphaIterCV' in kwargs['clf_name']:
    #     clf = logitnetAlphaIterCV(lbd=kwargs['lbd'], n_alpha=kwargs['n_alpha'], n_lambda=kwargs['n_lambda'], 
    #                           n_splits=kwargs['n_in'], fold_type=kwargs['in_fold'], scoring=kwargs['inner_score'],
    #                           standardize= False, fit_intercept=kwargs['fit_intercept'], prescreen=kwargs['prescreen'], 
    #                           thresh=kwargs['tol'], maxit=kwargs['max_iter'], n_jobs=None, verbose=kwargs['verbose']) 
    
    # elif 'logitnetAlphaCV' in kwargs['clf_name']:
    #     clf = logitnetAlphaCV(lbd=kwargs['lbd'], n_alpha=kwargs['n_alpha'], n_lambda=kwargs['n_lambda'], 
    #                           n_splits=kwargs['n_in'], fold_type=kwargs['in_fold'], scoring=kwargs['inner_score'],
    #                           standardize= False, fit_intercept=kwargs['fit_intercept'], prescreen=kwargs['prescreen'], 
    #                           thresh=kwargs['tol'], maxit=kwargs['max_iter'], n_jobs=None, verbose=kwargs['verbose']) 
    
    # elif 'logitnetIterCV' in kwargs['clf_name']:
    #     clf = logitnetIterCV(lbd=kwargs['lbd'], alpha=kwargs['alpha'], n_lambda=kwargs['n_lambda'], 
    #                          n_splits=kwargs['n_in'], fold_type=kwargs['fold_type'], scoring=kwargs['inner_score'],
    #                          standardize=True, fit_intercept=kwargs['fit_intercept'], prescreen=kwargs['prescreen'], 
    #                          thresh=kwargs['tol'], maxit=kwargs['max_iter'], n_jobs=None, verbose=kwargs['verbose']) 
    
    # elif 'logitnetCV' in kwargs['clf_name']:
    #     clf = logitnetCV(lbd=kwargs['lbd'], alpha=kwargs['alpha'], n_lambda=kwargs['n_lambda'], n_splits=kwargs['n_in'],
    #                      standardize=kwargs['standardize'], fit_intercept=kwargs['fit_intercept'],
    #                      prescreen=kwargs['prescreen'], confidence=kwargs['pval'],
    #                      fold_type=kwargs['in_fold'], scoring=kwargs['inner_score'],
    #                      thresh=kwargs['tol'] , maxit=kwargs['max_iter'], n_jobs=None) 
    
    # elif 'logitnet' in kwargs['clf_name']: 
    #     clf = logitnet(lbd=kwargs['lbd'], alpha=kwargs['alpha'], n_lambda=kwargs['n_lambda'],
    #                    prescreen=kwargs['prescreen'], 
    #                    standardize=False, fit_intercept=kwargs['fit_intercept'], 
    #                    scoring=kwargs['inner_score'],
    #                    thresh=kwargs['tol'], maxit=kwargs['max_iter'], verbose=0) 
    
    if 'bolasso' in kwargs['clf_name']:
        clf = bolasso(lbd=kwargs['lbd'], alpha=kwargs['alpha'], n_lambda=kwargs['n_lambda'], n_splits=kwargs['n_in'],
                      standardize=kwargs['standardize'], fit_intercept=kwargs['fit_intercept'], prescreen=kwargs['prescreen'],
                      fold_type=kwargs['in_fold'], scoring=kwargs['inner_score'], confidence=kwargs['pval'],
                      thresh=kwargs['tol'] , maxit=kwargs['max_iter'], n_jobs=-1, verbose=0) 
    
    # if 'keras_logreg' in kwargs['clf_name']:
    #     clf = keras_logreg(alpha=kwargs['alpha'], input_dim=kwargs['input_dim'])
    
    return clf

def set_scaler(clf, scaling):
    
    if scaling=='standard':
        pipe = Pipeline([('scaler', StandardScaler()), ('clf', clf)]) 
    elif scaling=='center':
        pipe = Pipeline([('scaler', StandardScaler(with_std=False)), ('clf', clf)])
    elif scaling=='robust':
        pipe = Pipeline([('scaler', RobustScaler()), ('clf', clf)]) 
    else:
        pipe = Pipeline([('clf', clf)])
    
    return pipe

def nested_cv(clf, X, y, param_grid, scaling='standard',
              in_fold='stratified', n_in=10, in_score='accuracy',
              out_fold='stratified', n_out=10, out_score='accuracy',
              random_state = None,
              fix_hyper=1, n_jobs=-1): 

    if fix_hyper:
        dum = n_jobs
    else:
        dum = None 
    
    grid = inner_cv(clf, param_grid, scaling=scaling,
                    n_in=n_in, fold_type=in_fold, scoring=in_score,
                    random_state=random_state, n_jobs=dum) # n_jobs=2 
    
    if fix_hyper==1: 
        grid.fit(X, y) 
        grid_cv = deepcopy(grid.best_estimator_) 
    else: 
        grid_cv = deepcopy(grid) 
    
    if out_fold == 'stratified': 
        cv_outer = StratifiedKFold(n_splits=n_out, shuffle=True, random_state=random_state) # outer cv loop for scoring 
    else: 
        cv_outer = LeaveOneOut() 
    
    cv_scores = cross_val_score(grid_cv, X, y, cv=cv_outer, scoring=out_score, n_jobs=n_jobs, error_score='raise') 
    
    mean_score = np.nanmean(cv_scores) 
    
    # if fix_hyper:
    #     return  mean_score , grid 
    # else:
    return  mean_score
    
def temp_nested_cv(clf, X_t_train, X_t_test, y, param_grid, scaling='standard',
                   in_fold='stratified', n_in=10, in_score='accuracy',
                   out_fold='stratified', n_out=10, out_score='accuracy',
                   fix_hyper=1, n_jobs=-1):
    
    if fix_hyper==1: 
        grid = inner_cv(clf, param_grid, scaling=scaling, n_in=n_in, fold_type=in_fold, scoring=in_score, n_jobs=n_jobs) # n_jobs=2 
        grid.fit(X_t_train, y) 
        grid_cv = deepcopy(grid.best_estimator_) 
    else: 
        grid_cv = deepcopy( inner_cv(clf, param_grid, scaling=scaling, n_in=n_in, fold_type=in_fold, scoring=in_score, n_jobs=None) ) 
    
    if out_fold == 'stratified': 
        cv_outer = StratifiedKFold(n_splits=n_out, shuffle=True, random_state=None) # outer cv loop for scoring 
    else: 
        cv_outer = LeaveOneOut() 
    
    cv_scores = temp_cross_val_score(grid_cv, X_t_train, X_t_test, y, cv=cv_outer, scoring=out_score, n_jobs=n_jobs, error_score='raise', verbose=0) 
    
    mean_score = np.nanmean(cv_scores) 
    
    return  mean_score 

def outer_cv(clf, X, y, scaling='standard', n_out=10, folds='stratified', outer_score='accuracy', inner_score='deviance', random_state=None, n_jobs=-1, return_clf=0): 
    
    # pipe = set_scaler(clf, scaling) 
    if folds == 'stratified': 
        cv_outer = StratifiedKFold(n_splits=n_out, shuffle=True, random_state=random_state) # outer cv loop for scoring 
    if folds == 'loo': 
        cv_outer = LeaveOneOut() 
    if folds == 'repeated':
        cv_outer = RepeatedStratifiedKFold(n_splits=n_out, n_repeats=100, random_state=random_state)
    
    cv_scores = cross_val_score(clf, X, y, cv=cv_outer, scoring=outer_score, n_jobs=n_jobs, verbose=1) 
    
    # cv_scores = glmnet_cv_loop(deepcopy(clf), X, y, cv=cv_outer, inner_score=inner_score, outer_score=outer_score, return_clf=0) 
    
    return np.nanmean(cv_scores) 

def outer_temp_cv(clf, X_t_train, X_t_test, y, scaling='standard', n_out=10, folds='stratified', outer_score='accuracy', inner_score='deviance', random_state=None, n_jobs=-1, return_clf=0): 
    
    # pipe = set_scaler(clf, scaling) 
    if folds == 'stratified': 
        cv_outer = StratifiedKFold(n_splits=n_out, shuffle=True, random_state=random_state) # outer cv loop for scoring 
    if folds == 'loo': 
        cv_outer = LeaveOneOut() 
    if folds == 'repeated':
        cv_outer = RepeatedStratifiedKFold(n_splits=n_out, n_repeats=100, random_state=random_state)
        
    cv_scores = temp_cross_val_score(clf, X_t_train, X_t_test, y, cv=cv_outer, scoring=outer_score, n_jobs=n_jobs) 
    
    # cv_scores = glmnet_temp_cv_loop(clf, X_t_train, X_t_test, y, cv=cv_outer, inner_score=inner_score, outer_score=outer_score, return_clf=0) 
    
    return np.nanmean(cv_scores) 

def inner_cv(clf, param_grid, scaling='standard', fold_type='stratified', n_in=10, scoring='accuracy', random_state=None, n_jobs=None): 
    
    pipe = clf 
    # pipe = set_scaler(clf, scaling) 
    
    if fold_type == 'stratified': 
        cv_inner = StratifiedKFold(n_splits=n_in, shuffle=True, random_state=random_state) 
    else: 
        cv_inner = LeaveOneOut() 
    
    grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv_inner, scoring=scoring, n_jobs=n_jobs) 
    
    return grid 

def glmnet_cv_loop(clf, X, y, cv, inner_score='deviance', outer_score='accuracy', return_clf=0):
    
    cv_scores = [] 
    cv_clfs = []
    
    for train_ix, test_ix in cv.split(X, y): 
        # split data 
        X_train, X_test = X[train_ix], X[test_ix] 
        y_train, y_test = y[train_ix], y[test_ix] 
        
        # fit the model 
        clf.scoring = inner_score 
        clf.fit(X_train, y_train) 
        # evaluate the model 
        clf.scoring = outer_score 
        cv_score = clf.score(X_test, y_test) 
        # store the result 
        cv_scores.append(cv_score) 
        
        if return_clf: 
            cv_clfs.append(deepcopy(clf)) 
        
    if return_clf: 
        return cv_scores, cv_clfs 
    else: 
        return cv_scores 

def glmnet_temp_cv_loop(clf, X_t_train, X_t_test, y, cv, inner_score='deviance', outer_score='accuracy', return_clf=0):
    
    cv_scores = [] 
    cv_clfs = []
    
    for train_ix, test_ix in cv.split(X_t_train, y): 
        # split data 
        X_train, X_test = X_t_train[train_ix], X_t_test[test_ix] 
        y_train, y_test = y[train_ix], y[test_ix] 
        
        # fit the model 
        clf.scoring = inner_score 
        clf.fit(X_train, y_train) 
        # evaluate the model 
        clf.scoring = outer_score 
        cv_score = clf.score(X_test, y_test) 
        # store the result 
        cv_scores.append(cv_score) 
        
        if return_clf: 
            cv_clfs.append(deepcopy(clf)) 
    
    if return_clf: 
        return cv_scores, cv_clfs 
    else: 
        return cv_scores 

def my_cv_loop(clf, X, y, cv): 
    
    cv_scores = [] 
    models = [] 
    Cs = [] 
    
    for train_ix, test_ix in cv.split(X, y):
        # split data 
        X_train, X_test = X[train_ix], X[test_ix] 
        y_train, y_test = y[train_ix], y[test_ix] 
        
        clf.fit(X_train, y_train)
        # get the best performing model fit on the whole training set 
        # best_model = clf.best_estimator_ 
        # evaluate the model 
        # cv_score = best_model.cv_score(X_test, y_test) 
        cv_score = clf.score(X_test, y_test) 
        # store the result 
        cv_scores.append(cv_score) 
        # models.append(best_model) 
        # Cs.append(best_model.named_steps['clf'].C) 
    
    # print('Cs', Cs) 
    
    return cv_scores 
    
def bootstrap_parloop(clf, X_t_train, X_t_test, y):
    
    np.random.seed(None) 
    clf_copy = deepcopy(clf) 
    
    n_0 = np.sum(y==0) 
    n_1 = np.sum(y==1) 
    
    idx_trials = np.hstack( ( np.random.randint(0, n_0, n_0), 
                              np.random.randint(n_0, X_t_train.shape[0], n_1) ) ) 
    
    X_sample = X_t_train[idx_trials] 
    y_sample = y[idx_trials] 
    
    array = np.arange(0, X_t_train.shape[0]) 
    idx_oob = np.delete(array, idx_trials) 
    
    X_oob = X_t_test[idx_oob] 
    y_oob = y[idx_oob] 
    
    clf_copy.fit(X_sample, y_sample) 
    oob_score = clf_copy.score(X_oob, y_oob) 
    
    return oob_score 

def bootstrap_score(clf, X_t_train, X_t_test, y, n_boots=1000, n_jobs=-1): 
    
    with pg.tqdm_joblib(pg.tqdm(desc='bootstrap', total=n_boots)) as progress_bar: 
        boots_score = Parallel(n_jobs=n_jobs)(delayed(bootstrap_parloop)(clf, X_t_train, X_t_test, y) for _ in range(n_boots) ) 
    boots_scores = np.array(boots_score) 
    
    return boots_scores

def bootstrap_score_ci(clf, X_t_train, X_t_test, y, n_boots=1000, n_jobs=-1):
    
    # grid = inner_cv(clf, X_t_train, y, param_grid, scaling=scaling, n_in=n_in, n_jobs=n_jobs) 
    # # tune hyperparams globally 
    # grid.fit(X, y)
    # grid = grid.best_estimator_
    
    boots_scores = bootstrap_score(clf, X_t_train, X_t_test, y, n_boots, n_jobs) 
    lower, upper = my_conf_int(boots_scores) 
    
    return np.mean(boots_scores), lower, upper 

def my_conf_int(stats, alpha=0.95):
    ostats = np.sort(stats, axis=0)
    mean_stat = np.mean(ostats, axis=0)
    
    p = ((1.0-alpha)/2.0) * 100    
    lower = mean_stat - max(0.0, np.percentile(ostats, p))
    p = (alpha+((1.0-alpha)/2.0)) * 100
    upper = -mean_stat + min(1.0, np.percentile(ostats, p)) 
    
    return lower, upper

def boots_coefs_parloop(clf, X, y):
    
    np.random.seed(None) 
    clf_copy = deepcopy(clf) 
    
    n_0 = np.sum(y==0) 
    n_1 = np.sum(y==1) 
    
    idx_trials = np.hstack( ( np.random.randint(0, n_0, n_0), 
                              np.random.randint(n_0, X.shape[0], n_1) ) ) 
    
    X_sample = X[idx_trials] 
    y_sample = y[idx_trials] 
    
    clf_copy.fit(X_sample, y_sample)
    coefs = clf_copy.coef_ 
    
    return coefs 

def bootstrap_coefs(clf, X, y, n_boots=1000, n_jobs=-1): 
    
    with pg.tqdm_joblib(pg.tqdm(desc='bootstrap coefs', total=n_boots)) as progress_bar: 
        boots_coefs = Parallel(n_jobs=n_jobs)(delayed(boots_coefs_parloop)(clf, X, y) for _ in range(n_boots) ) 
    boots_coefs = np.array(boots_coefs) 
    
    return boots_coefs 
    
    # return boots_coefs_parloop(clf, X, y)
    
