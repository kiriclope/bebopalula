import warnings 
warnings.filterwarnings('ignore') 

from copy import deepcopy
from joblib import Parallel, delayed, parallel_backend

import numpy as np
from scipy.stats import ttest_1samp, ttest_ind
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.linear_model import LogisticRegressionCV

import utils.progressbar as pg 
from .glmnet_wrapper import logitnetCV, logitnetAlphaCV, get_score 
# from .glms import bootstrap_coefs

def boots_coefs_parloop(clf, X, y, confidence=0.05):
    
    np.random.seed(None) 
    clf_copy = deepcopy(clf) 
    
    n_0 = np.sum(y==0) 
    n_1 = np.sum(y==1) 
    
    idx_trials = np.hstack( ( np.random.randint(0, n_0, n_0), 
                              np.random.randint(n_0, X.shape[0], n_1) ) ) 
    
    X_sample = X[idx_trials] 
    y_sample = y[idx_trials] 
    
    clf_copy.fit(X_sample, y_sample)
    coefs = np.zeros(X.shape[-1]) 
    
    try:
        pval = clf_copy.named_steps['filter'].pvalues_ 
        coef = clf_copy.named_steps['clf'].coef_ 
        idx = pval<=confidence
        if np.sum(idx)==0:
            coefs = coef
        else:
            print('non zero', np.sum(idx), np.sum(coef!=0) ) 
            coefs[idx] = coef 
    except:
        try:
            coef = clf_copy.named_steps['clf'].coef_ 
            coefs = coef 
        except:
            coefs = clf_copy.coef_ 
        
    return coefs 

def bootstrap_coefs(clf, X, y, n_boots=1000, n_jobs=-1): 
    
    with pg.tqdm_joblib(pg.tqdm(desc='bootstrap coefs', total=n_boots)) as progress_bar: 
        boots_coefs = Parallel(n_jobs=n_jobs)(delayed(boots_coefs_parloop)(clf, X, y) for _ in range(n_boots) ) 
    boots_coefs = np.array(boots_coefs) 
    
    return boots_coefs 
    
    # return boots_coefs_parloop(clf, X, y)

class bolasso(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, clf=None, n_boots=1000, confidence=0.05, alpha=0.5, n_lambda=10, lbd ='lambda_min', 
                 standardize=False, fit_intercept=True,
                 scoring='accuracy', fold_type='stratified',
                 n_splits=10, random_state=None, 
                 prescreen=False, f_screen='f_classif', 
                 thresh=1e-4 , maxit=1e3,
                 n_jobs=1, verbose=0):

        self.clf = clf 
        self.confidence = confidence 
        self.n_boots = n_boots
        
        self.alpha = alpha
        self.n_lambda = n_lambda
        self.lbd = lbd 
        
        self.n_splits = n_splits
        
        self.standardize = standardize 
        self.fit_intercept = fit_intercept 
        
        self.scoring = scoring 
        
        self.fold_type = fold_type 
        self.random_state = random_state 
        
        self.prescreen = prescreen 
        self.f_screen = f_screen 
        
        self.thresh = thresh
        self.maxit = maxit 
        
        self.n_jobs = n_jobs 
        self.verbose = verbose  
            
    def fit(self, X, y): 

        if self.clf is None:
            self.model_ = logitnetCV(alpha=self.alpha, n_lambda=self.n_lambda, lbd=self.lbd, 
                                     standardize=self.standardize, fit_intercept=self.fit_intercept, 
                                     scoring=self.scoring, fold_type=self.fold_type, 
                                     n_splits=self.n_splits, random_state=self.random_state, 
                                     prescreen=self.prescreen, f_screen=self.f_screen, 
                                     thresh=self.thresh, maxit=self.maxit, 
                                     n_jobs=None, verbose=0) 
        else: 
            self.model_ = self.clf 
        # self.model_ = logitnetAlphaCV(n_alpha=10, n_lambda=self.n_lambda, lbd=self.lbd, 
        #                               standardize=self.standardize, fit_intercept=self.fit_intercept, 
        #                               scoring=self.scoring, fold_type=self.fold_type, 
        #                               n_splits=self.n_splits, random_state=self.random_state, 
        #                               prescreen=self.prescreen, f_screen=self.f_screen, 
        #                               thresh=self.thresh, maxit=self.maxit, 
        #                               n_jobs=None, verbose=0) 
        
        self.boots_coef_ = bootstrap_coefs(self.model_, X, y, n_boots=self.n_boots, n_jobs=self.n_jobs) 
        mean_coefs = np.mean(self.boots_coef_, axis=0) 
        
        if self.verbose:
            print('boots_coefs', self.boots_coef_.shape, mean_coefs[:10]) 
        
        # _, self.p_val_ = ttest_1samp(self.boots_coef_, 0, axis=0) 
        _, self.p_val_ = ttest_ind(self.boots_coef_, np.zeros(self.boots_coef_.shape), axis=0, equal_var=0, nan_policy='omit') 
        
        if self.verbose: 
            print('p_val', self.p_val_.shape, self.p_val_[:5]) 
            print('coef', self.coef_.shape, self.coef_[:5]) 
                
        self.fs_idx_ = self.p_val_<=self.confidence 
        self.mean_coef_ = mean_coefs[self.fs_idx_] 
        
        X_fs = X[:, self.fs_idx_]
        
        print('X_fs', X_fs.shape, 'significant', np.sum(self.fs_idx_)) 
        
        self.model_ = LogisticRegressionCV(Cs=np.logspace(-4, 4, self.n_lambda), solver='liblinear',
                                           penalty='l2', l1_ratios=None,
                                           tol=self.thresh, max_iter=self.maxit, scoring=self.scoring, 
                                           fit_intercept=self.fit_intercept, intercept_scaling=1.0, 
                                           cv=self.n_splits, n_jobs=self.n_jobs, verbose=0) 
        
        # self.model_.alpha = 0 # ridge regression 
        # self.lbd = 0 
        self.model_.n_jobs = self.n_jobs 
        self.model_.fit(X_fs, y) 

        print('coefs',  self.model_.coef_.shape)
        
        self.coef_ = np.zeros(X.shape[-1]) 
        self.coef_[self.fs_idx_] = self.model_.coef_[0]
        # self.coef_[self.fs_idx_] = self.mean_coef_ 
        self.intercept_ = self.model_.intercept_ 
        
        return self 
    
    def predict(self, X):
        X_fs = X[:, self.fs_idx_] 
        return cvglmnetPredict(self.model_, newx=X_fs, ptype='class', s=self.lbd) 
    
    def predict_proba(self, X): 
        X_fs = X[:, self.fs_idx_] 
        return cvglmnetPredict(self.model_, newx=X_fs, ptype='response', s=self.lbd) 
    
    def score(self, X, y):
        
        if self.scoring =='accuracy' or self.scoring=='mean_squared_error' or self.scoring=='neg_mean_squared_error':
            y_pred = self.predict(X) 
        else: 
            y_pred = self.predict_proba(X) 
        
        return get_score(y, y_pred, self.scoring) 
