import warnings 
warnings.filterwarnings('ignore') 

from copy import deepcopy
from joblib import Parallel, delayed, parallel_backend

import numpy as np
from scipy.stats import ttest_1samp, ttest_ind, ttest_ind_from_stats 

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 
from sklearn.pipeline import Pipeline 
from sklearn.preprocessing import StandardScaler, RobustScaler 
from sklearn.model_selection import LeaveOneOut

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.base import BaseEstimator, ClassifierMixin

import utils.progressbar as pg 

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
    # coefs = clf_copy.coef_ 
    coefs = np.zeros(X.shape[-1]) 
    
    try:
        pval = clf_copy.named_steps['filter'].pvalues_ 
        coef = clf_copy.named_steps['clf'].coef_[0] 
        idx = pval<=confidence 
        # print('coef', coef.shape, 'non zero', np.sum(idx), np.sum(coef!=0) ) 
        coefs[idx] = coef 
    except:
        try:
            coef = clf_copy.named_steps['clf'].coef_ 
            coefs = coef 
        except:
            coefs = clf_copy.coef_ 
    
    return coefs 

def bootstrap_coefs(clf, X, y, n_boots=1000, confidence=0.05, n_jobs=-1): 
    
    with pg.tqdm_joblib(pg.tqdm(desc='bootstrap coefs', total=n_boots)) as progress_bar: 
        dum = Parallel(n_jobs=n_jobs)(delayed(boots_coefs_parloop)(clf, X, y, confidence=confidence) for _ in range(n_boots) )
    
    boots_coefs = np.array(dum) 
    
    return boots_coefs 
    
class bolasso(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, clf=None, n_boots=1000, confidence=0.05, alpha=0.5,
                 n_lambda=10, lbd=1, 
                 standardize='standard', fit_intercept=True,
                 scoring='neg_log_loss', fold_type='stratified',
                 n_splits=10, random_state=None, 
                 prescreen=True, f_screen='f_classif', 
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
            if self.fold_type == 'stratified':
                self.clf = LogisticRegressionCV(Cs=np.logspace(-4, 4, self.n_lambda), solver='liblinear',
                                                penalty='l1', l1_ratios=None,
                                                tol=self.thresh, max_iter=self.maxit, scoring=self.scoring, 
                                                fit_intercept=self.fit_intercept, intercept_scaling=1.0, 
                                                cv=self.n_splits, n_jobs=self.n_jobs, verbose=0) 
            else:
                self.clf = LogisticRegressionCV(Cs=np.logspace(-4, 4, self.n_lambda), solver='liblinear',
                                                penalty='l1', l1_ratios=None, 
                                                tol=self.thresh, max_iter=self.maxit, scoring=self.scoring, 
                                                fit_intercept=self.fit_intercept, intercept_scaling=1.0, 
                                                cv=LeaveOneOut(), n_jobs=self.n_jobs, verbose=0) 
        
        pipe = [] 
        if self.standardize=='standard':
            pipe.append( ('scaler', StandardScaler()) )
        if self.standardize=='center':
            pipe.append( ('scaler', StandardScaler(with_std=False)) )
        if self.standardize=='robust':
            pipe.append( ('scaler', RobustScaler()) )
        if self.prescreen:
            pipe.append( ('filter', SelectFpr(f_classif, alpha=self.confidence) ) )
        
        pipe.append(('clf', self.clf)) 
        self.model_ = Pipeline(pipe) 
        
        self.boots_coef_ = bootstrap_coefs(self.model_, X, y, n_boots=self.n_boots, confidence=self.confidence, n_jobs=self.n_jobs) 
        mean_coefs = np.nanmean(self.boots_coef_, axis=0) 
        self.mean_coef_ = mean_coefs 
        # std_coefs = np.nanstd(self.boots_coef_, axis=0) 
        
        # if self.verbose:
        #     print('boots_coefs', self.boots_coef_.shape, mean_coefs[:10]) 
        
        # _, self.p_val_ = ttest_1samp(self.boots_coef_, 0, axis=0, nan_policy='omit') 
        _, self.p_val_ = ttest_ind(self.boots_coef_, np.zeros(self.boots_coef_.shape), axis=0, equal_var=False, nan_policy='omit') 
        
        # # if self.verbose: 
        # print('p_val', self.p_val_.shape, self.p_val_[:5]) 
        # print('coefs', mean_coefs[:5]) 
        
        self.fs_idx_ = self.p_val_<=self.confidence 
        print('significant', np.sum(self.fs_idx_)) 
        # self.mean_coef_ = mean_coefs[self.fs_idx_] 
        
        # self.coef_ = self.mean_coef_ 
        
        self.coef_ = np.zeros(X.shape[-1]) 
        X_fs = X[:, self.fs_idx_] 
        
        # if self.verbose: 
        print('X_fs', X_fs.shape) 
        
        pipe = [] 
        if self.standardize: 
            pipe.append( ('scaler', StandardScaler()) )
        
        pipe.append(('clf', self.clf)) 
        self.model_ = Pipeline(pipe) 
        
        self.model_.named_steps['clf'].penalty = 'l2' 
        self.model_.named_steps['clf'].n_jobs = -1 
        
        self.model_.fit(X_fs, y) 
        
        coef = self.model_.named_steps['clf'].coef_[0] 
        print('coef', coef.shape, self.model_.named_steps['clf'].C_) 
        
        # if self.standardize:
        #     print('mean', self.model_.named_steps['scaler'].mean_.shape , 'var', self.model_.named_steps['scaler'].var_.shape) 
        #     coef /= self.model_.named_steps['scaler'].var_ 
        
        self.coef_[self.fs_idx_] = coef 
        self.intercept_ = self.model_.named_steps['clf'].intercept_ 
        
        # if self.standardize: 
        #     self.intercept_ -= np.dot(coef , self.model_.named_steps['scaler'].mean_) 
        
        return self 
    
    
