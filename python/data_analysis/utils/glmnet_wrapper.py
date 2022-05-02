import warnings 
warnings.filterwarnings('ignore') 

from copy import deepcopy

import pandas as pd

import scipy
import scipy.stats as stats
import numpy as np

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.linear_model._base import LinearClassifierMixin, SparseCoefMixin
from sklearn.metrics import roc_auc_score, accuracy_score, log_loss, mean_squared_error
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut
from sklearn.feature_selection import SelectKBest, f_classif, SelectFpr 
from sklearn.utils.extmath import safe_sparse_dot, softmax 

from glmnet_python.glmnetSet import glmnetSet 

from glmnet_python.glmnet import glmnet 
from glmnet_python.glmnetCoef import glmnetCoef 
from glmnet_python.glmnetPredict import glmnetPredict 
from glmnet_python.glmnetPlot import glmnetPlot

from glmnet_python.cvglmnet import cvglmnet 
from glmnet_python.cvglmnetCoef import cvglmnetCoef 
from glmnet_python.cvglmnetPredict import cvglmnetPredict 

from glmnet_python.cvglmnetPlot import cvglmnetPlot

from joblib import Parallel, delayed

from . import progressbar as pg

def set_lbd(lbd): 
    if isinstance(lbd, str): 
        lbd = scipy.float64([1]) 
    else: 
        lbd = scipy.float64([lbd])
        if lbd.shape == (1,1): 
            lbd = lbd[0]
    
    return lbd

def set_score(scoring):
                    
    if 'accuracy' in scoring: 
        scoring = 'class' 
    if 'roc_auc' in scoring: 
        scoring = 'auc' 
    if ('log_loss' or 'neg_log_loss') in scoring: 
        scoring = 'deviance' 
    if ('mean_squared_error' or 'neg_mean_squared_error') in scoring: 
        scoring = 'mse' 
    
    return scoring 

def get_score(y, y_pred, scoring):
    
    if scoring=='accuracy': 
        score =  accuracy_score(y, y_pred) 
    if scoring=='roc_auc': 
        score =  roc_auc_score(y, y_pred) 
    if scoring=='log_loss': 
        score = log_loss(y, y_pred) 
    if scoring=='neg_log_loss': 
        score = -log_loss(y, y_pred) 
    if scoring=='mean_squared_error': 
        score = mean_squared_error(y, y_pred) 
    if scoring=='neg_mean_squared_error': 
        score = -mean_squared_error(y, y_pred)  
    
    return score 
    
def pre_screen_fold(X, y, p_alpha=0.05, scoring=f_classif): 
    ''' X is trials x neurons 
    alpha is the level of significance 
    scoring is the statistics, use f_classif or mutual_info_classif 
    '''    
    model = SelectKBest(score_func=scoring, k=X.shape[1])    
    model.fit(X,y) 
    pval = model.pvalues_.flatten() 
    idx_out = np.argwhere(pval>p_alpha) 
    X_screen = np.delete(X, idx_out, axis=1) 
        
    return X_screen

def createFolds(X, y, fold_type='stratified', n_splits=10, random_state=None):

    if random_state is None:
        random_state = int(np.random.rand()*1000)         
    
    if fold_type=='stratified': 
        folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state) 
    elif fold_type == 'kfold': 
        folds = KFold(n_splits=n_splits, shuffle=True, random_state=random_state) 
    elif fold_type=='loo':
        folds = KFold(n_splits=X.shape[0], shuffle=True, random_state=random_state) 
        
    foldid = np.empty(y.shape[0]) 
    i_split = -1 
    for idx_train, idx_test in folds.split(X,y): 
        i_split = i_split+1 
        foldid[idx_test] = scipy.int32(i_split) 
            
    foldid = scipy.int32(foldid) 
    
    return foldid

def set_options(alpha=1, n_lambda=10, standardize=False, fit_intercept=True, thresh=1e-4, maxit=10000, prescreen=False, f_screen='f_classif', confidence=0.05):
    
    opts = dict() 
    opts['alpha'] = scipy.float64(alpha) 
    opts['nlambda'] = scipy.int32(n_lambda) 
    # opts['lambdau'] = scipy.float64( np.exp(np.linspace(-10, 0, opts['nlambda']) ) ) 
    # opts['lambdau'] = scipy.float64( np.logspace(-4, 4, n_lambda) ) 
    # opts['lambdau'] = scipy.array( -np.sort(-np.logspace(np.log(0.5), np.log(0.01), n_lambda, base=np.exp(1)) ) ) 
    
    opts['standardize'] = standardize 
    opts['intr'] = fit_intercept 
    
    opts['thresh'] = scipy.float64(thresh) 
    opts['maxit'] = scipy.int32(maxit) 
    
    opts['prescreen'] = prescreen 
    opts['f_screen'] = f_screen 
    opts['confidence'] = confidence
    
    options = glmnetSet(opts) 
    
    return options 


def pre_screen_fold(X, y, p_alpha=0.05, scoring=f_classif): 
    
    # model = SelectKBest(score_func=scoring, k=X.shape[1]) 
    model = SelectFpr(score_func=scoring, alpha=p_alpha) 
    model.fit(X,y) 
    pval = model.pvalues_ 
    
    idx_screen = pval<=p_alpha 
    X_screen = X[:, idx_screen] 
    
    return X_screen, idx_screen 

class logitnet(LinearClassifierMixin, SparseCoefMixin, BaseEstimator): 
    
    def __init__(self, alpha=1, n_lambda=10, lbd=1, scoring='accuracy',
                 prescreen=False, f_screen='f_classif',
                 standardize=True, fit_intercept=True,
                 thresh=1e-4 , maxit=1e6, verbose=0): 
        
        self.alpha = alpha 
        self.n_lambda = n_lambda 
        self.lbd = lbd 
        
        self.scoring = scoring 
        self.prescreen = prescreen 
        self.f_screen = f_screen 
        self.standardize = standardize 
        
        self.fit_intercept = fit_intercept 
        self.thresh = thresh 
        self.maxit = maxit 
        self.classes_ = np.array([-1, 1]) # useful if using cross_val_score 
        self.verbose = verbose 
        
    def fit(self, X, y):
        
        self.lbd = set_lbd(self.lbd) 
        
        self.options = set_options(self.alpha, self.n_lambda, self.standardize, self.fit_intercept,
                                   self.thresh, self.maxit, self.prescreen, self.f_screen) 

        if self.options['prescreen']:
            X_screen, self.idx_screen = pre_screen_fold(X, y, p_alpha=0.05, scoring=f_classif)             
            self.model_ = glmnet(x = X_screen.copy(), y = y.copy(), family = 'binomial', **self.options) 
        else:
            self.model_ = glmnet(x = X.copy(), y = y.copy(), family = 'binomial', **self.options)
            # self.model_ = model_ # for some reason I have to pass it like that a = funcnet() then self.a = a 
        
        coefs = self.get_coefs() 
        
        if self.verbose:
            print('alpha', self.alpha, 'lambda', self.lbd, 'coefs', coefs.shape) 
        
        return self 
    
    def get_coefs(self):
        coefs = glmnetCoef(self.model_, s=self.lbd, exact=False) 
        
        if self.fit_intercept: 
            self.intercept_ = coefs[0] 
            self.coef_ = coefs[1:] 
        else: 
            self.intercept_ = None 
            self.coef_ = coefs 
        
        if self.coef_.shape[-1] == 1:
            self.coef_ = self.coef_[:,0]
        
        return coefs 
    
    def predict(self, X): 
        if self.options['prescreen']:
            X = X[:, self.idx_screen]
            
        y_pred = glmnetPredict(self.model_, newx=X, ptype='class', s=self.lbd )
        # print(y_pred.shape) 
        return y_pred
        # return y_pred[..., 0] 
    
        # scores = self.decision_function(X)
        # if len(scores.shape) == 1:
        #     indices = (scores > 0).astype(int)
        # else:
        #     indices = scores.argmax(axis=1)
        # return self.classes_[indices]
    
    def predict_proba(self, X):  
        if self.options['prescreen']:
            X = X[:, self.idx_screen]
        
        y_pred = glmnetPredict(self.model_, newx=X, ptype='response', s=self.lbd )
        # print(y_pred.shape) 
        return y_pred
        # return y_pred[..., 0] 
        
        # decision = self.decision_function(X) 
        # if decision.ndim == 1: 
        #         decision_2d = np.c_[-decision, decision] 
        # else: 
        #     decision_2d = decision 
        # return softmax(decision_2d, copy=False)
    
    def score(self, X, y):
        
        if self.scoring =='accuracy' or self.scoring=='mean_squared_error' or self.scoring=='neg_mean_squared_error' :
            y_pred = self.predict(X) 
        else: 
            y_pred = self.predict_proba(X) 
        
        if self.verbose:
            print('y', y.shape, 'y_pred', y_pred.shape) 
        
        return get_score(y, y_pred, self.scoring) 
    
    def decision_function(self, X):
        if self.verbose:
            print('X', X.shape, 'coef_', self.coef_.shape, 'intercept_', self.intercept_.shape) 
        scores = safe_sparse_dot(X, self.coef_.T, dense_output=True) + self.intercept_ 
        return scores.ravel() if scores.shape[-1] == 1 else scores 
    
class logitnetCV(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, alpha=0.5, n_lambda=10, lbd ='lambda_min', 
                 standardize=False, fit_intercept=True,
                 scoring='accuracy', fold_type='stratified',
                 n_splits=10, random_state=None, 
                 prescreen=False, f_screen='f_classif', confidence=0.05,
                 thresh=1e-4 , maxit=1e6,
                 n_jobs=1, verbose=0): 
        
        self.alpha = alpha
        self.n_lambda = n_lambda
        self.lbd = lbd 
        
        self.n_splits = n_splits
        
        self.standardize = standardize 
        self.fit_intercept = fit_intercept 
        
        self.scoring = scoring # 'deviance', 'class', 'auc', 'mse' or 'mae' 
        
        self.fold_type = fold_type 
        self.random_state = random_state 
        
        self.prescreen = prescreen 
        self.f_screen = f_screen 
        self.confidence = confidence
        
        self.thresh = thresh
        self.maxit = maxit 
        
        self.n_jobs = n_jobs 
        self.verbose = verbose 

    def fit_func(self, X, y):
        
        self.model_ = cvglmnet(x = X.copy(), y = y.copy(), family = 'binomial',
                               ptype = self.glmnet_scoring, foldid=self.foldid,
                               n_jobs=self.n_jobs, **self.options) 
        
        self.lbd_min_ = self.model_['lambda_min']  
        self.lbd_1se_ = self.model_['lambda_1se'] 

        if self.verbose :
            print('alpha', self.alpha, 'lbd min', self.lbd_min_, 'lbd 1se', self.lbd_1se_) 
        
        self.non_zero_min_ = self.model_['nzero'][ self.model_['lambdau'] == self.lbd_min_ ] 
        self.non_zero_1se_ = self.model_['nzero'][ self.model_['lambdau'] == self.lbd_1se_ ] 
        
        if self.verbose : 
            print('non zero min', self.non_zero_min_, 'non zero 1se', self.non_zero_1se_)
            
        return self
    
    def fit(self, X, y): 
        self.glmnet_scoring = set_score(self.scoring) 
        
        self.options = set_options(self.alpha, self.n_lambda, self.standardize, self.fit_intercept,
                                   self.thresh, self.maxit, self.prescreen, self.f_screen, self.confidence) 
        
        self.foldid = createFolds(X, y, self.fold_type, self.n_splits, self.random_state) 
        
        self.fit_func(X, y)
        
        if self.lbd == 'lambda_min':
            i=0
            
            # if self.lbd_min_ == self.model_['lambdau'][0] :                
            #     print('warning lbd min is on the edge of the lambda path, log(lbd)=', np.log(self.lbd_min_) , 'non zero min', self.non_zero_min_ ) 
            #     self.options['lambdau'] = scipy.float64( np.logspace(-4-i, -3, self.n_lambda) ) 
            #     i+=1 
            #     self.fit_func(X, y) 
            
            # if self.lbd_min_ == self.model_['lambdau'][-1] : 
            #     print('warning lbd min is on the edge of the lambda path, log(lbd)=', np.log(self.lbd_min_) , 'non zero min', self.non_zero_min_ ) 
            #     self.options['lambdau'] = scipy.float64( np.logspace(1, 2+i, self.n_lambda) ) 
            #     i+=1 
            #     self.fit_func(X, y) 
        
        # if self.lbd == 'lambda_1se': 
        #     if self.lbd_1se_ == self.model_['lambdau'][0] or self.lbd_1se_ == self.model_['lambdau'][-1] : 
        #         print('warning lbd 1se is on the edge of the lambda path , log(lbd)=', np.log(self.lbd_1se_), 'non zero 1se', self.non_zero_1se_ ) 
        
        self.get_coefs() 
        
        return self
    
    def get_coefs(self):
        coefs = cvglmnetCoef(self.model_, s=self.lbd) 
        
        if self.fit_intercept: 
            self.coef_ = coefs[1:] 
            self.intercept_ = coefs[0] 
        else:
            self.coef_ = coefs 
            self.intercept_ = None 
            
        if self.coef_.shape[-1] == 1:
            self.coef_ = self.coef_[:,0]
        
        return self 
    
    def lasso_path(self): 
        cvglmnetPlot(self.model_) 
        return self 
    
    def predict(self, X): 
        return cvglmnetPredict(self.model_, newx=X, ptype='class', s=self.lbd) 
    
    def predict_proba(self, X): 
        return cvglmnetPredict(self.model_, newx=X, ptype='response', s=self.lbd) 
    
    def score(self, X, y):
        
        if self.scoring =='accuracy' or self.scoring=='mean_squared_error' or self.scoring=='neg_mean_squared_error':
            y_pred = self.predict(X) 
        else: 
            y_pred = self.predict_proba(X) 
        
        return get_score(y, y_pred, self.scoring) 

class logitnetAlphaCV(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, n_alpha=10, n_lambda=10, lbd='lbd_min', 
                 n_splits=10, fold_type='stratified', scoring='accuracy', 
                 standardize=True, fit_intercept=True, random_state=None,
                 prescreen=False, f_screen='f_classif', 
                 thresh=1e-4 , maxit=1e6, n_jobs=None, verbose=True): 

        self.alpha = None 
        self.n_alpha = n_alpha         
        self.n_lambda = n_lambda 
        self.lbd = lbd 
        
        self.n_splits = n_splits
        self.fold_type = fold_type 
        
        self.scoring = scoring 
        
        self.standardize = standardize 
        self.fit_intercept = fit_intercept 
        self.random_state = random_state
        
        self.prescreen = prescreen 
        self.f_screen = f_screen 
        
        self.thresh = thresh
        self.maxit = maxit
        
        self.n_jobs = n_jobs         
        self.verbose= verbose 
        
    def fit(self, X, y): 

        self.alpha_path = np.linspace(0, 1, self.n_alpha) 
        
        glmnet_scoring = set_score(self.scoring) 
        
        self.options = set_options(self.alpha, self.n_lambda, self.standardize, self.fit_intercept,
                                   self.thresh, self.maxit, self.prescreen, self.f_screen)
        
        foldid = createFolds(X, y, self.fold_type, self.n_splits, self.random_state) 
        
        # fit each model along the alpha path 
        with pg.tqdm_joblib(pg.tqdm(desc='alpha path', total= int(self.n_alpha * (self.n_splits+1) ), disable=True ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.fitFixedAlpha)(X, y, glmnet_scoring, i_alpha, self.alpha_path,
                                                                           foldid, self.options) 
                                               for i_alpha in range(self.n_alpha) ) 
        
        self.models_ = scipy.array(dum) 
        # if self.verbose: 
        #     print('models', self.models_.shape) 
        
        # compute min score 
        with pg.tqdm_joblib(pg.tqdm(desc='cvm min', total= int(self.n_alpha) , disable=True ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.minScoreAlpha)(self.models_[i_alpha]) for i_alpha in range(self.n_alpha) ) 
        self.cvms_min = scipy.array(dum) 
        
        self.idx_alpha_min_ = np.argmin(self.cvms_min) 
        self.alpha_ = self.alpha_path[self.idx_alpha_min_] 
        
        self.model_= self.models_[self.idx_alpha_min_] 
        self.lbd_1se_ = self.model_['lambda_1se'] 
        self.lbd_min_ = self.model_['lambda_min'] 
        
        self.non_zero_1se = self.model_['nzero'][ self.model_['lambdau'] == self.model_['lambda_1se'] ] 
        self.non_zero_min = self.model_['nzero'][ self.model_['lambdau'] == self.model_['lambda_min'] ] 
        
        if self.verbose:
            print('best alpha', self.alpha_, 'lambda min', self.lbd_min_, 'lambda 1se', self.lbd_1se_) 
            print('non zero min', self.non_zero_min, 'non zero 1se', self.non_zero_1se) 

        self.get_coefs() 
        
        if self.verbose:
            print(self.coef_.shape)
        
        return self 
    
    def get_coefs(self):
        coefs = cvglmnetCoef(self.model_, s=self.lbd) 
        if self.fit_intercept: 
            self.coef_ = coefs[1:] 
            self.intercept_ = coefs[0] 
        else:
            self.coef_ = coefs
            self.intercept_ = None 
        
        if self.coef_.shape[-1] == 1:
            self.coef_ = self.coef_[:,0]
        
        return self
    
    def fitFixedAlpha(self, X, y, scoring, i_alpha, alpha_path, foldid, options): 
        
        opts = options.copy() 
        opts['alpha'] = alpha_path[i_alpha] 
        
        model_ = cvglmnet(x = X.copy(), y = y.copy(), family = 'binomial',
                          ptype = scoring, foldid=foldid,
                          n_jobs=self.n_jobs, **opts) 
        return model_ 
    
    def minScoreAlpha(self, model_): 
        
        idx_lbd_min = model_['lambdau'] == model_['lambda_min'] 
        cvm_min = model_['cvm'][idx_lbd_min] 
        
        return cvm_min 
    
    def lasso_path(self): 
        cvglmnetPlot(self.model_) 
    
        return self
    
    def predict(self, X): 
        return cvglmnetPredict(self.model_, newx=X, ptype='class', s=self.lbd ) 
            
    def predict_proba(self, X): 
        return cvglmnetPredict(self.model_, newx=X, ptype='response', s=self.lbd) 

    def score(self, X, y):
        self.scoring = set_score(self.scoring) 
        
        if self.scoring =='accuracy' or self.scoring=='mean_squared_error' or self.scoring=='neg_mean_squared_error' :
            y_pred = self.predict(X) 
        else:
            y_pred = self.predict_proba(X) 
        
        return get_score(y, y_pred, self.scoring) 
    
class logitnetIterCV(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, n_iter=100, alpha=1, n_lambda=100, lbd='lbd_1se', 
                 n_splits=10, fold_type='kfold', scoring='accuracy',
                 standardize=True, fit_intercept=True,
                 prescreen=True, f_screen='f_classif',
                 thresh=1e-4 , maxit=1e6, n_jobs=1, verbose=False):
        
        opts = dict() 
        opts['nlambda'] = scipy.int32(n_lambda) 
        opts['lambdau'] = scipy.array( -np.sort(-np.logspace(np.log(0.5), np.log(0.01), opts['nlambda'], base=np.exp(1)) ) ) 
        
        opts['alpha'] = scipy.float64(alpha) 
        opts['standardize'] = standardize 
        opts['intr'] = fit_intercept 
        opts['prescreen'] = prescreen 
        opts['f_screen'] = f_screen 
        
        opts['thresh'] = scipy.float64(thresh) 
        opts['maxit'] = scipy.int32(maxit) 
        
        self.alpha_ = alpha
        self.n_iter = n_iter
        self.lbd = lbd
        self.options = glmnetSet(opts) 
        
        self.scoring = scoring # 'deviance', 'class', 'auc', 'mse' or 'mae'         
        if 'accuracy' in scoring: 
            self.scoring = 'class'  
        if 'roc_auc' in scoring: 
            self.scoring = 'auc' 
        if 'log_loss' in scoring: 
            self.scoring = 'deviance' 
            
        self.n_lambda = n_lambda
        
        self.n_splits = scipy.int32(n_splits) 
        self.fold_type = fold_type
        
        self.n_jobs = n_jobs 
        self.verbose = verbose 

    def foldsLoop(self, X, y): 
        self.random_state = np.random.randint(1e6) 
        
        if self.fold_type=='stratified':
            folds = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state) 
        else: 
            folds = KFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state) 
            
        foldid = np.empty(y.shape[0]) 
        i_split = -1 
        for idx_train, idx_test in folds.split(X,y): 
            i_split = i_split+1 
            foldid[idx_test] = i_split
            
        return scipy.int32(foldid)
        
    def fit(self, X, y):
        
        # create foldids for each iteration 
        with pg.tqdm_joblib(pg.tqdm(desc='create foldids', total= int(self.n_iter), disable=not self.verbose ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.foldsLoop)(X, y) for _ in range(self.n_iter) ) 
        
        foldid = scipy.array(dum) 
        if self.verbose: 
            print(foldid.shape) 
            
        # fit cvglmnet for each iteration and create dataframe 
        with pg.tqdm_joblib(pg.tqdm(desc='iter', total= int(self.n_iter * (self.n_splits+1) ),
                                    disable=not self.verbose ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.fit_one_iter)(X, y, i_iter, foldid[i_iter], self.options) 
                                               for i_iter in range(self.n_iter) )
        
        df = pd.concat(dum) 
        df = df.groupby(['lambdau'])['cvm', 'cvsd'].mean() 
        
        if self.verbose: 
            print(df)
        
        min_cvm = df['cvm'].min() 
        idx_min_cvm = df['cvm'].idxmin()
        
        if self.verbose: 
            print(idx_min_cvm, min_cvm) 
        
        self.lbd_min_ = scipy.float64( [idx_min_cvm] ) 

        # fit on the entire data 
        model_ = glmnet(x = X.copy(), y = y.copy(), family = 'binomial', **self.options) 
        self.model_ = model_ 
        
        # get coefs for lambda = lambda_min 
        coefs = glmnetCoef(self.model_, s = self.lbd_min_, exact = False) 
        self.intercept_ = coefs[0] 
        self.coef_ = coefs[1:] 
        
        self.non_zero_min = self.model_['df'][ self.model_['lambdau'] == self.lbd_min_ ] 
        
        if self.verbose:
            print('non zero min', self.non_zero_min) 
        
        return self 
    
    def fit_one_iter(self, X, y, i_iter, foldid, options): 
        
        opts = options.copy()         
        model_ = cvglmnet(x = X.copy(), y = y.copy(), family = 'binomial',ptype = self.scoring, foldid=foldid, 
                          nfolds=self.n_splits, n_jobs=self.n_jobs, **opts)
        
        df = pd.DataFrame({ 'i_iter': i_iter,
                            'lambdau': model_['lambdau'],
                            'cvm': model_['cvm'], 
                            'cvsd': model_['cvsd'] } )
        return df 
    
    def lasso_path(self):
        glmnetPlot(self.model_, xvar = 'lambda', label = False);
        
    def predict(self, X):
        return glmnetPredict(self.model_, newx=X, ptype='class', s = self.lbd_min_ ) 
            
    def predict_proba(self, X): 
        return glmnetPredict(self.model_, newx=X, ptype='response', s = self.lbd_min_ ) 
    
    def score(self, X, y): 
        if self.scoring=='class' or self.scoring=='accuracy': 
            y_pred = self.predict(X) 
            return accuracy_score(y, y_pred) 
        if self.scoring=='auc' or self.scoring=='roc_auc' : 
            y_pred = self.predict_proba(X) 
            return roc_auc_score(y, y_pred) 
        if self.scoring=='deviance' or self.scoring=='log_loss': 
            y_pred = self.predict_proba(X) 
            return log_loss(y, y_pred) 
        if self.scoring=='mse': 
            y_pred = self.predict(X) 
            return mean_squared_error(y, y_pred) 
        
class logitnetAlphaIterCV(BaseEstimator, ClassifierMixin): 
    
    def __init__(self, n_iter=100, n_alpha=10, n_lambda=100, lbd='lbd_1se', 
                 n_splits=10, fold_type='kfold', scoring='accuracy',
                 standardize=True, fit_intercept=True,
                 prescreen=True, f_screen='f_classif',
                 thresh=1e-4 , maxit=1e6, n_jobs=1, verbose=False):
        
        opts = dict() 
        opts['nlambda'] = scipy.int32(n_lambda) 
        # opts['lambdau']= scipy.array( -np.sort(-np.logspace(-4, -1, opts['nlambda'])) )  
        opts['lambdau'] = scipy.array( -np.sort(-np.logspace(np.log(0.5), np.log(0.01), opts['nlambda'], base=np.exp(1)) ) ) 
        
        opts['standardize'] = standardize 
        opts['intr'] = fit_intercept 
        opts['prescreen'] = prescreen
        opts['f_screen'] = f_screen 
        
        opts['thresh'] = scipy.float64(thresh) 
        opts['maxit'] = scipy.int32(maxit) 
        
        self.n_iter = n_iter
        self.lbd = lbd
        self.options = glmnetSet(opts)                
        
        self.scoring = scoring # 'deviance', 'class', 'auc', 'mse' or 'mae'         
        if 'accuracy' in scoring: 
            self.scoring = 'class'  
        if 'roc_auc' in scoring: 
            self.scoring = 'auc' 
        if 'log_loss' in scoring: 
            self.scoring = 'deviance' 
            
        self.n_alpha = n_alpha 
        self.alpha_path = np.linspace(.1, 1, n_alpha) 
        self.n_lambda = n_lambda
        
        self.n_splits = scipy.int32(n_splits) 
        self.fold_type = fold_type
        
        self.n_jobs = n_jobs 
        self.verbose = verbose 

    def foldsLoop(self, X, y):
        # fixed seed accross alphas 
        self.random_state = np.random.randint(1e6) 
        
        if self.fold_type=='stratified':
            folds = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state) 
        else: 
            folds = KFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state) 
            
        foldid = np.empty(y.shape[0]) 
        i_split = -1 
        for idx_train, idx_test in folds.split(X,y): 
            i_split = i_split+1 
            foldid[idx_test] = i_split
            
        return scipy.int32(foldid)
        
    def fit(self, X, y):
        
        # create foldids for each iteration 
        with pg.tqdm_joblib(pg.tqdm(desc='create foldids', total= int(self.n_iter), disable=not self.verbose ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.foldsLoop)(X, y) for _ in range(self.n_iter) ) 
        
        foldid = scipy.array(dum) 
        if self.verbose: 
            print(foldid.shape) 
            
        # fit cvglmnet for each iteration for each alpha and create dataframe
        with pg.tqdm_joblib(pg.tqdm(desc='alpha path', total= int(self.n_alpha * self.n_iter * (self.n_splits+1) ),
                                    disable=not self.verbose ) ) as progress_bar: 
            dum = Parallel(n_jobs=self.n_jobs)(delayed(self.fit_fixed_alpha)(X, y, i_alpha, i_iter, self.alpha_path, foldid[i_iter], self.options) 
                                               for i_alpha in range(self.n_alpha) for i_iter in range(self.n_iter) ) 
        df = pd.concat(dum) 
        df = df.groupby(['i_alpha','lambdau'])['cvm', 'cvsd'].mean() 
        if self.verbose: 
            print(df)
        
        min_cvm = df['cvm'].min() 
        idx_min_cvm = df['cvm'].idxmin() 
        if self.verbose: 
            print(idx_min_cvm, min_cvm) 
        
        self.idx_alpha_min_ = idx_min_cvm[0]
        self.lbd_min_ = scipy.float64( [idx_min_cvm[1]] )
        
        # fit on the entire data for alpha = alpha_min
        self.alpha_ =  self.alpha_path[self.idx_alpha_min_] 
        self.options['alpha'] = self.alpha_path[self.idx_alpha_min_] 
        model_ = glmnet(x = X.copy(), y = y.copy(), family = 'binomial', **self.options) 
        self.model_ = model_ 
        
        # get coefs for lambda = lambda_min
        coefs = glmnetCoef(self.model_, s = self.lbd_min_, exact = False) 
        self.intercept_ = coefs[0] 
        self.coef_ = coefs[1:] 
        
        self.non_zero_min = self.model_['df'][ self.model_['lambdau'] == self.lbd_min_ ] 
        
        if self.verbose:
            print('non zero min', self.non_zero_min, 'non zero 1se', self.non_zero_1se) 
        
        return self 
    
    def fit_fixed_alpha(self, X, y, i_alpha, i_iter, alpha_path, foldid, options): 
        
        opts = options.copy() 
        opts['alpha'] = alpha_path[i_alpha] 
        
        model_ = cvglmnet(x = X.copy(), y = y.copy(), family = 'binomial', ptype = self.scoring, foldid=foldid, grouped=True,
                          nfolds=self.n_splits, n_jobs=self.n_jobs, **opts)
        
        df = pd.DataFrame({ 'i_alpha': i_alpha, 
                            'i_iter': i_iter,
                            'lambdau': model_['lambdau'],
                            'cvm': model_['cvm'], 
                            'cvsd': model_['cvsd'] } )
        return df         
    
    def lasso_path(self):
        glmnetPlot(self.model_, xvar = 'lambda', label = False);
        
    def predict(self, X):
        return glmnetPredict(self.model_, newx=X, ptype='class', s = self.lbd_min_ ) 
            
    def predict_proba(self, X): 
        return glmnetPredict(self.model_, newx=X, ptype='response', s = self.lbd_min_ ) 
    
    def score(self, X, y): 
        if self.scoring=='class' or self.scoring=='accuracy': 
            y_pred = self.predict(X) 
            return accuracy_score(y, y_pred) 
        if self.scoring=='auc' or self.scoring=='roc_auc' : 
            y_pred = self.predict_proba(X) 
            return roc_auc_score(y, y_pred) 
        if self.scoring=='deviance' or self.scoring=='log_loss': 
            y_pred = self.predict_proba(X) 
            return log_loss(y, y_pred) 
        if self.scoring=='mse': 
            y_pred = self.predict(X) 
            return mean_squared_error(y, y_pred) 
    
