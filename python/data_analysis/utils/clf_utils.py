#!/usr/bin/env python3
import sys
import numpy as np

from sklearn.datasets import load_iris
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import KFold, StratifiedKFold 

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

def nested_cv(clf, X, y, param_grid, scaling='standard', n_in=10, n_out=10, n_jobs=-1):

    # pipeline with scaling 
    pipe = set_scaler(clf, scaling) 
    cv_inner = StratifiedKFold(n_splits=n_in, shuffle=True, random_state=None) # outer cv loop for scoring 
    grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv_inner, n_jobs=n_jobs ) # gridsearchcv for hyperparam tuning (inner cv loop) 
    
    # grid.fit(X, y) 
    # print(grid.best_params_) 
    
    cv_outer = StratifiedKFold(n_splits=n_out, shuffle=True, random_state=None) # outer cv loop for scoring 
    scores = cross_val_score(grid, X, y, cv=cv_outer, n_jobs=n_jobs)
    
    return np.mean(scores)

def outer_cv(clf, X, y, scaling='standard', n_splits=10, n_jobs=-1):
    
    pipe = set_scaler(clf, scaling) 
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=None) # outer cv loop for scoring 
    scores = cross_val_score(pipe, X, y, cv=folds) 
    
    return np.mean(scores)

def inner_cv(clf, X, y, param_grid, scaling='standard', n_splits=10, n_jobs=-1):
    
    # pipeline with scaling 
    pipe = set_scaler(clf, scaling) 
    cv_inner = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=None) # outer cv loop for scoring 
    grid = GridSearchCV(pipe, param_grid=param_grid, cv=cv_inner, n_jobs=n_jobs ) # gridsearchcv for hyperparam tuning (inner cv loop) 
    
    grid.fit(X, y) 
    
    # return grid.best_params_

if __name__ == '__main__':
    
    if(len(sys.argv)>1): 
        n_alphas = int(sys.argv[1]) 
    else:
        n_alphas = 10
    
    X, y = load_iris(return_X_y=True)
    X = X[:100] 
    y = y[:100] 
    
    print(X.shape, y.shape) 
        
    clf = KerasClassifier(build_fn=get_model, epochs=10, batch_size=1, verbose=0)
    # clf.fit(X,y)
    
    # # outer cv loop
    score = outer_cv(clf, X, y) 
    print('<outer_cv_score>', score) 
    
    alphas = np.linspace(0, 1, n_alphas) 
    param_grid = dict(clf__alpha=alphas) # this is important if using pipeline
    
    # inner cv loop
    score = inner_cv(clf, X, y, param_grid) 
    print('<inner_cv_score>', score) 
    
    # # # nested cv loop with hyperparam tuning 
    # score = nested_cv(clf, X, y, param_grid) 
    # print('<nested_cv_score>', score) 
    
