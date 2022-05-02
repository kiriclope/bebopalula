import numpy as np 

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr 
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler

def set_scaler(clf, scaling):
    
    if scaling=='standard':
        pipe = Pipeline([('scaler', StandardScaler()), ('filter', SelectFpr(f_classif, alpha=.05) ), ('clf', clf)]) 
    elif scaling=='center':
        pipe = Pipeline([('scaler', StandardScaler(with_std=False)), ('filter', SelectFpr(f_classif, alpha=.05) ), ('clf', clf)]) 
    elif scaling=='robust':
        pipe = Pipeline([('scaler', RobustScaler()), ('filter', SelectFpr(f_classif, alpha=.05) ), ('clf', clf)]) 
    else:
        pipe = Pipeline([('filter', SelectFpr(f_classif, alpha=.05) ), ('clf', clf)]) 
    
    return pipe

class logregFS():

    def __init__(Cs=10, fit_intercept=True, cv=None, dual=False, penalty='l2', scoring=None, solver='lbfgs', tol=0.0001, max_iter=100, class_weight=None, n_jobs=None, verbose=0, refit=True, intercept_scaling=1.0, multi_class='auto', random_state=None, l1_ratios=None):
                
    def fit(self, X, y): 

        clf = LogisticRegressionCV(solver=kwargs['solver'], penalty=kwargs['penalty'], l1_ratios=None, 
                                   tol=kwargs['tol'], max_iter=int(kwargs['max_iter']), scoring=kwargs['scoring'], 
                                   fit_intercept=kwargs['fit_intercept'], intercept_scaling=kwargs['intercept_scaling'], 
                                   cv=kwargs['n_splits'], n_jobs=None) 
        
