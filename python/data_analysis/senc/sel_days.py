from importlib import reload
import inspect, sys
import gc 
import numpy as np 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.feature_selection import f_classif, SelectPercentile, SelectFpr, SelectFdr, SelectFwe
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler

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

import senc.utils
reload(senc.utils)
from senc.utils import * 
from senc.plot_utils import * 
from senc.statistics import * 

def set_scaler(clf, scaling):
    
    if scaling=='standard':
        pipe = Pipeline([('scaler', StandardScaler()), ('filter', SelectFpr(f_classif, alpha=.01) ), ('clf', clf)]) 
    elif scaling=='center':
        pipe = Pipeline([('scaler', StandardScaler(with_std=False)), ('clf', clf)]) 
    elif scaling=='robust':
        pipe = Pipeline([('scaler', RobustScaler()), ('clf', clf)]) 
    else:
        pipe = Pipeline([('clf', clf)])
    
    return pipe

def sel_time(**options): 
    
    data.get_days() # do not delete that !!
    
    X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'], stimulus=options['stimulus'], task=options['task'], trials=options['trials']) 
    # that must be before bins 
    X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                         scaler=options['scaler_BL'],
                                         center=options['center_BL'],
                                         scale=options['scale_BL'],
                                         avg_mean=options['avg_mean_BL'],
                                         avg_noise=options['avg_noise_BL'],
                                         unit_var=options['unit_var']) 
    
    # print(options['epochs']) 
    # X_S1 = pp.avg_epochs(X_S1, options['epochs']) 
    # X_S2 = pp.avg_epochs(X_S2, options['epochs']) 
    
    X_S1_S2 = np.vstack( ( X_S1, X_S2 ) )
    y = np.hstack((np.zeros(X_S1.shape[0]), np.ones(X_S2.shape[0]) ))

    print('task', options['task'], 'day', options['day'], 'X_S1', X_S1.shape, 'X_S2', X_S2.shape, 'X_S1_S2', X_S1_S2.shape) 
    
    pipe = set_scaler(options['clf'], scaling='standard') 
    
    coefs = np.zeros( (X_S1.shape[-1], X_S1.shape[-2]) ) 
    scores = []
    non_zeros = []
    
    for i_epoch in range(X_S1.shape[-1]): 
        X = X_S1_S2[...,i_epoch] 
        pipe.fit(X, y) 

        print(X.shape)
        coef = pipe.named_steps['clf'].coef_[0] 
        print('coef', coef.shape) 
        # coefs[i_epoch, pval<=0.05] = coef 
        
        pval = pipe.named_steps['filter'].pvalues_ 
        non_zero = np.sum(pval<=0.001) / X_S1.shape[-2] * 100.0 
        print('pval', pval.shape, 'non_zero', non_zero ) 
        
        non_zeros.append(non_zero) 
        
        score = outer_cv(pipe, X, y, n_out=10, folds='stratified',
                         outer_score=options['outer_score'], 
                         inner_score=options['inner_score'], 
                         random_state=None) 
        
        scores.append(score) 
    
    # print(coefs.shape)
    scores = np.array(scores) 
    # sel = cos_between(coefs[0], coefs[1])
    
    sel = scores 
    print(sel) 
    
    # model = SelectPercentile(f_classif, percentile=10).fit(X_S1_S2, y) 
    # p_value = model.pvalues_ 
    
    # # _, p_value = f_classif(X_S1_S2, y)
    # print('pval', p_value.shape, np.sum(p_value<=0.05) )

    # idx_sel = p_value<=0.05
    
    # X_S1 = X_S1[:, idx_sel, :] 
    # X_S2 = X_S2[:, idx_sel, :] 
    
    # df = pd.DataFrame(X_S1_S2)
    # print(df.shape)
    # corr = df.corr()
    # print(corr.shape)
    
    # # sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool),
    # #             cmap=sns.diverging_palette(220, 10, as_cmap=True),
    # #             square=True)
    
    # columns = np.full((corr.shape[0],), True, dtype=bool)
    # for i in range(corr.shape[0]):
    #     for j in range(i+1, corr.shape[0]):
    #         if corr.iloc[i,j] >= 0.9:
    #             if columns[j]:
    #                 columns[j] = False
    # selected_columns = df.columns[columns] 
    # df_sel = df[selected_columns]
    
    # print(df_sel.shape)
    
    # print('task', options['task'], 'day', options['day'], 'X_S1', X_S1.shape, 'X_S2', X_S2.shape, 'X_S1_S2', X_S1_S2.shape) 
    
    # options['bins'] = None 
    # sel = get_sel(X_S1, X_S2, **options)[-1] 
    # print('sel', sel.shape, sel) 
    
    return sel 

if __name__ == '__main__':
    
    kwargs = dict() 
    
    kwargs['T_WINDOW'] = 0.5 
    kwargs['bins'] = 'ED' 
    kwargs['sample'] = 'S1' 
    
    kwargs['ci'] = 0
    kwargs['n_samples'] = 1000 
    kwargs['shuffle'] = 0 
    kwargs['n_shuffles'] = 1000 
    
    kwargs['pval']= .05 # .05, .01, .001 
    
    kwargs['scaler'] = 'standard' #'standard' # if not standardized gives strange results for norm 
    kwargs['scaler_BL'] = 'robust' 
    kwargs['avg_mean_BL'] = 0 
    kwargs['avg_noise_BL'] = 1
    
    kwargs['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo', 'Dual']) 
    # kwargs['tasks'] = ['DPA', 'Dual'] 
    
    kwargs['fold_type'] = 'stratified' 
    
    if(len(sys.argv)>1): 
        kwargs['i_mice'] = int(sys.argv[1]) 
        kwargs['task'] = sys.argv[2] 
        kwargs['day'] = sys.argv[3] 
        kwargs['trials'] = sys.argv[4] 
        kwargs['obj'] = sys.argv[5] 
        kwargs['stimulus'] = sys.argv[6] 
        kwargs['sample'] = sys.argv[7] 
    
    kwargs['n_alpha'] = 20 
    kwargs['alpha'] = 0.5 
    kwargs['n_lambda'] = 20 
    kwargs['alphas'] = np.linspace(0, 1, kwargs['n_alpha']) 
    kwargs['lbds'] = np.exp(np.linspace(-4, 4, kwargs['n_lambda']) ) 
    
    kwargs['out_fold'] = 'stratified' 
    kwargs['n_out'] = 10 
    kwargs['in_fold'] = 'stratified' 
    kwargs['n_in'] = 10 
    
    kwargs['clf_name'] = 'LogisticRegressionCV' 
    
    kwargs['param_grid'] = dict(lbd=kwargs['lbds']) # this is important if using pipeline 
    
    if kwargs['clf_name']=='logitnetCV':
        kwargs['param_grid'] = dict(clf__alpha=kwargs['alphas']) # this is important if using pipeline 
    
    kwargs['lbd'] = 'lambda_min' 
    kwargs['fold_type'] = 'stratified' 
    kwargs['inner_score']= 'neg_log_loss' 
    kwargs['outer_score']= 'roc_auc' 
    
    kwargs['n_splits'] = 10 
    kwargs['scoring'] = 'accuracy' 
    
    options = set_options(**kwargs) 
    set_globals(**options) 
    
    options['clf'] = get_clf(**options) 
    print('clf', options['clf']) 
    
    print(options['task'])
    
    options['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo', 'Dual']) 
    # options['tasks'] = np.array(['DPA','Dual']) 
    options['bins'] = 'ED' 
    options['stimulus']= 'sample' 
    options['i_task'] = np.argwhere(options['tasks']==options['task'])[0][0] 
    print(options['tasks'], options['task'], options['i_task']) 

    options['epochs'] = ['ED', 'MD', 'LD'] 
    sel = [] 
    # for i_day in range(2):
    
    options['day'] = 'all' 
    options['trials'] = 'correct' 
    sel = sel_time(**options) 
    
    # options['day'] = 'all' 
    # options['trials'] = 'incorrect' 
    # sel.append(sel_time(**options)) 
    
    plt.plot(sel, '-') 
    plt.ylabel('cosine')
    
    
