import multiprocessing 
from . import constants as gv 
import numpy as np

def set_globals(**opts):

    gv.code = opts['code']

    gv.laser_on = opts['laser'] 
    
    gv.num_cores = opts['n_jobs']
    gv.IF_SAVE = opts['IF_SAVE']
    gv.SYNTHETIC = opts['SYNTHETIC'] 
    gv.data_type = opts['type']

    gv.first_days = opts['firstDays']
    gv.last_days = opts['lastDays']
    gv.all_days = opts['allDays']
    
    gv.inner_scoring = opts['inner_score']
    
    # parameters 
    gv.mouse = gv.mice[opts['i_mice']] 
    
    # if opts['task']=='all': 
    #     gv.tasks = ['all'] 
    #     gv.task = ['all'] 
    # elif opts['task']=='Dual': 
    #     gv.tasks = ['Dual'] 
    #     gv.task = ['Dual']   
    # else:
    #     gv.tasks = opts['tasks'] 
    #     gv.task = gv.tasks[opts['i_task']] 
    
    gv.day = gv.days[opts['i_day']] 
    
    # if gv.code=='memory':
    gv.epochs = [ 'ED', 'MD', 'LD'] 
    # elif gv.code=='sensory':
    # gv.epochs = ['STIM', 'DIST', 'TEST'] 
    # else:
    gv.epochs = opts['epochs']
    
    gv.epoch = gv.epochs[opts['i_epoch']] 

    gv.n_days = opts['n_days'] 
    
    gv.SAME_DAYS = opts['same_days']
    # gv.cos_trials = opts['cos_trials']
    # gv.scores_trials = opts['scores_trials']
    # gv.inter_trials = opts['inter_trials'] 
    
    # gv.pal = ['#ff00ff','#ffff00','#00ffff'] 
    # gv.pal = ['black', 'dimgray', 'lightgray'] 
    
    # preprocessing
    gv.standardize = opts['scaler']
    
    gv.T_WINDOW = opts['T_WINDOW'] 

    # feature selection 
    gv.FEATURE_SELECTION = 0 
    gv.LASSOCV = 0 
        
    # bootstrap 
    gv.n_boots = opts['n_boots'] 
    gv.bootstrap_method = opts['bootstrap_method'] 
    
    # temporal decoder
    gv.fold_type = opts['fold_type']
    gv.n_iter = opts['n_iter']
    
    # classification parameters 
    gv.clf_name = opts['clf'] 
    # gv.scoring = opts['scoring'] 
    
    # dimensionality reduction 
    
    # PCA parameters
    gv.AVG_BEFORE_PCA = 1 
    gv.pca_model = opts['pca_model'] # PCA, sparsePCA, supervisedPCA or None
    gv.explained_variance = opts['exp_var']
    gv.n_components = opts['n_comp']
    gv.list_n_components = None 
    gv.inflection = opts['inflection']
    
    gv.sparse_alpha = 1 
    gv.ridge_alpha = .01
    
    gv.pca_method = opts['pca_method'] # 'hybrid', 'concatenated', 'averaged' or None
        
    gv.fix_alpha_lbd = opts['fix_alpha_lbd']
    
def set_options(**kwargs): 
    
    opts = dict()
    
    opts['obj'] = 'frac' # 'cos', 'norm'     
    opts['trial_type'] = 'correct'
    opts['trials'] = 'correct' 
    opts['stimulus'] = 'sample'
    opts['sample'] = 'S1' # S1, S2, or S1_S2  
    opts['n_samples'] = 1000 # for permutation test 
    opts['n_shuffles'] = 1000 # for permutation test 
    
    opts['ci']=1 
    opts['shuffle']=1 
    opts['perm_test']=1
    
    opts['mouse_name'] = ['Dumb', 'Alice', 'Bob', 'Charly', 'Dave', 'Eve', 'mPFC', 'ACC'] 
    opts['tasks'] = np.array(['DPA', 'DualGo', 'DualNoGo'])
    opts['code'] = 'memory' 
    opts['pval'] = .05 
    opts['verbose'] = 0 
    opts['type'] = 'raw' 
    opts['n_jobs'] = int(0.9*multiprocessing.cpu_count()) 
    opts['IF_SAVE'] = 1
    opts['add_vlines']=0
    opts['SYNTHETIC'] = 0
    
    opts['fix_alpha_lbd'] = 0
    opts['bins'] = 'ED'
    opts['Delta0']= None 
    
    opts['firstDays'] = 0 
    opts['lastDays'] = 0 
    opts['allDays'] = 0 
    
    opts['stim'] = 0
    opts['day'] = 'all'
    
    # globals 
    opts['i_mice'] = 1
    opts['i_day'] = -1
    opts['i_trial'] = 0  
    opts['i_epoch'] = 0
    opts['i_task'] = 0 
    opts['task'] = 'DPA' # DPA, DualGo, DualNoGo, Dual, or all 
    opts['n_days'] = 9
    
    opts['same_days'] = 1 
    opts['laser']=0 

    opts['feature_sel'] = 'lasso' # 'ttest_ind' or 'lasso' 
    # bootstrap
    opts['boots'] = False 
    opts['n_boots'] = int(1e3) 
    opts['bootstrap_method'] = 'block' # 'bayes', 'bagging', 'standard', 'block' or 'hierarchical' 
    opts['boot_cos'] = 0 
    opts['n_cos_boots'] = int(1e3) 
    
    # temporal decoder 
    opts['n_iter'] = 1 
    opts['fold_type'] = 'stratified' 
    
    # preprocessing parameters 
    opts['T_WINDOW'] = 0.5 
    
    opts['epochs'] = ['ED', 'MD', 'LD']
    # opts['epochs'] = ['STIM', 'DIST', 'LD'] 
    
    opts['savgol'] = 0 # sav_gol filter             
    opts['detrend'] = 0 # detrend the data 
    opts['order'] = 3
    
    opts['scaler_BL']= 'standard' # standard, robust, center
    opts['center_BL']= None 
    opts['scale_BL']= None 
    opts['avg_mean_BL']=0
    opts['avg_noise_BL']=1 
    
    opts['scaler']= 'standard' # standard, robust, center 
    opts['center']= None 
    opts['scale']= None 
    opts['return_center_scale'] = 0 
    
    opts['avg_mean']=0
    opts['avg_noise']=1 
    opts['unit_var']=0 
    
    # PCA parameters 
    opts['pca_model'] = None # PCA, sparsePCA, supervisedPCA or None
    opts['pca_method'] = 'hybrid' # 'hybrid', 'concatenated', 'averaged' or None
    opts['exp_var'] = 0.90 
    opts['n_comp'] = None
    opts['inflection'] = False 
    
    # classification parameters 
    opts['clf_name'] = 'logitnetCV' 
    opts['clf'] = None 
    
    # sklearn LogisticRegression, LogisticRegressionCV
    opts['random_state'] = None 
    opts['tol']=1e-4 
    opts['max_iter']= int(1e2) 
    
    opts['fit_intercept'] = True 
    opts['intercept_scaling'] = 1

    opts['off_diag'] = False 
    opts['t_train'] = 'ED' 
    
    opts['C']=1e2 
    opts['Cs'] = np.logspace(-4, 4, 10) 
    opts['penalty']='l2' 
    opts['solver']='liblinear' # liblinear or saga 
    opts['l1_ratios'] = np.linspace(0, 1, 10)         
    opts['l1_ratio'] = None 
    opts['param_grid'] = dict(clf__C=opts['Cs']) # this is important if using pipeline 
    # opts['param_grid'] = dict(clf__C=opts['Cs'], clf__l1_ratio=opts['l1_ratios']) # this is important if using pipeline                 
    
    opts['alpha'] = .5 
    opts['n_alpha'] = 10 
    opts['alphas'] = np.linspace(0, 1, opts['n_alpha']) 
    
    opts['lbd'] = 'lambda_min' 
    opts['n_lambda'] = 10 
    opts['lbds'] = np.exp(np.linspace(-10, 0, opts['n_lambda']) ) 
    opts['min_lambda_ratio'] = 1e-4 
    
    opts['standardize'] = True 
    opts['prescreen'] = True 
    
    if opts['clf_name']=='logitnet': 
        opts['param_grid'] = dict(lbd=opts['lbds'], alpha=opts['alphas']) 
    
    if opts['clf_name']=='logitnetCV': 
        opts['param_grid'] = dict(alpha=opts['alphas']) 
    
    opts['out_fold'] = 'repeated' 
    opts['n_out'] = 10 
    opts['outer_score']= 'roc_auc' 
    
    opts['in_fold'] = 'stratified' 
    opts['n_in'] = 10 
    opts['inner_score']= 'neg_log_loss' 
    
    opts.update(kwargs) 
    
    return opts 
