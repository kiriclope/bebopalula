from .libs import * 
from . import constants as gv 
from . import preprocessing as pp 
from . import featureSel as fs 
from . import progressbar as pg

from joblib import Parallel, delayed
from meegkit.detrend import detrend

def get_epochs():
    if gv.EDvsLD: 
        gv.epochs = ['ED', 'MD', 'LD'] 
        print('angle btw ED and other epochs') 
    else: 
        gv.epochs = ['Stim', 'ED', 'MD', 'LD'] 
        print('angle btw STIM and other epochs') 

def get_n_trials():
    if gv.mouse in [gv.mice[0]]: 
        gv.n_trials = 40 
    elif gv.mouse in [gv.mice[4]] or gv.mouse in [gv.mice[5]]: 
        gv.n_trials = 64*5  
    else:
        gv.n_trials = 32
            
def get_days():
    if (gv.mouse=='JawsM15') | (gv.mouse=='ChRM04') | (gv.mouse=='JawsM18'): 
        gv.days = np.arange(1, gv.n_days+1) 
    else:
        gv.days = [1,2,3,4,5]
        
def get_delays_times():
    if((gv.mouse=='ChRM04') | ('Jaws' in gv.mouse) | ('ACC' in gv.mouse) ): 
        gv.t_ED = [3, 4.5] 
        gv.t_MD = [5.5, 6.5] 
        gv.t_LD = [7.5, 9] 
    else: 
        gv.t_ED = [3, 6]
        gv.t_MD = [7, 8]
        gv.t_LD = [9, 12]

def get_stimuli_times():
    if gv.mouse=='C57_2_DualTask' : 
        gv.t_DIST = [6, 7] 
        gv.t_cue = [8, 8.5] 
        gv.t_DRT_reward = [8.5, 9] 
        gv.t_test = [12, 13] 
    else:
        gv.t_DIST = [4.5, 5.5] 
        gv.t_cue = [6.5, 7] 
        gv.t_DRT_reward = [7, 7.5] 
        gv.t_test = [9, 10] 
    
def get_sessions_mouse(): 
    if gv.mouse=='C57_2_DualTask' : 
        gv.sessions = list( map( str, np.arange(20200116, 20200121) ) ) 
    else : 
        gv.sessions = list( map( str, np.arange(20200521, 20200527) ) ) 
        
    gv.session=gv.sessions[gv.day-1] 
    
# def get_fluo_data(): 
                
#     if gv.mouse=='C57_2_DualTask' : 
#         data = loadmat(gv.path +  '/data/' + gv.mouse +  '/' + gv.session + '-C57-2-DualTaskAcrossDaySameROITrace' + '.mat') 
#         data_labels = loadmat(gv.path +  '/data/' + gv.mouse + '/' + gv.session + '-C57-2-DualTask-SumFluoTraceFile' + '.mat') 
        
#         X_data = np.rollaxis(data['SameAllCdf'],2,0)
#         # print('same neurons accross days')
#         # X_data = np.rollaxis(data['SamedFF0'],2,0) 
        
#         y_labels= data_labels['AllFileEvents'+gv.session][0][0][0].transpose() 
#         gv.frame_rate = 7.5
        
#     else: 
#         if gv.SAME_DAYS: 
#             print('same neurons accross days') 
#             if 'ACC' in gv.mouse:
#                 print(gv.path + '/data/' + gv.mouse + '/SamedROI/' + gv.mouse + '_day_' + str(gv.day) + '.mat' )
#                 data = loadmat(gv.path + '/data/' + gv.mouse + '/SamedROI/' + gv.mouse + '_day_' + str(gv.day) + '.mat' ) 
#             else:
#                 data = loadmat(gv.path + '/data/' + gv.mouse + '/SamedROI_0%dDays/' % len(gv.days) + gv.mouse + '_day_' + str(gv.day) + '.mat' ) 
#         else: 
#             data = loadmat(gv.path + '/data/' + gv.mouse + '/' + gv.mouse +'_day_' + str(gv.day) + '.mat') 
            
#         if 'raw' in gv.data_type:
#             # print('raw') 
#             X_data = np.rollaxis(data['Cdf_Mice'],1,0) 
#         else: 
#             # print('dF') 
#             X_data = np.rollaxis(data['dff_Mice'],1,0) 
            
#         y_labels = data['Events'].transpose() 
#         gv.frame_rate = 6 
        
#     # print(data.keys())    
#     # print(y_labels[4])
    
#     gv.duration = X_data.shape[2]/gv.frame_rate 
#     gv.time = np.linspace(0, gv.duration, X_data.shape[2]); 
#     gv.bins = np.arange(0, len(gv.time)) 
#     gv.n_neurons = X_data.shape[1] 
#     gv.trial_size = X_data.shape[2] 
    
#     get_stimuli_times() 
#     get_delays_times() 
#     get_n_trials() 
#     get_bins() 
    
#     print('mouse', gv.mouse, 'day', gv.day, 'type', gv.data_type, 'all data: X', X_data.shape,'y', y_labels.shape) 
#     # print(gv.n_neurons)
    
#     return X_data, y_labels

def get_fluo_data():
    
    get_days()
    print('days', gv.days)
    
    if gv.SAME_DAYS: 
        # print('same neurons accross days') 
        if 'ACC' in gv.mouse:
            # print(gv.path + '/data/' + gv.mouse + '/SamedROI/' + gv.mouse + '_day_' + str(gv.day) + '.mat' )
            data = loadmat(gv.path + '/data/' + gv.mouse + '/SamedROI/' + gv.mouse + '_day_' + str(gv.day) + '.mat' ) 
        else:
            data = loadmat(gv.path + '/data/' + gv.mouse + '/SamedROI_0%dDays/' % gv.n_days + gv.mouse + '_day_' + str(gv.day) + '.mat' )
                
    else: 
        data = loadmat(gv.path + '/data/' + gv.mouse + '/' + gv.mouse +'_day_' + str(gv.day) + '.mat') 
            
    if 'raw' in gv.data_type:
        # print('raw') 
        X_data = np.rollaxis(data['Cdf_Mice'],1,0) 
    else: 
        # print('dF') 
        X_data = np.rollaxis(data['dff_Mice'],1,0) 
            
    y_labels = data['Events'].transpose()
    
    if 'ACC' in gv.mouse:
        # print('mouse', gv.mouse, 'days', gv.days, 'type', gv.data_type, 'all data: X', X_data.shape,'y', y_labels.shape)
        
        X_data = X_data.reshape( (5, int(X_data.shape[0]/5), X_data.shape[1], X_data.shape[2]) )
        y_labels = y_labels.T.reshape( (5, int(y_labels.T.shape[0]/5), y_labels.T.shape[1]) ) 

        X_data = X_data[gv.day-1]
        y_labels = y_labels[gv.day-1].T 
        
    gv.frame_rate = 6 
    
    # print(data.keys()) 
    # print(y_labels[4])
    
    gv.duration = X_data.shape[2]/gv.frame_rate 
    gv.time = np.linspace(0, gv.duration, X_data.shape[2]); 
    gv.bins = np.arange(0, len(gv.time)) 
    gv.n_neurons = X_data.shape[1] 
    gv.trial_size = X_data.shape[2] 
    
    get_stimuli_times() 
    get_delays_times() 
    get_n_trials() 
    get_bins() 
    
    print('mouse', gv.mouse, 'day', gv.day, 'type', gv.data_type, 'all data: X', X_data.shape,'y', y_labels.shape) 
    # print(gv.n_neurons)
    # print(y_labels[4])
    return X_data, y_labels

def correct_incorrect_labels(y_labels):
    bool_correct = ( y_labels[2]==1 ) | ( y_labels[2]==4 ) 
    bool_incorrect = ( y_labels[2]==2 ) | ( y_labels[2]==3 ) 
    
    y_correct = np.argwhere( bool_correct ).flatten()
    y_incorrect = np.argwhere( bool_incorrect ).flatten()
    
    return y_correct, y_incorrect

def which_trials(y_labels): 
    y_trials = []
    
    if gv.laser_on:
        bool_ND = (y_labels[4]==0) & (y_labels[8]!=0) 
        bool_D1 = (y_labels[4]==13) & (y_labels[8]!=0) 
        bool_D2 = (y_labels[4]==14) & (y_labels[8]!=0) 
    else: 
        bool_ND = (y_labels[4]==0) & (y_labels[8]==0) 
        bool_D1 = (y_labels[4]==13) & (y_labels[8]==0) 
        bool_D2 = (y_labels[4]==14) & (y_labels[8]==0) 
 
    bool_S1 = (y_labels[0]==17) 
    bool_S2 = (y_labels[0]==18) 
    
    bool_T1 = (y_labels[1]==11) 
    bool_T2 = (y_labels[1]==12) 
    
    bool_pair = ( (y_labels[0]==17) & (y_labels[1]==11) ) | ( (y_labels[0]==18) & (y_labels[1]==12) ) 
    bool_unpair = ( (y_labels[0]==17) & (y_labels[1]==12) ) | ( (y_labels[0]==18) & (y_labels[1]==11) ) 
    
    if 'all' in gv.trial: 
        bool_trial = bool_ND & bool_D1 & bool_D2 
    if 'ND' in gv.trial: 
        bool_trial = bool_ND 
    if 'D1' in gv.trial: 
        bool_trial = bool_D1 
    if 'D2' in gv.trial: 
        bool_trial = bool_D2 
        
    if 'incorrect' in gv.trial:
        bool_incorrect = ( y_labels[2]==2 ) | ( y_labels[2]==3 ) 
        bool_trial = bool_trial & bool_incorrect 
        
    elif 'correct' in gv.trial: 
        bool_correct = ( y_labels[2]==1 ) | ( y_labels[2]==4 ) 
        bool_trial = bool_trial & bool_correct 
        
    if 'S1' in gv.trial: 
        y_trials = np.argwhere( bool_trial & bool_S1 ).flatten() 
    elif 'S2' in gv.trial: 
        y_trials = np.argwhere( bool_trial & bool_S2 ).flatten() 
        
    elif 'pair' in gv.trial: 
        y_trials = np.argwhere( bool_trial & bool_pair ).flatten() 
    elif 'unpair' in gv.trial:
        y_trials = np.argwhere( bool_trial & bool_unpair ).flatten() 
        
    # elif 'D1' in gv.trial:
    #     y_trials = np.argwhere( bool_D1 ).flatten() 
    # elif 'D2' in gv.trial: 
    #     y_trials = np.argwhere( bool_D2 ).flatten()
    
    elif 'T1' in gv.trial: 
        y_trials = np.argwhere( bool_trial & bool_T1 ).flatten() 
    elif 'T2' in gv.trial: 
        y_trials = np.argwhere( bool_trial & bool_T2 ).flatten() 
    else:
        y_trials = np.argwhere( bool_trial & bool_S1 & bool_S2 ).flatten()
    
    # print(gv.trial, y_trials) 
    
    return y_trials 

def get_pair_trials(X_data, y_labels): 
    
    trial = gv.trial 
    gv.trial = trial + "_pair" 
    y_pair = which_trials(y_labels) 
    
    trial = gv.trial 
    gv.trial = trial + "_unpair" 
    y_unpair = which_trials(y_labels) 
    
    gv.trial = trial 
    X_pair = X_data[y_pair] 
    X_unpair = X_data[y_unpair] 
    
    print('X_pair', X_pair.shape, 'X_unpair', X_unpair.shape) 
    return X_pair, X_unpair 
    
def get_S1_S2_trials(X_data, y_labels): 
    
    trial = gv.trial
    gv.trial = trial + "_S1" 
    y_S1 = which_trials(y_labels) 
    
    gv.trial = trial + "_S2" 
    y_S2 = which_trials(y_labels) 
    
    gv.trial = trial 
    X_S1 = X_data[y_S1] 
    X_S2 = X_data[y_S2] 
    
    # print('X_S1', X_S1.shape, 'X_S2', X_S2.shape)
    
    return X_S1, X_S2 

def get_X_S1_S2(X_data, y_labels, sample='sample'):
    
    X_S1_S2 = np.zeros( (len(gv.trials), len(gv.samples), int(gv.n_trials/len(gv.samples)), gv.n_neurons, gv.trial_size) ) 
    y_trials = np.zeros( (len(gv.trials), 3, len(gv.samples), int(gv.n_trials/len(gv.samples)) ) )
    
    if sample=='sample':
        dum_str = 'S'
    else:
        dum_str = 'T'
        
    _trial = gv.trial
    print(X_S1_S2.shape, X_data.shape)
    
    for i_trial, gv.trial in enumerate(gv.trials):
        
        trial = gv.trial
        gv.trial = trial + "_%s1" % dum_str
        y_S1 = which_trials(y_labels) 
        
        gv.trial = trial + "_%s2" % dum_str
        y_S2 = which_trials(y_labels) 

        gv.trial = trial + "_%s1_correct" % dum_str
        y_S1_correct = which_trials(y_labels) 

        gv.trial = trial + "_%s1_incorrect" % dum_str
        y_S1_incorrect = which_trials(y_labels) 
        
        gv.trial = trial + "_%s2_correct" % dum_str
        y_S2_correct = which_trials(y_labels) 

        gv.trial = trial + "_%s2_incorrect" % dum_str
        y_S2_incorrect = which_trials(y_labels) 
        
        y_S1_S2 = np.stack((y_S1, y_S2))
        
        y_correct = np.empty((2, y_S1_S2.shape[1])) * np.nan 
        y_correct[0,0:y_S1_correct.shape[0]] = y_S1_correct            
        y_correct[1,0:y_S2_correct.shape[0]] = y_S2_correct
        
        y_incorrect = np.empty((2, y_S1_S2.shape[1])) * np.nan 
        y_incorrect[0,0:y_S1_incorrect.shape[0]] = y_S1_incorrect 
        y_incorrect[1,0:y_S2_incorrect.shape[0]] = y_S2_incorrect 
        
        y_trials[i_trial] = np.stack( (y_S1_S2, y_correct, y_incorrect) ) 
        
        gv.trial = trial 
        X_S1 = X_data[y_S1] 
        X_S2 = X_data[y_S2] 
        
        X_S1_S2[i_trial, 0, 0:X_S1.shape[0]] = X_S1 
        X_S1_S2[i_trial, 1, 0:X_S2.shape[0]] = X_S2 
        
    # print('X_S1', X_S1.shape, 'X_S2', X_S2.shape, 'X_S1_S2', X_S1_S2.shape) 
    gv.trial = _trial 
    
    return X_S1_S2, y_trials

def get_trial_types(X_S1_trials): 
    gv.n_trials = 2*X_S1_trials.shape[0] 
    gv.trial_type = ['ND'] * gv.n_trials + ['D1'] * gv.n_trials + ['D2'] * gv.n_trials 

def get_bins():

    if(gv.T_WINDOW==0): 
        gv.bins_BL = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_BL[0]) and (gv.time[bin]<=gv.t_BL[1])] 
    
        gv.bins_STIM = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_STIM[0]) and (gv.time[bin]<=gv.t_STIM[1]) ] 
    
        gv.bins_ED = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_ED[0]) and (gv.time[bin]<=gv.t_ED[1]) ] 
        
        gv.bins_DIST = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_DIST[0]) and (gv.time[bin]<=gv.t_DIST[1]) ]
        
        gv.bins_MD = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_MD[0]) and (gv.time[bin]<=gv.t_MD[1]) ]
        
        gv.bins_LD = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_LD[0]) and (gv.time[bin]<=gv.t_LD[1]) ] 
        
        gv.bins_cue = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_cue[0]) and (gv.time[bin]<=gv.t_cue[1]) ] 

        gv.bins_DRT_rwd = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_DRT_reward[0]) and (gv.time[bin]<=gv.t_DRT_reward[1]) ] 

        gv.bins_test = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_test[0]) and (gv.time[bin]<=gv.t_test[1]) ]
        
    else:
        gv.bins_BL = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_BL[0]) and (gv.time[bin]<=gv.t_BL[1]) ] 
    
        gv.bins_STIM = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_STIM[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_STIM[1]+gv.T_WINDOW) ] 
    
        gv.bins_ED = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_ED[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_ED[1]+gv.T_WINDOW) ] 
        
        gv.bins_DIST = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_DIST[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_DIST[1]) ] 
        
        gv.bins_MD = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_MD[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_MD[1]+gv.T_WINDOW) ]
        
        gv.bins_LD = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_LD[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_LD[1]+gv.T_WINDOW) ] 
        
        gv.bins_cue = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_cue[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_cue[1]+gv.T_WINDOW) ] 

        gv.bins_DRT_rwd = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_DRT_reward[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_DRT_reward[1]+gv.T_WINDOW) ] 

        gv.bins_test = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_test[1]-gv.T_WINDOW) and (gv.time[bin]<=gv.t_test[1]+gv.T_WINDOW) ] 
        

    gv.bins_delay = [ bin for bin in gv.bins if (gv.time[bin]>=gv.t_ED[0]) and (gv.time[bin]<=gv.t_LD[1]) ]     
    gv.t_delay = [ gv.time[bin] for bin in gv.bins if (gv.time[bin]>=gv.t_ED[0]) and (gv.time[bin]<=gv.t_LD[1]) ] 
    
    # gv.bins_STIM_ED_MD_LD =  np.hstack( (gv.bins_STIM, gv.bins_ED, gv.bins_MD, gv.bins_LD) ) 
    # gv.t__STIM_ED_MD_LD = gv.time[gv.bins_STIM_ED_MD_LD] 
    
    # gv.bins_ED_MD_LD = np.hstack( (gv.bins_ED, gv.bins_MD, gv.bins_LD) )
    # gv.t_ED_MD_LD = gv.time[gv.bins_ED_MD_LD] 
    gv.bins_epochs = np.array( [gv.bins_ED, gv.bins_MD, gv.bins_LD] ) 
    
def get_X_y_epochs(X_S1_trials, X_S2_trials): 

    X_S1 = [] 
    X_S2 = [] 
    
    if 'all' in gv.epochs :
        X_S1=X_S1_trials
        X_S2=X_S2_trials

    if 'BL' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_BL-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_BL-gv.bin_start],axis=2)) 

    if 'STIM' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_STIM-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_STIM-gv.bin_start],axis=2)) 

    if 'ED' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_ED-gv.bin_start],axis=2))
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_ED-gv.bin_start],axis=2))

    if 'DIST' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_DIST-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_DIST-gv.bin_start],axis=2)) 

    if 'MD' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_MD-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_MD-gv.bin_start],axis=2)) 
        
    if 'LD' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_LD-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_LD-gv.bin_start],axis=2)) 

    if 'Cue' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_cue-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_cue-gv.bin_start],axis=2))
        
    if 'DRT_rwd' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_DRT_rwd-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_DRT_rwd-gv.bin_start],axis=2))

    if 'Test' in gv.epochs: 
        X_S1.append(np.mean(X_S1_trials[:,:,gv.bins_test-gv.bin_start],axis=2)) 
        X_S2.append(np.mean(X_S2_trials[:,:,gv.bins_test-gv.bin_start],axis=2))
        
    X_S1 = np.asarray(X_S1)  
    X_S2 = np.asarray(X_S2) 
    
    if 'all' in gv.epochs :
        X = np.concatenate([X_S1, X_S2], axis=0) 
    else: 
        X = np.concatenate([X_S1, X_S2], axis=1) 
        X = np.rollaxis(X,2,1).transpospe() 
        
    y_S1 = np.repeat(0, int(X_S1_trials.shape[0])) 
    y_S2 = np.repeat(1, int(X_S2_trials.shape[0]))
    
    y = np.concatenate((y_S1, y_S2)) 
    
    return X, y 

def get_X_y_trials(X_data, y_labels):
    X_S1_trials, X_S2_trials = get_S1_S2_trials(X_data, y_labels) 
    # print(y_S1_trials[0]) 
    
    X_S1 = bin_data(X_S1_trials, gv.n_bin, gv.n_bin) 
    X_S2 = bin_data(X_S2_trials, gv.n_bin, gv.n_bin) 
    
    # print(X_S1[0]) 
    X_trials = np.concatenate([X_S1,X_S2],axis=0) 
    # print(X_trials[0])
    y_S1 = np.repeat(0, int(X_S1_trials.shape[0]))
    y_S2 = np.repeat(1, int(X_S2_trials.shape[0]))

    y_trials = np.concatenate((y_S1, y_S2))
    
    return X_trials, y_trials

def get_X_y_mouse_session(synthetic=False):

    X, y = get_fluo_data() 
    print('mouse', gv.mouse, 'session', gv.session, 'data X', X.shape,'y', y.shape) 
    
    # compute DF over F0 
    if gv.F0_THRESHOLD is not None: 
        X = pp.dFF0_remove_silent(X) 
        gv.n_neurons = X.shape[1] 

    if gv.DECONVOLVE: 
        X = pp.deconvolveFluo(X) 
        
    if gv.mouse in [gv.mice[0]]: 
        if 'ND_D1' in gv.trials: 
            gv.n_trials = 40*2 
        else: 
            gv.n_trials = 40 
    else: 
        if 'ND_D1' in gv.trials: 
            gv.n_trials = 32*2 
        else: 
            gv.n_trials = 32
            
    X_trials = np.zeros( (len(gv.trials), len(gv.samples), int(gv.n_trials/len(gv.samples)), gv.n_neurons, gv.trial_size) ) 
        
    mins = [] 
    for n_trial, gv.trial in enumerate(gv.trials): 
        
        if gv.pair_trials:
            X_S1, X_S2 = get_pair_trials(X, y) 
        else:
            X_S1, X_S2 = get_S1_S2_trials(X, y)

        get_trial_types(X_S1) 
        
        if not gv.DECONVOLVE:
            # compute DF over F0 
            # X_S1 = pp.dFF0(X_S1) 
            # X_S2 = pp.dFF0(X_S2) 
            
            if gv.SAVGOL: 
                for trial in range(X_S2.shape[0]): 
                    X_S1[trial] = savgol_filter(X_S1[trial], int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = 5, deriv=0) 
                    X_S2[trial] = savgol_filter(X_S2[trial], int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = 5, deriv=0) 
                    
            if gv.DETREND :
                X_S1 = pp.detrend_X(X_S1,7)
                X_S2 = pp.detrend_X(X_S2,7)
            
            if gv.Z_SCORE : 
                for trial in range(min_S1_S2): 
                    X_S1[trial] = pp.z_score(X_S1[trial]) 
                    X_S2[trial] = pp.z_score(X_S2[trial]) 

        if X_S1.shape[0]!=X_trials.shape[3] or X_S2.shape[0]!=X_trials.shape[2]: 
            print('X_trials', X_trials.shape[2], 'X_S1', X_S1.shape[0], 'X_S2', X_S2.shape[0]) 
            min_S1_S2 = np.amin([X_S1.shape[0], X_S2.shape[0]]) 
        else: 
            min_S1_S2 = X_trials.shape[2] 
            
        mins.append(min_S1_S2)
        
        # if gv.DECONVOLVE:
        #     X_S1 = pp.deconvolveFluo(X_S1) 
        #     X_S2 = pp.deconvolveFluo(X_S2) 
            
        X_trials[n_trial, 0, 0:X_S1.shape[0], 0:X_S1.shape[1]] = X_S1 
        X_trials[n_trial, 1, 0:X_S2.shape[0], 0:X_S2.shape[1]] = X_S2 
        print('X_trials', X_trials.shape[3], 'X_S1', X_S1.shape[1], 'X_S2', X_S2.shape[1]) 
        
            
    mins = np.amin(mins) 
    y = np.array([np.zeros(mins), np.ones(mins)]).flatten() 
    X_trials = X_trials[:,:,0:mins]
    
    
    print('X_trials', X_trials.shape, 'y', y.shape) 
    
    return X_trials, y 

