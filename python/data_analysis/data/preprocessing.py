from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from scipy.signal import savgol_filter, butter, filtfilt
from scipy.ndimage import uniform_filter1d

from sklearn.feature_selection import SelectKBest, chi2, VarianceThreshold, f_regression, mutual_info_classif, f_classif

from joblib import Parallel, delayed, parallel_backend
from meegkit.detrend import detrend
from oasis.functions import deconvolve

from .libs import * 
from . import constants as gv
from . import progressbar as pg
from . import featureSel as fs

from dim_red.pca.pca_decomposition import pca_methods 

def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y 

def get_X_pca(X, **options):
    
    my_pca = pca_methods(pca_model = options['pca_model'],
                         pca_method = options['pca_method'],
                         n_components = options['n_comp'],
                         total_explained_variance = options['exp_var'],
                         inflection = options['inflection'],
                         verbose=options['verbose'])  
    
    X_pca = my_pca.fit_transform(X) 
    gv.list_n_components = my_pca.list_n_components 

    return X_pca 

def center(X): 
    scaler = StandardScaler(with_mean=True, with_std=False) 
    Xc = scaler.fit_transform(X.T).T 
    return Xc 

def z_score(X): 
    scaler = StandardScaler()
    
    if X.ndim>3: # task x sample x trial x neurons x time 
        Xz = X
        for i_task in range(3): 
            for i_sample in range(2): 
                for i in range(X.shape[2]): 
                    Xt = X[i_task, i_sample, i] 
                    if gv.Z_SCORE: 
                        Xz[i_task, i_sample, i] = scaler.fit_transform(Xt.T).T 
                    elif gv.Z_SCORE_BL:
                        scaler.fit(Xt[...,gv.bins_BL].T) 
                        Xz[i_task, i_sample, i] = scaler.transform(Xt.T).T                 
                        
    elif X.ndim>2: # trial x neurons x time 
        Xz = X 
        for i in range(X.shape[0]): 
            Xt = X[i] 
            if gv.Z_SCORE: 
                Xz[i] = scaler.fit_transform(Xt.T).T 
            elif gv.Z_SCORE_BL:
                scaler.fit(Xt[...,gv.bins_BL].T)  
                Xz[i] = scaler.transform(Xt.T).T 
                
                # mean = np.mean(Xt[:,gv.bins_BL],axis=1) # mean over baseline bins 
                
                # std = np.std(Xt[:,gv.bins_BL],axis=1) # std over baseline bins 
                # Xz[i] = ( Xt-mean[:,np.newaxis] ) /std[:,np.newaxis]
                
                # Xz[i] = ( Xt-mean[:,np.newaxis] ) 
                         
                # for i_neuron in range(X.shape[1]): 
                #     Xn = X[:,i_neuron] 
                #     std = np.std(Xn[:,gv.bins_BL]) # std over baseline in all trials 
                #     Xz[i, i_neuron] = Xz[i, i_neuron]/std
    else: # trial x neurons 
        
        if gv.Z_SCORE: 
            Xz = scaler.fit_transform(X.T).T 
        elif gv.Z_SCORE_BL : 
            scaler.fit(X[...,gv.bins_BL].T) 
            Xz = scaler.transform(X.T).T
            
    return Xz 

# def z_score_trials(X):
#     scaler = StandardScaler() 
    
#     # average over trials
#     avg_trials = np.mean(X, axis=0)
#     avg_trials = avg_trials[np.newaxis, ...]
        
#     # average and std over baseline bins 
#     avg_bl = np.mean( avg_trials[..., gv.bins_BL], axis=-1 ) 
#     std_bl = np.std( avg_trials[..., gv.bins_BL], axis=-1 )
        
#     avg_bl = avg_bl[..., np.newaxis]
#     std_bl = std_bl[..., np.newaxis]
    
#     # standardize 
#     z_trials = ( X - avg_bl ) / std_bl
        
#     return z_trials 

def z_score_trials(X):

    print('X', X.shape) 
    # average over trials
    if X.ndim>3: # task x sample x trial x neurons x time 
        z_trials = np.zeros(X.shape) 
        for i_task in range(3):
            X_BL = X[i_task, ..., gv.bins_BL]
            X_BL = np.moveaxis(X_BL, 0, -1)
            print('X_BL', X_BL.shape)
            m = np.nanmean( X_BL, axis=-1) 

            X_task = np.vstack(X[i_task])          
            # X_task = X_task[..., gv.bins_BL] 
            X_task = np.hstack(X_task) 
            print('X_task', X_task.shape)             
            s = np.nanstd( X_task, axis=-1 ) 
            
            print('m', m.shape, np.mean(m), 's', s.shape, np.std(s)) 
            z_trials[i_task] = ( X[i_task] - m[..., np.newaxis] ) / ( s[np.newaxis, np.newaxis, :, np.newaxis] + gv.eps) 
        
    else: # trials x neurons x time 
        X_BL = X[..., gv.bins_BL]
        
        print('X_BL', X_BL.shape) 
        
        X_trials = np.hstack(X) 
        print('X_trials', X_trials.shape) 
        
        m = np.nanmean( X_BL, axis=-1) 
        s = np.nanstd( X_trials, axis=-1) 
        
        print('m', m.shape, np.mean(m), 's', s.shape, np.std(s)) 
        z_trials = ( X - m[..., np.newaxis] ) / ( s[np.newaxis, :, np.newaxis] + gv.eps) 
        # z_trials = ( X - m[np.newaxis, ..., np.newaxis] ) / ( s[np.newaxis, :, np.newaxis] + gv.eps) 
            
    return z_trials

def normalize(X):
    # Xmin = np.amin(X, axis=-1) 
    # Xmax = np.amax(X, axis=-1) 
    
    Xmin = np.amin(X[..., gv.bins_BL], axis=-1) 
    Xmax = np.amax(X[..., gv.bins_BL], axis=-1) 
    
    Xmin = Xmin[..., np.newaxis] 
    Xmax = Xmax[..., np.newaxis] 
    
    return (X-Xmin)/(Xmax-Xmin+gv.eps) 

def normalize_trials(X):
    avg_trials = np.mean(X, axis=0)
    avg_trials = avg_trials[np.newaxis, ...]
    
    Xmin = np.amin(avg_trials[..., gv.bins_BL], axis=-1) 
    Xmax = np.amax(avg_trials[..., gv.bins_BL], axis=-1) 
    
    Xmin = Xmin[..., np.newaxis] 
    Xmax = Xmax[..., np.newaxis] 
    
    return (X-Xmin)/(Xmax-Xmin+gv.eps) 

def conf_inter(y): 
    ci = []
    for i in range(y.shape[0]):
        ci.append( stats.t.interval(0.95, y.shape[1]-1, loc=np.mean(y[i,:]), scale=stats.sem(y[i,:])) )
    ci = np.array(ci).T

    return ci

def dFF0_remove_silent(X): 
    ''' N_trials, N_neurons, N_times '''
    X = normalize(X) 
    
    if gv.AVG_F0_TRIALS: 
        F0 = np.mean( np.mean(X[...,gv.bins_BL],axis=-1), axis=0) 
        F0 = F0[np.newaxis,:, np.newaxis] 
    else:
        F0 = np.mean(X[...,gv.bins_BL],axis=-1) 
        F0 = F0[..., np.newaxis] 
        
    # print('X', X.shape,  X[0,0, 0:3])
    # print('F0', F0.shape, F0[0,0,0:3]) 
    
    if gv.F0_THRESHOLD is not None: 
        # removing silent neurons 

        # idx = np.where(F0<=gv.F0_THRESHOLD)        
        # F0 = np.delete(F0, idx, axis=-2) 
        # X = np.delete(X, idx, axis=-2)

        idx = np.where(F0<=gv.F0_THRESHOLD)
        F0[idx] = np.nan 
        X[idx] = np.nan 
        
    # print('X', X.shape, 'F0', F0.shape)
    
    return (X-F0) / (F0 + gv.eps) 
    
def dFF0(X): 
    if not gv.AVG_F0_TRIALS: 
        F0 = np.mean(X[...,gv.bins_BL],axis=-1)        
        # F0 = np.percentile(X, 15, axis=-1) 
        F0 = F0[..., np.newaxis] 
    else: 
        F0 = np.mean( np.mean(X[...,gv.bins_BL],axis=-1), axis=0) 
        F0 = F0[np.newaxis,:, np.newaxis]
        
    return (X-F0) / (F0 + gv.eps) 

def dF(X): 
    if not gv.AVG_F0_TRIALS: 
        F0 = np.mean(X[...,gv.bins_BL],axis=-1)        
        # F0 = np.percentile(X, 15, axis=-1) 
        F0 = F0[..., np.newaxis] 
    else: 
        F0 = np.mean( np.mean(X[...,gv.bins_BL],axis=-1), axis=0) 
        F0 = F0[np.newaxis,:, np.newaxis]
        
    return (X-F0)

def findBaselineF0(rawF, fs, axis=0, keepdims=False): 
    """Find the baseline for a fluorescence imaging trace line.
    
    The baseline, F0, is the 5th-percentile of the 1Hz
    lowpass filtered signal.
    
    Parameters
    ----------
    rawF : array_like
        Raw fluorescence signal.
    fs : float
        Sampling frequency of rawF, in Hz.
    axis : int, optional
        Dimension which contains the time series. Default is 0.
    keepdims : bool, optional
        Whether to preserve the dimensionality of the input. Default is
        `False`.
    
    Returns
    -------
    baselineF0 : numpy.ndarray
        The baseline fluorescence of each recording, as an array.
    
    Note
    ----
    In typical usage, the input rawF is expected to be sized
    `(numROI, numTimePoints, numRecs)`
    and the output will then be sized `(numROI, 1, numRecs)`
    if `keepdims` is `True`.
    """
    
    rawF = np.moveaxis(rawF.T,0,1)
    print('#neurons x #time x #trials', rawF.shape)
    
    # Parameters --------------------------------------------------------------
    nfilt = 30  # Number of taps to use in FIR filter
    fw_base = 1  # Cut-off frequency for lowpass filter, in Hz
    base_pctle = 5  # Percentile to take as baseline value
    
    # Main --------------------------------------------------------------------
    # Ensure array_like input is a numpy.ndarray
    rawF = np.asarray(rawF)
    
    # Remove the first datapoint, because it can be an erroneous sample
    rawF = np.split(rawF, [1], axis)[1]
    
    if fs <= fw_base:
        # If our sampling frequency is less than our goal with the smoothing
        # (sampling at less than 1Hz) we don't need to apply the filter.
        filtered_f = rawF
        
    else:
        # The Nyquist rate of the signal is half the sampling frequency
        nyq_rate = fs / 2.0
        
        # Cut-off needs to be relative to the nyquist rate. For sampling
        # frequencies in the range from our target lowpass filter, to
        # twice our target (i.e. the 1Hz to 2Hz range) we instead filter
        # at the Nyquist rate, which is the highest possible frequency to
        # filter 
        cutoff = min(1.0, fw_base / nyq_rate) 
        
        # Make a set of weights to use with our taps.
        # We use an FIR filter with a Hamming window.
        b = scipy.signal.firwin(nfilt, cutoff=cutoff, window='hamming')
        
        # The default padlen for filtfilt is 3 * nfilt, but in case our
        # dataset is small, we need to make sure padlen is not too big
        padlen = min(3 * nfilt, rawF.shape[axis] - 1)
        
        # Use filtfilt to filter with the FIR filter, both forwards and
        # backwards. 
        filtered_f = scipy.signal.filtfilt(b, [1.0], rawF, axis=axis, padlen=padlen) 
        
    # Take a percentile of the filtered signal
    baselineF0 = np.percentile(filtered_f, base_pctle, axis=axis, keepdims=keepdims)

    baselineF0 = baselineF0.T
    baselineF0 = baselineF0[:,np.newaxis,:]
    return baselineF0


def bin_data(data, bin_step, bin_size):
    # bin_step number of pts btw bins, bin_size number of size in each bin
    bin_array = [np.mean(np.take(data,np.arange(int(i*bin_step),int(i*bin_step+bin_size)), axis=2), axis=2) for i in np.arange(data.shape[2]//bin_step-1)]
    bin_array = np.array(bin_array)
    bin_array = np.rollaxis(bin_array,0,3)
    return bin_array

def detrend_loop(X, trial, neuron, order):
    X_det, _, _ = detrend(X[trial, neuron], order)
    return X_det

def detrend_X(X, order=3):
    with pg.tqdm_joblib(pg.tqdm(desc='detrend X', total=int(X.shape[0]*X.shape[1]) ) ) as progress_bar: 
        dum = Parallel(n_jobs=gv.num_cores)(delayed(detrend_loop)(X, trial, neuron, order) 
                                            for trial in range(X.shape[0]) 
                                            for neuron in range(X.shape[1]) )
               
        X = np.asarray(dum).reshape(X.shape[0], X.shape[1], X.shape[2])
    return X

def detrend_data(X_trial, poly_fit=1, degree=7): 
    """ Detrending of the data, if poly_fit=1 uses polynomial fit else linear fit. """
    # X_trial : # neurons, # times 
    
    model = LinearRegression()
    fit_values_trial = []

    indexes = range(0, X_trial.shape[1]) # neuron index 
    values = np.mean(X_trial,axis=0) # mean fluo value 
    
    indexes = np.reshape(indexes, (len(indexes), 1))
    
    if poly_fit:
        poly = PolynomialFeatures(degree=degree) 
        indexes = poly.fit_transform(indexes) 
            
    model.fit(indexes, values)
    fit_values = model.predict(indexes) 
    fit_values_trial = np.array(fit_values)
    
    # for i in range(0, X_trial.shape[0]): # neurons 
    #     indexes = range(0, X_trial.shape[1]) # neuron index 
    #     values = X_trial[i] # fluo value 
                
    #     indexes = np.reshape(indexes, (len(indexes), 1))

    #     if poly_fit:
    #         poly = PolynomialFeatures(degree=degree) 
    #         indexes = poly.fit_transform(indexes) 

    #     model.fit(indexes, values)
    #     fit_values = model.predict(indexes) 
        
    #     fit_values_trial.append(fit_values) 
        
    # fit_values_trial = np.array(fit_values_trial)
    return fit_values_trial

def feature_selection(X, method='variance'):

    X_avg = np.mean(X[:,:,gv.bins_ED_MD_LD],axis=-1) 

    if 'variance' in method :
        idx = fs.featSel.var_fit_transform(X_avg, threshold=threshold) 
        X_avg = np.delete(X_avg, idx, axis=1) 
        X = np.delete(X, idx, axis=1)
            
    if 'mutual' in method:
        idx = fs.featSel.select_best(X_avg, y, percentage=1-threshold) 
        X_avg = np.delete(X_avg, idx, axis=1) 
        X = np.delete(X, idx, axis=1)
            
    if 'correlation' in method:
        idx = fs.featSel.select_indep(X_avg, threshold=threshold) 
        X_avg = np.delete(X_avg, idx, axis=1) 
        X = np.delete(X, idx, axis=1) 

def avg_epochs(X, y=None): 

    print('average over epochs', gv.epochs) 
    X_avg = np.mean(X, axis=-1) 
    X_epochs = np.empty( tuple([len(gv.epochs)])+ X_avg.shape ) 
    print('X', X_epochs.shape, 'X_avg', X_avg.shape)

    print('start', gv.bin_start, 'epochs', gv.epochs) 
    
    for i_epoch, epoch in enumerate(gv.epochs):
        if epoch=='BL':
            X_BL = np.mean(X[...,gv.bins_BL[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_BL
        elif epoch == 'Stim':
            X_STIM = np.mean(X[...,gv.bins_STIM[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_STIM
        elif epoch == 'ED':
            X_ED = np.mean(X[...,gv.bins_ED[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_ED
        elif epoch == 'Dist':
            X_DIST = np.mean(X[...,gv.bins_DIST[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_DIST
        elif epoch == 'MD':
            X_MD = np.mean(X[...,gv.bins_MD[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_MD
        elif epoch == 'LD':
            X_LD = np.mean(X[...,gv.bins_LD[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_LD
        elif epoch == 'Test':
            X_TEST = np.mean(X[...,gv.bins_test[:]-gv.bin_start],axis=-1) 
            X_epochs[i_epoch] = X_TEST 
            
    X_epochs = np.moveaxis(X_epochs,0,-1)  
    
    return X_epochs 

def selectiveNeurons(X_S1, X_S2, Threshold=.01):
    X_S1n = X_S1 
    X_S2n = X_S2 
        
    for i in range(X_S1n.shape[0]): 
        # X_S1n[i] = normalize(X_S1n[i]) 
        # X_S2n[i] = normalize(X_S2n[i]) 
        X_S1n[i] = z_score(X_S1n[i]) 
        X_S2n[i] = z_score(X_S2n[i]) 
        
    sel_idx = (X_S1n - X_S2n)/(X_S1n + X_S2n + gv.eps) 
    
    if gv.ED_MD_LD: 
        sel_idx = np.mean(sel_idx, axis=-1) 
    else:
        sel_idx = np.mean(sel_idx[:,:,gv.bins_STIM-gv.bin_start], axis=-1) 
        
    idx = np.where(abs(sel_idx)<=Threshold) 
    X_S1 = np.delete(X_S1, idx, axis=1) 
    X_S2 = np.delete(X_S2, idx, axis=1) 
    
    return X_S1, X_S2, idx

def deconvolveFluo(X):

    # F0 = np.empty( (X.shape[0], X.shape[1]) ) 
    # F0[:] = np.mean( np.mean(X[...,gv.bins_BL],axis=-1), axis=0 ) 
    F0 = np.mean(X[..., gv.bins_BL], axis=-1)
    
    # F0 = np.percentile(X, 15, axis=-1) 
    
    # def F0_loop(X, n_trial, n_neuron, bins): 
    #     X_ij = X[n_trial, n_neuron]        
    #     c, s, b, g, lam = deconvolve(X_ij, penalty=1) 
    #     return b
    
    # # loop over trials and neurons 
    # with pg.tqdm_joblib(pg.tqdm(desc='F0', total=X.shape[0]*X.shape[1])) as progress_bar: 
    #     F0 = Parallel(n_jobs=gv.num_cores)(delayed(F0_loop)(X, n_trial, n_neuron, gv.bins_BL) 
    #                                         for n_trial in range(X.shape[0]) 
    #                                         for n_neuron in range(X.shape[1]) )
        
    # F0 = np.array(F0).reshape( (X.shape[0], X.shape[1]) ) 
    
    # def X_loop(X, F0, n_trial, n_neuron):
    #     X_ij = X[n_trial, n_neuron]
    #     F0_ij = F0[n_trial, n_neuron]
    #     c, s, b, g, lam = deconvolve(X_ij, penalty=1, b=F0_ij) 
    #     return c 
    
    def S_loop(X, F0, n_trial, n_neuron):
        X_ij = X[n_trial, n_neuron]
        F0_ij = F0[n_trial, n_neuron]
        c, s, b, g, lam = deconvolve(X_ij, penalty=1) 
        # c, s, b, g, lam = deconvolve(X_ij, penalty=1, b=F0_ij) 
        return s 
    
    # # loop over trials and neurons 
    # with pg.tqdm_joblib(pg.tqdm(desc='denoise', total=X.shape[0]*X.shape[1])) as progress_bar: 
    #     X_dcv = Parallel(n_jobs=gv.num_cores)(delayed(X_loop)(X, F0, n_trial, n_neuron) 
    #                                         for n_trial in range(X.shape[0]) 
    #                                         for n_neuron in range(X.shape[1]) ) 
    # X_dcv = np.array(X_dcv).reshape(X.shape) 
    
    with pg.tqdm_joblib(pg.tqdm(desc='deconvolve', total=X.shape[0]*X.shape[1])) as progress_bar: 
        S_dcv = Parallel(n_jobs=gv.num_cores)(delayed(S_loop)(X, F0, n_trial, n_neuron) 
                                              for n_trial in range(X.shape[0]) 
                                              for n_neuron in range(X.shape[1]) ) 
        
    S_dcv = np.array(S_dcv).reshape(X.shape)    
    # S_flt = savgol_filter(S_dcv, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = 5, deriv=0)
    
    def threshold_spikes(S_dcv, threshold): 
        # S_dcv[S_dcv<=threshold] = 0 
        # S_dcv[S_dcv>threshold] = 1 
        # S_dcv = uniform_filter1d( S_dcv, int(gv.frame_rate/2) ) 
        return S_dcv*1000
    
    S_th = threshold_spikes(S_dcv, gv.DCV_THRESHOLD)  
    S_avg = np.mean(S_th[...,gv.bins_BL],axis=-1) 
    S_avg = S_avg[..., np.newaxis]
    
    print('X_avg', np.mean(S_avg)) 
    # removing silent neurons 
    # idx = np.argwhere(S_avg<=0) 
    # S_th = np.delete(S_th, idx, axis=1)
    
    # print('X_dcv', S_th.shape[1]) 
    # gv.n_neurons = S_th.shape[1] 
    
    if gv.Z_SCORE | gv.Z_SCORE_BL: 
        
        if gv.Z_SCORE_BL: 
            gv.bins_z_score = gv.bins_BL 
        else: 
            gv.bins_z_score = gv.bins 
            
        def scaler_loop(S, n_trial, bins): 
            S_i = S[n_trial] 
            scaler = StandardScaler() 
            scaler.fit(S_i[:,bins].T) 
            return scaler.transform(S_i.T).T 
        
        with pg.tqdm_joblib(pg.tqdm(desc='standardize', total=X.shape[0])) as progress_bar: 
            S_scaled = Parallel(n_jobs=gv.num_cores)(delayed(scaler_loop)(S_th, n_trial, gv.bins_z_score) 
                                                     for n_trial in range(X.shape[0]) ) 
        
        # def scaler_loop(S, i_neuron, bins): 
        #     S_neuron = S[:, i_neuron] 
            
        #     avg_trial = np.mean(S_neuron[:, bins], axis=-1) 
        #     avg_trial = avg_trial[:, np.newaxis] 
            
        #     std_trial = np.std(S_neuron[:, bins], axis=-1) 
        #     std_trial = std_trial[:, np.newaxis] 
            
        #     return (S_neuron - avg_trial) / (std_trial) 
            
        #     # avg_all = np.mean(S_neuron[:, bins], axis=-1) 
        #     # std_all = np.std(S_neuron[:, bins]) 
        #     # return (S_neuron - avg_trial) / (std_all) 
        
        # with pg.tqdm_joblib(pg.tqdm(desc='standardize', total=X.shape[0])) as progress_bar: 
        #     S_scaled = Parallel(n_jobs=gv.num_cores)(delayed(scaler_loop)(S_th, i_neuron, gv.bins_z_score) 
        #                                              for i_neuron in range(gv.n_neurons) ) 
            
        S_scaled = np.asarray(S_scaled) 
        # S_scaled = np.swapaxes(S_scaled, 0, 1) 
        
        return S_scaled 
    
    return S_th 

def soft_thresholding():    
    ''' see Diagnosis of multiple cancer types by shrunken centroids of gene expression, Tibshirani et al. , 2002, PNAS

    Xij, i features, j samples 

    Xik = sum_j_Ck Xij/nk, sum on j in class Ck

    Xi = sum_i Xij /n, sum on the n samples , mean over samples 

    dik = Xik - Xi / mk (si +s0) where si^2= 1/(n-K) sum_k sum_j_Ck (Xij -Xik)^2 pooled within class standard deviation 
                                       s0 = median(si), guard against the possibility of large dik  
                                       mk = sqrt(1/nk + 1/n), so that mk*si is the std of the numerator in dik
    We rewrite as

    Xik = Xi + mk (si+s0) dik 

    and shrink the dik with soft thresholding defined as: 
    d'ik = sign(dik) TL(|dik|-D) where TL is t->t if t>0 (strictly), else t->0 

    This method has the desirable property that many of the features are eliminated from the class prediction as 
    the shrinkage parameter, D,  is increased.

    
    '''
    
    return 0 

def prescreening(X, y, alpha=0.05, scoring=f_classif): 
    ''' X is trials x neurons 
    alpha is the level of significance 
    scoring is the statistics, use f_classif or mutual_info_classif 
    '''
    
    model = SelectKBest(score_func=scoring, k=X.shape[1])    
    model.fit(X,y) 
    pval = model.pvalues_.flatten() 
    non_sel = np.argwhere(pval>alpha) 
    X_sel = np.delete(X, non_selected, axis=1) 
    return X_sel 

def preprocess_X(X):
    if gv.SAVGOL: 
        X = savgol_filter(X, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = gv.SAVGOL_ORDER, deriv=0, axis=-1) 
    
    if gv.F0_THRESHOLD is not None: 
        X = dFF0_remove_silent(X) 
        print(X.shape) 
        gv.n_neurons = X.shape[1] 
        
    if gv.DECONVOLVE: 
        if X.ndim>3: # task x sample x trial x neurons x time 
            for i_task in range(3): 
                for i_sample in range(2): 
                    X[i_task, i_sample] = deconvolveFluo(X[i_task, i_sample]) 
        else:
            X = deconvolveFluo(X) 
    else:            
        # if gv.SAVGOL: 
        #     X = savgol_filter(X, int(np.ceil(gv.frame_rate / 2.) * 2 + 1), polyorder = gv.SAVGOL_ORDER, deriv=0, axis=-1) 
            
        if gv.Z_SCORE | gv.Z_SCORE_BL :
            # print('z_score')
            X = z_score(X)
        if gv.Z_SCORE_TRIALS:
            X = z_score_trials(X) 
            
        if gv.NORMALIZE:
            # print('normalize') 
            X = normalize(X)
            
        if gv.NORMALIZE_TRIALS: 
            # print('normalize') 
            X = normalize_trials(X) 
            
        if gv.DETREND:
            X = detrend_X(X, order=gv.DETREND_ORDER) 

    if gv.standardize: 
        if X.ndim>3: # task x sample x trial x neurons x time 
            for i_task in range(3):
                X_task = np.vstack(X[i_task]) 
                m = np.nanmean( X_task, axis=0 ) 
                s = np.nanstd( X_task, axis=0 )
                
                # print('m', m.shape, np.mean(m), 's', s.shape, np.std(s)) 
                X[i_task] = ( X[i_task] - m ) / (s + gv.eps) 
    
    return X 
