import scipy as sp
import numpy as np
import params as gv 
import matplotlib.pyplot as plt

# def decode_m1_phi(signal):
    
#     dphi = np.pi / (double) n_per_pop[i_pop] ; 

#     x_coord = np.sum( rates * np.cos(2.0j * dphi) ) / .25
#     y_coord = np.sum( rates * np.sin(2.0j * dphi) ) / .25
    
    
#     m1 = ( 1.0 / (double) n_per_pop[i_pop]) * np.sqrt(x_coord * x_coord + y_coord * y_coord) 
#     phi = 0.5 * atan2(y_coord, x_coord) 

#     return m1, phi

def compute_phi(signal, axis=-1):
    # signal
    
    signal_copy = signal.copy()
    if axis!=-1 and signal.ndim!=1 :
        signal_copy = np.swapaxes(signal_copy, axis, -1) 
    
    length = signal_copy.shape[axis]
    dPhi = np.pi/length 
    
    dft = np.dot(signal_copy, np.exp(-2.0j * np.arange(length) * dPhi) )
    phi = 0.5 * ( np.arctan2(dft.imag, dft.real) % ( 2.0 * np.pi ) ) 
    
    if axis!=1 and signal.ndim!=1 :
        phi = np.swapaxes(phi, axis, -1) 
    
    return phi
    
def compute_m1(signal, axis=-1): 
    # if signal is nd array, need last dimension to be the dimension of interest for np.dot to work 
        
    signal_copy = signal.copy() 
    length = signal_copy.shape[axis]
    dPhi = np.pi/length
    
    if axis!=-1 and signal.ndim!=1 : 
        signal_copy = np.swapaxes(signal_copy, axis, -1) 
        
    # m1 = np.absolute(np.dot(signal_copy, np.exp(-2.0j * np.arange(length) * dPhi) ) ) / length 
    m1 = 2.0 * np.absolute(np.dot(signal_copy, np.exp(-2.0j * np.arange(length) * dPhi) ) ) / length 
    
    if axis!=1 and signal.ndim!=1 : 
        m1 = np.swapaxes(m1, axis, -1) 
        
    return m1 

def circular_convolution(signal, windowSize, axis=-1): 
    signal_copy = signal
    
    if axis!=-1 and signal.ndim!=1 : 
        signal_copy = np.swapaxes(signal_copy, axis, -1) 
    
    ker = np.concatenate((np.ones((windowSize, )), np.zeros((signal_copy.shape[-1] - windowSize, )))) 
    smooth_signal = np.real(np.fft.ifft( np.fft.fft(signal_copy, axis=-1) * np.fft.fft(ker, axis=-1) , axis=-1) ) * ( 1.0 / float(windowSize) ) 
    
    if axis!=1 and signal.ndim!=1 : 
        smooth_signal = np.swapaxes(smooth_signal, axis, -1) 
    
    return smooth_signal

def abs2(x):
    return x.real**2 + x.imag**2

def fourierfreq(signal):
    signal_copy = signal.copy()
    signal_copy -= np.mean(signal_copy)
    signal_copy *= sp.signal.windows.hann(len(signal_copy))
    
    fft = np.fft.rfft(signal_copy, norm="ortho")
    selfconvol=np.fft.irfft(abs2(fft), norm="ortho")

    return np.max(selfconvol) 
    

def get_m1_phi_smooth_rates(rates, n_size=gv.n_size, osc=0): 
    
    m1 = np.zeros( (rates.shape[1], rates.shape[0]) ) * np.nan 
    phi = np.zeros( (rates.shape[1], rates.shape[0]) ) * np.nan 
    smooth_rates = np.zeros( (rates.shape[1], rates.shape[0], rates.shape[-1]) ) * np.nan 
    m1_osc = np.zeros(rates.shape[1]) * np.nan 
    
    for i_pop in range(rates.shape[1]) : 
        
        pop_rates = rates[:, i_pop, : n_size[i_pop]].copy() 
        pop_rates  =  np.flip(pop_rates, axis=-1) 
        smooth_rates_pop = circular_convolution(pop_rates, int(pop_rates.shape[-1]*.001) ) # over neurons 
        
        print('smooth', smooth_rates_pop.shape )
        # smooth_rates[i_pop, :, :n_size[i_pop]] = np.flip(smooth_rates_pop.copy(), axis=-1) 
        m1[i_pop] = compute_m1(pop_rates.copy()) 
        phi[i_pop] = compute_phi(pop_rates.copy()) 
        # phi[i_pop] = circular_convolution(phi[i_pop], int(pop_rates.shape[0]*.1), axis=0 ) 
        # smooth_rates[i_pop, :, :n_size[i_pop]] = np.flip(smooth_rates_pop.copy(), axis=-1) 
        smooth_rates[i_pop, :, :n_size[i_pop]] = smooth_rates_pop.copy() 
        
        # phi_fft = abs(np.fft.fft(phi[i_pop]-np.mean(phi[i_pop]) ) ) / phi[i_pop].shape[0] 
        # m1_osc[i_pop] = np.max(phi_fft) * np.mean(m1[i_pop]) 
        # m1_osc[i_pop] = np.sqrt(np.mean( (phi[i_pop]-np.mean(phi[i_pop])**2 ) ) )
    
    if(osc):
        return m1, phi, smooth_rates, m1_osc 
    else:
        return m1, phi, smooth_rates 
    
def get_phi(rates, n_size=gv.n_size): 
    
    # phi = np.zeros( (rates.shape[1], rates.shape[0]) ) * np.nan 
    # for i_pop in range(rates.shape[1]):
        
    # phi = np.zeros( (rates.shape[1], rates.shape[0]) ) * np.nan
    i_pop = 0 
    pop_rates = rates[:, i_pop, : n_size[i_pop]]
    # smooth_rates_pop = circular_convolution(pop_rates, int(pop_rates.shape[-1]*.001) ) # over neurons 
    # smooth_rates_pop =  np.flip(smooth_rates_pop, axis=-1)
    smooth_rates_pop =  np.flip(pop_rates, axis=-1)
    
    # phi[i_pop] = compute_phi(smooth_rates_pop)
    phi = compute_phi(smooth_rates_pop)
    
    return phi

def get_avg_m1_phi_smooth_rates(rates, n_size=gv.n_size): 
    
    m1 = np.zeros(2) 
    phi = np.zeros(2) 
    smooth_rates = np.zeros( (2, rates.shape[-1]) ) 

    for i_pop in range(2) : 
        pop_rates = rates[i_pop, : n_size[i_pop]].copy()
        
        smooth_rates_pop = circular_convolution(pop_rates, int(pop_rates.shape[-1]*.001) ) # over neurons 
        
        m1[i_pop] = compute_m1(smooth_rates_pop.copy()) 
        phi[i_pop] = compute_phi(smooth_rates_pop.copy()) 
        smooth_rates[i_pop, :gv.n_size[i_pop]] = smooth_rates_pop.copy()
        # smooth_rates[i_pop, :gv.n_size[i_pop]] = np.flip(smooth_rates_pop.copy(), axis=-1) 
    
    return m1, phi, smooth_rates

def phitoPi(time, phi):

    new_time = []
    new_phi = [] 

    new_time.append(time[0]) 
    new_phi.append(phi[0]) 
    
    for i in range(phi.shape[0]-1):
        
        # if np.abs(phi[i]-phi[i+1])>=np.pi/4 :
        #     new_time.append((time[i+1]+time[i])/2.0) 
        #     new_phi.append(np.nan) 
            
        # elif phi[i]<phi[i+1] and phi[i+1]>=7*np.pi/8 :
        #     new_time.append((time[i+1]+time[i])/2.0) 
        #     new_phi.append(0) 
        # else:
        new_time.append(time[i+1]) 
        new_phi.append(phi[i+1]) 
    
    return new_time, new_phi
