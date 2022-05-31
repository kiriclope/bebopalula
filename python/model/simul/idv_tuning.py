import sys, os, importlib
from importlib import reload

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import params as gv 
importlib.reload(sys.modules['params']) 

gv.init_param()

from balance_inputs_dist import inputs_dist, vec_Phi
from utils import *
from get_m1 import compute_m1

rates_trial = []
n_trials = 25

for trial in range(1, n_trials+1):
    gv.TRIAL_ID = trial 
    gv.init_param() 
    time, rates = get_time_rates(path=gv.path) 
    rates_trial.append(np.mean(rates[20:40], axis=0))

rates_trial = np.asarray(rates_trial)
rates_trial = np.swapaxes(rates_trial, 0, -1)
print(rates_trial.shape)

# phi = compute_phi(rates_trial)

m0 = np.mean(rates_trial, axis=-1)
m1 = compute_m1(rates_trial)

print(m0.shape, m1.shape)
m1_m0_off = m1[:,0]/m0[:,0]

rates_trial = []
gv.folder = 'albert'
for trial in range(1, n_trials+1):
    gv.TRIAL_ID = trial 
    gv.init_param() 
    time, rates = get_time_rates(path=gv.path) 
    rates_trial.append(np.mean(rates[20:40], axis=0))

rates_trial = np.asarray(rates_trial)
rates_trial = np.swapaxes(rates_trial, 0, -1)
print(rates_trial.shape)

m0_on = np.mean(rates_trial, axis=0)
m1_on = compute_m1(rates_trial)

m1_m0_on = m1_on[:,0]/m0_on[:,0]

plt.scatter(m1_m0_off, m1_m0_on)
