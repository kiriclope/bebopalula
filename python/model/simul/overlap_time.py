import sys, os, importlib
from importlib import reload

import numpy as np
import matplotlib.pyplot as plt

import params as gv
importlib.reload(sys.modules['params'])
from utils import *
from plot_settings import *
SetPlotParams()

gv.init_param()

time, rates = get_time_rates(path=gv.path) 
print('rates',rates.shape) 

figtitle = 'rates_overlap_time_' 
fig = plt.figure(figtitle, figsize=(1.25*1.618*1.5*gv.RANK, 1.618*1.5)) 
ax = fig.add_subplot(1,gv.RANK,1)

overlap = get_overlap(rates, n_size=gv.n_size, ksi_path=gv.ksi_path, MAP=0) 
plt.plot(time, overlap[0], 'r-' ) 
plt.xlabel('Time (s)') 
plt.ylabel('Sample Overlap (a.u.)') 
add_vlines()
plt.xticks([0, 2, 4, 6, 8, 10, 12, 14]) 
plt.xlim([0, 14]) 

if(gv.RANK==2):
    overlap_1 = get_overlap(rates, n_size=gv.n_size, ksi_path=gv.ksi_path, MAP=1) 
    ax = fig.add_subplot(1,gv.RANK,2)
    plt.plot(time, overlap_1[0], 'r-' ) 
    plt.xlabel('Time (s)') 
    plt.ylabel('Dist. Overlap (a.u.)') 
    add_vlines()
    plt.xticks([0, 2, 4, 6, 8, 10, 12, 14]) 
    plt.xlim([0, 14]) 
    
