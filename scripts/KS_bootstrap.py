#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A program to determine the vailidity of power-law expoments
by the Kolmogorov - Smirnov goodness-of-fit test

@author: jdelv
"""

#Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from functions import log_power_law, KS_threshold, plot_KS_test

#Parameters and Variables
dataset = 'arson.csv'
data = pd.read_csv(dataset, sep = ',')

n_runs = 1000 #number of runs
KS_D = [] #array to hold KS distances

for _ in range(n_runs):
    
    #obtain sampled data
    smpl = data['Ranked km2'].sample(frac = 1, replace = True).values
    #turn emperical data to array
    real = data['Ranked km2'].values
    
    #find shared fire sizes (for discete KS testing)
    shared_area = np.intersect1d(real, smpl)
    
    #sort data to find cumulative distribution
    smpl = sorted(smpl)
    real = sorted(real)
    
    #find index of shared areas
    smpl_idx = []
    real_idx = []
    
    for burn_area in shared_area:
        for j in range(len(real)):
            if real[j] == burn_area:
                real_idx.append(j+1)
                break
        for k in range(len(smpl)):
            if smpl[k] == burn_area:
                smpl_idx.append(k+1)
                break
    
    #normalize index
    smpl_idx = np.multiply(smpl_idx, float(1.0/len(smpl)))
    real_idx = np.multiply(real_idx, float(1.0/len(real)))
    
    #Model: LLS Estimation
    y_r = np.polyfit( np.log10(shared_area), np.log10(real_idx), 1 )
    y_s = np.polyfit( np.log10(shared_area), np.log10(smpl_idx), 1 )

    #apply function to data
    LLS_fxn_r = []
    LLS_fxn_s = []
    for area in shared_area:
        LLS_fxn_r.append(log_power_law(10**y_r[1],-y_r[0],area))
        LLS_fxn_s.append(log_power_law(10**y_s[1],-y_s[0],area))
        
    #determine distance threshold for sample
    KS_D.append( KS_threshold(np.log10(real_idx), np.log10(smpl_idx), LLS_fxn_r, LLS_fxn_s, shared_area) )
    
plot_KS_test(KS_D)
    