#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program to determine power-law fits for several methods
with confidence intervals

@author: jdelv
"""

#Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functions import mle_beta, mle_alpha

#Load data and setup variables
dataset = 'child.csv'
data = pd.read_csv(dataset, sep = ',')

n_runs = 1000

params = {
    'MLE_exp' : [],
    'MLE_cof' : [],
    'LLS_exp' : [],
    'LLS_cof' : [],
    }

for _ in range(n_runs):
    
    #Sample the data
    smpl = pd.DataFrame()
    smpl['area'] = data['Ranked km2'].sample(frac = 0.8, replace = True)
    
    #Bin data logarithmically
    smpl['logArea'] = np.log10(smpl['area'])
    _, bins = pd.qcut(smpl['logArea'], 15, duplicates = 'drop', retbins = True)
    #Unlog bins
    bins = np.power(10, bins)
    #Note: must extend bin widths slightly for proper sorting
    bins[0] -= 1e-6
    bins[-1] += 1e-6
    #Sort data into bins
    smpl['bin'] = pd.cut(smpl['area'], bins, right = False)
    #Get counts / bin
    bin_counts = smpl['bin'].value_counts()
    bin_counts = bin_counts.sort_index()
    
    #Determine Non-normalized frequency ( num in bin / bin width )
    freq = []
    std_dev = []
    bin_midpoints = (bins[1:] + bins[:-1]) / 2.
    for i in bin_counts.index:
        std_dev.append(2.0*np.sqrt(bin_counts[i]))
        freq.append( bin_counts[i] / ( i.right - i.left ) )
    
    #Normalize frequency by years of data (30 years)
    freq = np.multiply(freq, 1./30.)
    
    #Linear Least Squares Estimation
    y = np.polyfit(np.log10(bin_midpoints), np.log10(freq), 1)
    
    #MLE
    a = smpl['area'].min(0)
    b = smpl['area'].max(0)
    size = smpl['area'].size
    g = np.power(10, np.average(np.log10(smpl['area'])))
    beta = mle_beta(g,a)
    alpha = mle_alpha(beta, a, b)
    
    #Store Values
    params['MLE_exp'].append(-beta)
    params['MLE_cof'].append(alpha)
    params['LLS_exp'].append(y[0])
    params['LLS_cof'].append(10**y[1])
    
#Print or store results
print('MLE_exp average: ', np.average(params['MLE_exp']))
print('MLE_exp max: ', max(params['MLE_exp']))
print('MLE_exp min: ', min(params['MLE_exp']))
print('MLE_beta average: ', np.average(params['MLE_cof']))
print('MLE_beta max: ', max(params['MLE_cof']))
print('MLE_beta min: ', min(params['MLE_cof']))
print ('LLS_exp average: ', np.nanmean(params['LLS_exp']))
print ('LLS_exp max: ', max(params['LLS_exp']))
print ('LLS_exp min: ', min(params['LLS_exp']))
print ('LLS_beta average: ', np.nanmean(params['LLS_cof']))
print ('LLS_beta max: ', max(params['LLS_cof']))
print ('LLS_beta min: ', min(params['LLS_cof']))
    