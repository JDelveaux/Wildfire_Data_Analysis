#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for data analysis on wildfires dataset

@author: jdelv
"""

#Imports
import numpy as np
import matplotlib.pyplot as plt

#General functions
def power_law(a,b,x):
    """
    Power law function

    Parameters
    ----------
    a : float
        Power law coefficient.
    b : float
        Power law exponent.

    """
    return a*x**(b)

def log_power_law(a,b,x):
    """
    Logarithmic power law function

    Parameters
    ----------
    a : float
        Power law coefficient.
    b : float
        Power law exponent.

    """
    return np.log10(a) - b * np.log10(x)

#Test statistics
def r2(log_data, log_mean, log_fxn):
    """
    Calculate the coefficient of determination
    (R^2) value for a test function

    Parameters
    ----------
    log_data : array
        data points to be used.
    log_mean : float
        mean of data.
    log_fxn : array
        test function applied to data.

    Returns
    -------
    float
        R^2 value.

    """
    ss_tot = 0
    ss_res = 0
    cnt = 0
    for i in log_data:
        ss_tot += np.power( (i-log_mean) , 2.0 )
        ss_res += np.power( (i-log_fxn[cnt]) , 2.0 )
        cnt += 1
    return 1.0 - (ss_res/ss_tot)

def KS_distance(log_data, log_fxn):
    """
    

    Parameters
    ----------
    log_data : TYPE
        data points to be used.
    log_fxn : array
        test function applied to data.

    Returns
    -------
    float
        Kolomogorov - Smirnof test statistic.

    """
    KS = []
    cnt = -1
    for i in log_data:
        cnt += 1
        if cnt<(0.1*log_data.size): continue
        #if cnt>(0.8*log_data.size): break
        KS.append(np.abs(log_fxn[cnt]-i))
    return max(KS)


#MLE Parameters
def mle_beta(g,a):
    """
    Returns emperical beta exponent for 
    MLE power law fit. See Deluca and Corral 
    (2013)

    Parameters
    ----------
    g : float
        geometric mean of data.
    a : float
        minimum value of data.

    Returns
    -------
    float.
        Beta exponent

    """
    
    return 1.0 + ( 1.0 / np.log(g/a))

def mle_alpha(beta, a, b):
    """
    Returns emperical alpha coefficient for 
    MLE power law fit. See Deluca and Corral 
    (2013)

    Parameters
    ----------
    beta : float
        MLE beta.
    a : flaot
        minimum value of data.
    b : float
        maximum value of data.

    Returns
    -------
    float
        MLE alpha.

    """
    return b * ((beta - 1.0) / ((np.power(a, (1.0 - beta))) - (np.power(b, (1.0 - beta)))))


#Kolmogorov - Smirnov Bootstrapping
def KS_threshold(log_data_r, log_data_s, log_fxn_r, log_fxn_s, pnts):
    """
    Kolmogorov - Smirnov threshold as determined in Babu et, al. (2004)

    Parameters
    ----------
    log_data_r : float
        emperical cumulative output with original data.
    log_data_s : float
        emperical cumulative output with sample data.
    log_fxn_r : float
        fitted cumulative output with original data.
    log_fxn_s : float
        fitted cululative output with resampled data.
    pnts : array
        list of data points.

    Returns
    -------
    float
        Threshold distance for Kolmogorov - Smirnov test.

    """
    KS = []
    cnt = -1
    for i in pnts:
        cnt += 1
        if cnt<(0.1*pnts.size): continue
        KS.append( np.abs( log_data_s[cnt] - log_fxn_s[cnt] - log_data_r[cnt] + log_fxn_r[cnt] ) )
    return max(KS)


#Plotting functions
def plot_KS_test(KS):
    """
    

    Parameters
    ----------
    KS : array
        list of KS distances between emperical and
        resampled data.

    Returns
    -------
    None.

    """
    #set x-axis
    n_points = len(KS)
    x = np.linspace(0, n_points, n_points)
    KS = sorted(KS)
    fiveP = round(n_points*0.05)
    ninetyfiveP = n_points - fiveP
    
    plt.figure(dpi=1200)
    
    #scatter plot of results
    plt.scatter(x,KS, c = 'k')
    
    #fill 95 percent blue
    plt.fill_between( np.linspace(fiveP,n_points,ninetyfiveP),
                     np.zeros(ninetyfiveP),
                     KS[fiveP:], color = 'b' )
    
    #fill 5 percent red
    plt.fill_between( np.linspace(0,fiveP,fiveP),
                     np.zeros(fiveP),
                     KS[:fiveP], color = 'r' )
    
    
    #plot parameters
    plt.ylabel('KS Distance')
    plt.xlabel('Trials')
    plt.xlim(0,n_points)
    plt.ylim(0)
    plt.xticks([50, n_points], [50, n_points])
    
    #display results
    plt.show()