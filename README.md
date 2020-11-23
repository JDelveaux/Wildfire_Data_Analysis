# Wildfire Data Analysis
This repo features an example of data analysis on a wildfires dataset.

## Data Description
A complete description of the data, including aquisition, can be found in Brown, TJ (2004) or in the 'Data' section of 'WILDFIRES.pdf'.

Two datasets were used: fires started by arsonists from 1970-2000 and fires started by children in the same time period -- both for fires within the contiguous United States.

## Method of analysis
An overview of the data, including a spatial graph, can be found in section 2 of 'WILDFIRES.pdf'. The desired outcome was to determine power-law parameters for a frequency-size analysis of wildfires started by the two groups. Several methods were used to determine the coefficients, including Maximum Likelihood Estimation (MLE), Maximum Goodness-of-Fit (MGF), and variations of logarithmic Linear Least Squares (LLS).

Dependencies: numpy, matplotlib, pandas

## Method of verification
The power-law parameters were verified by the R^2 coefficient of determination and the Kolomogorov-Smirnov goodness of fit using bootstrapping methods on the dataset.
