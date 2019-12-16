#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 11 October 2012
# Last update: 3 May 2019
# Location: /data/local/hadkw/ISTI/PROGS/	
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# All code in this file calculates linear trends
#
# MedianPairwise:
#     Median of Pairwise slopes estimation of a linear trend
#     95 percent confidence intervals based on weighted percentiles linked to degrees of freedom
#     Copes with missing data
#     After Sen, 1961
#     Based on IDL code - Mark McCarthy
# 
# Sen, P. K.: Estimates of the regression coef?cient based on Kendall?s tau, J. Am. Stat. Assoc., 63, 1379?1389, 1968.
# Helsel and Hirsch 1995, Statistical Methods in Water Resources, page 273 http://pubs.usgs.gov/twri/twri4a3/pdf/twri4a3-new.pdf
#
#
# OLS_AR1Corr
#	Ordinary Least Squares fit of a linear trends
# 	90 pct (2sigma) confidence intervals from Standard Error of the slope of the trend adjusted for AR(1) correlation
#	Following the method of Santer et al., 2008
#       For this we use statsmodels and pandas
#
# Santer, B. D., P. W. Thorne,  L. Haimberger,  K. E. Taylor,  T. M. L. Wigley,
# J. R. Lanzante,  S. Solomon,  M. Free,  P. J. Gleckler,  P. D. Jones,  T. R. Karl,
# S. A. Klein,  C. Mears,  D. Nychka,  G. A. Schmidt,  S. C. Sherwood, and  F. J. Wentz, 2008:
# Consistency of modelled and observed temperature trends in the tropcial troposphere. 
# International Journal of Climatology, 28 (13), 1703-1722.
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt: (may not all be required actually)
# import numpy as np
# import matplotlib.pyplot as plt
# import sys, os
# from scipy.optimize import curve_fit,fsolve,leastsq
# from scipy import pi,sqrt,exp
# from scipy.special import erf
# import scipy.stats
# from math import sqrt,pi
# import struct
#
# import statsmodels.api as sm
# import statsmodels.formula as smf
# import pandas as pd
# 
# -----------------------
# DATA
# -----------------------
# MedianPairwise:
#    TheData - a numpy array of data which can contain missing data
#    TheMDI - the missing data indicator
#    TheSlope - a 3 element array or list to contain [trend per time unit, lower bound confidence range, upper bound confidence rnage]
# 
# OLS_AR1Corr
#    TheData - a numpy array of data which can contain missing data
#    TheMDI - the missing data indicator
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# MedianPairwise:
#     from LinearTrends import MedianPairwise
#     TheSlope=[0.,0.,0.]
#     TheSlope=MedianPairwise(TheData,TheMDI,TheSlope)
#
#
# OLS_AR1Corr:
#     from LinearTrends import OLS_AR1Corr
#     TheSlope=OLS_AR1Corr(TheData,TheMDI)
#
# -----------------------
# OUTPUT
# -----------------------
# MedianPairwise:
# 	Outputs a 3 element list containing: 
#		the trend per 1 time unit, 
#		the lower bound of 95th conf, 
#		the upper bouund of 95th conf
#
# OLS_AR1Corr:
# 	Outputs a 3 element list containing: 
#		the trend per 1 time unit, 
#		the lower bound of 90th conf, (2 sigma = 2 * standard error) 
#		the upper bouund of 90th conf (2 sigma = 2 * standard error)
#
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 3 16th Dec 2019
# ---------
#  
# Enhancements
# NOw added OLS_AR1Corr which utilises statsmodels.formula and pandas
#  
# Changes
#  
# Bug fixes
#
# Version 2 3rd May 2019
# ---------
#  
# Enhancements
#  
# Changes
# Now python 3
#  
# Bug fixes
#
# 
# Version 1 (8th October 2015)
# ---------
#  
# Enhancements
# Added catch and change to make sure ranks don't go below zero or exceed N
#  
# Changes
#  
# Bug fixes
# No -1 for rankL needed
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
# This code was originally coded in IDL by Mark McCarthy
# An IDL version exists called median_pairwise.pro
#
# Just playing with slopes
# moo = np.random.rand(100) # 100 random numbers between 0 and 1
# moo = np.sin(np.arange(100) # 100 numbers between 0 and 1 with a sine curve
# moo = np.sin(np.random.normal(np.zeros(100))) # 100 normally distributed numbers with a mean of 0 and a sine curve
# moo = [(i*0.1)+ii for i,ii in enumerate(moo)]
#
#************************************************************************
# Set up python imports
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
from math import sqrt,pi
import struct
import pdb

import statsmodels.api as sm
import statsmodels.formula.api as smf
# >dir(smf) prints list of available models. Lower case ones accept formula 
# and df (for pandas data frame) arguments. Upper case take endog and exog design matrices. 
import pandas as pd

#************************************************************************
# Subroutines
#************************************************************************
# OLS_AR1Corr
def OLS_AR1Corr(TheData,TheMDI): # ,Lowee=Lowee,Highee=Highee):
    ''' Calculates the linear trend using Ordinary Least Squares regression '''
    ''' Can cope with specified missing data indicator TheMDI '''
    ''' Outputs the slope at a rate of unit per time step '''
    ''' Outputs the 90th pct standard error confidence intervals (2sigma) around the slope corrected for AR(1) correlation '''
    ''' Santer et al., 2008 - methodology '''
    ''' If Lowee and/or Highee are set they will come out as changed values '''

    # Set up empty list for returning the output values of the slope per unit time, lower (10th pct) CI, upper (90th pct) CI]
    TheSlope=[0.,0.,0.]		# median, 5th adn 95th percentile rate of change per time step

    # Convert the data to a pandas dataframe
    # First set any missing values to NaNs
    TheDataNANs = np.copy(TheData) # copying just to be safe
    TheDataNANs[np.where(TheDataNANs == TheMDI)] = float('NaN')
    DataDF = pd.DataFrame(TheDataNANs, index=np.arange(len(TheDataNANs)),columns=['variable'],dtype=float,copy=True)
    # This needs to be a copy otherwise changes to TheData lead to changes to DataDF.variable[]
    
    # As we are using formula= then we don't need to specify a column of 1s and an intercept is calculated
    olsmod = smf.ols(formula='variable ~ np.arange(len(TheDataNANs))',data=DataDF,missing='drop') # drop the NaNs
    olsres = olsmod.fit()
    # olsmod.summary() prints all
    # olsmod.params prints the slopes for each item 0 = intercept, 1 = variable
    # olsmod.predict() prints the predicted values using the model fit, including intercept
    # olsmod.bse prints the standard errors (1 sigma using n-2 degrees of freedom) 0 = intercept, 1 = variable
    TheSlope[0] = olsres.params[1]
    print('Trend: ',TheSlope[0])
    #pdb.set_trace()
    
    # Now get the original 1 sigma standard error which uses n-2 degrees of freedom and then calculate the corrected one
    # First we need to use masked arrays to make sure we account for missing data
    TheDataMSK = np.ma.masked_equal(TheData,TheMDI) # better hope that this captures them all and nothing silly with floats

# We need the AR(1) correlation of the regression residuals, not the actual datapoints
#    # We need the AR(1) valuee and shoudl probably test to make sure the data are autocorrelated
#    # Using the np.ma.corrcoef works even if all data are present
#    Lag1AR = np.ma.corrcoef(TheDataMSK[0:-1],TheDataMSK[1:])[0][1]
#    print('Autocorrelation at lag 1: ',Lag1AR) 
#    #pdb.set_trace()
    
    # Now get the time series of regression residuals for each data point
    # First create a masked array of missing data
    TheResids = np.ma.masked_equal(np.repeat(TheMDI,len(TheDataMSK)),TheMDI)
    # Get a pointer array to the non-missing data
    MaskofPoints = np.where(np.ma.getmask(TheDataMSK) == False)
    # Subtract the predicted values for each time point from the actual values
    TheResids[MaskofPoints] = TheDataMSK[MaskofPoints] - olsres.predict() # if there are missing values then these won't be predicted so we need to fill back in

    # We need the AR(1) values of the regression residuals and shoudl probably test to make sure the data are autocorrelated
    # Using the np.ma.corrcoef works even if all data are present
    Lag1AR = np.ma.corrcoef(TheResids[0:-1],TheResids[1:])[0][1]
    print('Autocorrelation at lag 1: ',Lag1AR) 
    #pdb.set_trace()

    # This is the original number of samples of time NOT INCLUDING MISSING DATA POINTS
    nORIG = np.ma.count(TheDataMSK) # number of data points
    
    # Now get the effective number of samples dependent on the degree of autocorrelation at lag 1
    nEFF = nORIG * ((1 - Lag1AR) / (1 + Lag1AR))
    print('Original no. time points: ',nORIG)
    print('Effective no. time points: ',nEFF)
    #pdb.set_trace()
    
    # Now get the variance of the regression residuals s_e^2
    s_eSQ = (1 / (nEFF - 2)) * np.ma.sum(TheResids**2)
    # and just for comparison get it for the original number of samples
    s_eSQORIG = (1 / (nORIG - 2)) * np.ma.sum(TheResids**2)

    print('Original variance of regression residuals: ',s_eSQORIG)
    print('Effective variance of regression residuals: ',s_eSQ)
    #pdb.set_trace()
    
    # Now calculate the 1 sigma standard error
    s_1sig = (s_eSQ / np.sum((MaskofPoints - np.mean(MaskofPoints))**2))**0.5
    # and just for comparison get it for the original number of samples
    s_1sigORIG = (s_eSQORIG / np.sum((MaskofPoints - np.mean(MaskofPoints))**2))**0.5

    print('Original 1 sigma standard error: ',s_1sigORIG)
    print('Effective 1 sigma standard error: ',s_1sig)
    #pdb.set_trace()
    
    # Now populate TheSlope array with the lower and upper bound 
    # I assume that this is slope - 2*s_1sig and slope + 2*s_1sig
    # When later the slope may be multiplied to get decadal trend the standard errors should be multiplied likewise
    TheSlope[1] = TheSlope[0] - (2 * s_1sig) 
    TheSlope[2] = TheSlope[0] + (2 * s_1sig) 

    print('AR(1) corrected 2 sigma standard error confidence intervals: ',TheSlope[1],TheSlope[2])
    #pdb.set_trace()
    
    return  TheSlope # ReadData

#***********************************************************************
# MedianPairwise
def MedianPairwise(TheData,TheMDI,TheSlope): # ,Lowee=Lowee,Highee=Highee):
    ''' Calculates slope from every point to every other point '''
    ''' Outputs the median of those slopes at a rate of unit per time step '''
    ''' Optionally outputs the 5th and 95th percentiles of those slopes as uncertainty ranges '''
    ''' If Lowee and/or Highee are set they will come out as changed values '''
    TheSlope=[0.,0.,0.]		# median, 5th adn 95th percentile rate of change per time step
    PairwiseSlopes=[]
    for i,ii in enumerate(TheData):
#        print i,ii
        if ii > TheMDI:
            for r in range(len(TheData)-1,i,-1):
                if TheData[r] > TheMDI:
                    PairwiseSlopes=np.append(PairwiseSlopes,((TheData[r]-ii)/(r-i)))
#		    print i,r,ii,TheData[r],(TheData[r]-ii)/(r-i)
    
#    print(len(PairwiseSlopes))
#    print "%5.2f"*len(PairwiseSlopes) % tuple(PairwiseSlopes)
    if len(PairwiseSlopes) > 10:	# check there are sufficient data
        TheSlope[0]=np.median(PairwiseSlopes)
        nData=len(list(TheData[np.where(TheData > TheMDI)]))
        DegofFree=(nData * (nData-1))/2.
        weight=np.sqrt(nData * (nData-1) * ((2*nData)+5) / 18.)
#	print "No. PW, No. Data Present, Deg Free ",len(PairwiseSlopes),nData, DegofFree
#	print "WEIGHT ", weight
        rankL=int(((DegofFree-1.96*weight)/2.))
        rankU=int(((DegofFree+1.96*weight)/2.)+1)

# Checks to make sure the ranks are actually sensible
        if (rankU >= len(PairwiseSlopes)):
            rankU=len(PairwiseSlopes)-1
        if (rankU < 0):
            rankU=0
        if (rankL < 0):
            rankL=0

#	print "RANKS ",rankL,rankU
        PairwiseSlopes.sort()
        TheSlope[1]=PairwiseSlopes[rankL]
        TheSlope[2]=PairwiseSlopes[rankU]

         	
#        print(TheSlope)
#        num_bins = 50
#        # the histogram of the data
#        n, bins, patches = plt.hist(PairwiseSlopes, num_bins, normed=1, facecolor='green', alpha=0.5)
#        plt.show()
    
    ''' At present I think the 5th and 95th are too sensitive to variability '''
    ''' This is a quicklook Kate version of the code and needs to be statistically verified '''
    return  TheSlope # ReadData

#************************************************************************
