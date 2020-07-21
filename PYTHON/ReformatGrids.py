#!/usr/local/sci/bin/python
# PYTHON3
# 
# Author: Kate Willett
# Created: 16 October 2015
# Last update: 20 July 2020
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# GetAnomalies:
#    This renormalises monthly mean (anomaly or actual) data to the desired climatology period from actuals or anomalies
#    Can cope with missing data
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt: (may not all be required actually)
# import numpy as np
# import scipy.stats
# import pdb # pdb.set_trace() or c
#
# Kate's:
#
# -----------------------
# DATA
# -----------------------
# GetAnomalies:
#    The code requires a 3D monthly resolution gridded dataset as time, lat, lon (anomalies or monthly means) 
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# GetAnomalies:
#    INPUTS:
#    TheData: A 3D numpy array (months, latitude, longitude) of monthly means or anomalies 
#    TheStYr = 1973 The start year of the provided data 
#    TheEdYr = 2014 The end year of the provided data 
#    TheClimSt = 1981 The start year of the new climatology period 
#    TheClimEd = 2010 The end year of the new climatology period 
#    TheMDI = -1e30 The missing data indicator 
#        TheMDI=-1e30 # DEFAULT
#    TheMDITol = 0.6 The proportion of data required for a gridbox climatology to be calculated from 0 to 1 
#        TheMDITol=0.6 # DEFAULT 
#
#    python2.7
#    from ReformatGrids import GetAnomalies
#    TmpAnoms = GetAnomalies(TmpData,StYr,EdYr,ClimSt,ClimEd)
#
# -----------------------
# OUTPUT
# -----------------------
# GetAnomalies:
#    OUTPUTS:
#    TheAnoms: a 3D array identical in shape to TheData for the output, utilises missing data indicator 	    
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 2 (20th July 2020)
# ---------
#  
# Enhancements
# Now python 3 rather than 2.7
#  
# Changes
#  
# Bug fixes

# Version 1 (16th October 2015)
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
#
#************************************************************************
# Functions:
#########################################################################
# Get Anomalies
def GetAnomalies(TheData,
		 TheStYr,
		 TheEdYr,
		 TheClimSt,
		 TheClimEd,
                 TheMDI=-1e30,
		 TheMDITol=0.6):
    ''' This code takes in a 3D gridded field of monthly mean/anomalies (time,lat,lon)
        and anomalises/renormalises (climate anomalies) to a new climatology period
	which is supplied along side start and end years of the data.
	It assumes all data go from Jan to Dec 
        INPUTS:
	TheData: A 3D numpy array (months, latitude, longitude) of monthly means or anomalies 
        TheStYr: The start year of the provided data 
        TheEdYr: The end year of the provided data 
        TheClimSt: The start year of the new climatology period 
        TheClimEd: The end year of the new climatology period 
        TheMDI: The missing data indicator 
        TheMDI=-1e30 # DEFAULT
        TheMDITol: The proportion of data required for a gridbox climatology to be calculated from 0 to 1 
        TheMDITol=0.6 # DEFAULT 
	OUTPUTS:
        TheAnoms: a 3D array identical in shape to TheData for the output, utilises missing data indicator '''	    
    
    # Set up python imports
    import numpy as np
    import scipy.stats
    import pdb # pdb.set_trace() or c
    
    # Work out how many years of data in total and for the climatology period
    TheYCount = (TheEdYr - TheStYr) + 1
    TheCCount = (TheClimEd - TheClimSt) + 1
    
    # make empty array for anomalies
    TheAnoms = np.reshape(np.repeat(TheMDI,TheYCount*12*len(TheData[0,:,0])*len(TheData[0,0,:])),(TheYCount*12,len(TheData[0,:,0]),len(TheData[0,0,:])))
    
    # Loop through each latitude and longitude to work on a gridbox at a time
    for lnn in range(len(TheData[0,0,:])):
        for ltt in range(len(TheData[0,:,0])):
	    # Have to use copies!
            subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(TheYCount,12)))
            newsubarr = np.empty_like(subarr)
            newsubarr.fill(TheMDI)
            # loop through each month
            for mm in range(12):
                subclimarr = np.copy(subarr[(TheClimSt-TheStYr):(TheClimEd-TheStYr)+1,mm])	# always need +1 for working with subranges
                subsubarr = np.copy(subarr[:,mm])
                subclimarr[subclimarr == TheMDI] = np.nan
                subsubarr[subsubarr == TheMDI] = np.nan
                #print(np.float(len(subclimarr[np.isfinite(subclimarr)]))/np.float(TheCCount))
                if (np.float(len(subclimarr[np.isfinite(subclimarr)]))/np.float(TheCCount) >= TheMDITol):
                    subsubarr[np.isfinite(subsubarr)] = subsubarr[np.isfinite(subsubarr)] - np.nanmean(subclimarr)
                    subsubarr[np.isnan(subsubarr)] = TheMDI
                    newsubarr[:,mm] = np.copy(subsubarr)
            TheAnoms[:,ltt,lnn] = np.copy(np.reshape(newsubarr,TheYCount*12))
	    
    return TheAnoms # GetAnomalies

##########################################################################
## TESTING CODE ##########################################################
##########################################################################
## Check if GetAnomalies works - tested for 1981-2010
## create a data array with an identical field for each month within year but increments annually
#TmpCandFields =  np.reshape(np.array(np.repeat(range(NYrs),12*3*7),dtype=float),(NMons,3,7))
## Set up input vars:
#StYr = 1973
#EdYr = 2014
#ChooseClimSt = 1981 
#ChooseClimEd = 2010 
#
#TmpData = GetAnomalies(TmpCandFields,StYr,EdYr,ChooseClimSt,ChooseClimEd) # leave others as default
#
#moo=np.array(range(42),dtype=float)
#print(np.mean(moo[(ChooseClimSt-StYr):(ChooseClimEd-StYr)+1])) # need +1 for ranges
## 22.5
#print(moo-np.mean(moo[(ChooseClimSt-StYr):(ChooseClimEd-StYr)+1])) # need +1 for ranges
##array([-22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5
## -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5
##   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5
##  13.5  14.5  15.5  16.5  17.5  18.5]
#
#foo=np.reshape(CandFields[:,0,0],(NYrs,12))
#print(foo[:,0])
##array([-22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5
##       -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5
##         1.5	2.5   3.5   4.5   5.5	6.5   7.5   8.5   9.5  10.5  11.5  12.5
##        13.5  14.5  15.5  16.5  17.5  18.5]
#
#pdb.set_trace()#
#
## GetAnomalies works!
##########################################################################
    
##################################################################################################
