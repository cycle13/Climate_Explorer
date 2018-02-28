#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 16 October 2015
# Last update: 16 October 2015
# Location: /data/local/hadkw/HADCRUH2/UPDATE2014/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code can make a lat, lon gridded field of:
#	a single month, 
#	an average of months within a year (or adjacent for DJF) up to annual - set minimum data presence
#	an average of single months across a period of years (climatology) - set minimum data presence
# 	an average of several months across a period of years (climatology) up to annual - set minimum data presence
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
# The code requires a 3D monthly resolution gridded dataset as time, lat, lon (anomalies or monthly means) 
# It also needs to know about the years/months contained
# It assumes data from January to December for each year
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# INPUTS:
# TheData: A 3D numpy array (months, latitude, longitude) of monthly means or anomalies 
# TheStYr = 1973 The start year of the provided data 
# TheEdYr = 2014 The end year of the provided data 
# TheChosenMonth = 1981 The start year of the new climatology period 
# TheChosenYear = 2010 The end year of the new climatology period 
# TheTheMDI = -1e30 The missing data indicator 
#     TheTheMDI=-1e30 # DEFAULT
# TheTheMDITol = 0.6 The proportion of data required for a gridbox climatology to be calculated from 0 to 1 
#     TheTheMDITol=0.6 # DEFAULT 
#
# python2.7
# from SelectSlice import SelectSlice
# TmpData = SelectSlice(TheData,TheStYr,TheEdYr,TheChosenMonth,TheChosenYear,TheTheMDI,TheTheMDITol)
#
# -----------------------
# OUTPUT
# -----------------------
# OUTPUTS:
# TmpData: a 3D array identical in lat, long shape to TheData for the output, utilises missing data indicator		
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
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
# SelectSlice
def SelectSlice(TheData,
		 TheStYr,
		 TheEdYr,
		 TheChosenMonth,
		 TheChosenYear,
                 TheMDI=-1e30,
		 TheMDITol=0.6):
    ''' This code takes in a 3D gridded field of monthly mean/anomalies (time,lat,lon)
        and anomalises/renormalises (climate anomalies) to a new climatology period
	which is supplied along side start and end years of the data.
	It assumes all data go from Jan to Dec 
        INPUTS:
	TheData: A 3D numpy array (months, latitude, longitude) of monthly means or anomalies 
        TheTheStYr: The start year of the provided data 
        TheTheEdYr: The end year of the provided data 
        TheChosenMonth: a LIST of month or months to select or average over
            Select your month of choice, or a range for an average
            0...11 represent Jan..Dec, [2,4] for Mar-Apr-May average, [0,11] for annual average, [11,1] for Dec-Jan-Feb average
            For month ranges that span 11 to 0, December will be taken from the first year of ChooseYr - will NOT work for last year!
            TheChosenMonth = [11] # [11,1]
        TheChosenYear: a LIST of year or years to select or average over
            Select your year of choice, or a range for an average
            1973...2014 for individual years, [1973,1982] for decadal average etc
            TheChosenYear = [2014] # [1981, 2010] 
        TheTheMDI: The missing data indicator 
        TheTheMDI=-1e30 # DEFAULT
        TheTheMDITol: The proportion of data required for a gridbox climatology to be calculated from 0 to 1 
        TheTheMDITol=0.6 # DEFAULT 
	OUTPUTS:
        TheField: a 3D array identical in shape to TheData for the output, utilises missing data indicator '''	    
    
    # Set up python imports
    import numpy as np
    import scipy.stats
    import pdb # pdb.set_trace() or c

    # Make empty array for the derived field filled with missing data
    TheField = np.reshape(np.repeat(TheMDI,len(TheData[0,:,0])*len(TheData[0,0,:])),(len(TheData[0,:,0]),len(TheData[0,0,:])))

    # Set up time pointers
    NYrs = (TheEdYr - TheStYr) + 1
    if (len(TheChosenMonth) > 1) | (len(TheChosenYear) > 1): # Only need to point to relevant years, this works if BOTH are > 1
        StTimePointer = (TheChosenYear[0]-TheStYr)
        if (len(TheChosenYear) > 1):
	    EdTimePointer = (TheChosenYear[1]-TheStYr)
        else:
	    EdTimePointer = StTimePointer
    else:
        StTimePointer = (TheChosenYear[0]-TheStYr)*12+TheChosenMonth[0]
        EdTimePointer = StTimePointer

    # Extract chosen month/year average
    
    # Easy case of single month from single year
    if (len(TheChosenMonth) == 1) & (len(TheChosenYear) == 1):
        print("One Month One Year")
        TheField = TheData[StTimePointer,:,:]
    
    # Easy-ish case of single month from multiple years - requires ?% presence for each gridbox
    elif (len(TheChosenMonth) == 1) & (len(TheChosenYear) > 1):
        print("One Month X Year")
        for lnn in range(len(TheData[0,0,:])):
            for ltt in range(len(TheData[0,:,0])):
	        subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(NYrs,12))) # fills columns first
                subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,TheChosenMonth[0]]) # need +1 for ranges
		subsubarr[subsubarr == TheMDI] = np.nan # set TheMDI to NaN
	        if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= TheMDITol):
	            TheField[ltt,lnn] = np.nanmean(subsubarr)

    # Slightly harder: multiple months from a single year (unless crossing DEC-JAN)
    elif (len(TheChosenMonth) > 1) & (len(TheChosenYear) == 1):
        print("X Month One Year")
        if (TheChosenMonth[1] > TheChosenMonth[0]):	# simple run of months within a year
            for lnn in range(len(TheData[0,0,:])):
                for ltt in range(len(TheData[0,:,0])):
	            subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(NYrs,12))) 
	            subsubarr = np.copy(subarr[StTimePointer,TheChosenMonth[0]:TheChosenMonth[1]+1]) # need +1 for ranges
		    subsubarr[subsubarr == TheMDI] = np.nan # set TheMDI to NaN
	            if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= TheMDITol):
	                TheField[ltt,lnn] = np.nanmean(subsubarr)
        else: # more complex as need to pull out from two years
            for lnn in range(len(TheData[0,0,:])):
                for ltt in range(len(TheData[0,:,0])):
	            subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(NYrs,12)))
	            subsubarr = np.copy(subarr[StTimePointer,TheChosenMonth[0]:12]) # need +1 for ranges
		    subsubarr = np.append(subsubarr,subarr[StTimePointer+1,0:TheChosenMonth[1]+1]) # need +1 for ranges
		    subsubarr[subsubarr == TheMDI] = np.nan # set TheMDI to NaN
	            if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= TheMDITol):
	                TheField[ltt,lnn] = np.nanmean(subsubarr)

    # Hardest: multiple months and multiple years
    else: # now we're dealing with seasonal/annual average climatology
        print("X Month X Year")
        if (TheChosenMonth[1] > TheChosenMonth[0]):	# simple run of months and run of years
            for lnn in range(len(TheData[0,0,:])):
                for ltt in range(len(TheData[0,:,0])):
	            subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(NYrs,12)))
	            subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,TheChosenMonth[0]:TheChosenMonth[1]+1]) # need +1 for ranges
		    subsubarr[subsubarr == TheMDI] = np.nan # set TheMDI to NaN
	            if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= TheMDITol):
	                TheField[ltt,lnn] = np.nanmean(subsubarr)
        else: # more complex as need to pull out month runs across years 
            if (EdTimePointer < TheEdYr): # then we can go to the next year to get the extra months
                ExtraPointer=EdTimePointer+1
            else:
	        ExtraPointer=EdTimePointer
	    for lnn in range(len(TheData[0,0,:])):
                for ltt in range(len(TheData[0,:,0])):
	            subarr = np.copy(np.reshape(TheData[:,ltt,lnn],(NYrs,12)))
	            subsubarr = np.copy(subarr[StTimePointer:EdTimePointer+1,TheChosenMonth[0]:12,]) # need +1 for ranges
		    subsubarr = np.append(subsubarr,subarr[StTimePointer+1:ExtraPointer,0:TheChosenMonth[1]+1])
		    subsubarr[subsubarr == TheMDI] = np.nan # set TheMDI to NaN
	            if (np.float(len(subsubarr[np.isfinite(subsubarr)]))/np.float(len(subsubarr)) >= TheMDITol):
	                TheField[ltt,lnn] = np.nanmean(subsubarr)
	    
    return TheField # SelectSlice

##########################################################################
## TESTING CODE ##########################################################
##########################################################################
## Check if SelectSlice works 
## create a data array with an identical field for each month within year but increments annually
#TmpCandFields =  np.reshape(np.array(np.repeat(range(NYrs),12*3*7),dtype=float),(NMons,3,7))
#
## Check the selection output works on actual values - all should be ltt,lnn arrays of identical numbers
# SelectSlice(TmpCandFields,1973,2014,[6],[1980],-1e30,0.6)
## One month, one year: tested for June, 1980 = 7
## This works!
## One month, multiple years: tested October, 2000-2010 = mean of 27:37 = 32
## This works!
## Multiple months, one year: tested MAM, 1991 = mean of [18,18,18] = 18, tested DJF, 1992 = mean of [19,20,20] = 19.66666.
## This works for both!
## Multiple months, multiple years: tested SON, 1973-1982 = mean of 0:9,0:9,0:9 = 4.5, tested JAN-DEC, 1981-2010 mean of 8:37, 30 times = 22.5
## This works for both!
##########################################################################
#
## GetAnomalies works!
##########################################################################
    
##################################################################################################
