#!/usr/local/sci/bin/python

#***************************************
# 21st July 2015
# This version reads in from /data/local/hadkw/HADCRUH2
# Reads in a gridded netCDF file of some variable
# Computes various statistics as desired:

# GLOBMEAN 21st July 2015
# Calculate the global (or other region) mean value based on data present
# Weight by the cosine of the latitude
# Combine the uncertainties (in quadrature) to give an uncertainty based on the stations and within gridbox sampling
# OPTION: read in ERA-Interim masked and non-masked, use the difference to add uncertainty due to coverage?
#	On a year by year basis - does greater coverage provide a bias/systematic difference?
#	If so then this isn't really equivalent to a 2 sigma uncertainty!!!
# OPTION: provide % of gridboxes (land only) with data present
# 	This isn't perfect because the land/sea mask may not contain all small islands. 
#	It does contain some at least. 
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetMean_JUL2015.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from scipy.io import netcdf
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm

# Set up initial run choices
param='RH' # 'RH','q','e','T','Tw','Td','DPD'
param2='rh' # 'rh','q','e','t','tw','td','dpd'
runtype='IDPHA'	#'IDPHA','PHA','PHADPD','RAW@
stattype='globmean'	#'globmean','nhrmmean','shemmean','tropmean','fullglobmean'
StYear=1973	# Actual start year
EdYear=2014	# Actual end year
ClimSt=1976	# Start of climatology period
ClimEd=2005	# End of climatology period
PeriodSt=1979	# Chosen start year for statistic
PeriodEd=2003	# Chosen end year for statistic

ERA=True	# True,False - use ERA to look at coverage uncertainty
LAND=True	# True,False - use land cover mask to output % of land gridboxes with data in region

INDIRH='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
INFIL='HadISDH.land'+param+'.2.0.1.2014p_FLATgrid'+runtype+'5by5_JAN2015_cf'

INDIRO='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
INERA=param2+'2m_monthly_5by5_ERA-Interim_data_19792014'

#incover='new_coverpercentjul08'
INCOV='HadCRUT.4.3.0.0.land_fraction'

nowmon='JUL'
nowyear='2015'

# Set up variables
nlats=36		#set once file read in
nlons=72		#set once file read in
CandData=[]
ERAData=[]
LandCover=[]
LatList=np.arange(-87.5,92.5,5.)
LonList=np.arange(-177.5,182.5,5.)

# Missing data
MDI=-1e30 # may set up as masked arrays later

#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFGRID
def ReadNetCDFGrid(FileName,Var_Name,Unc_Name):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data '''

    ncf=netcdf.netcdf_file(FileName,'r')
    # ncf.variables this lists the variable names
#    var=ncf.variables['q_MPtrend']
#    var=ncf.variables['RH_MPtrend']
    if (Unc_Name):	# if false then do nothing
        var=ncf.variables[Unc_Name]
        TheUnc=np.array(var.data)
    else:
        TheUnc=[]	
    var=ncf.variables[Var_Name]
    TheData=np.array(var.data)
#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheUnc # ReadNetCDFGrid

#************************************************************************
# GlobMean
def GlobMean(TheData,TheUnc,TheLandCover,TheERAData,TheLatList,TheLonList,TheClimSt,TheClimEd,ThePeriodSt,ThePeriodEd,TheStYear,TheEdYear,TheMDI):
    ''' Create Period Mean for Each Gridbox - maximises coverage - must have 50% of years with 100% of months   '''
    ''' Create Region Mean from each Period Mean gridbox - weight by cosine of latitude '''
    ''' Create Period Mean uncertainty for each Gridbox - combine in quadrature which assumes these are NOT correlated errors '''
    ''' Create REgion mean uncertainty - combine in quadrature, weighting by cosine of latitude, which assumes these are NOT correlated errors '''
    ''' IF ERA switch is on then look at uncertainty due to missing gridboxes '''
    ''' This will mean renormalising to climatology, masking and generally involved stuff '''
    ''' uses ERACoverageUnceratinty  - which uses CalcAnomalies'''
    ''' IF LAND switch is on then provide % coverage of region in terms of gridboxes over land '''
    ''' uses GetPercentCoverage ''' 
    ''' Print out results '''
    
    # Set up the working arrays
    PeriodMeans=np.empty((len(TheLatList),len(TheLonList)))
    PeriodMeans.fill(TheMDI)
    PeriodMeanUncs=np.empty((len(TheLatList),len(TheLonList)))
    PeriodMeanUncs.fill(TheMDI)
    GlobMeanVal=TheMDI
    GlobMeanUnc=TheMDI  
    
    # Set up the working variables
    NAYrs=(TheStYear-TheEdYear)+1
    NPYrs=(ThePeriodEd-ThePeriodSt)+1
    if (NAYrs != NPYrs):	# if these aren't the same then you'll need to subset the data to the desired length
        StMon=(ThePeriodSt-TheStYear)*12
	EdMon=((ThePeriodEd-TheStYear)*12)+12  
	TheData=TheData[StMon:EdMon,:,:]
        TheUnc=TheUnc[StMon:EdMon,:,:]
    
    # Fill the Period Means (just wondering whether we should use climatology added back to anomalies - ah well...)
    # Simulataneously combine in quadrature the PeriodMeanUncs - assume these are independent
    for ltt in range(len(TheLatList)):
        for lnn in range(len(TheLonList)):
            # Are there enough years with all 12 months present?
	    Subarr=np.reshape(TheData[:,ltt,lnn],(NPYrs,12))
	    Subunc=np.reshape(TheUnc[:,ltt,lnn],(NPYrs,12))
	    Counts=0
	    TotalVals=0.
	    TotalUncs=0.
	    for yy in range(NPYrs): # v slow way of working through but checking there is enough data
	        gots=np.where(Subarr[yy,:] > TheMDI)[0]
	        gotsU=np.where(Subunc[yy,:] > TheMDI)[0]
	        gotsZ=np.where(Subarr[yy,:] == 0.0)[0]
	        gotsZU=np.where(Subarr[yy,:] == 0.0)[0]
		#if (len(gots) != len(gotsU)): # NO CASES FOUND
		#    print('non-match')
		#    stop
		if ((len(gotsZ) > 0) | (len(gotsZU) > 0)):
		    print(ltt,lnn,yy,'Odd zeros') # THERE ARE SOME CASES (5 for q) - FORCE TO BE 0.000001
		    Subarr[yy,np.where(Subarr[yy,:] == 0.0)[0]]=0.01	# Not ideal but stops silly numbers
		    Subunc[yy,np.where(Subunc[yy,:] == 0.0)[0]]=0.01
		if (len(gots) == 12):	# there are 12 months with data present so add
		    Counts=Counts+12
		    TotalVals=TotalVals+np.sum(Subarr[yy,:])
		    TotalUncs=TotalUncs+np.sum((Subunc[yy,:]/Subarr[yy,:])**2) # must normalise to get relative uncertainties in order to combine over time???
#	    print(ltt,lnn)
#	    if ((ltt == 31) & (lnn == 62)):
#	        stop
		    #print('GOT ONE!',yy,Subarr[yy,:])
		    #stop
	    
	    if (Counts > (NPYrs*12)*0.5):	# must be at least 50% of complete years present
	        PeriodMeans[ltt,lnn]=TotalVals/Counts
	        PeriodMeanUncs[ltt,lnn]=np.sqrt(TotalUncs)
	    	

    #stop
    # if there is LandCover data then calculate the % gridbox coverage over land	
    if (len(TheLandCover) > 0):
        TheLandCover=GetPercentCoverage(PeriodMeans,TheLandCover,TheLatList,TheLonList,TheMDI,[-90,90])	

    # if there is ERAData then calculate the coverage uncertainty component using ERA
    # THERE MUST BE LAND COVER FOR THIS!!!
    if (len(TheERAData) > 0):
        # do something
	print('Working with ERA')
	ERACoverageUncertainty(TheData,TheLandCover,TheERAData,TheLatList,TheLonList,TheClimSt,TheClimEd,ThePeriodSt,ThePeriodEd,TheStYear,TheEdYear,TheMDI)
		
    # Get Cosine of Latitude Weighted Region Mean
    # Simulataneously get combined uncertainty in quadrature - also weighted????
    ## TEST CODE GETTING GLOBAL MEAN:
    # PeriodMeans.fill(TheMDI)
    # PeriodMeans[1:10,5:15]=1.5
    # PeriodMeanUncs.fill(TheMDI)
    # PeriodMeanUncs[1:10,5:15]=3.5
    ## END TEST CODE - gave expected answer of 1.5 for PeriodMeans, tested PeriodMeanUncs
    TotWeights=0.
    TotVals=0.
    TotUncs=0.
    for ltt in range(len(TheLatList)):
    	TotVals=TotVals+np.sum(np.cos(np.radians(TheLatList[ltt]))*PeriodMeans[ltt,np.where(PeriodMeans[ltt,:] > TheMDI)[0]])
	TotWeights=TotWeights+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(PeriodMeans[ltt,:] > TheMDI)[0]))
    	TotUncs=TotUncs+np.sum(np.cos(np.radians(TheLatList[ltt]))*((PeriodMeanUncs[ltt,np.where(PeriodMeans[ltt,:] > TheMDI)[0]]/PeriodMeans[ltt,np.where(PeriodMeans[ltt,:] > TheMDI)[0]])**2))
#    	TotUncs=TotUncs+np.sum(np.cos(np.radians(TheLatList[ltt]))*((PeriodMeanUncs[ltt,np.where(PeriodMeans[ltt,:] > TheMDI)[0]])**2))
        #print(ltt,TotUncs)
	
    GlobMeanVal=TotVals/TotWeights
    GlobMeanUnc=np.sqrt(TotUncs) # guessing at the methodology here!	
    #stop		
    return GlobMeanVal,GlobMeanUnc # GlobMean
#*************************************************************************
# ERACoverageUncertainty
def ERACoverageUncertainty(TheData,TheLandCover,TheERAData,TheLatList,TheLonList,TheClimSt,TheClimEd,ThePeriodSt,ThePeriodEd,TheStYear,TheEdYear,TheMDI):
    ''' Print out average for same period for ERA full (or 1979+ at least) '''
    ''' Print out average for same period for ERA masked (or 1979+ at least) '''
    ''' Try to get a handle on the coverage uncertainty by looking at difference for each year on masked vs full '''
    
    # Set up the working arrays
    ERAPeriodMeansALL=np.empty((len(TheLatList),len(TheLonList)))
    ERAPeriodMeansALL.fill(TheMDI)
    ERAPeriodMeans=np.empty((len(TheLatList),len(TheLonList)))
    ERAPeriodMeans.fill(TheMDI)
    ERAPeriodMeansOCN=np.empty((len(TheLatList),len(TheLonList)))
    ERAPeriodMeansOCN.fill(TheMDI)
    ERAPeriodMeansMSK=np.empty((len(TheLatList),len(TheLonList)))
    ERAPeriodMeansMSK.fill(TheMDI)
    ERAPeriodMeanUncs=np.empty((len(TheLatList),len(TheLonList)))
    ERAPeriodMeanUncs.fill(TheMDI)
    ERAGlobMeanVal=TheMDI
    ERAGlobMeanValMSK=TheMDI
    ERAGlobMeanUnc=TheMDI  
    
    # Set up the working variables
    NAYrs=(TheStYear-TheEdYear)+1
    NPYrs=(ThePeriodEd-ThePeriodSt)+1
    if (NAYrs != NPYrs):	# if these aren't the same then you'll need to subset the data to the desired length
        StMon=(ThePeriodSt-1979)*12
	EdMon=((ThePeriodEd-1979)*12)+12  
	TheERAData=TheERAData[StMon:EdMon,:,:]
    
    # Fill the Period Means (just wondering whether we should use climatology added back to anomalies - ah well...)
    # Simulataneously combine in quadrature the PeriodMeanUncs - assume these are independent
    for ltt in range(len(TheLatList)):
        for lnn in range(len(TheLonList)):
            # Are there enough years with all 12 months present?
	    SubERAarr=np.reshape(TheERAData[:,ltt,lnn],(NPYrs,12))
	    Subarr=np.reshape(TheData[:,ltt,lnn],(NPYrs,12))
	    CountsALL=0
	    TotalValsALL=0.
	    CountsOCN=0
	    TotalValsOCN=0.
	    Counts=0
	    TotalVals=0.
	    CountsMSK=0
	    TotalValsMSK=0.
	    for yy in range(NPYrs): # v slow way of working through but checking there is enough data
	        gots=np.where(SubERAarr[yy,:] > TheMDI)[0]
	        gotsH=np.where(Subarr[yy,:] > TheMDI)[0]
		if (len(gots) == 12):	# there are 12 months with data present so add
		    CountsALL=CountsALL+12
		    TotalValsALL=TotalValsALL+np.sum(SubERAarr[yy,:])
		if ((len(gots) == 12) & (TheLandCover[ltt,lnn] < 0.5)):	# there are 12 months with data present so add
		    CountsOCN=CountsOCN+12
		    TotalValsOCN=TotalValsOCN+np.sum(SubERAarr[yy,:])
		if ((len(gots) == 12) & (TheLandCover[ltt,lnn] >= 0.5)):	# there are 12 months with data present so add
		    Counts=Counts+12
		    TotalVals=TotalVals+np.sum(SubERAarr[yy,:])
		if (len(gotsH) == 12):	# there are 12 months with data present so add
		    CountsMSK=CountsMSK+12
		    TotalValsMSK=TotalValsMSK+np.sum(SubERAarr[yy,:])

#	    print(ltt,lnn)
#	    if ((ltt == 31) & (lnn == 62)):
#	        stop
		    #print('GOT ONE!',yy,Subarr[yy,:])
		    #stop
	    
	    if (CountsALL > (NPYrs*12)*0.5):	# must be at least 50% of complete years present
	        ERAPeriodMeansALL[ltt,lnn]=TotalValsALL/CountsALL
	    if (CountsOCN > (NPYrs*12)*0.5):	# must be at least 50% of complete years present
	        ERAPeriodMeansOCN[ltt,lnn]=TotalValsOCN/CountsOCN
	    if (Counts > (NPYrs*12)*0.5):	# must be at least 50% of complete years present
	        ERAPeriodMeans[ltt,lnn]=TotalVals/Counts
	    if (CountsMSK > (NPYrs*12)*0.5):	# must be at least 50% of complete years present
	        ERAPeriodMeansMSK[ltt,lnn]=TotalValsMSK/CountsMSK
	    	
	
    # Get Cosine of Latitude Weighted Region Mean
    TotWeightsALL=0.
    TotValsALL=0.
    TotWeightsOCN=0.
    TotValsOCN=0.
    TotWeights=0.
    TotVals=0.
    TotWeightsMSK=0.
    TotValsMSK=0.
    for ltt in range(len(TheLatList)):
    	TotValsALL=TotValsALL+np.sum(np.cos(np.radians(TheLatList[ltt]))*ERAPeriodMeansALL[ltt,np.where(ERAPeriodMeansALL[ltt,:] > TheMDI)[0]])
	TotWeightsALL=TotWeightsALL+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(ERAPeriodMeansALL[ltt,:] > TheMDI)[0]))
    	TotValsOCN=TotValsOCN+np.sum(np.cos(np.radians(TheLatList[ltt]))*ERAPeriodMeansOCN[ltt,np.where(ERAPeriodMeansOCN[ltt,:] > TheMDI)[0]])
	TotWeightsOCN=TotWeightsOCN+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(ERAPeriodMeansOCN[ltt,:] > TheMDI)[0]))
    	TotVals=TotVals+np.sum(np.cos(np.radians(TheLatList[ltt]))*ERAPeriodMeans[ltt,np.where(ERAPeriodMeans[ltt,:] > TheMDI)[0]])
	TotWeights=TotWeights+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(ERAPeriodMeans[ltt,:] > TheMDI)[0]))
    	TotValsMSK=TotValsMSK+np.sum(np.cos(np.radians(TheLatList[ltt]))*ERAPeriodMeansMSK[ltt,np.where(ERAPeriodMeansMSK[ltt,:] > TheMDI)[0]])
	TotWeightsMSK=TotWeightsMSK+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(ERAPeriodMeansMSK[ltt,:] > TheMDI)[0]))
	
    GlobMeanValALL=TotValsALL/TotWeightsALL
    GlobMeanValOCN=TotValsOCN/TotWeightsOCN
    GlobMeanVal=TotVals/TotWeights
    GlobMeanValMSK=TotValsMSK/TotWeightsMSK
    
    print('ERA GLOBMEAN ALL: = ',GlobMeanValALL)	
    print('ERA GLOBMEAN OCEAN: = ',GlobMeanValOCN)	
    print('ERA GLOBMEAN LAND: = ',GlobMeanVal)	
    print('ERA GLOBMEAN MASKED : = ',GlobMeanValMSK)	
    
    
    return # ERACoverageUncertainty
#*************************************************************************
# GetPercentCoverage
def GetPercentCoverage(TheMeans,TheLandCover,TheLatList,TheLonList,TheMDI,TheLimits):
    # output pct of land boxes represented per region
    for ltt in range(len(LatList)):
        # First merge LandCover with Data coverage so that at least some small islands are counted initially
        TheLandCover[ltt,np.where(TheMeans[ltt,:] > TheMDI)[0]]=0.5
        
    StLat=np.flipud(np.where(TheLatList-2.5 <= TheLimits[0])[0])[0]
    EdLat=np.flipud(np.where(TheLatList+2.5 <= TheLimits[1])[0])[0]
    TotLats=(EdLat-StLat)+1
    print(TheLimits,StLat,EdLat)
    
    GlobCount=len(np.where(TheLandCover[StLat:EdLat+1,:] > 0)[0])

    ActGlobCount=len(np.where(TheMeans[StLat:EdLat+1,:] != TheMDI)[0])
	
    print('REGION: ',72*TotLats,' boxes in total of which ',GlobCount,' are land and ',ActGlobCount,' have obs')
    print('Percent land gridboxes in region: ',float(GlobCount)/(72*TotLats))
    print('Percent of region land gridboxes with obs: ',float(ActGlobCount)/GlobCount)	

    return TheLandCover # GetPercentCoverage
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in grids
MyFile=INDIRH+INFIL+'.nc'
VarName=param2+'_abs'
UncName=param2+'_combinederr'
CandData,CandUnc=ReadNetCDFGrid(MyFile,VarName,UncName)

if (ERA): # read in the ERA data
    MyFile=INDIRO+INERA+'.nc'
    VarName='actuals'
    UncName=False	# can then be tested
    ERAData,ERAUnc=ReadNetCDFGrid(MyFile,VarName,UncName) # ERAUnc is null

if (LAND): # read in the ERA data
    MyFile=INDIRO+INCOV+'.nc'
    VarName='land_area_fraction'
    UncName=False	# can then be tested
    LandCover,LandUnc=ReadNetCDFGrid(MyFile,VarName,UncName) # LandUnc is null
    LandCover=np.flipud(LandCover)
    LandCover=LandCover[0,:,:]

# TEST THESE ARE ALL THE RIGHT WAY UP
    # YES THEY ARE!!!
    #stop

# Now run which ever statistic you have chosen

if (stattype == 'globmean'):
    GlobMeanVal,GlobMeanUnc=GlobMean(CandData,CandUnc,LandCover,ERAData,LatList,LonList,ClimSt,ClimEd,PeriodSt,PeriodEd,StYear,EdYear,MDI)
    print('GLOBMEAN '+param+': = ',GlobMeanVal)	
    print('GLOBMEAN UNCERTAINTY '+param+': = ',GlobMeanUnc)	
		
#    stop()

print("And, we are done!")

