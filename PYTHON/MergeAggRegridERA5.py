#!/usr/local/sci/bin/python
# PYTHON3.6.1
# 
# Author: Kate Willett
# Created: 18 Mar 2020
# Last update: 18 Mar 2020
# Location: /data/users/hadkw/WORKING_HADISDH/UPDATE2019/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/HadISDH_Build					
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# THIS CODE DOES MANY THINGS 
# NOTE THAT FOR ANY 1BY1 OUTPUT IT REGRIDS TO BE 89.5 to -89.5 rather than 90 - -90 (180 boxes rather than 181!!!)
# AND ROLLS LONGITUDE TO -179.5 to 179.5
#
# AT THE MOMENT THIS ASSUMES COMPLETE FIELDS SO WON'T WORK FOR SST!!! 
#
# ANOTHER ISSUE IS LAND / SEA MASKING - TOO MUCH LAND COVER, TOO MUCH SEA COVER - SO THERE WILL BE CONTAMINATION!
# I COMPUTE ANOMALIES AT 1by1 RES BEFORE REGRIDDING TO 5by5 TO MINIMISE THIS>
#
#
# This code reads in the ERA5 months of 1by1 daily or monthly variables
# (e.g., T, Td and Surface Pressure etc) for the full time period
#
# It averages to pentad means, anoms, clims and stdevs and saves to netCDF:
#	days since 19790101 (float), 180 lats 89.5 to -89.5, 360 lons -179.5 to 179.5, <var>2m
# It averages to monthly means, anoms, clims and stdevs and saves to netCDF:
#	days since 19790101 (float), 180 lats 89.5 to -89.5, 360 lons -179.5 to 179.5, <var>2m
# It regrids to 5by5 monthly means, anoms, clims and stdevs and saves to netCDF
#	days since 19790101 (int), 36 lats -87.5 to 87.5, 72 lons -177.5 to 177.5, actuals
# For anomalies it also creates a land only and ocean only set of grids (based on ERA5 land sea mask at 1by1 degree) to save along side the complete 
# fields
#	anomalies_land, anomalies_sea
# Climatology is 1981-2010
# 
# This can cope with an annual update and just append to existing files
#
# use submit_spice_ERA5download.bash and get_era5.py to download ERA5 months
# This may require a key to be set up in .ecmwfapirc annually - obtained from logging in to ECMWF
# https://confluence.ecmwf.int/display/WEBAPI/How+to+retrieve+ECMWF+Public+Datasets
# It also requires ecmwfapi to be downloaded and in the directory as you are running to code from
#
# Copy previous years of monthly ERA5 data from the previous 
# UPDATE<yyyy>/OTHERDATA/<var>2m_daily_1by1_ERA5_data_1979<yyyy>.nc
# UPDATE<yyyy>/OTHERDATA/<var>2m_pentad_1by1_ERA5_data_1979<yyyy>.nc
# UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_1by1_ERA5_data_1979<yyyy>.nc
# UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_5by5_ERA5_data_1979<yyyy>.nc
# to OTHERDATA/
#
# <references to related published material, e.g. that describes data set>
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# inbuilt:
# from datetime import datetime
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.dates import date2num,num2date
# import sys, os
# from scipy.optimize import curve_fit,fsolve,leastsq
# from scipy import pi,sqrt,exp
# from scipy.special import erf
# import scipy.stats
# from math import sqrt,pi
# import struct
# from netCDF4 import Dataset
# from netCDF4 import stringtoarr # for putting strings in as netCDF variables
# import pdb
#
# Kates:
# import CalcHums - written by kate Willett to calculate humidity variables
# import TestLeap - written by Kate Willett to identify leap years
# from ReadNetCDF import GetGrid4 - written by Kate Willett to pull out netCDF data
# from ReadNetCDF import GetGrid4Slice - written by Kate Willett to pull out a slice of netCDF data
# from GetNiceTimes import make_days_since
#
#-------------------------------------------------------------------
# DATA
# -----------------------
# ERA5 1by1 daily gridded data for each year and month
# THIS CODE ALWAYS WORKS FROM DAILY 1BY1 TO START
# EITHER TO BUILD ALL OUTPUTS FROM SCRATCH FROM INDIVIDUAL MONTHS:
# ERA5 = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<yyyy><mm>_daily_<variable>.nc
# OR TO APPEND NEW INDIVIDUAL MONTHS TO EXISTING OUTPUTS
# ERA5 = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<yyyy><mm>_daily_<variable>.nc
# ERA5OLD = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_daily_1by1_ERA5_1979<yyyy-1>.nc
# ERA5OLD = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_pentad_1by1_ERA5_1979<yyyy-1>.nc
# ERA5OLD = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_1by1_ERA5_1979<yyyy-1>.nc
# ERA5OLD = /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_5by5_ERA5_1979<yyyy-1>.nc
#
# It also needs a land -sea mask so use the ERA5 one ( created by get_era5_lsm.py)
# LSM = /data/users//hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/197901_hourly_land_sea_mask_ERA5.nc
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# First make sure the New ERA-Interim data are in the right place.
# Also check all editables in this file are as you wish
# module load scitools/default-current
# python MergeAggRegridERA5.py
#
# -----------------------
# OUTPUT
# -----------------------
# New ERA5 data for 1979 to present
# NewERA<var>:
# /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_daily_1by1_ERA5_1979<yyyy>.nc
# /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_pentad_1by1_ERA5_1979<yyyy>.nc
# /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_1by1_ERA5_1979<yyyy>.nc
# /data/users/hadkw/WORKING_HADISDH/UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_5by5_ERA5_1979<yyyy>.nc
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 1 (18 Mar 2020)
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
#                                 START
#************************************************************************
# inbuilt:
import datetime as dt
from datetime import datetime
#import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.dates import date2num,num2date
#import sys, os
#from scipy.optimize import curve_fit,fsolve,leastsq
#from scipy import pi,sqrt,exp
#from scipy.special import erf
#import scipy.stats
#from math import sqrt,pi
#import struct
from netCDF4 import Dataset
from netCDF4 import stringtoarr # for putting strings in as netCDF variables
import pdb

# Kates:
import CalcHums 
import TestLeap 
from ReadNetCDF import GetGrid4
from ReadNetCDF import GetGrid4Slice
from GetNiceTimes import MakeDaysSince

### START OF EDITABLES ###############################

# Set up initial run choices
# Start and end years
actstyr    = 1979 # start year of dataset
styr       = 1979 # start year of data to work with e.g., same as edyr
edyr       = 2019
thisedyr   = 2019
stmon      = 1
edmon      = 12

## Set up output variables - for q, e, RH, dpd, Tw we will need to read in multiple input files
#Var = 'q' # this can be 't','td','q','rh','e','dpd','tw','ws','slp','sp','uv','sst'

# Is this a new run or an update?
ThisProg = 'Build' 
# Update for updating an existing file ( from new 1by1 dailies and preious outputs daily monthly or pentad)
# Build for building from scratch (from 1by1 dailies)
	# THIS AUTOMATICALLY REGRIDS LATS TO BE 180 RATHER THAN 181!!!

# Is this ERA-Interim or ERA5?
ThisRean = 'ERA5' # 'ERA5' or 'ERA-Interim'

# Do you want to create anomalies and if so, what climatology period? We will output absolutes anyway
#MakeAnoms = 1 # 1 for create anomalies (and clim and stdev fields), 0 for do NOT create anomalies
ClimStart = 1981 # any year but generally 1981
ClimEnd = 2010 # any year but generally 2010

### END OF EDITABLES ################

# Set up file paths and other necessary things
#if (MakeAnoms == 1): # set the filename string for anomalies
#    AnomsStr = 'anoms'+str(ClimStart)+'-'+str(ClimEnd)+'_'
#else:
#    AnomsStr = ''

# Set up file locations
updateyy  = str(thisedyr)[2:4]
updateyyyy  = str(thisedyr)
workingdir  = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+updateyyyy
## RANDOM BIT FOR TESTING 
#edyr=1980
# dictionary: chosen and output var, actual var, ERA5 read in name
VarDict = dict([('q',['q','specific_humidity']),
                ('rh',['RH','relative_humidity']),
		('e',['e','vapour_pressure']),
		('t',['T','2m_temperature']),
		('tw',['Tw','wetbulb_temperature']),
		('td',['Td','2m_dewpoint_temperature']),
		('dpd',['DPD','dewpoint_depression']),
		('p',['P','surface_pressure']),
		('slp',['SLP','sea_level_pressure'])])

# Dictionary for looking up variable names for netCDF read in of variables
NameDict = dict([('q','q2m'),
                 ('rh','rh2m'),
	         ('e','e2m'),
	         ('tw','tw2m'),
	         ('t','t2m'),
	         ('td','td2m'),
	         ('dpd','dpd2m'),
	         ('slp','msl'),
	         ('p','p2m'), # does this matter that it is p not sp?
	         ('uv',['u10','v10']), # this one might not work
	         ('ws','si10'),
		 ('sst','sst')])

LandMask = workingdir+'/OTHERDATA/197901_hourly_land_sea_mask_ERA5.nc' # 0 = 100% sea, 1 = 100% land - no islands!, latitude, longitude, land_area_fraction, -87.5 to 87.5, -177.5 to 177.5

# Set up variables
mdi = -1e30

# Required variable names for reading in from ERA-Interim
LatInfo = ['latitude'] 
LonInfo = ['longitude'] 


# Dictionary for looking up variable standard (not actually always standard!!!) names for netCDF output of variables
StandardNameDict = dict([('q','specific_humidity'),
             ('rh','relative_humidity'),
	     ('e','vapour_pressure'),
	     ('tw','wetbulb_temperature'),
	     ('t','drybulb_temperature'),
	     ('td','dewpoint_temperature'),
	     ('dpd','dewpoint_depression'),
	     ('slp','mean_sea_level_pressure'),
	     ('p','surface_pressure'), # does this matter that its p and not sp?
	     ('uv',['10_metre_U_wind_component','10_metre_V_wind_component']), # this one might not work
	     ('ws','10_metre_windspeed'),
	     ('sst','sea_surface_temperature')])

# Dictionary for looking up variable long names for netCDF output of variables
LongNameDict = dict([('q','2m specific humidity from 1by1 hrly T, Td and p '+ThisRean),
             ('rh','2m relative humidity from 1by1 hrly T, Td and p '+ThisRean),
	     ('e','2m vapour pressure from 1by1 hrly T, Td and p '+ThisRean),
	     ('tw','2m wetbulb temperature from 1by1 hrly T, Td and p '+ThisRean),
	     ('t','2m drybulb temperature from 1by1 hrly T '+ThisRean),
	     ('td','2m dewpoint temperature from 1by1 hrly Td '+ThisRean),
	     ('dpd','2m dewpoint depression from 1by1 hrly T, Td and p '+ThisRean),
	     ('slp','2m mean_sea level pressure from 1by1 hrly msl '+ThisRean),
	     ('p','2m surface pressure from 1by1 hrly sp '+ThisRean), # does this matter that its p not sp
	     ('uv',['10 metre U wind component from 1by1 hrly '+ThisRean,'10 metre V wind component from 1by1 6hrly'+ThisRean]), # this one might not work
	     ('ws','10 metre windspeed from 1by1 hrly'+ThisRean),
	     ('sst','sea surface temperature from 1by1 hrly'+ThisRean)])

# Dictionary for looking up unit of variables
UnitDict = dict([('q','g/kg'),
             ('rh','%'), # made this % rather than %rh as iris doesn't seem to like it
	     ('e','hPa'),
	     ('tw','deg C'),
	     ('t','deg C'),
	     ('td','deg C'),
	     ('dpd','deg C'),
	     ('slp','hPa'),
	     ('p','hPa'), # does this matter that it is p not sp
	     ('uv','m/s'),
	     ('ws','m/s'),
	     ('sst','deg C')])
	    

nyrs = (edyr+1)-styr
nmons = nyrs*12
npts = nyrs*73
ndays = (dt.date(edyr+1,1,1).toordinal() - dt.date(styr,1,1).toordinal())

# set up nlons and nlats depending on what we are reading in and out
nlonsIn = 360
nlatsIn= 18 # ERA style to have grids over the poles rather than up to the poles
nlonsOut1 = 360
nlatsOut1 = 180 # ERA style to have grids over the poles rather than up to the poles but this will be changed here with Build or Regrid
nlonsOut5 = 72 # assuming this is correct
nlatsOut5 = 36 # assuming this is correct

#************************************************************
# SUBROUTINES
#************************************************************
# RegridField
def RegridField(TheOldData,TheMDI):
    '''
    This function does a simple regridding of data by averaging over the larger gridboxes
    NO COSINE WEIGHTING FOR LATITUDE!!!!
    
    NOTE: 
    FOR OutputGrid = 5by5 THIS AUTOMATICALLY FLIPS LATITUDE AND ROLLS LONGITUDE TO BE -87.5 to 87.5 and -177,5 to 177.5
    FOR OutputGrid = 1by1 THIS JUST REGRIDS LATITUDE FROM 181 boxes 90 to -90 TO 180 boxes 89.5 to -89.5 and rolls longitude to -179.5 to 179.5
    
    Assumes input grid is always 1by1
    
    INPUTS:
    TheOldData[:,:,:] - time, lat, long numpy array of complete field in original grid resolution
    OUTPUTS:
    TheNewData[:,:,:] - time, lat, long numpy array of complete field in new grid resolution
    
    I'm hoping that things set above are seen by the function rather than being passed explicitly
    
    '''
    
    # Set up the desired output array
    TheNewData = np.empty((len(TheOldData[:,0,0]),36,72),dtype = float)
    TheNewData.fill(TheMDI)
    
    # Then we know we're reading in my converted ERA5 data which has 180 lats and already has lons rolled 180 degrees.
	
    # flip lats to go south to north
    # regrid to 5by5

    # Regrid to 5by5 by simple averaging
    # Data input here should already be 89.5 to -89.5 lat and -179.5 to 179.5 long!!!
    StLt = 0
    EdLt = 0
    # Loop through the OutputGrid (5by5) lats and lons
    for ltt in range(36):
	
    	# create pointers to the five lats to average over
        StLt = np.copy(EdLt)
        EdLt = EdLt + 5

        StLn = 0
        EdLn = 0

        for lnn in range(72):
	
    	    # create pointers to the five lons to average over
            StLn = np.copy(EdLn)
            EdLn = EdLn + 5
	    #print(ltt,lnn,StLt,EdLt,StLn,EdLn)
		    
    	    # Loop over each time point
            for mm in range(len(TheNewData[:,0,0])):
    	    
		# Create a subarr first so that we can deal with missing data
                subarr = TheOldData[mm,StLt:EdLt,StLn:EdLn]			
                gots = np.where(subarr > TheMDI)
                if (len(gots[0]) > 0):
    		
    		    # FILL THE LATITUDES BACKWARDS SO THAT THIS REVERSES THEM!!!
                    TheNewData[mm,35-ltt,lnn] = np.mean(subarr[gots])		    
    		    
    #pdb.set_trace()			
		        
    return TheNewData
    
#************************************************************
# BuildFieldOLD - this is only done at the daily level 
def BuildFieldOLD(BuildType, TheVar, DirStr, InFileStr, TheStYr, TheEdYr,InFileOLD = ''):
    ''' function for building complete reanalyses files over the period specified 
        this can be very computationally expensive so do it by year
	This requires initial reanalysis data to be read in in chunks of 1 month previously

        INPUTS:
	BuildType - 'Update' or 'Build'
	TheVar - string lower case character of q, rh, t, td, dpd, tw, e, msl, sp, ws +2m or appropriate
	InFileStr - string of dir+file string to read in 
	TheStYr = integer start year of data - assume Jan 1st (0101) start 
	TheEdYr = integer end year of data - assume Dec 31st (1231) end	
	InFileOLD = optional string for old data file to be read in
	OUTPUTS:
	TheNewData[:,:,:] - time, lat, long numpy array of complete field in new time resolution
	
	'''
    # Should just be able to read in the data and then append to build complete file
    nyrs = (TheEdYr - TheStYr) + 1

    if (BuildType == 'Build'):
    
        FuncStYr = TheStYr
        NewDataArr = np.array(()) # This will be set up on first read in - this has len() of 0!!!


    elif (BuildType == 'Update'):
    
        FuncStYr = TheEdYr
	
	# Now read in the old data to start array to append to
        NewDataArr,Latitudes,Longitudes = GetGrid4(InFileOLD,[TheVar],['latitude'],['longitude'])
    
    for y in np.arange(FuncStYr,TheEdYr+1):
     
    	# Get actual year we're working on
        print('Working Year: ',y)
    	
        # Loop through each month or pentad depending on thing
        for m in range(12):
		
    	    ## string for file name
            mm = '%02i' % (m+1)

            # Now read in the old data to start array to append to
            TmpDataArr,Latitudes,Longitudes = GetGrid4(DirStr+str(y)+mm+InFileStr,[TheVar],['latitude'],['longitude'])
	    
            if (len(NewDataArr) == 0):
	    
                # this is the start of the build
                NewDataArr = np.copy(TmpDataArr)
		
            else:
	    
                NewDataArr = np.append(NewDataArr,np.copy(TmpDataArr),0)

    #print('Check built array')
    #pdb.set_trace()
    
    return NewDataArr

#************************************************************
# GetDaily - this is only done at the daily level 
def GetDaily(TheVar, DirStr, InFileStr, TheYr):
    ''' function for building complete reanalyses files over the period specified 
        this can be very computationally expensive so do it by year
	This requires initial reanalysis data to be read in in chunks of 1 month previously

        INPUTS:
	TheVar - string lower case character of q, rh, t, td, dpd, tw, e, msl, sp, ws +2m or appropriate
	InFileStr - string of dir+file string to read in 
	TheYr = integer start year of data - assume Jan 1st (0101) start 
	OUTPUTS:
	TheNewData[:,:,:] - time, lat, long numpy array of complete field in new time resolution
	
	'''
    NewDataArr = np.array(()) # This will be set up on first read in - this has len() of 0!!!

    # Loop through each month or pentad depending on thing
    for m in range(12):
		
        ## string for file name
        mm = '%02i' % (m+1)

        # Now read in the old data to start array to append to
        TmpDataArr,Latitudes,Longitudes = GetGrid4(DirStr+str(TheYr)+mm+InFileStr,[TheVar],['latitude'],['longitude'])
	    
        if (len(NewDataArr) == 0):
	    
            # this is the start of the build
            NewDataArr = np.copy(TmpDataArr)
		
        else:
	    
            NewDataArr = np.append(NewDataArr,np.copy(TmpDataArr),0)
    
    return NewDataArr

#************************************************************
# MaskLandS
def MaskLandS(TheDataArray,LandMaskFile,TheMDI):
    ''' This function takes in any field (anoms, clims, actuals etc)
    and returns a land-masked (lsm < 0.5 - do not want inland lakes!!!) for sea and a sea-masked (lsm > 0) for land ''' 

    # Read in land mask file
    LSMArr,Latitudes,Longitudes = GetGrid4(LandMaskFile,['lsm'],['latitude'],['longitude'])
    LSMArr = np.reshape(LSMArr,(1,len(LSMArr[:,0]),len(LSMArr[0,:])))
      
    # loop through each time step and mask
    for tt in range(len(TheDataArray[:,0,0])):
    
        # create a temporary array for this time step with an added time dimension of 1
        TmpArr = np.reshape(TheDataArray[tt,:,:],(1,len(TheDataArray[0,:,0]),len(TheDataArray[0,0,:])))

        # Set up temporary MDI land and sea arrays
        LandTmp = np.empty_like(TmpArr,dtype = float)
        LandTmp[:,:,:] = TheMDI
    
        SeaTmp = np.empty_like(TmpArr,dtype = float)
        SeaTmp[:,:,:] = TheMDI
   
        # Fill the land and sea arrays - trying more inclusive approach
        LandTmp[np.where(LSMArr > 0.)] = TmpArr[np.where(LSMArr > 0.)]
        SeaTmp[np.where(LSMArr < 0.5)] = TmpArr[np.where(LSMArr < 0.5)]
#        LandTmp[np.where(LSMArr >= 0.25)] = TmpArr[np.where(LSMArr >= 0.25)]
#        SeaTmp[np.where(LSMArr <= 0.75)] = TmpArr[np.where(LSMArr <= 0.75)]
	
	# Now build complete array
        if (tt == 0):
	
            LandData = LandTmp
            SeaData = SeaTmp
	    
        else:
	
            LandData = np.append(LandData,np.copy(LandTmp),0)
            SeaData = np.append(SeaData,np.copy(SeaTmp),0)
	    
        #print('Test land sea mask')
        #pdb.set_trace()
	
    LSMArr = 0
	    
    return LandData, SeaData
   
#************************************************************
# MakeMonths
def MakeMonths(TheDataArr,TheStYr,TheEdYr):
    ''' THis function takes a 3D array of daily data and
    averages it up to monthly means '''
    
    # Set up date bits and bobs
    nyrs = TheEdYr - TheStYr + 1
    nmns = nyrs*12
    
    # Set up empty arrays
    MonthlyData = np.empty((nmns,len(TheDataArr[0,:,0]),len(TheDataArr[0,0,:])))

    mncounter = 0 # counter for months    
    mnst = 0 # counter for days in months
    mned = 0#

    # Loop through each year
    for y in np.arange(TheStYr, TheEdYr+1):	
	    
        for m in np.arange(1,13):	    
            
            if (m < 12):
                    
                mndays = (dt.date(y,m+1,1).toordinal() - dt.date(y,m,1).toordinal())

            else:
                    
                mndays = (dt.date(y+1,1,1).toordinal() - dt.date(y,m,1).toordinal())

            mnst = mned
            mned = mned + mndays

            MonthlyData[mncounter,:,:] = np.mean(TheDataArr[mnst:mned,:,:],0)
            mncounter += 1
		    
#                    if (mncounter == 1):
#                    
#                        print('Test Months: ')
#                        pdb.set_trace()



    
#    # Loop through each gridbox
#    for ln in range(len(TheDataArr[0,0,:])):
#    
#        for lt in range(len(TheDataArr[0,:,0])):
#	
#            mncounter = 0 # counter for months    
#            mnst = 0 # counter for days in months
#            mned = 0#
#
#	    # Loop through each year
#            for y in np.arange(TheStYr, TheEdYr+1):	
#	    
#                for m in np.arange(1,13):	    
#            
#                    if (m < 12):
#                    
#                        mndays = (dt.date(y,m+1,1).toordinal() - dt.date(y,m,1).toordinal())#
#
#                    else:
#                    
#                        mndays = (dt.date(y+1,1,1).toordinal() - dt.date(y,m,1).toordinal())#
#
#                    mnst = mned
#                    mned = mned + mndays#
#
#                    MonthlyData[mncounter,lt,ln] = np.mean(TheDataArr[mnst:mned,lt,ln])
#                    mncounter += 1
#		    
##                    if (mncounter == 1):
##                    
##                        print('Test Months: ')
##                        pdb.set_trace()
   
    return MonthlyData
   
#************************************************************
# MakePentads
def MakePentads(TheDataArr,TheStYr,TheEdYr):
    ''' This function takes a 3D array of daily data and 
    averages it up to pentad means '''   

    # Set up date bits and bobs
    nyrs = TheEdYr - TheStYr + 1
    npts = nyrs*73
    
    # Set up empty arrays
    PentadData = np.empty((npts,len(TheDataArr[0,:,0]),len(TheDataArr[0,0,:])))
    
    # Loop through each gridbox
    for ln in range(len(TheDataArr[0,0,:])):
    
        for lt in range(len(TheDataArr[0,:,0])):
	
            ptst = 0 # counter for pentads
            pted = 0
            yrst = 0 # pointer for start of year
            yred = 0 # pointer for end of year
	    
	    # Loop through each year
            for y in np.arange(TheStYr, TheEdYr+1):	
            
                yrdays = (dt.date(y+1,1,1).toordinal() - dt.date(y,1,1).toordinal())
                # these don't change when yred changes so do not need to be copies
                yrst = yred
                yred = yrst + yrdays
		
                ptst = pted
                pted = pted + 73
	
	        # Easy if its not a leap year
                if (yrdays == 365):	
		
                    tmpdata = np.reshape(TheDataArr[yrst:yred,lt,ln],(73,5))
		    
		    # get means over each 5 day period    
                    PentadData[ptst:pted,lt,ln] = np.mean(tmpdata,1)
		    
#                    if (y == TheStYr):
#                        
#                        print('Test year pentads: ')
#                        pdb.set_trace()
		    
                elif (yrdays == 366): # a LEAP year
		
                    tmpdata = TheDataArr[yrst:yred,lt,ln]
		    
		    
		    # get means over each 5 day period    
                    PentadData[ptst:ptst+11,lt,ln] = np.mean(np.reshape(tmpdata[0:55],(11,5)),1)
                    PentadData[ptst+11,lt,ln] = np.mean(tmpdata[55:61])
                    PentadData[ptst+12:pted,lt,ln] = np.mean(np.reshape(tmpdata[61:366],(61,5)),1)
                    
                    #print('Test LEAP year pentads: ',y,yrdays)
                    #pdb.set_trace()
		
    return PentadData
        
#************************************************************
# CreateAnoms
def CreateAnoms(TheClimSt,TheClimEd,TheStYr,TheEdYr,TheInData,TheMDI):
    '''
    This function takes any grid, any var, any time resolution and computes climatologies/stdevs over given period and then anomalies
    
    INPUTS:
    TheClimSt - interger start year of climatology Always Jan start
    TheClimEd - integer end  year of climatology Always Dec end
    TheStYr - integer start year of data to find climatology
    TheEdYr - integer end year of data to find climatology
    TheInData[:,:,:] - time, lat, lon array of actual values
    
    OUTPUTS:
    AllAnomsArr[:,:,:] - time, lat, lon array of anomalies
    ClimsArr[:,:,:] - time, lat, lon array of climatologies
    StDevsArr[:,:,:] - time, lat, lon array of stdeviations
    
    '''

    # Set up for time
    nyrs = TheEdYr - TheStYr + 1
    ntims = len(TheInData[:,0,0])

    if (ntims/nyrs == 12): # rough check - should be 12 obvs
     
        # Monthly data
        nclims = 12
	
    elif (ntims/nyrs == 73):
    
        # Pentad data
        nclims = 73
	
    elif (ntims/nyrs > 300): # rough check as leap yrs make it screwy		
    
        # Daily data
        nclims = 365 # Feb 29th will just be ignored to create a clim, use Feb 28th for anomaly	     
    	    
    # first create empty arrays
    AllAnomsArr = np.empty_like(TheInData)
    AllAnomsArr.fill(mdi)
    ClimsArr = np.copy(AllAnomsArr[0:nclims,:,:])
    StDevsArr = np.copy(AllAnomsArr[0:nclims,:,:])
        
    # loop through gridboxes
    for lt in range(len(TheInData[0,:,0])):

        for ln in range(len(TheInData[0,0,:])):
	
	    # if monthly or pentad then easy peasy
            if (nclims < 300):

	        # pull out gridbox and reform to years by nclims (months or pentads)
                SingleSeries = np.reshape(TheInData[:,lt,ln],(nyrs,nclims)) # nyrs rows, nclims columns
            
                ClimsArr[:,lt,ln] = np.mean(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,:],0)
                StDevsArr[:,lt,ln] = np.std(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,:],0)
	    
                SingleSeriesAnoms = SingleSeries - ClimsArr[:,lt,ln]
	    
                #print('Test pentad and monthlyanomalies, clims and stdevs')
                #pdb.set_trace()

                AllAnomsArr[:,lt,ln] = np.copy(np.reshape(SingleSeriesAnoms,ntims))	    

            # dailies more tricky
            else:

	        # pull out gridbox and reform to years by nclims (months or pentads)
                TmpSeries = TheInData[:,lt,ln]
		
		# Make an nyrs (rows) by 366 days array
                SingleSeries = np.empty((nyrs,366),dtype = float)
                SingleSeries.fill(TheMDI)
		
		# Now populate year by year
		# Non leap years have a missing 29th Feb
                StDy = 0
                EdDy = 0
		
                for y in range(nyrs):
		
                    # Is this a leap year?
                    TotalDays = dt.date(y+TheStYr+1,1,1).toordinal() - dt.date(y+TheStYr,1,1).toordinal()     
                    
                    #print('Getting daily clims: ',y+TheStYr, TotalDays)
                    StDy = EdDy
                    EdDy = StDy + TotalDays
		    
                    if (TotalDays == 366):
		
                        SingleSeries[y,:] = TmpSeries[StDy:EdDy]

                    else:
		    
                        SingleSeries[y,0:59] = TmpSeries[StDy:StDy+59]
                        SingleSeries[y,60:] = TmpSeries[StDy+59:EdDy]
                        	 				    					      
                # Calculate clims and stdevs for each of the 365 days ignoring Feb 29th
                ClimsArr[0:59,lt,ln] = np.mean(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,0:59],0)
                StDevsArr[0:59,lt,ln] = np.std(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,0:59],0)
                ClimsArr[59:,lt,ln] = np.mean(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,60:],0)
                StDevsArr[59:,lt,ln] = np.std(SingleSeries[TheClimSt-TheStYr:(TheClimEd-TheStYr)+1,60:],0)
	    
                #print('Got the clims, stdevs')
                #pdb.set_trace()
	    
		# Now repopulate year by year
		# Non leap years have a missing 29th Feb
                StDy = 0
                EdDy = 0
                SingleSeriesAnoms = np.empty(ntims,dtype = float)
		
                for y in range(nyrs):
		
                    # Is this a leap year?
                    TotalDays = dt.date(y+TheStYr+1,1,1).toordinal() - dt.date(y+TheStYr,1,1).toordinal()     
                    StDy = EdDy
                    EdDy = StDy + TotalDays
                    #print('Getting daily anoms: ',y+TheStYr, TotalDays,StDy,EdDy)
		    
                    if (TotalDays == 366):
		
                        SingleSeriesAnoms[StDy:StDy+59] = SingleSeries[y,0:59] - ClimsArr[0:59,lt,ln]
                        # Use Feb 28th as climatology for 29th
                        SingleSeriesAnoms[59] = SingleSeries[y,59] - ClimsArr[58,lt,ln]
                        SingleSeriesAnoms[StDy+60:EdDy] = SingleSeries[y,60:] - ClimsArr[59:,lt,ln]

                    else:
		    
                        SingleSeriesAnoms[StDy:StDy+59] = SingleSeries[y,0:59] - ClimsArr[0:59,lt,ln]
                        SingleSeriesAnoms[StDy+59:EdDy] = SingleSeries[y,60:] - ClimsArr[59:,lt,ln]
	    
                    #print('Test daily anomalies, clims and stdevs',y+TheStYr, TotalDays)
                    #pdb.set_trace()
	    
                AllAnomsArr[:,lt,ln] = np.copy(SingleSeriesAnoms)	    
	    	    
    return AllAnomsArr, ClimsArr, StDevsArr
    
#************************************************************
# WriteNetCDF
def WriteNetCDF(Filename,TheOutputTime, TheOutputGrid, TheOutputVar, TheFullArray, TheFullArrayAnoms, TheLandArrayAnoms, TheOceanArrayAnoms, TheClimsArray, TheStDevsArray,
            TheStYr, TheEdYr, TheClimStart, TheClimEnd, TheName, TheStandardName, TheLongName, TheUnit, TheMDI):
    '''
    This function writes out a NetCDF 4 file
    
    NOTE: 
    All 1by1 outputs will have lats 89.5 to -89.5 and lons -179.5 to 179.5
    All 5by5 outputs will have lats -87.5 to 87.5 and lons -177.5 to 177.5
    
    INPUTS:
    FileOut - string file name
    TheOutputTime - string monthly or pentad or daily
    TheOutputGrid - string 1by1 or 5by5
    TheOutputVar - string lower case variable name
    TheFullArray[:,:,:] - time, lat, lon array of actual values
    TheFullArrayAnoms[:,:,:] - time, lat, lon array of anomalies
    TheLandArrayAnoms[:,:,:] - time, lat, lon array of land anomalies
    TheOceanArrayAnoms[:,:,:] - time, lat, lon array of ocean anomalies 
    TheClimsArray[:,:,:] - time(12 or 73), lat, lon array of climatology 
    TheStDevsArray[:,:,:] - time(12 or 73, lat, lon array of st devs
    TheStYr - integer start year assumes Jan start 
    TheEdYr - integer end year assumes Dec start 
    TheClimStart - integer start of clim Jan start 
    TheClimEnd - integer end of clim Dec start
    TheName - string short name of var q2m
    TheStandardName - string standard name of variable
    TheUnit - string unit of variable
    TheMDI - missing data indicator
    OUTPUTS:
    None
    
    '''
    
    # Sort out times in days since 1979-01-01    
    # Sort out climatology time
    if (TheOutputTime == 'monthly'):
        nClims = 12
        TimPoints = MakeDaysSince(TheStYr,1,TheEdYr,12,'month') # use 'day','month','year'
    elif (TheOutputTime == 'pentad'):
        nClims = 73
        TimPoints = MakeDaysSince(TheStYr,1,TheEdYr,73,'pentad') # use 'day','month','year'
    elif (TheOutputTime == 'daily'):
        nClims = 365
        TimPoints = MakeDaysSince(TheStYr,1,TheEdYr,365,'day') # will work even if its a leap year
    nTims = len(TimPoints)
    		    
    # Sort out Lats, Lons and LatBounds and LonBounds
    if (TheOutputGrid == '1by1'):
        LatList = np.flip(np.arange(180)-89.5)
        LonList = np.arange(360)-179.5

        LatBounds = np.empty((len(LatList),2),dtype='float')
        LonBounds = np.empty((len(LonList),2),dtype='float')

        LatBounds[:,0] = LatList + ((LatList[0]-LatList[1])/2.)
        LatBounds[:,1] = LatList - ((LatList[0]-LatList[1])/2.)

        LonBounds[:,0] = LonList - ((LonList[1]-LonList[0])/2.)
        LonBounds[:,1] = LonList + ((LonList[1]-LonList[0])/2.)    

        nlatsOut = 180
        nlonsOut = 360

    elif (TheOutputGrid == '5by5'):
        LatList = (np.arange(36)*5)-87.5
        LonList = (np.arange(72)*5)-177.5
    
        LatBounds = np.empty((len(LatList),2),dtype='float')
        LonBounds = np.empty((len(LonList),2),dtype='float')

        LatBounds[:,0] = LatList - ((LatList[1]-LatList[0])/2.)
        LatBounds[:,1] = LatList + ((LatList[1]-LatList[0])/2.)

        LonBounds[:,0] = LonList - ((LonList[1]-LonList[0])/2.)
        LonBounds[:,1] = LonList + ((LonList[1]-LonList[0])/2.)    
	
        nlatsOut = 36
        nlonsOut = 72

    # No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
    # Not sure what this does to float precision though...

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    ncfw = Dataset(Filename,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first

    # Set up the dimension names and quantities
    ncfw.createDimension('time',nTims)
    ncfw.createDimension('latitude',nlatsOut)
    ncfw.createDimension('longitude',nlonsOut)
    
    # If there are climatologies to be written then also set up clim dimension
    if (len(np.shape(TheClimsArray)) > 1):
        
        if (TheOutputTime == 'monthly'):
            ncfw.createDimension('month_time',nClims)
        elif (TheOutputTime == 'pentad'):
            ncfw.createDimension('pentad_time',nClims)
        elif (TheOutputTime == 'daily'):
            ncfw.createDimension('day_time',nClims)

    # Go through each dimension and set up the variable and attributes for that dimension if needed
    MyVarT = ncfw.createVariable('time','f4',('time',))
    MyVarT.standard_name = 'time'
    MyVarT.long_name = 'time'
    MyVarT.units = 'days since 1979-1-1 00:00:00'
    MyVarT.start_year = str(TheStYr)
    MyVarT.end_year = str(TheEdYr)
    MyVarT[:] = TimPoints

    MyVarLt = ncfw.createVariable('latitude','f4',('latitude',))
    MyVarLt.standard_name = 'latitude'
    MyVarLt.long_name = 'gridbox centre latitude'
    MyVarLt.units = 'degrees_north'
    MyVarLt[:] = LatList

    MyVarLn = ncfw.createVariable('longitude','f4',('longitude',))
    MyVarLn.standard_name = 'longitude'
    MyVarLn.long_name = 'gridbox centre longitude'
    MyVarLn.units = 'degrees_east'
    MyVarLn[:] = LonList

    # If there are climatologies to be written then also set up clim dimension
    if (len(np.shape(TheClimsArray)) > 1):
        
        if (TheOutputTime == 'monthly'):
            MyVarM = ncfw.createVariable('month_time','i4',('month_time',))
            MyVarM.long_name = 'months of the year'
            MyVarM.units = 'months'
            MyVarM[:] = np.arange(nClims)
        elif (TheOutputTime == 'pentad'):
            MyVarM = ncfw.createVariable('pentad_time','i4',('pentad_time',))
            MyVarM.long_name = 'pentads of the year'
            MyVarM.units = 'pentads'
            MyVarM[:] = np.arange(nClims)
        elif (TheOutputTime == 'daily'):
            MyVarM = ncfw.createVariable('day_time','i4',('day_time',))
            MyVarM.long_name = 'days of the year'
            MyVarM.units = 'days'
            MyVarM[:] = np.arange(nClims)

    # Go through each variable and set up the variable attributes
    # I've added zlib=True so that the file is in compressed form
    # I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
    MyVarD = ncfw.createVariable(TheName,'f4',('time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
    MyVarD.standard_name = TheStandardName
    MyVarD.long_name = TheLongName
    MyVarD.units = TheUnit
    MyVarD.valid_min = np.min(TheFullArray)
    MyVarD.valid_max = np.max(TheFullArray)
    MyVarD.missing_value = TheMDI
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarD[:,:,:] = TheFullArray[:,:,:]  	

    # If there are climatologies etc to be written then also set them up
    if (len(np.shape(TheClimsArray)) > 1):
    
        MyVarA = ncfw.createVariable(TheName+'_anoms','f4',('time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        MyVarA.standard_name = TheStandardName+'_anomalies'
        MyVarA.long_name = TheLongName+' anomalies from 1981-2010'
        MyVarA.units = TheUnit
        MyVarA.valid_min = np.min(TheFullArrayAnoms)
        MyVarA.valid_max = np.max(TheFullArrayAnoms)
        MyVarA.missing_value = TheMDI
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarA[:,:,:] = TheFullArrayAnoms[:,:,:]  	
    
        MyVarAL = ncfw.createVariable(TheName+'_anoms_land','f4',('time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        MyVarAL.standard_name = TheStandardName+'_anomalies'
        MyVarAL.long_name = TheLongName+' land anomalies from 1981-2010'
        MyVarAL.units = TheUnit
        MyVarAL.valid_min = np.min(TheLandArrayAnoms)
        MyVarAL.valid_max = np.max(TheLandArrayAnoms)
        MyVarAL.missing_value = TheMDI
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarAL[:,:,:] = TheLandArrayAnoms[:,:,:]  	

        MyVarAO = ncfw.createVariable(TheName+'_anoms_ocean','f4',('time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        MyVarAO.standard_name = TheStandardName+'_anomalies'
        MyVarAO.long_name = TheLongName+' ocean anomalies from 1981-2010'
        MyVarAO.units = TheUnit
        MyVarAO.valid_min = np.min(TheOceanArrayAnoms)
        MyVarAO.valid_max = np.max(TheOceanArrayAnoms)
        MyVarAO.missing_value = TheMDI
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarAO[:,:,:] = TheOceanArrayAnoms[:,:,:]  	

        if (TheOutputTime == 'monthly'):
            MyVarC = ncfw.createVariable(TheName+'_clims','f4',('month_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'pentad'):
            MyVarC = ncfw.createVariable(TheName+'_clims','f4',('pentad_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'daily'):
            MyVarC = ncfw.createVariable(TheName+'_clims','f4',('day_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        MyVarC.standard_name = TheStandardName+'_climatologies'
        MyVarC.long_name = TheLongName+' climatology over 1981-2010'
        MyVarC.units = TheUnit
        MyVarC.valid_min = np.min(TheClimsArray)
        MyVarC.valid_max = np.max(TheClimsArray)
        MyVarC.missing_value = TheMDI
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarC[:,:,:] = TheClimsArray[:,:,:]  	

        if (TheOutputTime == 'monthly'):
            MyVarS = ncfw.createVariable(TheName+'_stdevs','f4',('month_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'pentad'):
            MyVarS = ncfw.createVariable(TheName+'_stdevs','f4',('pentad_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'daily'):
            MyVarS = ncfw.createVariable(TheName+'_stdevs','f4',('day_time','latitude','longitude',),fill_value = TheMDI,zlib=True,least_significant_digit=4)
        MyVarS.standard_name = TheStandardName+'_climatological_standard_deviations'
        MyVarS.long_name = TheLongName+' climatological standard deviation over 1981-2010'
        MyVarS.units = TheUnit
        MyVarS.valid_min = np.min(TheStDevsArray)
        MyVarS.valid_max = np.max(TheStDevsArray)
        MyVarS.missing_value = TheMDI
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarS[:,:,:] = TheStDevsArray[:,:,:]  	

    ncfw.close()

    return
    
#************************************************************
# MAIN
#************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--var', dest='var', action='store', default='q', type=str,
                        help='Variable [q]')
## Set up output variables - for q, e, RH, dpd, Tw we will need to read in multiple input files
#Var = 'q' # this can be 't','td','q','rh','e','dpd','tw','ws','slp','sp','uv','sst'
    parser.add_argument('--freq', dest='freq', action='store', default='monthly', type=str,
                        help='Time Frequency [monthly]')
 
    args = parser.parse_args()         

    print(args)

    Var = args.var
    Freq = args.freq

    InputERA = '_daily_'+VarDict[Var][1]+'.nc'

    # only needed if we're updating...
    OldERAStrPT = workingdir+'/OTHERDATA/'+NameDict[Var]+'_pentad_1by1_'+ThisRean+'_1979'+str(edyr-1)+'.nc'
    OldERAStrMN = workingdir+'/OTHERDATA/'+NameDict[Var]+'_monthly_1by1_'+ThisRean+'_1979'+str(edyr-1)+'.nc'

    if (Freq == 'pentad'):

#    NewERAStrD1   = workingdir+'/OTHERDATA/'+NameDict[Var]+'_daily_1by1_'+ThisRean+'_1979'+updateyyyy+'.nc'
        NewERAStrP1   = workingdir+'/OTHERDATA/'+NameDict[Var]+'_pentad_1by1_'+ThisRean+'_1979'+updateyyyy+'.nc'

    else:
        
        NewERAStrM1   = workingdir+'/OTHERDATA/'+NameDict[Var]+'_monthly_1by1_'+ThisRean+'_1979'+updateyyyy+'.nc'
        NewERAStrM5   = workingdir+'/OTHERDATA/'+NameDict[Var]+'_monthly_5by5_'+ThisRean+'_1979'+updateyyyy+'.nc'
    
    # What are we working on?
    print('Working variable: ',Var)
    print('Working frequency: ',Freq)
    print('Type of run: ',ThisProg, styr, edyr)
    print('Reanalysis: ',ThisRean)

    # If its an update read in old pentad fields and old monthly fields to append to,
    # then read in year of daily data, process to pentad and monthly and then append
    # If its a total build then read in year by year, process to pentads and monthlies and append

    if (Freq == 'pentad'):

        if (ThisProg == 'Build'):
    
            PentadDataArr = np.array(()) # This will be set up on first read in - this has len() of 0!!!

        elif (ThisProg == 'Update'):
    
	    # Now read in the old data to start array to append to
            PentadDataArr,Latitudes,Longitudes = GetGrid4(InFileOLDPT,[TheVar],['latitude'],['longitude'])

    else:
       
        if (ThisProg == 'Build'):
    
            MonthDataArr = np.array(()) # This will be set up on first read in - this has len() of 0!!!

        elif (ThisProg == 'Update'):
    
            # Now read in the old data to start array to append to
            MonthDataArr,Latitudes,Longitudes = GetGrid4(InFileOLDMN,[TheVar],['latitude'],['longitude'])
    
    
    # Loop through the years
    for years in range(styr,edyr+1):
    
    	# Get actual year we're working on
        print('Working Year: ',years)
    	
        # read in a year of data
        # First make up the full daily field for one year - build or update
        DailyData = GetDaily(NameDict[Var],workingdir+'/OTHERDATA/ERA5/',InputERA,years)
        print('Data all read in')

# TOO TIME/MEMORY HUNGRY
#    ## Now make daily anoms, clims, stdevs
#    DailyAnoms,DailyClims,DailyStdevs = CreateAnoms(ClimStart,ClimEnd,styr,edyr,DailyData,mdi)
#
#    print('Created daily anomalies')
#
#    # Now get land and sea masked pentad anomalies, actuals, clims and stdevs
#    DailyAnomsLand, DailyAnomsSea = MaskLandS(DailyAnoms,LandMask,mdi)
#    #DailyLand, DailySea = MaskLandS(DailyData,LandMask,mdi)
#    #DailyClimsLand, DailyClimsSea = MaskLandS(DailyClims,LandMask,mdi) 
#    #DailyStdevsLand, DailyStdevsSea = MaskLandS(DailyStdevs,LandMask,mdi)
#
#    print('Created daily land and sea masks')
#
#    # Now save to nc
#    WriteNetCDF(NewERAStrD1,'daily', '1by1', Var, DailyData, DailyAnoms, DailyAnomsLand, DailyAnomsSea, DailyClims, DailyStdevs,
#            styr, edyr, ClimStart, ClimEnd, NameDict[Var], StandardNameDict[Var], LongNameDict[Var], UnitDict[Var],mdi) 
#    print('Written out dailies')
#    # clean up
#    DailyAnoms = 0
#    DailyAnomsLand = 0
#    DailyAnomsSea = 0
#    DailyClims = 0
#    DailyStdev = 0

        if (Freq == 'pentad'):

            # Now make pentads
            PentadData = MakePentads(DailyData,years,years)
            # clean up
            DailyData = 0
            print('Created pentads')

            if (len(PentadDataArr) == 0):
	    
                # this is the start of the build
                PentadDataArr = np.copy(PentadData)
		
            else:
	    
                PentadDataArr = np.append(PentadDataArr,np.copy(PentadData),0)

        else:
	
            # Now make monthlies
            MonthData = MakeMonths(DailyData,years,years)
            # clean up
            DailyData = 0
            print('Created months')
            #pdb.set_trace()
	 
            if (len(MonthDataArr) == 0):
	    
                # this is the start of the build
                MonthDataArr = np.copy(MonthData)
		
            else:
	    
                MonthDataArr = np.append(MonthDataArr,np.copy(MonthData),0)

    if (Freq == 'pentad'):

        # Now make pentad anoms, clims, stdevs
        PentadAnoms,PentadClims,PentadStdevs = CreateAnoms(ClimStart,ClimEnd,actstyr,edyr,PentadDataArr,mdi)
        print('Created pentad anomalies')
    #    pdb.set_trace()

        # Now get land and sea masked pentad anomalies, actuals, clims and stdevs
        PentadAnomsLand, PentadAnomsSea = MaskLandS(PentadAnoms,LandMask,mdi)
        #PentadLand, PentadSea = MaskLandS(PentadData,LandMask,mdi)
        #PentadClimsLand, PentadClimsSea = MaskLandS(PentadClims,LandMask,mdi)
        #PentadStdevsLand, PentadStdevsSea = MaskLandS(PentadStdevs,LandMask,mdi)

        print('Created pentad land and sea masks')

        # Now save to nc
        WriteNetCDF(NewERAStrP1,'pentad', '1by1', Var, PentadDataArr, PentadAnoms, PentadAnomsLand, PentadAnomsSea, PentadClims, PentadStdevs,
            actstyr, edyr, ClimStart, ClimEnd, NameDict[Var], StandardNameDict[Var], LongNameDict[Var], UnitDict[Var],mdi) 

        # clean up
        PentadAnoms = 0
        PentadAnomsLand = 0
        PentadAnomsSea = 0
        PentadClims = 0
        PentadStdev = 0

        print('Written out pentads')

    else:
 
        # Now make monthly anoms, clims, stdevs
        MonthAnoms,MonthClims,MonthStdevs = CreateAnoms(ClimStart,ClimEnd,actstyr,edyr,MonthDataArr,mdi)

        print('Created month anomalies')

        # Now get land and sea masked monthly anomalies, actuals, clims and stdevs
        MonthAnomsLand, MonthAnomsSea = MaskLandS(MonthAnoms,LandMask,mdi)
        #MonthLand, MonthSea = MaskLandS(MonthData,LandMask,mdi)
        #MonthClimsLand, MonthClimsSea = MaskLandS(MonthClims,LandMask,mdi)
        #MonthStdevsLand, MonthStdevsSea = MaskLandS(MonthStdevs,LandMask,mdi)

        print('Created month land and sea masks')

        # Now save to nc
        WriteNetCDF(NewERAStrM1,'monthly', '1by1', Var, MonthDataArr, MonthAnoms, MonthAnomsLand, MonthAnomsSea, MonthClims, MonthStdevs,
            actstyr, edyr, ClimStart, ClimEnd, NameDict[Var], StandardNameDict[Var], LongNameDict[Var], UnitDict[Var],mdi) 

        print('Written out months')

        # Now make monthly 5by5s - regrid
        Month5by5Data = RegridField(MonthDataArr,mdi)
        MonthDataArr = 0
        print('Created month 5by5s')
        #pdb.set_trace()

        Month5by5Anoms = RegridField(MonthAnoms,mdi)
        MonthAnoms = 0
        Month5by5Clims = RegridField(MonthClims,mdi)
        MonthClims = 0       
        Month5by5Stdevs = RegridField(MonthStdevs,mdi)
        MonthStdevs = 0

        #Month5by5Land = RegridField(MonthLand,mdi)
        Month5by5AnomsLand = RegridField(MonthAnomsLand,mdi)
        MonthAnomsLand = 0
        #Month5by5ClimsLand = RegridField(MonthClimsLand,mdi)
        #Month5by5StdevsLand = RegridField(MonthStdevsLand,mdi)

        #Month5by5Sea = RegridField(MonthDataSea,mdi)
        Month5by5AnomsSea = RegridField(MonthAnomsSea,mdi)
        MonthAnomsSea = 0
        #Month5by5ClimsSea = RegridField(MonthClimsSea,mdi)
        #Month5by5StdevsSea = RegridField(MonthStdevsSea,mdi)

        # Now save to nc
        WriteNetCDF(NewERAStrM5,'monthly', '5by5', Var, Month5by5Data, Month5by5Anoms, Month5by5AnomsLand, Month5by5AnomsSea, Month5by5Clims, Month5by5Stdevs,
            actstyr, edyr, ClimStart, ClimEnd, NameDict[Var], StandardNameDict[Var], LongNameDict[Var], UnitDict[Var],mdi) 

        print('Written out month 5by5s')
    
    print('And we are done!')
