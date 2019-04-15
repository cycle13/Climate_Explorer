#!/usr/local/sci/bin/python
# PYTHON3.6.1
# 
# Author: Kate Willett
# Created: 18 Jul 2018
# Last update: 15 Apr 2019
# Location: /data/local/hadkw/HADCRUH2/UPDATE2017/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/HadISDH_Build					
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# THIS CODE DOES MANY THINGS BUT ONLY ONE THING AT A TIME! SO RE-RUN FOR MULTIPLE THINGS
# NOTE THAT FOR ANY 1BY1 OUTPUT IT REGRIDS TO BE 89.5 to -89.5 rather than 90 - -90 (180 boxes rather than 181!!!)
# AND ROLLS LONGITUDE TO -179.5 to 179.5
#
# AT THE MOMENT THIS ASSUMES COMPLETE FIELDS SO WON'T WORK FOR SST!!! 
#
# ANOTHER ISSUE IS LAND / SEA MASKING - TOO MUCH LAND COVER, TOO MUCH SEA COVER - SO THERE WILL BE CONTAMINATION!
# I COMPUTE ANOMALIES AT 1by1 RES BEFORE REGRIDDING TO 5by5 TO MINIMISE THIS>
#
#
# This code reads in the ERA-Interim months of 1by1 6 hourly or monthly variables
# (e.g., T, Td and Surface Pressure etc) for the full time period
#
# If desired it converts to humidity variables 
# If desired it averages to monthly means and saves to netCDF:
#	days since 19790101 (float), 181 lats 90 to -90, 360 lons 0 to 359, <var>2m
# If desired it regrids to 5by5 (monthly means only) and saves to netCDF
#	days since 19790101 (int), 36 lats -87.5 to 87.5, 72 lons -177.5 to 177.5, actuals
# If desired it calculates anomalies over a climatological references period given (default 1981-2010) 
# and saves to netCDF
# For anomalies it also creates a land only and ocean only set of grids to save along side the complete 
# fields
#	days since 19790101 (int), 36 lats -87.5 to 87.5, 72 lons -177.5 to 177.5, anomalies, 
#	anomalies_land, anomalies_sea
#
# The ERA-Interim updates have to be downloaded from ERADownload.py code
# This requires a key to be set up in .ecmwfapirc annually - obtained from logging in to ECMWF
# https://confluence.ecmwf.int/display/WEBAPI/How+to+retrieve+ECMWF+Public+Datasets
# It also requires ecmwfapi to be downloaded and in the directory as you are running to code from
#
# The ERA5 updates have to be downloaded using ERA5Download.py which is in cdsapi-0.1.3/

# Each time you download change the filename to ERAINTERIM_6hr_1by1_MMYYYY.nc
# Save to /data/local/hadkw/HADCRUH2/UPDATE<yyyy>/OTHERDATA/
# Copy previous years of monthly ERAINTERIM data from the previous 
# UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_1by1_ERA-Interim_data_1979<yyyy>.nc
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
# ERA-Interim 1by1 6 hrly gridded data
# ERA<Mmm> = /data/local/hadkw/HADCRUH2/UPDATE<yyyy>/OTHERDATA/ERAINTERIM_<var>_6hr_1by1_<MMYYYY>.nc
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# First make sure the New ERA-Interim data are in the right place.
# Also check all editables in this file are as you wish
# python2.7 ExtractMergeRegridERA_JUL2018.py
#
# -----------------------
# OUTPUT
# -----------------------
# New ERA-Interim 1by1 monthly gridded data for 1979 to present
# NewERA<var> = /data/local/hadkw/HADCRUH2/UPDATE<yyyy>/OTHERDATA/<var>2m_monthly_1by1_ERA-Interim_data_1979<yyyy>.nc
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 1 (18 Jul 2018)
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
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import date2num,num2date
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
from math import sqrt,pi
import struct
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
styr       = 1979
edyr       = 2018
#edOLD      = (edyr-styr)*12
stmon      = 1
edmon      = 12

# Set up output variables - for q, e, RH, dpd, Tw we will need to read in multiple input files
OutputVar = 'dpd' # this can be 't','td','q','rh','e','dpd','tw','ws','slp','sp','uv','sst'

# Is this a new run or an update?
ThisProg = 'Regrid' 
# Update for updating an existing file (1by1 monthly or pentad)
# Build for building from scratch (1by1 6hr to 1by1 monthly or pentad)
	# THIS AUTOMATICALLY REGRIDS LATS TO BE 180 RATHER THAN 181!!!
# Regrid for changing spatial res from 1by1 to 6hr
	# IF OutputGrid = 1by1 then this just changes lats from 181 to 180
	# IF OutputGrid = 5by5 then this changes to 36 lats (-87.5 to 87.5) and 72 lons (-177.5 to 177.5)

# Is this ERA-Interim or ERA5?
ThisRean = 'ERA-Interim' # 'ERA5' or 'ERA-Interim'

# Are you reading in hourlies or monthlies?
ReadInTime = 'monthly' # this can be '1hr', '6hr' or 'month' or maybe 'day' later

# Are you converting to monthlies? We will output hourlies anyway if they are read in
OutputTime = 'monthly' # this could be 'monthly' or 'pentad'

# Are you reading in 1by1 or 5by5? We will output 1by1 anyway if they are read in.
ReadInGrid = '1by1' # this can be '1by1' or '5by5'

# Are you converting to 5by5?
OutputGrid = '5by5' # this can be '1by1' or '5by5'

# Do you want to create anomalies and if so, what climatology period? We will output absolutes anyway
MakeAnoms = 1 # 1 for create anomalies (and clim and stdev fields), 0 for do NOT create anomalies
ClimStart = 1981 # any year but generally 1981
ClimEnd = 2010 # any year but generally 2010

### END OF EDITABLES ################

# Set up file paths and other necessary things
if (MakeAnoms == 1): # set the filename string for anomalies
    AnomsStr = 'anoms'+str(ClimStart)+'-'+str(ClimEnd)+'_'
else:
    AnomsStr = ''

# Set up file locations
updateyy  = str(edyr)[2:4]
updateyyyy  = str(edyr)
workingdir  = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+updateyyyy

if (ReadInGrid == '5by5'):
    LandMask = workingdir+'/OTHERDATA/HadCRUT.4.3.0.0.land_fraction.nc' # 0 = 100% sea, 1 = 100% land - no islands!, latitude, longitude, land_area_fraction, -87.5 to 87.5, -177.5 to 177.5
elif (ReadInGrid == '1by1'):
    LandMask = workingdir+'/OTHERDATA/lsmask.nc' # 1 = sea, 0 = land - no islands! lat, lon, mask 89.5 to -89.5Lat, 0.5 to 359.5 long 

if (OutputVar in ['t','td']): # these are the simple ones that do not require conversion
    InputERA = ThisRean+'_'+ReadInGrid+'_'+ReadInTime+'_'+OutputVar+'2m_'
    if (ThisProg == 'Update'):
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+ThisRean+'_data_1979'+str(edyr-1)+'.nc'
    else:
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+AnomsStr+ThisRean+'_data_1979'+str(edyr)+'.nc'
    NewERAStr   = OutputVar+'2m_'+OutputGrid+'_'+OutputTime+'_'+AnomsStr+ThisRean+'_data_1979'+updateyyyy+'.nc'

elif (OutputVar in ['ws','uv']):
    InputERA = ThisRean+'_'+ReadInGrid+'_'+ReadInTime+'_'+OutputVar+'10m_'
    if (ThisProg == 'Update'):
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+ThisRean+'_data_1979'+str(edyr-1)+'.nc'
    else:
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+AnomsStr+ThisRean+'_data_1979'+str(edyr)+'.nc'
    NewERAStr   = OutputVar+'10m_'+OutputGrid+'_'+OutputTime+'_'+AnomsStr+ThisRean+'_data_1979'+updateyyyy+'.nc'

elif (OutputVar in ['slp','sp','sst']):
    InputERA = ThisRean+'_'+ReadInGrid+'_'+ReadInTime+'_'+OutputVar+'_'
    if (ThisProg == 'Update'):
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+ThisRean+'_data_1979'+str(edyr-1)+'.nc'
    else:
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+AnomsStr+ThisRean+'_data_1979'+str(edyr)+'.nc'
    NewERAStr   = OutputVar+'_'+OutputGrid+'_'+OutputTime+'_'+AnomsStr+ThisRean+'_data_1979'+updateyyyy+'.nc'

elif (OutputVar in ['tw','q','rh','e','dpd']): # these require T, Td and SLP
    InputERA = ThisRean+'_'+ReadInGrid+'_'+ReadInTime+'_'
    if (ThisProg == 'Update'):
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+ThisRean+'_data_1979'+str(edyr-1)+'.nc'
    else:
        OldERAStr   = OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'+AnomsStr+ThisRean+'_data_1979'+str(edyr)+'.nc'
    NewERAStr   = OutputVar+'2m_'+OutputGrid+'_'+OutputTime+'_'+AnomsStr+ThisRean+'_data_1979'+updateyyyy+'.nc'
	
    # Might have some other options

# Set up variables
mdi = -1e30

# Required variable names for reading in from ERA-Interim
LatInfo = ['latitude'] 
LonInfo = ['longitude'] 

# Dictionary for looking up variable names for netCDF read in of variables
NameDict = dict([('q','q2m'),
                 ('rh','rh2m'),
	         ('e','e2m'),
	         ('tw','tw2m'),
	         ('t','t2m'),
	         ('td','td2m'),
	         ('dpd','dpd2m'),
	         ('slp','msl'),
	         ('sp','sp'),
	         ('uv',['u10','v10']), # this one might not work
	         ('ws','si10'),
		 ('sst','sst')])

# Dictionary for looking up variable standard (not actually always standard!!!) names for netCDF output of variables
StandardNameDict = dict([('q','specific_humidity'),
             ('rh','relative_humidity'),
	     ('e','vapour_pressure'),
	     ('tw','wetbulb_temperature'),
	     ('t','drybulb_temperature'),
	     ('td','dewpoint_temperature'),
	     ('dpd','dewpoint depression'),
	     ('slp','mean_sea_level_pressure'),
	     ('sp','surface_pressure'),
	     ('uv',['10 metre U wind component','10 metre V wind component']), # this one might not work
	     ('ws','10 metre windspeed'),
	     ('sst','sea_surface_temperature')])

# Dictionary for looking up variable long names for netCDF output of variables
LongNameDict = dict([('q','specific_humidity'),
             ('rh','2m relative humidity from 1by1 6hrly T and Td '+ThisRean),
	     ('e','2m vapour_pressure from 1by1 6hrly T and Td '+ThisRean),
	     ('tw','2m wetbulb_temperature from 1by1 6hrly T and Td '+ThisRean),
	     ('t','2m drybulb_temperature from 1by1 6hrly T '+ThisRean),
	     ('td','2m dewpoint_temperature from 1by1 6hrly Td '+ThisRean),
	     ('dpd','2m dewpoint depression from 1by1 6hrly T and Td '+ThisRean),
	     ('slp','2m mean_sea_level_pressure from 1by1 6hrly msl '+ThisRean),
	     ('sp','2m surface_pressure from 1by1 6hrly sp '+ThisRean),
	     ('uv',['10 metre U wind component from 1by1 6hrly '+ThisRean,'10 metre V wind component from 1by1 6hrly'+ThisRean]), # this one might not work
	     ('ws','10 metre windspeed from 1by1 6hrly'+ThisRean),
	     ('sst','sea surface temperature from 1by1 6hrly'+ThisRean)])

# Dictionary for looking up unit of variables
UnitDict = dict([('q','g/kg'),
             ('rh','%rh'),
	     ('e','hPa'),
	     ('tw','deg C'),
	     ('t','deg C'),
	     ('td','deg C'),
	     ('dpd','deg C'),
	     ('slp','hPa'),
	     ('sp','hPa'),
	     ('uv','m/s'),
	     ('ws','m/s'),
	     ('sst','deg C')])
	    

nyrs = (edyr+1)-styr
nmons = nyrs*12
npts = nyrs*73
#ndays = 
#n6hrs = 
#n1hrs = 

# set up nlons and nlats depending on what we are reading in and out
if (ReadInGrid == '1by1'):
    nlonsIn = 360
    nlatsIn= 181 # ERA style to have grids over the poles rather than up to the poles
elif (ReadInGrid == '5by5'):
    nlonsIn = 72 # assuming this is correct
    nlatsIn = 36 # assuming this is correct

if (OutputGrid == '1by1'):
    nlonsOut = 360
    nlatsOut = 180 # ERA style to have grids over the poles rather than up to the poles but this will be changed here with Build or Regrid
elif (OutputGrid == '5by5'):
    nlonsOut = 72 # assuming this is correct
    nlatsOut = 36 # assuming this is correct


## Array for monthly mean data for q, RH, e, T, Tw, Td, DPD one at a time though
##FullMonthArray = np.empty((nmons,nlats,nlons,7),dtype = float)
#FullMonthArray = np.empty((nmons,nlats,nlons),dtype = float)
#FullMonthArray.fill(mdi)

#************************************************************
# SUBROUTINES
#************************************************************
# GetHumidity
def GetHumidity(TheTDat,TheTdDat,TheSPDat,TheVar):
    ''' Calculates the desired humidity variable if the code is set up to output humidity '''
    ''' REQUIRES: '''
    ''' CalcHums.py file to be in the same directory as this file '''
    
    if (TheVar == 't'):
        TheHumDat = TheTDat

    elif (TheVar == 'td'):
        TheHumDat = TheTdDat

    elif (TheVar == 'q'):
        TheHumDat = CalcHums.sh(TheTdDat,TheTDat,TheSPDat,roundit=False)
    
    elif (TheVar == 'e'):
        TheHumDat = CalcHums.vap(TheTdDat,TheTDat,TheSPDat,roundit=False)

    elif (TheVar == 'rh'):
        TheHumDat = CalcHums.rh(TheTdDat,TheTDat,TheSPDat,roundit=False)

    elif (TheVar == 'tw'):
        TheHumDat = CalcHums.wb(TheTdDat,TheTDat,TheSPDat,roundit=False)

    elif (TheVar == 'dpd'):
        TheHumDat = CalcHums.dpd(TheTdDat,TheTDat,roundit=False)
    
    return TheHumDat

#************************************************************
# RegridField
def RegridField(TheOutputGrid,TheOldData):
    '''
    This function does a simple regridding of data by averaging over the larger gridboxes
    NO COSINE WEIGHTING FOR LATITUDE!!!!
    
    NOTE: 
    FOR OutputGrid = 5by5 THIS AUTOMATICALLY FLIPS LATITUDE AND ROLLS LONGITUDE TO BE -87.5 to 87.5 and -177,5 to 177.5
    FOR OutputGrid = 1by1 THIS JUST REGRIDS LATITUDE FROM 181 boxes 90 to -90 TO 180 boxes 89.5 to -89.5 and rolls longitude to -179.5 to 179.5
    
    Assumes input grid is always 1by1
    
    INPUTS:
    TheOutputGrid - string of 1by1 or 5by5
    TheOldData[:,:,:] - time, lat, long numpy array of complete field in original grid resolution
    OUTPUTS:
    TheNewData[:,:,:] - time, lat, long numpy array of complete field in new grid resolution
    
    I'm hoping that things set above are seen by the function rather than being passed explicitly
    
    '''
    
    # Set up the desired output array
    TheNewData = np.empty((len(TheOldData[:,0,0]),nlatsOut,nlonsOut),dtype = float)
    TheNewData.fill(mdi)
    
    if (TheOutputGrid == '1by1'):

        # Then we know we're reading in original ERA-Interim or ERA5 data which has 181 lats
	# regrid to 0.5 by 0.5 degree gridboxes and then reaverage over 89.5 to -89.5 lats 
        # shift lons back to -179.5 to 179.5 from 0 to 359
        # regrid to 5by5
	
        # First sort out the latitudes
        for ln in range(nlonsIn):
            for tt in range(len(TheNewData[:,0,0])):
                subarr = np.repeat(TheOldData[tt,:,ln],2) 
        	# this creates 362 grid boxes where each is repeated: [0a, 0b, 1a, 1b ...180a, 180b]
                subarr = subarr[1:361]
                # This removes the superfluous 90-90.5 and -90 to -90.5 boxes
                subarr = np.reshape(subarr,(180,2))
                # This now reshapes to 180 rows and 2 columns so that we can average the gridboxes across the columns
                TheNewData[tt,:,ln] = np.mean(subarr,axis = 1) # hopefully this should work!
                #pdb.set_trace()
        # Then sort out the longitudes
        for tt in range(len(TheNewData[:,0,0])):
            TheNewData[tt,:,:] = np.roll(TheNewData[tt,:,:],180,axis = 1)
		
    if (TheOutputGrid == '5by5'):

        # Then we know we're reading in my converted ERA-Interim / ERA5 data which has 180 lats and already has lons rolled 180 degrees.
	
	# flip lats to go south to north
        # regrid to 5by5

        # Regrid to 5by5 by simple averaging
        # Data input here should already be 89.5 to -89.5 lat and -179.5 to 179.5 long!!!
        StLt = 0
        EdLt = 0
        # Loop through the OutputGrid (5by5) lats and lons
        for ltt in range(nlatsOut):
            
	    # create pointers to the five lats to average over
            StLt = np.copy(EdLt)
            EdLt = EdLt + 5

            StLn = 0
            EdLn = 0
	
            for lnn in range(nlonsOut):
            
	        # create pointers to the five lons to average over
                StLn = np.copy(EdLn)
                EdLn = EdLn + 5
                #print(ltt,lnn,StLt,EdLt,StLn,EdLn)
            	        
		# Loop over each time point
                for mm in range(len(TheNewData[:,0,0])):
	        
                    # Create a subarr first so that we can deal with missing data
                    subarr = TheOldData[mm,StLt:EdLt,StLn:EdLn]    		    
                    gots = np.where(subarr > mdi)
                    if (len(gots[0]) > 0):
		    
		        # FILL THE LATITUDES BACKWARDS SO THAT THIS REVERSES THEM!!!
                        TheNewData[mm,35-ltt,lnn] = np.mean(subarr[gots])    		
			
    #pdb.set_trace()			    
		        
    return TheNewData
    
#************************************************************
# BuildField
def BuildField(TheOutputVar, TheInputTime, TheOutputTime, InFileStr, TheStYr, TheEdYr):
    ''' function for building complete reanalyses files over the period specified 
        this can be very computationally expensive so do it by year
	This requires initial reanalysis data to be read in in chunks of 1 year
	I may change this to month later and will slice out 1 month at a time anyway
	For derived variables this will read in the source vars and compute
	
	NOTE: THIS AUTOMATICALLY REGRIDS LATITUDE TO BE 180 RATHER THAN 181 BOXES AND ROLLS LONGITUDE TO -179.5 to 179.5
	
        INPUTS:
	TheOutputVar - string lower case character of q, rh, t, td, dpd, tw, e, msl, sp, ws
	TheInputTime - string of 1hr or 6hr
	TheOutputTime - string of monthly or pentad
	#OutputGrid - string of 1by1 or 5by5 (WHICH SHOULD BE SAME AS INPUT GRID) - ASSUME THIS IS ALWAYS 1by1 FOR NOW
	InFileStr - string of dir+file string to read in 
	TheStYr = integer start year of data - assume Jan 1st (0101) start 
	TheEdYr = integer end year of data - assume Dec 31st (1231) end	
	OUTPUTS:
	TheNewData[:,:,:] - time, lat, long numpy array of complete field in new time resolution
	
	'''
    # Set up the desired output array 
    if (TheOutputTime == 'monthly'):
        TheNewData = np.empty((nmons,nlatsOut,nlonsOut),dtype = float)
    elif (TheOutputTime == 'pentad'):
        TheNewData = np.empty((npts,nlatsOut,nlonsOut),dtype = float)
    TheNewData.fill(mdi)

    # The input grids are different to the output grids (181 lat boxes rather than 180) so we need a TmpNewData first
    TmpNewData = np.empty((len(TheNewData[:,0,0]),nlatsIn,nlonsIn),dtype = float)
    TmpNewData.fill(mdi)

    nyrs = (TheEdYr - TheStYr) + 1

    # Begin the time counter for this dec - may need to do hourly in 5 year or 1 year chunks
    # 0 to ~87600 + leap days for 1 hourly data (24*365*10)
    # 0 to ~14600 + leap days for 6 hourly data (4*365*10)
    # 0 to 120 for monthly data
    HrStPoint = 0 # set as HrEdPoint which is actual ed point +1
    HrEdPoint = 0 # set as HrStPoint + MonthHours or Month(must be +1 to work in Python!!!)
    	    
    # Loop through the years 
    for y in range(nyrs):
     

    	# Get actual year we're working on
        yr = y + StYr
        print('Working Year: ',yr)

        # First work out the time pointers for the year we're working with
        if (TheOutputTime == 'monthly'):	
            mnarr    = [31,29,31,30,31,30,31,31,30,31,30,31]
            nbits = 12
        elif (TheOutputTime == 'pentad'):
            mnarr    = list(np.repeat(5,73))
            nbits = 73
        # Is it a leap year?
        if (TestLeap.TestLeap(yr) == 0.0):
            if (TheOutputTime == 'monthly'):	
                mnarr[1] = 29
            elif (TheOutputTime == 'pentad'):
                mnarr[11] = 6

        print('TestLeap (m, pt): ',mnarr[1],mnarr[11], yr)
    	
        # Loop through each month or pentad depending on thing
        for m in range(nbits):
		
    	   ## string for file name
           #mm = '%02i' % (m+1)

	   # Month pointer
           MonthPointer = (y * nbits)+m
           print('Month/Pentad Pointer: ',m)

    	   # Set the time counter for this dec in either 1hr or 6hrs
    	   # 0 to ~14600 + leap days
           HrStPoint = np.copy(HrEdPoint)  # set as HrEdPoint which is actual end point +1
           if (ReadInTime == '1hr'):
               HrEdPoint = HrStPoint + (mnarr[m]*24) # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
           elif (ReadInTime == '6hr'):
               HrEdPoint = HrStPoint + (mnarr[m]*4) # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
	   
           print('Hr Pointies for this month: ',HrStPoint,HrEdPoint)

    	   # Open and read in the reanalysis files for the month
	   # Sort out time pointers to pull out month
	   # This assumes we're always reading in 1by1!!!!
           SliceInfo = dict([('TimeSlice',[HrStPoint,HrEdPoint]),
           		 ('LatSlice',[0,181]),
           		 ('LonSlice',[0,360])]) 
	   		 
           # Are we working on a direct variable or do we need to read in lots of variables and convert (e.g., humidity)
	   # For humidity variables
           if (TheOutputVar in ['q','rh','e','tw','dpd']):

	       # DOES automatically unpack the scale and offset 
	       # However, SP is Pa and T and Td are Kelvin
	       # This kills memory so need to be tidy
               ReadInfo = ['t2m']
               FileName = InFileStr+'t2m_'+str(yr)+'0101'+str(yr)+'1231.nc'
               T_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	       # Unpack t
               T_Data  = T_Data-273.15

	       # DOES automatically unpack the scale and offset 
	       # However, SP is Pa and T and Td are Kelvin
	       # This kills memory so need to be tidy
               ReadInfo = ['d2m']
               FileName = InFileStr+'td2m_'+str(yr)+'0101'+str(yr)+'1231.nc'
               Td_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	       # Unpack td
               Td_Data = Td_Data-273.15

	       # DOES automatically unpack the scale and offset 
	       # However, SP is Pa and T and Td are Kelvin
	       # This kills memory so need to be tidy
               ReadInfo = ['sp']
               FileName = InFileStr+'sp_'+str(yr)+'0101'+str(yr)+'1231.nc'
               SP_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	       # Unpack sp
               SP_Data = SP_Data/100.
           
               # Convert to desired humidity variable
               TmpData = GetHumidity(T_Data,Td_Data,SP_Data,TheOutputVar)

               # Empty the SP_Data array
               SP_Data = 0
               T_Data = 0
               Td_Data = 0

           else:

	       # DOES automatically unpack the scale and offset 
	       # However, SP is Pa and T and Td are Kelvin
	       # This kills memory so need to be tidy
               ReadInfo = [NameDict[TheOutputVar]] # the variable name to read in
               #pdb.set_trace()
               FileName = InFileStr+str(yr)+'0101'+str(yr)+'1231.nc'
               TmpData,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)
               #pdb.set_trace()

	       # Is there an unpack thing like for T - -273.15?
               if (TheOutputVar in ['t','td','sst']): # t
                   if (TheOutputVar == 'sst'): # there are missing values over land.
                       TmpData[np.where(TmpData < 270.03)] = mdi # ERA Mdi is actually -32767 but in ncview it is 270.024
                   TmpData[np.where(TmpData > mdi)] = TmpData[np.where(TmpData > mdi)]-273.15
               elif (TheOutputVar in ['slp','sp']): # pressure are Pa so need to be converted to hPa
                   TmpData = TmpData/100.
    
	   # Create monthly or pentad means
           for ltt in range(nlatsIn):
               for lnn in range(nlonsIn):

                   TmpNewData[MonthPointer,ltt,lnn] = np.mean(TmpData[:,ltt,lnn])
	   	   
	   # Empty the data arrays
           TmpData = 0 
	   
    # Now regrid to 180 latitude boxes 89.5 to -89.5 and longitude from -179.5 to 179.5	
    TheNewData = RegridField('1by1',TmpNewData)  

    return TheNewData
        
#************************************************************
# CreateAnoms
def CreateAnoms(TheInputGrid,TheOutputTime,TheClimSt,TheClimEd,TheStYr,TheEdYr,TheInData):
    '''
    This function takes any grid and any var, computes climatologies/stdevs over given period and then anomalies
    It also outputs land only and ocean only anomalies dependning on the grid
    if (TheInputGrid == '5by5'):
    LandMask = workingdir+'/OTHERDATA/HadCRUT.4.3.0.0.land_fraction.nc' # 0 = 100% sea, 1 = 100% land - no islands!, latitude, longitude, land_area_fraction, -87.5 to 87.5, -177.5 to 177.5
    elif (TheInputGrid == '1by1'):
    LandMask = workingdir+'/OTHERDATA/lsmask.nc' # 1 = sea, 0 = land - no islands! lat, lon, mask 89.5 to -89.5Lat, 0.5 to 359.5 long 
    
    INPUTS:
    TheInputGrid - string of 1by1 or 5by5 to determine the land mask to use 
    TheOutputTime - string of monthly or pentad
    TheClimSt - interger start year of climatology Always Jan start
    TheClimEd - integer end  year of climatology Always Dec end
    TheStYr - integer start year of data to find climatology
    TheEdYr - integer end year of data to find climatology
    TheInData[:,:,:] - time, lat, lon array of actual values
    
    OUTPUTS:
    AllAnomsArr[:,:,:] - time, lat, lon array of anomalies
    LandAnomsArr[:,:,:] - time, lat, lon array of land anomalies
    OceanAnomsArr[:,:,:] - time, lat, lon array of ocean anomalies
    ClimsArr[:,:,:] - time, lat, lon array of climatologies
    StDevsArr[:,:,:] - time, lat, lon array of stdeviations
    
    '''

    # Set up for time
    if (TheOutputTime == 'monthly'):
        nclims = 12
    elif (TheOutputTime == 'pentad'):
        nclims = 73	
    nyrs = (TheEdYr - TheStYr) + 1	
        
    # Get land/sea mask and format accordingly
    if (TheInputGrid == '1by1'):
        MaskData,Lats,Longs = GetGrid4(LandMask,['mask'],['lat'],['lon']) 
	# Check shape and force to be 2d
        if (len(np.shape(MaskData)) == 3):
            MaskData = np.reshape(MaskData,(180,360))
	# roll the longitudes
        MaskData = np.roll(MaskData[:,:],180,axis = 1)
	# swap the land/sea so that land = 1
        land = np.where(MaskData == 0)
        MaskData[np.where(MaskData == 1)] = 0
        MaskData[land] = 1
    elif (TheInputGrid == '5by5'):
        MaskData,Lats,Longs = GetGrid4(LandMask,['land_area_fraction'],LatInfo,LonInfo) 
        if (len(np.shape(MaskData)) == 3):
            MaskData = np.reshape(MaskData,(36,72))
    
    # first create empty arrays
    AllAnomsArr = np.empty_like(TheInData)
    AllAnomsArr.fill(mdi)
    LandAnomsArr = np.copy(AllAnomsArr)
    OceanAnomsArr = np.copy(AllAnomsArr)
    ClimsArr = np.copy(AllAnomsArr[0:nclims,:,:])
    StDevsArr = np.copy(AllAnomsArr[0:nclims,:,:])
        
    # loop through gridboxes
    for lt in range(len(TheInData[0,:,0])):
        for ln in range(len(TheInData[0,0,:])):
	
	    # pull out gridbox and reform to years by nclims (months or pentads)
            SingleSeries = np.reshape(TheInData[:,lt,ln],(nyrs,nclims)) # nyrs rows, nclims columns
            
	    # create an empty array to fill with anomalies
            NewSingleSeries = np.empty_like(SingleSeries)
            NewSingleSeries.fill(mdi)
	    
	    # loop through clims 1 to 12 or 73
            for m in range(nclims):
	    
	        # create, save and subtract climatological mean 
		# THERE ARE NO MISSING DATA IN ERA INTERIM but sst is missing over land
		# test first time value only
                if (SingleSeries[0,m] > mdi):
                    ClimsArr[m,lt,ln] = np.mean(SingleSeries[(TheClimSt-TheStYr):((TheClimEd-TheStYr)+1),m])
                    StDevsArr[m,lt,ln] = np.std(SingleSeries[(TheClimSt-TheStYr):((TheClimEd-TheStYr)+1),m])
                    NewSingleSeries[:,m] = SingleSeries[:,m] - ClimsArr[m,lt,ln]
	    
	    # fill new arrays
            AllAnomsArr[:,lt,ln] = np.reshape(NewSingleSeries,nyrs*nclims)
	    
	    # is there any land?
            if (MaskData[lt,ln] > 0):
                LandAnomsArr[:,lt,ln] = np.reshape(NewSingleSeries,nyrs*nclims)
		
	    # is there any sea?
            if (MaskData[lt,ln] < 1):
                OceanAnomsArr[:,lt,ln] = np.reshape(NewSingleSeries,nyrs*nclims)
	    
    return AllAnomsArr, LandAnomsArr, OceanAnomsArr, ClimsArr, StDevsArr
    
#************************************************************
# WriteNetCDF
def WriteNetCDF(Filename,TheOutputTime, TheOutputGrid, TheOutputVar, TheFullArray, TheFullArrayAnoms, TheLandArrayAnoms, TheOceanArrayAnoms, TheClimsArray, TheStDevsArray,
            TheStYr, TheEdYr, TheClimStart, TheClimEnd, TheName, TheStandardName, TheLongName, TheUnit):
    '''
    This function writes out a NetCDF 4 file
    
    NOTE: 
    All 1by1 outputs will have lats 89.5 to -89.5 and lons -179.5 to 179.5
    All 5by5 outputs will have lats -87.5 to 87.5 and lons -177.5 to 177.5
    
    INPUTS:
    FileOut - string file name
    TheOutputTime - string monthly or pentad
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
    nTims = len(TimPoints)
    		    
    # Sort out Lats, Lons and LatBounds and LonBounds
    if (TheOutputGrid == '1by1'):
        LatList = np.flip(np.arange(180)-89.5)
        LonList = np.arange(360)-179.5

        LatBounds = np.empty((len(LatList),2),dtype='float')
        LonBounds = np.empty((len(LonList),2),dtype='float')

        LatBounds[:,0] = LatList + ((LatList[0]-LatList[1])/2.)
        LatBounds[:,1] = LatList - ((LatList[0]-Latitudes[1])/2.)

        LonBounds[:,0] = LonList - ((LonList[1]-LonList[0])/2.)
        LonBounds[:,1] = LonList + ((LonList[1]-LonList[0])/2.)    

    elif (TheOutputGrid == '5by5'):
        LatList = (np.arange(36)*5)-87.5
        LonList = (np.arange(72)*5)-177.5

    
        LatBounds = np.empty((len(LatList),2),dtype='float')
        LonBounds = np.empty((len(LonList),2),dtype='float')

        LatBounds[:,0] = LatList - ((LatList[1]-LatList[0])/2.)
        LatBounds[:,1] = LatList + ((LatList[1]-Latitudes[0])/2.)

        LonBounds[:,0] = LonList - ((LonList[1]-LonList[0])/2.)
        LonBounds[:,1] = LonList + ((LonList[1]-LonList[0])/2.)    

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

    # Go through each variable and set up the variable attributes
    # I've added zlib=True so that the file is in compressed form
    # I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
    MyVarD = ncfw.createVariable(TheName,'f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarD.standard_name = TheStandardName
    MyVarD.long_name = TheLongName
    MyVarD.units = TheUnit
    MyVarD.valid_min = np.min(TheFullArray)
    MyVarD.valid_max = np.max(TheFullArray)
    MyVarD.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarD[:,:,:] = TheFullArray[:,:,:]  	

    # If there are climatologies etc to be written then also set them up
    if (len(np.shape(TheClimsArray)) > 1):
    
        MyVarA = ncfw.createVariable(TheName+'_anoms','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        MyVarA.standard_name = TheStandardName+'_anomalies'
        MyVarA.long_name = TheLongName+' anomalies from 1981-2010'
        MyVarA.units = TheUnit
        MyVarA.valid_min = np.min(TheFullArrayAnoms)
        MyVarA.valid_max = np.max(TheFullArrayAnoms)
        MyVarA.missing_value = mdi
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarA[:,:,:] = TheFullArrayAnoms[:,:,:]  	
    
        MyVarAL = ncfw.createVariable(TheName+'_anoms_land','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        MyVarAL.standard_name = TheStandardName+'_anomalies'
        MyVarAL.long_name = TheLongName+' anomalies from 1981-2010'
        MyVarAL.units = TheUnit
        MyVarAL.valid_min = np.min(TheLandArrayAnoms)
        MyVarAL.valid_max = np.max(TheLandArrayAnoms)
        MyVarAL.missing_value = mdi
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarAL[:,:,:] = TheLandArrayAnoms[:,:,:]  	

        MyVarAO = ncfw.createVariable(TheName+'_anoms_ocean','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        MyVarAO.standard_name = TheStandardName+'_anomalies'
        MyVarAO.long_name = TheLongName+' anomalies from 1981-2010'
        MyVarAO.units = TheUnit
        MyVarAO.valid_min = np.min(TheOceanArrayAnoms)
        MyVarAO.valid_max = np.max(TheOceanArrayAnoms)
        MyVarAO.missing_value = mdi
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarAO[:,:,:] = TheOceanArrayAnoms[:,:,:]  	

        if (TheOutputTime == 'monthly'):
            MyVarC = ncfw.createVariable(TheName+'_clims','f4',('month_time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'pentad'):
            MyVarC = ncfw.createVariable(TheName+'_clims','f4',('pentad_time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        MyVarC.standard_name = TheStandardName+'_climatologies'
        MyVarC.long_name = TheLongName+' climatology over 1981-2010'
        MyVarC.units = TheUnit
        MyVarC.valid_min = np.min(TheClimsArray)
        MyVarC.valid_max = np.max(TheClimsArray)
        MyVarC.missing_value = mdi
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarC[:,:,:] = TheClimsArray[:,:,:]  	

        if (TheOutputTime == 'monthly'):
            MyVarS = ncfw.createVariable(TheName+'_stdevs','f4',('month_time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        elif (TheOutputTime == 'pentad'):
            MyVarS = ncfw.createVariable(TheName+'_stdevs','f4',('pentad_time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
        MyVarS.standard_name = TheStandardName+'_climatological_standard_deviations'
        MyVarS.long_name = TheLongName+' climatological standard deviation over 1981-2010'
        MyVarS.units = TheUnit
        MyVarS.valid_min = np.min(TheStDevsArray)
        MyVarS.valid_max = np.max(TheStDevsArray)
        MyVarS.missing_value = mdi
        # Provide the data to the variable - depending on howmany dimensions there are
        MyVarS[:,:,:] = TheStDevsArray[:,:,:]  	

    ncfw.close()

    return
    
#************************************************************
# MAIN
#************************************************************
# What are we working on?
print('Working variable: ',OutputVar)
print('Input Time and Grid: ',ReadInTime,ReadInGrid)
print('Output Time and Grid: ',OutputTime,OutputGrid)
print('Type of run: ',ThisProg, styr, edyr, MakeAnoms)
print('Reanalysis: ',ThisRean)

# For ThisProg = Convert or Update read in monthly or pentad 1by1 (to present or previous year)
if (ThisProg != 'Build'):

    ReadInfo = [OutputVar+'2m']
    FileName = workingdir+'/OTHERDATA/'+OldERAStr
    TheData,Latitudes,Longitudes = GetGrid4(FileName,ReadInfo,LatInfo,LonInfo)

    # For Update we also need to read in most recent year (or months) of data and convert to desired variable (BuildField)
    if (ThisProg == 'Update'):

        print('Creating Update')

        # Set up the desired output array 
        if (OutputTime == 'monthly'):
            FullArray = np.empty((nmons,nlatsOut,nlonsOut),dtype = float)
        elif (OutputTime == 'pentad'):
            FullArray = np.empty((npts,nlatsOut,nlonsOut),dtype = float)
        FullArray.fill(mdi)

        # Build the most recent year
        if (OutputTime == 'monthly'):
            RecentField = FullArray[0:nmons-12,:,:]
        elif (OutputTime == 'pentad'):
            RecentField = FullArray[0:npts-73,:,:]
        RecentField = BuildField(OutputVar, ReadInTime, OutputTime, workingdir+'/OTHERDATA/'+InputERA, edyr, edyr, RecentField)
	
	# Fill the full array with data
        if (OutputTime == 'monthly'):
            FullArray[0:nmons-12,:,:] = TheData
            FullArray[nmons-12:nmons,:,:] = RecentField
        elif (OutputTime == 'pentad'):
            FullArray[0:npts-73,:,:] = TheData
            FullArray[npts-73:nmons,:,:] = RecentField

        # Do we need to create anomalies?
        # Just in case we don't, create blank arrays for write out
        FullArrayAnoms = 0 
        LandArrayAnoms = 0 
        OceanArrayAnoms = 0 
        ClimsArray = 0
        StDevArray = 0
        if (MakeAnoms == 1):

            print('Creating anomalies')
    
            FullArrayAnoms, LandArrayAnoms, OceanArrayAnoms, ClimsArray, StDevArray = CreateAnoms(ReadInGrid,OutputTime,ClimStart,ClimEnd,styr,edyr,FullArray)

		
    elif (ThisProg == 'Regrid'):

        print('Creating Regrid')

        # Do we need to create anomalies?
        # Just in case we don't, create blank arrays for write out
        FullArrayAnoms = 0 
        LandArrayAnoms = 0 
        OceanArrayAnoms = 0 
        ClimsArray = 0
        StDevArray = 0
        if (MakeAnoms == 1):

            print('Creating anomalies')
    
            TheDataAnoms, LandDataAnoms, OceanDataAnoms, ClimsDataArray, StDevDataArray = CreateAnoms(ReadInGrid,OutputTime,ClimStart,ClimEnd,styr,edyr,TheData)
            #pdb.set_trace()

        # Regrid the fields to desired resolution
        FullArray = RegridField(OutputGrid,TheData)
        if (MakeAnoms == 1):
            FullArrayAnoms = RegridField(OutputGrid,TheDataAnoms)
            LandArrayAnoms = RegridField(OutputGrid,LandDataAnoms)
            OceanArrayAnoms = RegridField(OutputGrid,OceanDataAnoms)
            ClimsArray = RegridField(OutputGrid,ClimsDataArray)
            StDevArray = RegridField(OutputGrid,StDevDataArray)

	
# For ThisProg = Build then loop through the decs for Build only
elif (ThisProg == 'Build'):

    print('Creating Build')

    FullArray = BuildField(OutputVar, ReadInTime, OutputTime, workingdir+'/OTHERDATA/'+InputERA, styr, edyr)

    # Do we need to create anomalies?
    # Just in case we don't, create blank arrays for write out
    FullArrayAnoms = 0 
    LandArrayAnoms = 0 
    OceanArrayAnoms = 0 
    ClimsArray = 0
    StDevArray = 0
    if (MakeAnoms == 1):

        print('Creating anomalies')
    
        FullArrayAnoms, LandArrayAnoms, OceanArrayAnoms, ClimsArray, StDevArray = CreateAnoms(ReadInGrid,OutputTime,ClimStart,ClimEnd,styr,edyr,FullArray)


# I've moved the anomaly bit to before any regridding because its better to do this (and the land / sea masking 
# at the highest resolution stage

## Do we need to create anomalies?
## Just in case we don't, create blank arrays for write out
#FullArrayAnoms = 0 
#LandArrayAnoms = 0 
#OceanArrayAnoms = 0 
#ClimsArray = 0
#StDevArray = 0
#if (MakeAnoms == 1):
#
#    print('Creating anomalies')
#    
#    FullArrayAnoms, LandArrayAnoms, OceanArrayAnoms, ClimsArray, StDevArray = CreateAnoms(OutputGrid,OutputTime,ClimStart,ClimEnd,styr,edyr,FullArray)
        
# Write out
print('Writing out interim monthly array: ',OutputVar)
# Now write out netcdf of field    
WriteNetCDF(workingdir+'/OTHERDATA/'+NewERAStr,OutputTime, OutputGrid, OutputVar, FullArray, FullArrayAnoms, LandArrayAnoms, OceanArrayAnoms, ClimsArray, StDevArray,
            styr, edyr, ClimStart, ClimEnd, NameDict[OutputVar], StandardNameDict[OutputVar], LongNameDict[OutputVar], UnitDict[OutputVar]) 

print('And we are done!')
