#!/usr/local/sci/bin/python
# PYTHON3.6.1
# 
# Author: Kate Willett
# Created: 18 Jul 2018
# Last update: 19 Jul 2018
# Location: /data/local/hadkw/HADCRUH2/UPDATE2017/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/HadISDH_Build					
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
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


# I DONT THINK THE BELOW APPLIES ANYMORE!!!
# This can be done using as individual months  and their file
# names changed which is tedious.
# Go to: http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/
# Log in with ECMWF log in details
# Download month by month - all hours, 0 time step, 2mT, 2mTdewpoint, Surface Pressure
# Change to 1x1 degree (select on second page!!!) before downloading as netCDF
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
# USE python2.7
# python2.7 MakeERAMonthlies.py
#
# REQUIRES
# TestLeap.py
# CalcHums.py
# ReadNetCDF.py
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

### START OF EDITABLES ###############################

# Set up initial run choices
# Start and end years
styr       = 1979
edyr       = 2017
#edOLD      = (edyr-styr)*12
stmon      = 1
edmon      = 12

# Decs for working on monthly data
#Decs = ['197901197912']
Decs = ['197901197912',
        '198001198912',
	'199001199912',
	'200001200912',
	'201001201712'] # edit the last dec?

## Decs for working on 6hourly data?
#Decs = ['1979010119791231',
#        '1980010119891231',
#	'1990010119991231',
#	'2000010120091231'
#	'2010010120171231'] # edit the last dec?

# Starting year of Decs
StDec = [1979,1980,1990,2000,2010] # this will be edited annually

# Endind year of Decs
EdDec = [1979,1989,1999,2009,2017] # this will be edited annually

# Are you reading in hourlies or monthlies?
ReadInTime = 'monthly' # this can be '6hr' or 'month' or maybe 'day' later

# Are you converting to monthlies? We will output hourlies anyway if they are read in
OutputTime = 'monthly' # this could be 'monthly' or '6hr'

# Are you reading in 1by1 or 5by5? We will output 1by1 anyway if they are read in.
ReadInGrid = '1by1' # this can be '1by1' or '5by5'

# Are you converting to 5by5?
OutputGrid = '5by5' # this can be '1by1' or '5by5'

# Do you want to create anomalies and if so, what climatology period? We will output absolutes anyway
CreateAnoms = 1 # 1 for create anomalies, 0 for do NOT create anomalies
ClimStart = 1981 # any year but generally 1981
ClimEnd = 2010 # any year but generally 2010
if (CreateAnoms == 1): # set the filename string for anomalies
    AnomStr = 'anoms'+str(ClimStart)+'-'+str(ClimEnd)
else:
    AnomStr = ''

# Set up file locations
ThisYearDir = 2017
updateyy  = str(ThisYearDir)[2:4]
updateyyyy  = str(ThisYearDir)
workingdir  = '/data/local/hadkw/HADCRUH2/UPDATE'+updateyyyy

LandMask = workingdir+'/OTHERDATA/HadCRUT.4.3.0.0.land_fraction.nc' # 0 = 100% sea, 1 = 100% land - no islands!

# Set up output variables - for q, e, RH, dpd, Tw we will need to read in multiple input files
OutputVar = 'sst' # this can be 't','td','q','rh','e','dpd','tw','ws','slp','sp','uv','sst'

### END OF EDITABLES ################

if (OutputVar in ['t','td']): # these are the simple ones that do not require conversion
    InputERA = 'era_interim_'+OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'
elif (OutputVar in ['ws','uv']):
    InputERA = 'era_interim_'+OutputVar+'10m_'+ReadInGrid+'_'+ReadInTime+'_'
elif (OutputVar in ['slp','sp','sst']):
    InputERA = 'era_interim_'+OutputVar+'_'+ReadInGrid+'_'+ReadInTime+'_'
else:
    if (OutputVar in ['tw','q','rh','e','dpd']): # these require T, Td and SLP
        # This is phoney and never used but needs to be defined
        #InputERA = 'era_interim_'+OutputVar+'2m_'+ReadInGrid+'_'+ReadInTime+'_'
        InputERAT = 'era_interim_t2m_'+ReadInGrid+'_'+ReadInTime+'_'
        InputERATd = 'era_interim_td2m_'+ReadInGrid+'_'+ReadInTime+'_'
        InputERAsp = 'era_interim_sp2m_'+ReadInGrid+'_'+ReadInTime+'_'
	
    # Might have some other options

#NewERAStr   = OutputVar+'2m_'+OutputTime+'_'+OutputGrid+'_ERA-Interim_data_'+AnomStr+'_'+str(styr)+str(edyr)+'.nc'

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

# set up nlons and nlats depending on what we are reading in and out
if (ReadInGrid == '1by1'):
    nlonsIn = 360
    nlatsIn= 181 # ERA style to have grids over the poles rather than up to the poles
elif (ReadInGrid == '5by5'):
    nlonsIn = 72 # assuming this is correct
    nlatsIn = 36 # assuming this is correct

if (OutputGrid == '1by1'):
    nlonsOut = 360
    nlatsOut = 181 # ERA style to have grids over the poles rather than up to the poles
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
# MakeDaysSince
def MakeDaysSince(TheStYr,TheStMon,TheEdYr,TheEdMon):
    ''' This function makes lists of times for monthly data suitable for saving to a netCDF file '''
    ''' Take counts of months since styr, stmn (assume 15th day of month) '''
    ''' Work out counts of days since styr,stmn, January - incl leap days '''
    ''' Also work out time boundaries 1st and last day of month '''
    ''' This can cope with incomplete years or individual months '''
    ''' REQUIRES: 
        from datetime import datetime
	import numpy as np '''
    
    # set up arrays for month mid points and month bounds
    DaysArray = np.empty(((TheEdYr-TheStYr)+1)*((TheEdMon-TheStMon)+1))
    BoundsArray = np.empty((((TheEdYr-TheStYr)+1)*((TheEdMon-TheStMon)+1),2))
    
    # make a date object for each time point and subtract start date
    StartDate = datetime(TheStYr,TheStMon,1,0,0,0)	# January
    TheYear = TheStYr
    TheMonth = TheStMon
    for mm in range(len(DaysArray)):
        if (TheMonth < 12):
            DaysArray[mm] = (datetime(TheYear,TheMonth+1,1,0,0,0)-datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days
            BoundsArray[mm,0] = (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1] = (datetime(TheYear,TheMonth+1,1,0,0,0)-StartDate).days
        else:
            DaysArray[mm] = (datetime(TheYear+1,1,1,0,0,0)-datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days	
            BoundsArray[mm,0] = (datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1] = (datetime(TheYear+1,1,1,0,0,0)-StartDate).days
        TheMonth=TheMonth+1
        if (TheMonth == 13):
            TheMonth = 1
            TheYear = TheYear+1
	    
    return DaysArray,BoundsArray

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
# MAIN
#************************************************************
# What are we working on?
print('Working variable: ',OutputVar)
print('Input Time and Grid: ',ReadInTime,ReadInGrid)
print('Output Time and Grid: ',OutputTime,OutputGrid)

# Set up the interim output arrays - ReadInGrid monthly - interim arrays will be saved
FullInterimArray = np.empty((nmons,nlatsIn,nlonsIn),dtype = float)
FullInterimArray.fill(mdi)

# Set up the desired output array - OutputGrid monthly interim arrays will be saved
FullOutputArray = np.empty((nmons,nlatsOut,nlonsOut),dtype = float)
FullOutputArray.fill(mdi)

# Loop through the decs
for d,dec in enumerate(Decs):
        
    # Number of years and months in dec
    ndecyrs = (EdDec[d] - StDec[d]) + 1
    ndecmons = ndecyrs*12

    print('Working Decs: ',d,dec,ndecyrs)
    #pdb.set_trace()

    # Begin the time counter for this dec
    # 0 to ~14600 + leap days for 6 hourly data (4*365*10)
    # 0 to 120 for monthly data
    HrStPoint = 0 # set as HrEdPoint which is actual ed point +1
    HrEdPoint = 0 # set as HrStPoint + MonthHours or Month(must be +1 to work in Python!!!)
    	    
    # Loop through the years in the dec
    for y in range(ndecyrs):
     
        print('Working Year: ',StDec[d]+y)

    	# Get actual year we're working on
        yr = y + StDec[d]

        # First work out the time pointers for the year we're working with
        # ONLY USED IF WORKING ON 6 HOURLY
        # Is it a leap year or not?
        if (TestLeap.TestLeap(yr) == 0.0):
            mnarr    = [31,29,31,30,31,30,31,31,30,31,30,31]
        else:
            mnarr    = [31,28,31,30,31,30,31,31,30,31,30,31]

        print('TestLeap: ',mnarr[1], yr)
    	
        # Loop through each month
        # This is a long way of doing this but it allows the code to work on both 6hourly and monthly input data
        for m in range(12):
		
    	    # string for file name
            mm = '%02i' % (m+1)

	    # Month pointer
            MonthPointer = ((yr - styr) * 12)+m
            #print('Month Pointer: ',m)

    	    # Set the time counter for this dec
    	    # 0 to ~14600 + leap days
            HrStPoint = np.copy(HrEdPoint)  # set as HrEdPoint which is actual end point +1
            if (ReadInTime == '6hourly'):
                HrEdPoint = HrStPoint + (mnarr[m]*4) # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
            else:
                HrEdPoint = HrStPoint + 1 # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
	    
            print('Hr Pointies for this month: ',HrStPoint,HrEdPoint)

            #print('Reading in Month: ',mm)

    	    # Open and read in the ERAINTERIM files for the month
	    # Sort out time pointers to pull out month
            SliceInfo = dict([('TimeSlice',[HrStPoint,HrEdPoint]),
                              ('LatSlice',[0,181]),
                              ('LonSlice',[0,360])]) 
			      
            # Are we working on a direct variable or do we need to read in lots of variables and convert (e.g., humidity)
	    # For humidity variables
            if (OutputVar in ['q','rh','e','tw','dpd']):

                # Read in the data from the dec
                #print('Reading in T for dec: ',dec)

	        # DOES automatically unpack the scale and offset 
	        # However, SP is Pa and T and Td are Kelvin
	        # This kills memory so need to be tidy
                ReadInfo = ['t2m']
                FileName = workingdir+'/OTHERDATA/'+InputERAT+dec+'.nc'
                T_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	        # Unpack t
                T_Data  = T_Data-273.15

                #print('Reading in Td for dec: ',dec)

	        # DOES automatically unpack the scale and offset 
	        # However, SP is Pa and T and Td are Kelvin
	        # This kills memory so need to be tidy
                ReadInfo = ['d2m']
                FileName = workingdir+'/OTHERDATA/'+InputERATd+dec+'.nc'
                Td_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	        # Unpack td
                Td_Data = Td_Data-273.15

                #print('Reading in SP for dec: ',dec)

	        # DOES automatically unpack the scale and offset 
	        # However, SP is Pa and T and Td are Kelvin
	        # This kills memory so need to be tidy
                ReadInfo = ['sp']
                FileName = workingdir+'/OTHERDATA/'+InputERAsp+dec+'.nc'
                SP_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

	        # Unpack sp
                SP_Data = SP_Data/100.
            	
            	# Convert to desired humidity variable
                New_Data = GetHumidity(T_Data,Td_Data,SP_Data,OutputVar)

            	# Empty the SP_Data array
                SP_Data = 0
                T_Data = 0
                Td_Data = 0

            else:

        	# Read in the data from the dec
                #print('Reading in OutputVar for dec: ',dec)

	        # DOES automatically unpack the scale and offset 
	        # However, SP is Pa and T and Td are Kelvin
	        # This kills memory so need to be tidy
                ReadInfo = [NameDict[OutputVar]] # the variable name to read in
                #pdb.set_trace()
                FileName = workingdir+'/OTHERDATA/'+InputERA+dec+'.nc'
                New_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)
                #pdb.set_trace()

	        # Is there an unpack thing like for T - -273.15?
                if (OutputVar in ['t','td','sst']): # t
                    if (OutputVar == 'sst'): # there are missing values over land.
                        New_Data[np.where(New_Data < 270.03)] = mdi # ERA Mdi is actually -32767 but in ncview it is 270.024
                    New_Data = New_Data-273.15
                elif (OutputVar in ['slp','sp']): # pressure are Pa so need to be converted to hPa
                    New_Data = New_Data/100.

	    # If we're creating monthly means from 6hr then average and place in the full array
            if (ReadInTime == '6hr') & (OutputTime == 'monthly'):
                for ltt in range(nlatsIn):
                    for lnn in range(nlonsIn):

                        FullInterimArray[MonthPointer,ltt,lnn] = np.mean(New_Data[:,ltt,lnn])
			
            else:
                
                FullInterimArray[MonthPointer,:,:] = New_Data[:,:,:]	        		

	    # Empty the data arrays
            New_Data = 0 

# Now we should have a complete monthly mean dataset in the native array
# save as days since 19790101

# Write out
print('Writing out interim monthly array: ',OutputVar)

# Write out 
TimPoints,TimBounds = MakeDaysSince(styr,stmon,edyr,edmon)
nTims = len(TimPoints)
    		    
# Sort out LatBounds and LonBounds
LatBounds = np.empty((len(Latitudes),2),dtype='float')
LonBounds = np.empty((len(Longitudes),2),dtype='float')

LatBounds[:,0] = Latitudes - ((Latitudes[1]-Latitudes[0])/2.)
LatBounds[:,1] = Latitudes + ((Latitudes[1]-Latitudes[0])/2.)

LonBounds[:,0] = Longitudes - ((Longitudes[1]-Longitudes[0])/2.)
LonBounds[:,1] = Longitudes + ((Longitudes[1]-Longitudes[0])/2.)    

#pdb.set_trace()

# No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
# Not sure what this does to float precision though...

# Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
OutFile = workingdir+'/OTHERDATA/'+NameDict[OutputVar]+'_'+OutputTime+'_'+ReadInGrid+'_ERA-Interim_data_'+str(styr)+str(edyr)+'.nc'
ncfw = Dataset(OutFile,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first

# Set up the dimension names and quantities
ncfw.createDimension('time',nmons)
ncfw.createDimension('latitude',nlatsIn)
ncfw.createDimension('longitude',nlonsIn)

# Go through each dimension and set up the variable and attributes for that dimension if needed
MyVarT = ncfw.createVariable('time','f4',('time',))
MyVarT.standard_name = 'time'
MyVarT.long_name = 'time'
MyVarT.units = 'days since 1979-1-1 00:00:00'
MyVarT.start_year = str(styr)
MyVarT.end_year = str(edyr)
MyVarT[:] = TimPoints

MyVarLt = ncfw.createVariable('latitude','f4',('latitude',))
MyVarLt.standard_name = 'latitude'
MyVarLt.long_name = 'gridbox centre latitude'
MyVarLt.units = 'degrees_north'
MyVarLt[:] = Latitudes

MyVarLn = ncfw.createVariable('longitude','f4',('longitude',))
MyVarLn.standard_name = 'longitude'
MyVarLn.long_name = 'gridbox centre longitude'
MyVarLn.units = 'degrees_east'
MyVarLn[:] = Longitudes

# Go through each variable and set up the variable attributes
# I've added zlib=True so that the file is in compressed form
# I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
MyVarD = ncfw.createVariable(NameDict[OutputVar],'f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
MyVarD.standard_name = StandardNameDict[OutputVar]
MyVarD.units = UnitDict[OutputVar]
MyVarD.valid_min = np.min(FullInterimArray)
MyVarD.valid_max = np.max(FullInterimArray)
MyVarD.missing_value = mdi
# Provide the data to the variable - depending on howmany dimensions there are
MyVarD[:,:,:] = FullInterimArray[:,:,:]  	
#pdb.set_trace()    	
ncfw.close()
#pdb.set_trace()

# Now if desired convert the data to 5by5 grids
if (OutputGrid == '5by5'):

    # flip lats to go south to north
    # shift lons back to -180 to 180 from 0 to 359
    # save as days since 19790101
    # regrid to 5by5

    # Regrid to 5by5 by simple averaging
    # ERAI is 181 lats 90.5 to -90.5 (GB centres 90 to -90) therefore we need to regrid to 0.5by0.5 degrees first and then average over 5by5 box
    # ERAI is 360 longs -0.5 to 359.5 (GB centres 0 to 359) therefore we need to regrid to 0.5by0.5 degrees first and then average over 5by5 box
    hlt = 0
    hln = 0
    # Loop through the OutputGrid (5by5) lats and lons
    for ltt in range(nlatsOut):
        for lnn in range(nlonsOut):
            #print(ltt,lnn,hlt,hln)
            # Loop over each month - slow!
            for mm in range(nmons):
	        # create a temporary array of 12x12 0.5by0.5 boxes to be filled with 6by6 lats and lons
		# this is really confusing - the first 12 will cover 90.5 to 84.5 N and -0.5 W to 5.5 E
		# we then take an average only over 90.0-85.0 N and 0-5E.
                subhi = np.empty((12,12))
                subhi.fill(mdi) 
                for i in range(6):
                    for j in range(6):
                        icol = i*2
                        jrow = j*2
                        if (hln == 360): 
		            # Make sure it wraps around so reset 360 to be 0
                            hln = 0
                        subhi[jrow:jrow+2,icol:icol+2] = FullInterimArray[mm,hlt,hln] # +2 because python references arrays with 0-N rather 0-(N-1)
                        hlt = hlt+1
                    hlt = hlt-6
                    hln = hln+1
                hln = hln-6
                subsubhi = np.copy(subhi[1:11,1:11])
                gots = np.where(subsubhi > mdi)
                if (len(gots[0]) > 0):
                    FullOutputArray[mm,ltt,lnn] = np.mean(subsubhi[gots])
#                if (FullOutputArray[mm,ltt,lnn] < -10.):
#                    print('Got Low')
#                    pdb.set_trace()
            hln = hln+5
        hln = 0
        hlt = hlt+5
    #pdb.set_trace()
    
    # now flip and shift so its -90 to 90, -180 to 180
    for mm in range(nmons):
        FullOutputArray[mm,:,:] = np.flipud(FullOutputArray[mm,:,:]) # reverse latitudes
	# going from 2.5 to 357.5 to -177.5 to 177.5 means rolling forward such that the current 0 element lies at 36
        FullOutputArray[mm,:,:] = np.roll(FullOutputArray[mm,:,:],np.int(nlonsOut/2),axis=1)

    # Recreate Latitudes and Longitudes lists
    Latitudes = np.arange(-87.5,92.5,5.)
    Longitudes = np.arange(-177.5,182.5,5.)    

    # Now we should have a complete monthly mean dataset in the 5by5 array
    # Write out
    print('Writing out final monthly array: ',OutputVar)

    # Write out
    # Times already sorted from above 
    		    
    # Sort out LatBounds and LonBounds
    LatBounds = np.empty((len(Latitudes),2),dtype='float')
    LonBounds = np.empty((len(Longitudes),2),dtype='float')

    LatBounds[:,0] = Latitudes - ((Latitudes[1]-Latitudes[0])/2.)
    LatBounds[:,1] = Latitudes + ((Latitudes[1]-Latitudes[0])/2.)

    LonBounds[:,0] = Longitudes - ((Longitudes[1]-Longitudes[0])/2.)
    LonBounds[:,1] = Longitudes + ((Longitudes[1]-Longitudes[0])/2.)    

    #pdb.set_trace()

    # No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
    # Not sure what this does to float precision though...

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    OutFile = workingdir+'/OTHERDATA/'+NameDict[OutputVar]+'_'+OutputTime+'_'+OutputGrid+'_ERA-Interim_data_'+str(styr)+str(edyr)+'.nc'
    ncfw=Dataset(OutFile,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first

    # Set up the dimension names and quantities
    ncfw.createDimension('time',nmons)
    ncfw.createDimension('latitude',nlatsOut)
    ncfw.createDimension('longitude',nlonsOut)

    # Go through each dimension and set up the variable and attributes for that dimension if needed
    MyVarT = ncfw.createVariable('time','f4',('time',))
    MyVarT.standard_name = 'time'
    MyVarT.long_name = 'time'
    MyVarT.units = 'days since 1979-1-1 00:00:00'
    MyVarT.start_year = str(styr)
    MyVarT.end_year = str(edyr)
    MyVarT[:] = TimPoints

    MyVarLt = ncfw.createVariable('latitude','f4',('latitude',))
    MyVarLt.standard_name = 'latitude'
    MyVarLt.long_name = 'gridbox centre latitude'
    MyVarLt.units = 'degrees_north'
    MyVarLt[:] = Latitudes

    MyVarLn = ncfw.createVariable('longitude','f4',('longitude',))
    MyVarLn.standard_name = 'longitude'
    MyVarLn.long_name = 'gridbox centre longitude'
    MyVarLn.units = 'degrees_east'
    MyVarLn[:] = Longitudes

    # Go through each variable and set up the variable attributes
    # I've added zlib=True so that the file is in compressed form
    # I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
    MyVarD = ncfw.createVariable('actuals','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarD.standard_name = StandardNameDict[OutputVar]
    MyVarD.units = UnitDict[OutputVar]
    MyVarD.valid_min = np.min(FullOutputArray)
    MyVarD.valid_max = np.max(FullOutputArray)
    MyVarD.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarD[:,:,:] = FullOutputArray[:,:,:]  	
    	
    ncfw.close()


# Now if desired calculate anomalies and land/sea masked versions
if (CreateAnoms == 1):

    # Calculate anomalies over climatology period
    # at the same time mask out land and sea in respective fields
    # first create empty arrays
    FullOutputArrayAnoms = np.empty_like(FullOutputArray)
    FullOutputArrayAnoms.fill(mdi)
    FullOutputArrayAnomsLand = np.empty_like(FullOutputArray)
    FullOutputArrayAnomsLand.fill(mdi)
    FullOutputArrayAnomsSea = np.empty_like(FullOutputArray)
    FullOutputArrayAnomsSea.fill(mdi)
    
    # now read in the land sea mask
    MaskData,Lats,Longs = GetGrid4(LandMask,['land_area_fraction'],LatInfo,LonInfo)
    
    # loop through gridboxes
    for lt in range(nlatsOut):
        for ln in range(nlonsOut):
	
	    # pull out gridbox and reform to years by months
            SingleSeries = np.reshape(FullOutputArray[:,lt,ln],(nyrs,12)) # nyrs rows, 12 columns
            # create an empty array to fill with anomalies
            NewSingleSeries = np.empty_like(SingleSeries)
            NewSingleSeries.fill(mdi)
	    
	    # loop through months 1 to 12
            for m in range(12):
	    
	        # subtract climatological mean 
		# THERE ARE NO MISSING DATA IN ERA INTERIM but sst is missing over land
		# test first time value only
                if (SingleSeries[0,m] > mdi):
                    NewSingleSeries[:,m] = SingleSeries[:,m] - np.mean(SingleSeries[(ClimStart-styr):((ClimEnd-styr)+1),m])
	    
	    # fill new arrays
            FullOutputArrayAnoms[:,lt,ln] = np.reshape(NewSingleSeries,nmons)
	    
	    # is there any land?
            if (MaskData[0,lt,ln] > 0):
                FullOutputArrayAnomsLand[:,lt,ln] = np.reshape(NewSingleSeries,nmons)
		
	    # is there any sea?
            if (MaskData[0,lt,ln] < 1):
                FullOutputArrayAnomsSea[:,lt,ln] = np.reshape(NewSingleSeries,nmons)
	    
    # Now we should have a complete monthly mean dataset in the 5by5 array as anomalies and land/sea masked versions
    # Write out - already sorted times and lats and lons
    print('Writing out final monthly anomaly array: ',OutputVar)

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    OutFile = workingdir+'/OTHERDATA/'+NameDict[OutputVar]+'_'+OutputTime+'_'+OutputGrid+'_ERA-Interim_data_'+str(styr)+str(edyr)+'_'+AnomStr+'.nc'
    ncfw=Dataset(OutFile,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first

    # Set up the dimension names and quantities
    ncfw.createDimension('time',nmons)
    ncfw.createDimension('latitude',nlatsOut)
    ncfw.createDimension('longitude',nlonsOut)

    # Go through each dimension and set up the variable and attributes for that dimension if needed
    MyVarT = ncfw.createVariable('time','f4',('time',))
    MyVarT.standard_name = 'time'
    MyVarT.long_name = 'time'
    MyVarT.units = 'days since 1979-1-1 00:00:00'
    MyVarT.start_year = str(styr)
    MyVarT.end_year = str(edyr)
    MyVarT[:] = TimPoints

    MyVarLt = ncfw.createVariable('latitude','f4',('latitude',))
    MyVarLt.standard_name = 'latitude'
    MyVarLt.long_name = 'gridbox centre latitude'
    MyVarLt.units = 'degrees_north'
    MyVarLt[:] = Latitudes

    MyVarLn = ncfw.createVariable('longitude','f4',('longitude',))
    MyVarLn.standard_name = 'longitude'
    MyVarLn.long_name = 'gridbox centre longitude'
    MyVarLn.units = 'degrees_east'
    MyVarLn[:] = Longitudes

    # Go through each variable and set up the variable attributes
    # I've added zlib=True so that the file is in compressed form
    # I've added least_significant_digit=4 because we do not need to store information beyone 4 significant figures.
    MyVarD = ncfw.createVariable('anomalies','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarD.standard_name = StandardNameDict[OutputVar]
    MyVarD.units = UnitDict[OutputVar]
    MyVarD.valid_min = np.min(FullOutputArrayAnoms)
    MyVarD.valid_max = np.max(FullOutputArrayAnoms)
    MyVarD.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarD[:,:,:] = FullOutputArrayAnoms[:,:,:]  	

    MyVarDL = ncfw.createVariable('anomalies_land','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarDL.standard_name = StandardNameDict[OutputVar]
    MyVarDL.units = UnitDict[OutputVar]
    MyVarDL.valid_min = np.min(FullOutputArrayAnomsLand)
    MyVarDL.valid_max = np.max(FullOutputArrayAnomsLand)
    MyVarDL.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarDL[:,:,:] = FullOutputArrayAnomsLand[:,:,:]  	

    MyVarDS = ncfw.createVariable('anomalies_sea','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarDS.standard_name = StandardNameDict[OutputVar]
    MyVarDS.units = UnitDict[OutputVar]
    MyVarDS.valid_min = np.min(FullOutputArrayAnomsSea)
    MyVarDS.valid_max = np.max(FullOutputArrayAnomsSea)
    MyVarDS.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarDS[:,:,:] = FullOutputArrayAnomsSea[:,:,:]  	
    	
    ncfw.close()

print('And we are done!')
