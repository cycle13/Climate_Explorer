# REDUNDANT CODE - USE ExtractMergeRegridERA...

#!/usr/local/sci/bin/python
# PYTHON2.7
# 
# Author: Kate Willett
# Created: 8 Mar 2018
# Last update: 8 Mar 2018
# Location: /data/local/hadkw/HADCRUH2/UPDATE2017/PROGS/HADISDH_BUILD/	
# GitHub: https://github.com/Kate-Willett/HadISDH_Build					
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code reads in the ERA-Interim monthly 1by1 6 hourly T, Td and Surface Pressure for the full time period
#
# It converts to humidity variables and monthly means and saves.
#
# You can then run regridERA_HadISDH_MAY2015.pro (IDL) to create ERA-Interim
# anomalies (1981-2010) and regrid to 5by5 degree so that it can be compared
# with HadISDH. This code will eventually be converted to Python too.
#
# The ERA-Interim updates have to be downloaded as individual months and their file
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
# First make sure the New and Old ERA-Interim data are in the right place.
# python2.7 MakeERAMonthlies.py
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
# Version 1 (8 Feb 2018)
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

# Set up initial run choices
# Start and end years
styr       = 1979
edyr       = 2017
#edOLD      = (edyr-styr)*12
stmon      = 1
edmon      = 12

# Set up file locations
updateyy  = str(edyr)[2:4]
updateyyyy  = str(edyr)
workingdir  = '/data/local/hadkw/HADCRUH2/UPDATE'+updateyyyy

#OldERAStr   = '2m_monthly_1by1_ERA-Interim_data_1979'+str(edyr-1)+'.nc'
NewERAStr   = '2m_monthly_1by1_ERA-Interim_data_1979'+updateyyyy+'.nc'
# ERA updates from 2017 onwards
UpdateERAStr = 'ERAINTERIM_6hr_1by1_'
# ERA updates pre2017 
UpdateERAStrT = 'ERAINTERIM_drybulbT_6hr_1by1_'
UpdateERAStrTd = 'ERAINTERIM_dewpointT_6hr_1by1_'
UpdateERAStrSP = 'ERAINTERIM_sp_6hr_1by1_'
# DateString for pre-2017
Decs = ['1979010119881231','1989010119981231','1999010120081231','2009010120161231','2017']
StDec = [1979,1989,1999,2009,2017]
EdDec = [1988,1998,2008,2016,2017]

# Set up variables
mdi = -1e30

VarArray = ['q','rh','e','t','tw','td','dpd']
StNameArr = ['specific_humidity','relative_humidity','vapour_pressure','drybulb_temperature','wetbulb_temperature','dewpoint_temperature','dewpoint_depression']
UnitArr = ['g/kg','%rh','hPa','deg C','deg C','deg C','deg C']

nyrs = (edyr+1)-styr
nmons = nyrs*12

nlons = 360
nlats = 181

## Array for monthly mean data for q, RH, e, T, Tw, Td, DPD one at a time though
##FullMonthArray = np.empty((nmons,nlats,nlons,7),dtype = float)
#FullMonthArray = np.empty((nmons,nlats,nlons),dtype = float)
#FullMonthArray.fill(mdi)

#************************************************************
# SUBROUTINES
#************************************************************
# MakeDaysSince
def MakeDaysSince(TheStYr,TheStMon,TheEdYr,TheEdMon):
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
    ''' Calculates the desired humidity variable '''
    
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
# For memory we had better do one variable at a time or its going to get BIG!
LatInfo = ['latitude'] 
LonInfo = ['longitude'] 

for v,var in enumerate(VarArray):

    # If you want to skip a variable then do it here
    if ((v == 3) or (v == 5)):
        continue
    
    print(v,var)
    
    # Array for monthly mean data for q, RH, e, T, Tw, Td, DPD one at a time though
    FullMonthArray = np.empty((nmons,nlats,nlons),dtype = float)
    FullMonthArray.fill(mdi)

    # Loop through the decs
    for d,dec in enumerate(Decs):
        
	# Number of years and months in dec
	ndecyrs = (EdDec[d] - StDec[d]) + 1
	ndecmons = ndecyrs*12
	
	# Begin the hour counter for this dec
	# 0 to ~14600 + leap days
	HrStPoint = 0 # set as HrEdPoint which is actual ed point +1
	HrEdPoint = 0 # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
		
        # Loop through the years in the dec
	for y in range(ndecyrs):
	
	    # Get actual year we're working on
	    yr = y + StDec[d]

            # First work out the time pointers for the year we're working with
            # Is it a leap year or not?
            if (TestLeap.TestLeap(yr) == 0.0):
                mnarr    = [31,29,31,30,31,30,31,31,30,31,30,31]
            else:
                mnarr    = [31,28,31,30,31,30,31,31,30,31,30,31]
    
            print('TestLeap: ',mnarr[1], yr)
	    
            # Loop through each month
	    for m in range(12):
                    
		# string for file name
                mm = '%02i' % (m+1)

                # Month pointer
                MonthPointer = ((yr - styr) * 12)+m
                print('Month Pointer: ',m)

	        # Set the hour counter for this dec
	        # 0 to ~14600 + leap days
	        HrStPoint = np.copy(HrEdPoint)  # set as HrEdPoint which is actual end point +1
	        HrEdPoint = HrStPoint + (mnarr[m]*4) # set as HrStPoint + MonthHours (must be +1 to work in Python!!!)
                print('Hr Pointies for this month: ',HrStPoint,HrEdPoint)
    
                print('Reading in Month: ',mm)

	        # Open and read in the ERAINTERIM files for the month
	        # Until 2017 there is a seperate file for t, td and sp
	        # From 2017 onwards each month is seperate but contains all variables
	        if (yr < 2017): # test the string length to check its a pre-2017 dec

                    # Sort out time pointers to pull out month
		    SliceInfo = dict([('TimeSlice',[HrStPoint,HrEdPoint]),
		                      ('LatSlice',[0,181]),
				      ('LonSlice',[0,360])]) 
	
	            # Read in the data from the dec
                    print('Reading in T for dec: ',dec)

                    # DOES automatically unpack the scale and offset 
                    # However, SP is Pa and T and Td are Kelvin
                    # This kills memory so need to be tidy
                    ReadInfo = ['t2m']
                    FileName = workingdir+'/OTHERDATA/'+UpdateERAStrT+dec+'.nc'
                    T_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

                    # Unpack t
                    T_Data  = T_Data-273.15

                    print('Reading in Td for dec: ',dec)

                    # DOES automatically unpack the scale and offset 
                    # However, SP is Pa and T and Td are Kelvin
                    # This kills memory so need to be tidy
                    ReadInfo = ['d2m']
                    FileName = workingdir+'/OTHERDATA/'+UpdateERAStrTd+dec+'.nc'
                    Td_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

                    # Unpack td
                    Td_Data = Td_Data-273.15

                    print('Reading in SP for dec: ',dec)

                    # DOES automatically unpack the scale and offset 
                    # However, SP is Pa and T and Td are Kelvin
                    # This kills memory so need to be tidy
                    ReadInfo = ['sp']
                    FileName = workingdir+'/OTHERDATA/'+UpdateERAStrSP+dec+'.nc'
                    SP_Data,Latitudes,Longitudes = GetGrid4Slice(FileName,ReadInfo,SliceInfo,LatInfo,LonInfo)

                    # Unpack sp
                    SP_Data = SP_Data/100.

	        # IF its 2017+ then we'll need to loop through each month and read in the data
	        elif (yr >= 2017):

                    # Read in the New Files for each month
                    ReadInfo = ['sp','t2m','d2m']

                    # DOES automatically unpack the scale and offset 
                    # However, SP is Pa and T and Td are Kelvin
                    # This kills memory so need to be tidy
                    FileName = workingdir+'/OTHERDATA/'+UpdateERAStr+mm+str(yr)+'.nc'

                    TmpData,Latitudes,Longitudes = GetGrid4(FileName,ReadInfo,LatInfo,LonInfo)

                    # Unpack into month of t, td, and sp
                    T_Data  = np.copy(TmpData[1]-273.15)
                    Td_Data = np.copy(TmpData[2]-273.15)
                    SP_Data = np.copy(TmpData[0]/100.)
                    # Empty the TmpData array
                    TmpData = 0

	    
	        # Convert to desired humidity variable
 	        Hum_Data = GetHumidity(T_Data,Td_Data,SP_Data,var)

	        # Empty the SP_Data array
                SP_Data = 0
                T_Data = 0
                Td_Data = 0

                # Place the monthly means in the full array
                # Can we do this on mass or by lat/lon individually?
                # I think we have to loop over the gridboxes
                for ltt in range(nlats):
                    for lnn in range(nlons):
    
                        FullMonthArray[MonthPointer,ltt,lnn] = np.mean(Hum_Data[:,ltt,lnn])

                # Empty the data arrays
                Hum_Data = 0 

    # Now we should have a complete monthly mean dataset
    # Write out
    print('Writing out: ',var)

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
    OutFile = workingdir+'/OTHERDATA/'+var+NewERAStr
    ncfw=Dataset(OutFile,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first
    
    # Set up the dimension names and quantities
    ncfw.createDimension('time',nmons)
    ncfw.createDimension('latitude',nlats)
    ncfw.createDimension('longitude',nlons)
	
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
    MyVarD = ncfw.createVariable(var+'2m','f4',('time','latitude','longitude',),fill_value = mdi,zlib=True,least_significant_digit=4)
    MyVarD.standard_name = StNameArr[v]
    MyVarD.units = UnitArr[v]
    MyVarD.missing_value = mdi
    # Provide the data to the variable - depending on howmany dimensions there are
    MyVarD[:,:,:] = FullMonthArray	    
	    
    ncfw.close()

print('And we are done!')
