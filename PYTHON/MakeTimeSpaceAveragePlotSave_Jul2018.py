#!/usr/local/sci/bin/python
# PYTHON3
# Python3
# Author: Kate Willett
# Created: 26 July 2018
# Last update: 2 Aug 2018
# Location: /data/local/hadkw/HADCRUH2/UPDATE2017/PROGS/PYTHON/	# this will probably change
# GitHub: https://github.com/Kate-Willett/Climate_Explorer/tree/master/PYTHON/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in monthly grids from HadISDH or ERA-Interim and creates a time average
# e.g., across a number of months or years. It can create a month or seasonal 
# field or average a month or seasonal field over a specified number of years.
# It can plot a map of that time average
# It can save the field of the time average to netCDF
# It can work on a subset region only or a global field. 
# If its a region if can produce a plot of all gridbox time series for that region
# It can produce a global
# or regional average time series of the monthly or time averaged values.
# It can plot a time series of the spatial average and save it to netCDF
# 
# <references to related published material, e.g. that describes data set>
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# Inbuilt:
#import matplotlib.pyplot as plt
#import numpy as np
#import numpy.ma as ma
#import sys, os, getopt
#import scipy.stats
#import struct
#import cartopy.crs as ccrs
#import cartopy.feature as cpf
#import datetime as dt
#import netCDF4 as nc4 # To solve error "is not a valid NetCDF 3 file"
#from matplotlib.dates import date2num,num2date
#from matplotlib import ticker
#import matplotlib.cm as mpl_cm
#import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)
#
# Other:
# 
# -----------------------
# DATA
# -----------------------
# directory for trendmaps:
# /data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/
# directory for land cover map:
# /data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/
# land cover file used to calculate % land represented:
# HadCRUT.4.3.0.0.land_fraction.nc
# files currently worked on:
# HadISDH.landq.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc
# HadISDH.landRH.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc
# HadISDH.landT.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc
# q2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc
# rh2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc
# t2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
#
# # set up all editables as you wish within the script 
# > module load scitools/experimental-current
# > python MakeTimeSpaceAveragePlotSave_JUL2018.py
#
# or with some/all command line arguments
# > module load scitools/experimental-current
# >python MakeTimeSpaceAveragePlotSave_Jul2018.py --sy 1973 --ey 1973 --sm 1 --em 1 --ss 1 --sc JJA --gs 0 --slt -10.0 --sln -10.0 etc
#
# -----------------------
# OUTPUT
# -----------------------
# directory for output images:
# /data/local/hadkw/HADCRUH2/UPDATE2017/IMAGES/ANALYSIS/
# Output image files:
# TimeAverage_HadISDH.landq.4.0.0.2017f_'...',png/.eps
# RegionAverage_HadISDH.landq.4.0.0.2017f_'...',png/.eps
# RegionAverage_HadISDH.landq.4.0.0.2017f_'...'multi.png/.eps
# Output netCDF files:
# TimeAverage_HadISDH.landq.4.0.0.2017f_'...',png/.eps
# RegionAverage_HadISDH.landq.4.0.0.2017f_'...',png/.eps
#
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 1 26 July 2018)
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
# MAIN PROGRAM
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
# some jiggery pokery to make this work with scitools/experimental-current
plt.switch_backend('tkagg')
plt.clf()
# https://matplotlib.org/users/customizing.html 
# I don't have anything in my .config/matplotlib/ - should there be matplotlibrc? 
import numpy as np
import numpy.ma as ma
import sys, os, getopt
import scipy.stats
import struct
import cartopy.crs as ccrs
import cartopy.feature as cpf
import datetime as dt
import netCDF4 as nc4 # To solve error "is not a valid NetCDF 3 file"
from matplotlib.dates import date2num,num2date
#from scipy.io import netcdf
from matplotlib import ticker
#import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)


def main(argv):
    ' Some info about what this code does and any things we might want to specify in the command line when running the program '
    ' python PlotAverageMap_cartopy_Jul2018_timeplotaverage.py --sy 1979 --ey 2017 ... '
    ' You can specify all, some or no variables at the command line using the '"--XX"' codes below '
    ' Any command line specifications will overwrite the listed variable set ups within the code '
    ' --sy = styr start year of desired output '
    ' --ey = edyr end year of desired output '
    ' --sm = stmn start month of desired output should be 1-12 '
    ' --em = edmn end month of desired output should be 1-12 '
    ' --ss = SeasonSwitch should be 0 or 1 # 1 = create seasonal average, 0 = create entire time period average '
    ' --sc = SeasonChoice specify season to average over e.g., '"JJA"' see SeasonDict dictionary below/add new entries '
    ' --gs = GlobSwitch shoudl be 0 or 1 # 1 means global (all lat and lon), 0 means regions '
    ' --slt = StLt starting SOUTHERNMOST latitude (not gridbox centre) '
    ' --sln = StLn starting WESTERNMOST longitude (not gridbox centre) '
    ' --elt = EdLt ending NORTHERNMOST latitude (not gridbox centre) '
    ' --eln EdLn ending EASTERNMOST longitude (not gridbox centre) '
    ' --ts TimeSeriesSwitch should be 0 or 1 # 1 = time series for each grid will be plotted '
    ' --sls SliceSwitch should be 0 or 1 # 1 = time average will be created, plotted and saved as a netCDF file'
    ' --sps SpatialSwitch should be 0 or 1 # 1 = time series of space average will be created, plotted and saved as a netCDF file'
    ' --mb MyBundle dataset to work on: '"HadISDH.landq.ID"' '
    '                                    "HadISDH.landRH.ID"'
    '                                    "HadISDH.landT.ID"'
    '                                    "q2m_monthly_5by5_ERA-Interim_data_"' 
    '                                    "q2m_monthly_5by5_ERA-Interim_data_abs"'
    '                                    "rh2m_monthly_5by5_ERA-Interim_data_"' 
    '                                    "t2m_monthly_5by5_ERA-Interim_data_"' 
    '                                    "si10_monthly_5by5_ERA-Interim_data_"'
    '                                    "msl_monthly_5by5_ERA-Interim_data_"'     

    #************************************************************************
    # INPUT PARAMETERS AS STRINGS!!!!
    # Generic things that we can edit for each run of the code or that can be passed from the command line
    # Date selection for creating a time average slice
    styr = 1979
    edyr = 2017
    stmn = 1  # 1-12
    edmn = 12  # 1-12
    
    SeasonSwitch = 1 # 1 means over selected seasons; 0 means over one entire year (12 month)
    SeasonChoice = 'JJA' # This can be any word from the SeasonDict dictionary below or add new entries to the dictionary

    GlobSwitch = 1 # 1 means global (all lat and lon), 0 means chosing regions
    StLt = -20 # -90 - 90
    StLn = -65 # - 180 - 180
    EdLt = 0 # -90 - 90
    EdLn = -35 # - 180 - 180

    TimeSeriesSwitch = 0 # If 1, time series for each grid will be produced.
    SliceSwitch = 0 # If 1, time averages will be created
    SpatialSwitch = 0 # If 1, space average will be created

    # Run choice bundle for input/output files, units, names, letters, read in varnames, colourmap
    # CHOOSE/ADD A DICTIONARY BUNDLE!!!
    #MyBundle = 'HadISDH.landq.ID'
    #MyBundle = 'HadISDH.landRH.ID'
    #MyBundle = 'HadISDH.landT.ID'
    MyBundle = 'q2m_monthly_5by5_ERA-Interim_data_' # anoms1981-2010
    #MyBundle = 'q2m_monthly_5by5_ERA-Interim_data_abs' #absolute
    #MyBundle = 'rh2m_monthly_5by5_ERA-Interim_data_' # anoms1981-2010
    #MyBundle = 't2m_monthly_5by5_ERA-Interim_data_' # anoms1981-2010
    #MyBundle = 'si10_monthly_5by5_ERA-Interim_data_' # anoms1981-2010
    #MyBundle = 'msl_monthly_5by5_ERA-Interim_data_' # anoms1981-2010    

    # Seasonal Info Dictionary    
    SeasonDict = dict([('Jan',[0]), # Arrays start from 0.
                       ('Feb',[1]),
		       ('Mar',[2]),
		       ('Apr',[3]),
		       ('May',[4]),
		       ('Jun',[5]),
		       ('Jul',[6]),
		       ('Aug',[7]),
		       ('Sep',[8]),
		       ('Oct',[9]),
		       ('Nov',[10]),
		       ('Dec',[11]),
		       ('DJF',[11,0,1]),
		       ('MAM',[2,3,4]),
		       ('JJA',[5,6,7]),
		       ('SON',[8,9,10]),
                       ('Annual',[0,1,2,3,4,5,6,7,8,9,10,11])]) # add extended autum, winter; and wet/dry seasons!

#### Overwriting variables with any given values from the command line and testing their validity ####
### This may not test validity of variables specified within the code rather than at the command line - check ####
    # Look for any variables specified at the command line #
    try:
        opts, args = getopt.getopt(argv, "hi:",
	                           ["sy=","ey=","sm=","em=","ss=","sc=",
				    "gs=","slt=","sln=","elt=","eln=",
				    "ts=","sls=","sps=","mb="])
    # If there is some command line input but it doesn't look right then exit with an error statement
    except getopt.GetoptError:
        print('Usage (as strings) PlotAverageMap_cartopy_Jul2018_timeplotaverage.py --sy <1973> --ey <2017> ... ') 
        sys.exit(2)

    # Go through each arguement entered at the command line and use it to set up the variables if valid
    for opt, arg in opts:
        if opt == "--sy":
            try:
                styr = int(arg)
            except:
                sys.exit("Failed: styr not an integer")
        elif opt == "--ey":
            try:
                edyr = int(arg)
            except:
                sys.exit("Failed: edyr not an integer")
        elif opt == "--sm":
            try:
                stmn = int(arg)
            except:
                sys.exit("Failed: stmn not an integer")
            else:
                if (stmn < 1) | (stmn > 12):
                    sys.exit('stmn too high or too low, should be 1-12')	       
        elif opt == "--em":
            try:
                edmn = int(arg)
            except:
                sys.exit("Failed: edmn not an integer")
            else:
                if (edmn < 1) | (edmn > 12):
                    sys.exit('edmn too high or too low, should be 1-12')	       
        elif opt == "--ss":
            try:
                SeasonSwitch = int(arg)
            except:
                sys.exit("Failed: SeasonSwitch not an integer")
            else:
                if (SeasonSwitch  not in [0,1]):
                    sys.exit('SeasonSwitch invalid - should be 0 or 1')	       
        elif opt == "--sc":
            try:
                SeasonChoice = str(arg)
            except:
                sys.exit("Failed: SeasonChoice not a string")
            else:
                if (SeasonChoice  not in SeasonDict):
                    sys.exit('SeasonChoice invalid - should be in SeasonDict')	       
        elif opt == "--gs":
            try:
                GlobSwitch = int(arg)
            except:
                sys.exit("Failed: GlobSwitch not an integer")
            else:
                if (GlobSwitch  not in [0,1]):
                    sys.exit('GlobSwitch invalid - should be 0 or 1')	       
        elif opt == "--slt":
            try:
                StLt = float(arg)
            except:
                sys.exit("Failed: StLt not a float")
            else:
                if (StLt < -90.) | (StLt > 90.):
                    sys.exit('StLt invalid - should be -90.0 to 90.0')	       
        elif opt == "--sln":
            try:
                StLn = float(arg)
            except:
                sys.exit("Failed: StLn not a float")
            else:
                if (StLn < -180.) | (StLn > 180.):
                    sys.exit('StLn invalid - should be -180.0 to 180.0')	       
        elif opt == "--elt":
            try:
                EdLt = float(arg)
            except:
                sys.exit("Failed: EdLt not a float")
            else:
                if (EdLt < -90.) | (EdLt > 90.):
                    sys.exit('EdLt invalid - should be -90.0 to 90.0')	       
        elif opt == "--eln":
            try:
                EdLn = float(arg)
            except:
                sys.exit("Failed: EdLn not a float")
            else:
                if (EdLn < -180.) | (EdLn > 180.):
                    sys.exit('EdLn invalid - should be -180.0 to 180.0')	       
        elif opt == "--ts":
            try:
                TimeSeriesSwitch = int(arg)
            except:
                sys.exit("Failed: TimeSeriesSwitch not an integer")
            else:
                if (TimeSeriesSwitch  not in [0,1]):
                    sys.exit('TimeSeriesSwitch invalid - should be 0 or 1')	       
        elif opt == "--sls":
            try:
                SliceSwitch = int(arg)
            except:
                sys.exit("Failed: SliceSwitch not an integer")
            else:
                if (SliceSwitch  not in [0,1]):
                    sys.exit('SliceSwitch invalid - should be 0 or 1')	       
        elif opt == "--sps":
            try:
                SpatialSwitch = int(arg)
            except:
                sys.exit("Failed: SpatialSwitch not an integer")
            else:
                if (SpatialSwitch  not in [0,1]):
                    sys.exit('SpatialSwitch invalid - should be 0 or 1')	       
        elif opt == "--mb":
            try:
                MyBundle = str(arg)
            except:
                sys.exit("Failed: MyBundle not a string")

# A sanity check that something is as expected
    assert edyr >= styr, "End year earlier than start year!!!"

    print('Set up for this run: ')
    print('Start and End Year: ',styr,edyr)
    print('Start and End Month: ',stmn,edmn)
    if (SeasonSwitch == 1):
        print('Creating Seasonal Average over',SeasonChoice)
    if (GlobSwitch == 0):
        print('Working on the region: ',StLt,StLn,EdLt,EdLn)
    if (TimeSeriesSwitch == 1):
        print('Creating plots of time series each gridbox')
    if (SliceSwitch == 1):
        print('Creating a time average for each gridbox')
    if (SpatialSwitch == 1):
        print('Creating a spatial average for each timepoint')
    print('Working on dataset: ',MyBundle)
 
    #pdb.set_trace() 

    # Set up directories and files
    # working directory
    INDIR='/data/local/hadkw/HADCRUH2/UPDATE2017/'
    # Output directory for plots
    OUTDIRP=INDIR+'/IMAGES/ANALYSIS/'
    # Output directory for grids
    OUTDIRDG=INDIR+'/STATISTICS/GRIDS/'
    # Output directory for time series
    OUTDIRDT=INDIR+'/STATISTICS/TIMESERIES/'

    # Land cover file:
    incover = INDIR+'OTHERDATA/HadCRUT.4.3.0.0.land_fraction'

#**********************************************************
# Things that we shouldn't need to edit but need to be set for the program to run
    # Missing data
    mdi = -1e30 # may set up as masked arrays later

    if (MyBundle == 'HadISDH.landq.ID'):
        candidate = 'STATISTICS/GRIDS/HadISDH.landq.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc'
        OUTPLOT = '_HadISDH.landq.4.0.0.2017f_'
        Unit = 'g kg$^{-1}$'  
        StdName = 'specific_humidity'
        Namey = 'HadISDH.landq.4.0.0.2017f slice '
        nlats = 36	       #set once file read in
        nlons = 72	       #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['q_anoms'])
    #    ReadInfo = list(['q_abs'])
        ColourMapChoice = ('BrBG','noflip')
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1973 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'
#        DataRes = 'abs'

    if (MyBundle == 'HadISDH.landRH.ID'):
        candidate = 'STATISTICS/GRIDS/HadISDH.landRH.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc'
        OUTPLOT = '_HadISDH.landRH.4.0.0.2017f_'
        Unit = '%rh'  
        StdName = 'relative_humidity'
        Namey = 'HadISDH.landRH.4.0.0.2017f slice '
        nlats = 36	       #set once file read in
        nlons = 72	       #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['rh_anoms'])
    #    ReadInfo = list(['rh_abs'])
        ColourMapChoice = ('BrBG','noflip')
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1973 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms' 
#        DataRes = 'abs' 

    if (MyBundle == 'HadISDH.landT.ID'):
        candidate='STATISTICS/GRIDS/HadISDH.landT.4.0.0.2017f_FLATgridIDPHA5by5_anoms8110_MAR2018_cf.nc'
        OUTPLOT = '_HadISDH.landT.4.0.0.2017f_'
        Unit='$^{o}$C'  #'degrees C'
        StdName = 'air_temperature'
        Namey = 'HadISDH.landT.4.0.0.2017f slice '
        nlats=36	      #set once file read in
        nlons=72	      #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['t_anoms'])
        ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1973 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'
#        DataRes = 'abs'

    if (MyBundle == 'q2m_monthly_5by5_ERA-Interim_data_'):
        candidate = 'OTHERDATA/q2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
        OUTPLOT = '_q2m_monthly_5by5_ERA-Interim_'
        Unit = 'g kg$^{-1}$'  
        StdName = 'specific_humidity'
        Namey = 'q2m_monthly_5by5_ERA-Interim '
        nlats = 36	       #set once file read in
        nlons = 72	       #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['anomalies_land'])
        ColourMapChoice = ('BrBG','noflip')
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'

    if (MyBundle == 'q2m_monthly_5by5_ERA-Interim_data_abs'):
        candidate = 'OTHERDATA/q2m_monthly_5by5_ERA-Interim_data_19792017.nc'
        OUTPLOT = '_q2m_monthly_5by5_ERA-Interim_'
        Unit = 'g kg$^{-1}$'
        StdName = 'specific_humidity'
        Namey = 'q2m_monthly_5by5_ERA-Interim '
        nlats = 36           #set once file read in
        nlons = 72           #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['actuals'])
        ColourMapChoice = ('BrBG','noflip')
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'abs'

    if (MyBundle == 'rh2m_monthly_5by5_ERA-Interim_data_'):
        candidate = 'OTHERDATA/rh2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
        OUTPLOT = '_rh2m_monthly_5by5_ERA-Interim_'
        Unit = '%rh'  
        StdName = 'relative_humidity'
        Namey = 'rh2m_monthly_5by5_ERA-Interim slice '
        nlats = 36	       #set once file read in
        nlons = 72	       #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['anomalies_land'])
        ColourMapChoice = ('BrBG','noflip')
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'

    if (MyBundle == 't2m_monthly_5by5_ERA-Interim_data_'):
        candidate = 'OTHERDATA/t2m_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
        OUTPLOT = '_t2m_monthly_5by5_ERA-Interim_'
        Unit='$^{o}$C'  #'degrees C'
        StdName = 'air_temperature'
        Namey = 't2m_monthly_5by5_ERA-Interim slice '
        nlats=36	      #set once file read in
        nlons=72	      #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['anomalies_land'])
        ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'

    if (MyBundle == 'si10_monthly_5by5_ERA-Interim_data_'):
        candidate = 'OTHERDATA/si10_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
        OUTPLOT = '_si10_monthly_5by5_ERA-Interim_'
        Unit='m s$^{-1}$'
        StdName = 'wind_speed'
        Namey = 'si10_monthly_5by5_ERA-Interim slice '
        nlats=36	      #set once file read in
        nlons=72	      #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['anomalies_land'])
        ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'

    if (MyBundle == 'msl_monthly_5by5_ERA-Interim_data_'):
        candidate = 'OTHERDATA/msl_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
        OUTPLOT = '_msl_monthly_5by5_ERA-Interim_'
        Unit='hPa'
        StdName = 'mean_sea_level_pressure'
        Namey = 'msl_monthly_5by5_ERA-Interim slice '
        nlats=36	      #set once file read in
        nlons=72	      #set once file read in
        LatInfo = list(['latitude',nlats,-87.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        LonInfo = list(['longitude',nlons,-177.5,5.]) # name of latitudes in file, number of latitudes, starting latitude gribox centre, latitude box height
        ReadInfo = list(['anomalies_land'])
        ColourMapChoice=('coolwarm','noflip') # could also be 'bwr'
        IsLand = True # True for land, False for marine, None for blend
        DataStYr = 1979 # ERA starts in 1979
        DataStMn = 1 # 1 to 12
        DataRes = 'anoms'

    # Time elements
    NYears = (edyr-styr)+1
    NMonths = NYears*12

    # Season elements
    if SeasonSwitch == 1:
        Packages = 'years'
    else:
        Packages = 'month'
		       
    # Creates an array of month (better than a list here) that that could creates averages over repetitive periods, i.e. seasons and single months 
    SeasonInfo = np.array(SeasonDict[SeasonChoice])

    if (GlobSwitch == 1): # if global then we don't need lat lon info in title
        if (SeasonSwitch == 0): # if for all month we don't need season info in the title
            slicechoice = str(styr) + "%02d" % stmn + '_' + str(edyr) + "%02d" % edmn +'_'+DataRes
        else:
            slicechoice = str(styr) + "%02d" % stmn + '_' + str(edyr) + "%02d" % edmn + '_'+SeasonChoice+'_'+DataRes
    else:
        if (SeasonSwitch == 0): # if for all month we don't need season info in the title
            slicechoice = str(styr) + "%02d" % stmn + '_' + str(edyr) + "%02d" % edmn + '_'+'_'+DataRes+"%02d" % StLt + "%02d" % StLn + "%02d" % EdLt + "%02d" % EdLn
        else:
            slicechoice = str(styr) + "%02d" % stmn + '_' + str(edyr) + "%02d" % edmn + '_'+SeasonChoice+'_'+'_'+DataRes + "%02d" % StLt + "%02d" % StLn + "%02d" % EdLt + "%02d" % EdLn

    OUTPLOTMAP = 'TIMEAVERAGE_'+OUTPLOT+slicechoice
    OUTPLOTTS = 'REGIONAVERAGE_'+OUTPLOT+slicechoice
    Namey = Namey+slicechoice
    
    # Variables
    InputData = []
    LatList = []
    LonList = []

####### END OF DEFINITIONS ###########
####### START OF MAIN PROGRAM ########
    
    # Test for silly values
    if ((DataStMn < 1) | (DataStMn > 12)):
        print("Silly DataStMn - should be between 1 and 12")
        pdb.set_trace()

    if ((stmn < 1) | (stmn > 12)):
        print("Silly stsn - should be between 1 and 12")
        pdb.set_trace()

    if ((edmn < 1) | (edmn > 12)):
        print("Silly edsn - should be between 1 and 12")
        pdb.set_trace()

    if (styr < DataStYr):
        print("Silly styr or DataStYr - styr should be greater than or equal to DataStYr")
        pdb.set_trace()

    if (edyr < styr):
        print("Silly edyr - should be greater than styr")
        pdb.set_trace()
	
	
    ############################ Setting up file name for NetCDF file to read in
    
    MyFile = INDIR+candidate
  
    ########################### Select time - Note that two data sets can have two different periods of time.
    
    TimeInfo = WhichMonths(DataStYr,DataStMn,styr,edyr,stmn,edmn)
    # This is a stopper that stops the code and allows you to print and play with the variables
    # Type c <return> to carry on
    print('TimeInfo: ',TimeInfo)
    #pdb.set_trace()
    
    ########################### READ IN DATA FROM NETCDF - Select a region if GlobSwitch = 0
    ########################### Read data in from NetCDF for selections

    if GlobSwitch == 1: # means all data, no region choice
        # Read in Data and Pull out only the time slice of interest.
        CandData,LatList,LonList = ReadNetCDFGrid(MyFile,ReadInfo,LatInfo,LonInfo,TimeInfo)
        print('Read in netcdf for global')
        #pdb.set_trace()
	
    else:
        # Get the region info to pull out 
        RegionInfo = WhichRegions(LatInfo,LonInfo,StLt,EdLt,StLn,EdLn)
        print('RegionInfo: ',RegionInfo)
        #pdb.set_trace()

        # Read in Data and Pull out only the time and region slice of interest.
        CandData,LatList,LonList = ReadNetCDFGridRegion(MyFile,ReadInfo,LatInfo,LonInfo,TimeInfo,RegionInfo)
        print('Read in netcdf for region')
        # reset nlons and nlats to the region
        nlons = len(LonList)
        nlats = len(LatList)
        #pdb.set_trace()

    ########################### Select a season if SeasonSwitch = 1
    ########################### Read data in from NetCDF for selections

    if SeasonSwitch == 1: # 0 means all data, no season choice
        # Read in Data and Pull out only the time slice of interest.
        print('SeasonInfo: ',SeasonInfo)
        NewCandData = CreateSeasons(CandData,mdi,nlats,nlons,SeasonInfo) 
        print('Averaged to season')
        #pdb.set_trace()
	
    else:

        # Changing the name of CandData to NewCandData
        NewCandData = CandData #for seperate copy it would be: np.copy()
        print('Renamed data array for further processing')
        #pdb.set_trace()
    
    ############################ Create a spatial average of the data for each time step

    if (SpatialSwitch == 1):
        # Create Spatial average of slices
        MySpatials = CreateSpatials(NewCandData,mdi,LatList,LonList)
        print('Averaged over region for each time point')

        # pass to plotter
        MyFile = OUTDIRP+OUTPLOTTS    
        PlotTimeSeries(MyFile,MySpatials,Packages,Unit,Namey,NMonths,NYears,styr,edyr,mdi)
        print('Plotted single time series')

       # Write the time series of data out to a netCDF file
        # Dictionary for looking up variable names for netCDF read in of variables
        PassingOutDict = dict([('VarName','spatial_average'),
                 ('StdName',StdName),
                 ('LongName',slicechoice),
                 ('Unit',Unit),
                 ('vmin',np.min(MySpatials[np.where(MySpatials > mdi)])),
                 ('vmax',np.max(MySpatials)),
                 ('mdi',mdi)])

        OutFile = OUTDIRDT+OUTPLOTTS+'.nc'
	# pass 0s for LatList,LonList,nlats and nlons to tell the code its a time series
        WriteOut(OutFile,MySpatials,0,0,styr,stmn,edyr,edmn,NMonths,0,0,PassingOutDict)
        print('Written time series netCDF')

    ############################ Create a time average of the data for each gridbox

    if (SliceSwitch == 1):
        # Create Time average of slice
        # This has the option of plotting a time series on the fly for each gridbox up to 10 figures
        MySlice = CreateSlice(NewCandData,mdi,nlats,nlons,LatList,LonList,Packages,Unit,TimeSeriesSwitch)
        print('Averaged over all time points for each gridbox')

        # pass to plotter
        MyFile = OUTDIRP+OUTPLOTMAP    
        PlotMap(MyFile,LatList,LonList,MySlice,Unit,Namey,ColourMapChoice,IsLand,LatInfo,LonInfo,GlobSwitch)
        print('PLotted map of time average grids')
	
        # Write the slice of data out to a netCDF file
        # Dictionary for looking up variable names for netCDF read in of variables
        PassingOutDict = dict([('VarName','time_average'),
                 ('StdName',StdName),
                 ('LongName',slicechoice),
                 ('Unit',Unit),
                 ('vmin',np.min(MySlice[np.where(MySlice > mdi)])),
                 ('vmax',np.max(MySlice)),
                 ('mdi',mdi)])

        OutFile = OUTDIRDG+OUTPLOTMAP+'.nc'
	# pass 0s for NMonts to tell the code its a map with no time elements
        WriteOut(OutFile,MySlice,LatList,LonList,styr,stmn,edyr,edmn,0,nlats,nlons,PassingOutDict)
        print('Written out time average map to netCDF')

    ########################### Plot all gridbox time series if a region and TimeSeriesSwitch is 1

    if (GlobSwitch == 0) & (TimeSeriesSwitch == 1):
        MyFile = OUTDIRP+OUTPLOTTS+'_multi'
        PlotMultipleTimeSeries(MyFile,NewCandData,LatList,LonList,Packages,Unit,Namey,NMonths,NYears,styr,edyr,mdi)
        print('Plotted multiple time series')


    print("And, we are done!")

######### END OF MAIN PROGRAM ##############################
#************************************************************************
# Subroutines
#************************************************************************
# READNETCDFGRID # NETCDF4
def ReadNetCDFGrid(FileName,ReadInfo,LatInfo,LonInfo,TimeInfo):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data 
	FileName: stroing containing filepath/name
	ReadInfo: list of strings of variable names for the trend, 5th percentile, 95th percentile
	LatInfo: list of a string for the latitude variable name, an integer for number of lats, a float for the start latitude
	LonInfo: list of a string for the longitude variable name, an integer for number of lons, a float for the start longitude 
	RegionInfo: If chosen GlobSwitch = 0, then  list of a string for the chosen regions'''

    ncf = nc4.Dataset(FileName,'r')

    # ncf.variables this lists the variable names
    # var=f.variables['latitude']
    # TheLatList=var.data
    # lats currently screwy so make up
    if (LatInfo[2] < 0):
        TheLatList=np.arange(LatInfo[2], LatInfo[2]+180.,(180./LatInfo[1]))
    else:
        TheLatList=np.arange(LatInfo[2], LatInfo[2]-180.,-(180./LatInfo[1]))    
    # var=f.variables['longitude']
    # TheLonList=var.data
    # lons currently screwy so make up
    #TheLonList=np.arange(-177.5,182.5,LonInfo[3])
    if (LonInfo[2] < 10):
        TheLonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
    else:
        TheLonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    

    var=ncf.variables[ReadInfo[0]]
    
    # You should be able to slice the data here by selecting only a few time points
    # TheData = np.array( [TimeInfo[0]:TimeInfo[1],:,:] # select period already
    # You could do this by adding a list called TimeInfo to the things called into this function at the top
    # This would have the start and end years and months which you can then work out which months to pull out
    # You would then need to know the start year and month of the dataset you're working with - HadISDH is Jan 1973, ERA is Jan 1979
#    TheData = np.array(var.data)
    TheData = np.array(var[TimeInfo[0]:TimeInfo[1],:,:])

#    var=ncf.variables[ReadInfo[1]]
#    TheLower=np.array(var.data)
#    var=ncf.variables[ReadInfo[2]]
#    TheUpper=np.array(var.data)

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#***********************************************************************
# READNETCDFGRID WITH REGIONINFO
def ReadNetCDFGridRegion(FileName,ReadInfo,LatInfo,LonInfo,TimeInfo,RegionInfo):
    ''' Open the NetCDF File
        Get the list of latitudes
        Get the list of longitudes
        Get the data 
	FileName: stroing containing filepath/name
	ReadInfo: list of strings of variable names for the trend, 5th percentile, 95th percentile
	LatInfo: list of a string for the latitude variable name, an integer for number of lats, a float for the start latitude
	LonInfo: list of a string for the longitude variable name, an integer for number of lons, a float for the start longitude 
	RegionInfo: If chosen GlobSwitch = 0, then  list of a string for the chosen regions'''

    ncf = nc4.Dataset(FileName,'r')

    # ncf.variables this lists the variable names
    # var=f.variables['latitude']
    # TheLatList=var.data
    # lats currently screwy so make up
    if (LatInfo[2] < 0):
        TheLatList=np.arange(LatInfo[2], LatInfo[2]+180.,(180./LatInfo[1]))
    else:
        TheLatList=np.arange(LatInfo[2], LatInfo[2]-180.,-(180./LatInfo[1]))    
    # var=f.variables['longitude']
    # TheLonList=var.data
    # lons currently screwy so make up
    #TheLonList=np.arange(-177.5,182.5,LonInfo[3])
    if (LonInfo[2] < 10):
        TheLonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
    else:
        TheLonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    

    TheLatList = TheLatList[RegionInfo[0]:RegionInfo[1]]
    TheLonList = TheLonList[RegionInfo[2]:RegionInfo[3]]
   
    var=ncf.variables[ReadInfo[0]]
    
    # You should be able to slice the data here by selecting only a few time points
    # TheData = np.array( [TimeInfo[0]:TimeInfo[1],:,:] # select period already
    # You could do this by adding a list called TimeInfo to the things called into this function at the top
    # This would have the start and end years and months which you can then work out which months to pull out
    # You would then need to know the start year and month of the dataset you're working with - HadISDH is Jan 1973, ERA is Jan 1979
#    TheData = np.array(var) # old version (var.data), new version for NetCDF: var
    TheData = np.array(var[TimeInfo[0]:TimeInfo[1],RegionInfo[0]:RegionInfo[1],RegionInfo[2]:RegionInfo[3]])

#    var=ncf.variables[ReadInfo[1]]
#    TheLower=np.array(var.data)
#    var=ncf.variables[ReadInfo[2]]
#    TheUpper=np.array(var.data)

#    # Maybe I've done something wrong but its reading it transposed
#    TheData=np.transpose(TheData)
    ncf.close()
    
    return TheData,TheLatList,TheLonList # ReadNetCDFGrid

#***********************************************************************
# WhichMonths
def WhichMonths(DStYr,DStMn,StYr,EdYr,StMn,EdMn):
    ''' This function works out the months to point to when '''
    ''' selecting a slice based on the actual start year and month of the ''' 
    ''' data set '''
    ''' It returns a list called TimeInfo with the months to point to '''
    ' e.g., DStYr = 1979, DStMn = 1 '
    ' e.g., StYr = 1979, EdYr = 1979 '
    ' e.g., StMn = 1, EdMn = 3 '
    
    # Set up TimeInfo list
    TimeInfo = list([0,0])
    
    # Find the starting month
    # Start by assuming all datasets start in January and then have a check later
    YearCount = (StYr - DStYr) * 12
    
    MonthCount = StMn - 1
    
#    print(YearCount,MonthCount)

    if (DStMn != 1):
        YearCount = YearCount - (DStMn-1)
    
    TimeInfo[0] = YearCount + MonthCount
    
#    print(YearCount,MonthCount,TimeInfo[0])

    # Find the ending month
    # Start by assuming all datasets start in January and then have a check later
    YearCount = (EdYr - DStYr) * 12
    
    MonthCount = EdMn - 1

#    print(YearCount,MonthCount)
    
    if (DStMn != 1):
        YearCount = YearCount - (DStMn-1)
    
    TimeInfo[1] = YearCount + MonthCount + 1
    
#    print(YearCount,MonthCount,TimeInfo[1])

    return TimeInfo

#**********************************************************************
# WhichRegions
def WhichRegions(LatInfo,LonInfo,StLt,EdLt,StLn,EdLn):
    ''' This function works out the cells to point to when '''
    ''' selecting a slice based on the actual start lat and lon of the ''' 
    ''' data set '''
    ''' It returns a list called RegionInfo with the coordinated to point to '''
    ' e.g., DStYr = 1, DStMn = 1 '
    ' e.g., StLt = 1, EdLt = 1 '
    ' e.g., StLn = 1, EdLn = 3 '
    
    # Set up RegionInfo list
    RegionInfo = list([0,0,0,0])
       
    if (LatInfo[2] < 0): # second point of LatInfo = -85.7 (starting latitude gridbox)
        TheLatList=np.arange(LatInfo[2], LatInfo[2]+180.,(180./LatInfo[1]))
    else:
        TheLatList=np.arange(LatInfo[2], LatInfo[2]-180.,-(180./LatInfo[1]))    
    # var=f.variables['longitude']
    # TheLonList=var.data
    # lons currently screwy so make up
    if (LonInfo[2] < 10):
        TheLonList=np.arange(LonInfo[2], LonInfo[2]+360.,(360./LonInfo[1]))
    else:
        TheLonList=np.arange(LonInfo[2], LonInfo[2]-360.,-(360./LonInfo[1]))    

    LatPointers = np.where(((TheLatList+(LatInfo[3]/2.)) > StLt) & ((TheLatList-(LatInfo[3]/2.)) < EdLt))
    LonPointers = np.where(((TheLonList+(LonInfo[3]/2.)) > StLn) & ((TheLonList-(LonInfo[3]/2.)) < EdLn))
    
    RegionInfo = list([LatPointers[0][0],(LatPointers[0][-1]+1),LonPointers[0][0],(LonPointers[0][-1]+1)])

    return RegionInfo
    
#**********************************************************************
# CreateSeasons
def CreateSeasons(TheData,TheMDI,NLats,NLons,SeasonInfo):
    ''' This is a function to create a season, or a selection of month of all data '''
    ''' from monthly mean grids of latitudes and longitudes; time will always be 12 month * X '''
    
    # Create a for loop that loops through each year and averages all time points in
    # the selected period.
    
    # Make sure you check for missing data and only avereage over the data points that are non-missing
    # np.where! 
    
    # existing /  missing data threshold
    pTresh = 0.6 # MissThresh = 0.4 # If 2/3 month are available the average will still be calculated.
    
    # empty array for filling with time means - set up with missing data indicator
    #pdb.set_trace
    NYears = np.int(len(TheData[:,0,0])/12)
    OutputData = np.empty((NYears,NLats,NLons))
    OutputData.fill(TheMDI)
    
    for lt in range(NLats): # range nlats goes through 0...36
        for ln in range(NLons): # range nlons goes through 0...72

	    # CHECK TO MAKE SURE THERE ARE SOME DATA TO WORK ON
            if len(np.where(TheData[:,lt,ln] > TheMDI)[0]) > 0:

                MonthArray = np.reshape(TheData[:,lt,ln],(NYears,12)) # creates a new array with 12 colums and NYears rows
                #pdb.set_trace()
	   
	    
                for yr in range(NYears): # now on year-level; goes through all the years
                    
		    # If we are working on DJF or a complex season that crosses over the year end:
                    if (len(SeasonInfo) > 1) & (SeasonInfo[0] > SeasonInfo[-1]): # Seasonal average over winter, DJF: first given month (e.g. Dec, 12) is greater than the last given month (e.g. Feb, 2)
                        earlyMap = np.where(SeasonInfo < SeasonInfo[0]) # Creates two different maps. The early map would be for January, February...
                        lateMap = np.where(SeasonInfo > SeasonInfo[-1]) # The late map would be for December, also for November...
                        #pdb.set_trace()
			
			# Create a single array of data that contain the months we want to average over
                        SubArray = [] # set up an empty thing
                        if (yr > 0):
                            SubArray = np.concatenate((SubArray,MonthArray[yr-1,SeasonInfo[lateMap]]),0) # join the empty thing with the late year elements of MonthArray from the year before into a np.array 
			
                        SubArray = np.concatenate((SubArray,MonthArray[yr,SeasonInfo[earlyMap]]),0) # join the SubArray with the early year elements of MonthArray from the year of the loop
                        #pdb.set_trace()
			
			# test the SubArray to see if there is any/enough data in there to create a seasonal average
			
                        gotdata = np.where(SubArray > TheMDI)
			
                        if (len(gotdata[0]) / np.float(len(SeasonInfo)) > pTresh):
                            OutputData[yr,lt,ln] = np.mean(SubArray[gotdata])
                            #pdb.set_trace()
                    
		    # Here we are working on an easy season that is all within a year
                    else:
		    
	                # create a map pointing to the elements of the array that are non-missing
                        gotdata = np.where(MonthArray[yr,SeasonInfo] > TheMDI)
			#pdb.set_trace()
			
		    	   
	    # now test to check whether we have enough data - greater than the missing data threshold
                        if (len(gotdata[0]) / np.float(len(SeasonInfo)) > pTresh):
			    # This seems a very long winded way of using gotdata to point to elements of SeasonInfo that points to elements of MonthArray[yr,:]
			    # In IDL I would typ MonthArray[yr,SeasonInfo[gotdata]] but that doesn't work here because gotdata is an array not a list.
			    # Actually if we just set SeasonInfo up as a np.array in the first place then it works because we are then subsetting an array
			    # with an array and then a list rather than an array with a list and then an array
                            #SubArr = MonthArray[yr,SeasonInfo]
                            #OutputData[yr,lt,ln] = np.mean(SubArr[gotdata])
                            OutputData[yr,lt,ln] = np.mean(MonthArray[yr,SeasonInfo[gotdata]])
                            #pdb.set_trace()
			    
    return OutputData

#***********************************************************************    
# Create Spatial average of slices
def CreateSpatials(TheData,TheMDI,TheLatList,TheLonList):

    RegionMeanVal = np.empty(len(TheData[:,0,0]))
    RegionMeanVal.fill(TheMDI)
    
    for n in range(len(TheData[:,0,0])): # for lat and lons the time elements are the same, so 0,0
    
        TotWeights = 0
        TotVals = 0
	
        for ltt in range(len(TheLatList)):
            TotVals = TotVals+np.sum(np.cos(np.radians(TheLatList[ltt]))*TheData[n,ltt,np.where(TheData[n,ltt,:] > TheMDI)[0]])
            TotWeights = TotWeights+np.sum(np.cos(np.radians(TheLatList[ltt]))*len(np.where(TheData[n,ltt,:] > TheMDI)[0]))
    
        #pdb.set_trace()
        RegionMeanVal[n] = TotVals/TotWeights # 1 dim array
    
    return RegionMeanVal

#**********************************************************************
# CreateSlice
def CreateSlice(TheData,TheMDI,NLats,NLons,TheLatList,TheLonList,Packages,Unit,TimeSeriesSwitch):
    ''' This is a function to create a time averaged slice of data '''
    ''' from monthly mean grids of latitudes and longitudes '''
    ' ROWS THEN COLUMNS - LATS THEN LONS'
    
    # Create a for loop that loops through each grid box and averages all time points in
    # the slice.
    
    # Make sure you check for missing data and only avereage over the data points that are non-missing
    # np.where! 
    
    # existing /  missing data threshold
    pTresh = 0.8 # MissThresh = 0.2
    
    # empty array for filling with time means - set up with missing data indicator
    #pdb.set_trace()
    OutputData = np.empty((NLats,NLons))
    OutputData.fill(TheMDI)
    
    PlotCounter = 0
    
    avg = 0. # At the beginning the average is zero, then it will be calculated.

    for lt in range(NLats): # range nlats goes through 0...36
        for ln in range(NLons): # range nlons goes through 0...72
#	    p = 0. # p is a counter for existing months
#	    q = 0. # q is a counter for all months
				

## Long hand way of testing for missing data
#	    for month in range(Month):
#    		q = q+1 # Just counts that it\u2019s processing one month.					# avg = avg/p
#		if TheData[month,lt,ln] > TheMDI : # means they are valid
#                    avg = avg + TheData[month,lt,ln]
#                    p = p+1 # count up p
#		    if p/q >= pTresh:
#		        OutputData[lt,ln] = avg/p
					
# Short way
            # create a map pointing to the elements of the array that are non-missing
            gotdata = np.where(TheData[:,lt,ln] > TheMDI)
	    
	    # now test to check whether we have enough data - greater than the missing data threshold
            if (len(gotdata[0]) / len(TheData[:,lt,ln]) > pTresh):
                OutputData[lt,ln] = np.mean(TheData[gotdata,lt,ln])
		# OutputData[lt,ln] = np.corrcoef(TheData1[gotdata,lt,ln],TheData2[gotdata,lt,ln])[0,1] # This is for Corr. What does [0,1] mean?
                # checking this out in more detail by stopping for each gridbox and plotting
		# TheData1 in blue where it is present
		# Just checking it looks sensible
		# Could add an if loop and specify an exact lt and ln to look at a specific grid box of interest
        # Or at averages over many grid boxes
        #pdb.set_trace()
                if (TimeSeriesSwitch == 1)&(PlotCounter < 10): # this is a regional projection
                    plt.clf()
                    fig,ax = plt.subplots()
                    ax.plot(range(len(gotdata[0])),TheData[gotdata,lt,ln][0])
                    #pdb.set_trace()
                    plt.xlabel('Time ['+Packages+']')
                    plt.ylabel('Anomaly ['+Unit+']')
                    plt.figtext(0.5,0.8,'Data presence threshold '+"%4.3f" % (np.float(len(gotdata[0])) / len(TheData[:,lt,ln])),size=12,ha='left')
                    plt.figtext(0.5,0.75,'TheData',color='blue',size=10,ha='left')
                    plt.figtext(0.5,0.2,'Lat: '+"%6.3f" % TheLatList[lt]+', Lon: '+"%7.3f" % TheLonList[ln],color='black',size=10,ha='left')
                    plt.show()
                    PlotCounter = PlotCounter+1
                    #pdb.set_trace()
                    #else: # this is global projection
    return OutputData

#**********************************************************************
# PlotMultipleTimeSeries
def PlotMultipleTimeSeries(TheFile,TheData,TheLatList,TheLonList,ThePackages,TheUnitee,TheNamey,TheMCount,TheYCount,TheStYr,TheEdYr,TheMDI):
    ''' Plot a panel for each element of TheHvars '''
    ''' Add Coverage, Sampling and Station uncertainty ranges '''
    ''' Add lines for any extra estimates (TheVars) and HadISDH MASKED versions '''
    ''' Save as png and eps '''
    ''' TheHvars is a multi-row array: rows for vars, columns for months '''
    ''' Ditto TheHuncs C=coverage, Sp=sampling, St=station '''
    ''' TheLablees is the name list for all other vars '''
    ''' TheUnitees is the units name list '''
    ''' TheVars is a multirow array if there is 1+ var available - or [] '''
    ''' TheMASKVars - ditto above but masked to HadISDH coverage '''
    ''' TheMCount - number of months, TheStYr/EdYr - start and end years '''
    ''' TheMDI - missing data indicator for masking '''
    ''' TheColls - dictionary of colours for each dataset '''
    
    
    # set up number of panels and number of pages
    nplots=len(TheLatList)*len(TheLonList)
    print('PLOT NUMBERS: ',nplots)
    NRows = 5; NColls = 5
    NPageColls = np.ceil(len(TheLonList)/np.float(NColls))
    NPageRows = np.ceil(len(TheLatList)/np.float(NRows))
    TotPages = NPageColls*NPageRows
    print('NPageColls:',NPageColls)
    print('NPageRows:',NPageRows)
    print('TotPages:',TotPages)
    #pdb.set_trace()
    Letteree=[chr(i) for i in np.arange(97,97+np.int(NRows*NColls))] # a = 97 # Quick Loop to creat a range and give it a character

    # set up x axes
    if ThePackages == 'month':
        TheMonths=[]
        yr=TheStYr
        mon=1
        for m in range(TheMCount):
            TheMonths.append(dt.date(yr,mon,1))
            mon=mon+1
            if mon == 13:
                mon=1
                yr=yr+1
        TheMonths=np.array(TheMonths)
    else:
        TheMonths=[]
        yr=TheStYr
        mon=1
        for y in range(TheYCount):
            TheMonths.append(dt.date(yr,mon,1))
            yr=yr+1
        TheMonths=np.array(TheMonths)

#pdb.set_trace()

    xtitlee='Years'

    latpoints = 0 # in total, 0...NLats for the whole thing, does not care about pages
    lonpoints = 0

    for tp in range(np.int(TotPages)):
        # set up dimensions and plot
        print('Page:',tp)
        xpos=[]
        ypos=[]
        xfat=[]
        ytall=[]
        totalyspace=0.8    # start 0.15 end 0.95, 80/5 = 16 for each plot
        totalxspace=0.8    # start 0.15 end 0.95
        xwidth=0.8/NColls
        ywidth=0.8/NRows # height
        for c in range(NColls):
            for r in range(NRows):
                xpos.append(0.15+r*xwidth)
                ypos.append(0.15+c*ywidth)
                xfat.append(xwidth)
                ytall.append(ywidth)
        print('Set up plot positions')
        #pdb.set_trace()
	
	# Check where we are with latpoints and lonpoints
        latpoints = tp*NRows # e.g., 0, 5, 10 etc
        if (latpoints >= len(TheLatList)): # we've over shot so go back
            latpoints = NRows*(tp-1) # e.g., 0, 5, 10 etc
        lonpoints = tp*NColls
        if (lonpoints >= len(TheLonList)):
            lonpoints = NColls*(tp-1)

        #plt.clf()
        #    fig = plt.figure(1,figsize=(8,12))
        #    plt.axes([0.15,0.1,0.8,0.80])
        f,axarr=plt.subplots(NColls*NRows,figsize=(10,10)) #)    #6,18

        CollPoint = 0 # just for the pages, always 0...5
        RowPoint = 0

        for pp in range(NColls*NRows):
            axarr[pp].set_position([xpos[pp],ypos[pp],xfat[pp],ytall[pp]])
            print(xpos[pp],ypos[pp],xfat[pp],ytall[pp])
            
            print('Plot positions')
            #pdb.set_trace()

            if ThePackages == 'month':
                axarr[pp].set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
                lw = 0.5
            else:
                axarr[pp].set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
                axarr[pp].set_ylim([np.floor(np.min(TheData[np.where(TheData > TheMDI)])), np.ceil(np.max(TheData))])
                lw=2
      
            # Test whether we actually have data for this lat/lon col/row or whether we are plotting a null
            if (latpoints >= len(TheLatList)) | (lonpoints >= len(TheLonList)):
                print('Printing null timeseries')
                axarr[pp].plot(TheMonths,np.zeros(len(TheMonths)),c='black',linewidth=lw)
            else:
                DataMask = np.where(TheData[:,latpoints,lonpoints] > TheMDI)
                axarr[pp].plot(TheMonths[DataMask],TheData[DataMask[0],latpoints,lonpoints],c='black',linewidth=lw)
                axarr[pp].annotate('Lat: '+"%6.3f" % TheLatList[latpoints]+', Lon: '+"%7.3f" % TheLonList[lonpoints],xy=(0.03,0.8),xycoords='axes fraction',color='black',size=8)

            axarr[pp].annotate(Letteree[pp]+')',xy=(0.03,0.9),xycoords='axes fraction',size=10)

            
	    # Now we only want axis labels and ticks along the left and bottom
	    # need to do some work on the xticks/xticklabels because there are too many so it looks rubbish
            if RowPoint == 0:
                axarr[pp].set_xlabel(xtitlee,fontsize=10)
            else:
                axarr[pp].set_xticklabels([]) # set xticklabels to blank if not bottom row
            if CollPoint == 0:
                axarr[pp].set_ylabel(TheUnitee,fontsize=10)
            else:
                axarr[pp].set_yticklabels([]) # set ytickabels to black if not left most plots

            print(latpoints,lonpoints,CollPoint,RowPoint)
#pdb.set_trace()
            CollPoint = CollPoint+1
            lonpoints = lonpoints+1
#print((NColls-1),(len(TheLonList)-1),(NRows-1),(len(TheLatList)-1))
            if (CollPoint > (NColls-1)): # | (lonpoints == len(TheLonList)):
                print('ColumnLoop:')
                lonpoints = lonpoints-CollPoint
                CollPoint = 0
                RowPoint = RowPoint+1
                latpoints = latpoints+1
            if (RowPoint > (NRows-1)): # & (latpoints < len(TheLatList)):
                print('RowLoop:')
                latpoints = latpoints-RowPoint
                RowPoint = 0
                CollPoint = CollPoint+1
                lonpoints = lonpoints+1

# Figure Watermark and Labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


#plt.show()
        plt.savefig(TheFile+'_'+str(tp)+".eps")
        plt.savefig(TheFile+'_'+str(tp)+".png")
	
    return

#**********************************************************************
# PlotTimeSeries
def PlotTimeSeries(TheFile,TheData,ThePackages,TheUnitee,TheNamey,TheMCount,TheYCount,TheStYr,TheEdYr,TheMDI):
    ''' Plot a panel for each element of TheHvars '''
    ''' Add Coverage, Sampling and Station uncertainty ranges '''
    ''' Add lines for any extra estimates (TheVars) and HadISDH MASKED versions '''
    ''' Save as png and eps '''
    ''' TheHvars is a multi-row array: rows for vars, columns for months '''
    ''' Ditto TheHuncs C=coverage, Sp=sampling, St=station '''
    ''' TheLablees is the name list for all other vars '''
    ''' TheUnitees is the units name list '''
    ''' TheVars is a multirow array if there is 1+ var available - or [] '''
    ''' TheMASKVars - ditto above but masked to HadISDH coverage '''
    ''' TheMCount - number of months, TheStYr/EdYr - start and end years '''
    ''' TheMDI - missing data indicator for masking '''
    ''' TheColls - dictionary of colours for each dataset '''
    
    # set up x axes
    if ThePackages == 'month':
        TheMonths=[]
        yr=TheStYr
        mon=1
        for m in range(TheMCount):
            TheMonths.append(dt.date(yr,mon,1))
            mon=mon+1
# FLORENTINE WILL NEED TO MOVE THIS IN TOO
            if mon == 13:
                mon=1
                yr=yr+1
        TheMonths=np.array(TheMonths)
    else:
        TheMonths=[]
        yr=TheStYr
        mon=1
        for y in range(TheYCount):
            TheMonths.append(dt.date(yr,mon,1))
            yr=yr+1
        TheMonths=np.array(TheMonths)

    #pdb.set_trace()

    xtitlee='Years'

    plt.clf()
    f = plt.figure(1,figsize=(12,6))
    axarr = plt.axes([0.15,0.1,0.8,0.80])

    print('Plot positions')
#    pdb.set_trace()

    if ThePackages == 'month':
        axarr.set_xlim([TheMonths[0],TheMonths[TheMCount-1]])
        lw = 0.5
    else:
        axarr.set_xlim([TheMonths[0],TheMonths[TheYCount-1]])
        axarr.set_ylim([np.floor(np.min(TheData[np.where(TheData > TheMDI)])), np.ceil(np.max(TheData))])
        lw=2
      
#    pdb.set_trace()

    DataMask = np.where(TheData > TheMDI)
    axarr.plot(TheMonths[DataMask],TheData[DataMask[0]],c='black',linewidth=lw)

#    pdb.set_trace()
  
# Now we only want axis labels and ticks along the left and bottom
# need to do some work on the xticks/xticklabels because there are too many so it looks rubbish
    axarr.set_xlabel(xtitlee,fontsize=10)
    axarr.set_ylabel(TheUnitee,fontsize=10)

# Figure Watermark and Labels
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


#plt.show()
    #plt.savefig(TheFile+'_'+str(tp)+".eps")
    plt.savefig(TheFile+'_CreateSpatials'+".png")
    
    return
       
#**********************************************************************************************    
# PlotMap
def PlotMap(TheFile,TheLatList,TheLonList,TheCandData,TheUnitee,TheNamee,ColourMapChoice,IsLand,TheLatInfo,TheLonInfo,GlobSwitch):
    ''' Create a masked array of trends
        Create a masked array of significant trends
        Plot trends on map, bounded if significant
        Add vertical latitude/trend scatter with average trend overlaid
        Save as eps and png '''

    # Missing data
    mdi = -1e30 # may set up as masked arrays later
      
    # Create the masked array of trends
    MSKTheCandData=ma.masked_where(TheCandData == mdi,TheCandData)
    
    # make 2d arrays of lats and lons
    # nudge -2.5 degrees to make them south/west gridbox corners, not centres
    # add extra row/column to bound the data
    ArrLons,ArrLats = np.meshgrid(TheLonList,TheLatList)
    LngArrLons,LngArrLats = np.meshgrid(np.append(TheLonList-(TheLonInfo[3]/2.),(TheLonList[-1]+(TheLonInfo[3]/2.))),
                                        np.append(TheLatList-(TheLatInfo[3]/2.),(TheLatList[-1]+(TheLatInfo[3]/2.))))
    #LngArrLons,LngArrLats = np.meshgrid(np.append(TheLonList-2.5,180.),np.append(TheLatList-2.5,90.))  
      
    # set up plot
    fig = plt.figure(figsize=(7,5))
    plt.clf()
    ax = plt.axes([0.05,0.05,0.9,0.9],projection=ccrs.Robinson()) # left, bottom, width, height - with cartopy this plots actual axes which we do not want
    #plt.show()
    

    # THIS IS CARTOPY
    # plot map without continents and coastlines
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    #ax = plt.axes(projection=ccrs.Robinson())
    #if (GlobSwitch == 1): # this is a global projection
    #    ax.set_global()
    #    ax.set_extent([-180.,180.,-90.0,90.0],crs=ccrs.PlateCarree())
    #
    #else: # this is regional
    #    ax.set_extent(RegionExtents) # THIS WILL NEED TESTING AT SOME POINT
	
    ax.coastlines()
    #ax.set_boundary # not sure what this does? maybe useful for plotting regions?
    ax.gridlines(draw_labels=False) # probably a way to specify these exactly if we want something different to default
    # This line background fills the land with light grey
    #ax.add_feature(cpf.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k") # may or may not need this
    
    ext = ax.get_extent()
    
    #plt.show()
    # draw paralells and medians, no labels
    #    if (TheLatInfo[1] == len(TheLatList)) & (TheLonInfo[1] == len(TheLonList)):
    #        m.drawparallels(np.arange(-90,90.,30.))
    #        m.drawmeridians(np.arange(-180,180.,60.))
    # END OF CARTOPY STUFF - OR YOU CAN RUN USING python2.7 PlotTrendMap_Florentine_Jul2018.py
    
    # make up a blue to red (reverse) colour map
    cmap = plt.get_cmap(ColourMapChoice[0])
    
#    fullcols=list()
#    [fullcols.append(bluelist[i]) for i in range(len(bluelist))]
#    [fullcols.append(redlist[i]) for i in range(len(redlist))]
    
    # make a list of colours across the range
    cmaplist=[cmap(i) for i in range(cmap.N)]
    # Remove the middle of the range because the colours are very pale # CHECK IF LEAVING THE MIDDLE WHEN IT IS ABOUT CORRELATIONS!
    for loo in range(np.int((cmap.N/2)-10),np.int((cmap.N/2)+10)):
        cmaplist.remove(cmaplist[np.int((cmap.N/2)-10)]) # remove the very pale colours in the middle
    
    # Some variables/datasets require a reversed colour scheme so this tests whether 'flip' has been set and reverses
    if (ColourMapChoice[1] == 'flip'):	# then reverse the colours
        cmaplist.reverse()
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle
    cmap = cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    
    # work out best max and min values for colourbar and try to make them 'nice
    # must be an odd number of steps
    # for less than 0.1 must be 14 or fewer steps
    vmax = np.int(np.ceil(np.max(abs(MSKTheCandData))*10))/10.    
    vmin = -vmax
    nsteps = 9
    if (vmax <= 0.2):
        vmax = 0.2
        vmin = -0.2
    if (vmax <= 0.3):
        #        nsteps=np.int((vmax-vmin)/0.05)+1
        vmax = 0.32
        vmin = -0.32
    elif (vmax <= 0.4):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.06)*0.06
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.06)+1
        vmax = 0.4
        vmin = -0.4
    elif (vmax <= 0.6):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.08)*0.08
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.08)+1
        vmax = 0.6
        vmin = -0.6
    elif (vmax <= 0.8):
        vmax = 0.8
        vmin = -0.8
    elif (vmax <= 1.0):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.1)*0.1
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.1)+1
        vmax = 1.0
        vmin = -1.0
    elif (vmax <= 1.2):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.2)+1
        vmax = 1.2
        vmin = -1.2
    elif (vmax <= 1.6):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.2)*0.2
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.2)+1
        vmax = 1.6
        vmin = -1.6
    elif (vmax <= 2.0):
        #        vmax=np.ceil(np.max(abs(MSKTheCandData))/0.3)*0.3
        #        vmin=-vmax
        #        nsteps=np.int((vmax-vmin)/0.5)+1
        vmax = 2.0
        vmin = -2.0
    elif (vmax <= 3.0):
        vmax = 3.0
        vmin = -3.0
#    pdb.set_trace() # stop here and play
    
    bounds = np.linspace(vmin,vmax,nsteps)
    strbounds = ["%4.1f" % i for i in bounds]
    norm = mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    # USIGN m AGAIN SO WILL NEED CHANGING WHEN WE MOVE TO CARTOPY INSTEAD OF BASEMAP
    grids = ax.pcolormesh(LngArrLons,LngArrLats,MSKTheCandData,transform = ccrs.PlateCarree(),cmap=cmap,norm=norm) # ,latlon='TRUE')

    # adding the colourbar to the plot
    cbax = fig.add_axes([0.05,0.1,0.9,0.03])
    cb = plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=10) 
    plt.figtext(0.5,0.01,'Time Period Average Anomaly ('+TheUnitee+')',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.5,0.95,TheNamee,size=12,ha='center')
    
#    print(CountGoods,CountLargeNegs,CountSmallNegs,CountSmallPos,CountLargePos)
#    pctLNs=round((float(CountLargeNegs)/CountGoods)*100.,1)
#    pctSNs=round((float(CountSmallNegs)/CountGoods)*100.,1)
#    pctVSs=round((float(CountVSmalls)/CountGoods)*100.,1)
#    pctSPs=round((float(CountSmallPos)/CountGoods)*100.,1)
#    pctLPs=round((float(CountLargePos)/CountGoods)*100.,1)
    
#    ax1=plt.axes([0.01,0.01,0.64,0.09],frameon=False) # map only
#    ax1.set_ylim(0,1)
#    ax1.set_xlim(0,1)
#    ax1.axes.get_yaxis().set_visible(False)
#    ax1.axes.get_xaxis().set_visible(False)
#    plt.annotate(str(pctLNs)+"% ratio < -1",xy=(0.07,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.05,0.8,markersize=10,marker='s',color='Firebrick')
#    plt.annotate(str(pctSNs)+"% ratio -1 to -0.05",xy=(0.07,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.05,0.2,markersize=10,marker='s',color='LightCoral')
#    plt.annotate(str(pctVSs)+"% small trend",xy=(0.40,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.38,0.8,markersize=10,marker='s',color='Khaki')
#    plt.annotate(str(pctSPs)+"% ratio 0.05 to 1",xy=(0.73,0.7),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.71,0.8,markersize=10,marker='s',color='LightSkyBlue')
#    plt.annotate(str(pctLPs)+"% ratio > 1",xy=(0.73,0.1),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.plot(0.71,0.2,markersize=10,marker='s',color='DodgerBlue')

    
#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotMap

#**********************************************************************************************    
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
    StartDate = dt.datetime(TheStYr,TheStMon,1,0,0,0)	# January
    TheYear = TheStYr
    TheMonth = TheStMon
    for mm in range(len(DaysArray)):
        if (TheMonth < 12):
            DaysArray[mm] = (dt.datetime(TheYear,TheMonth+1,1,0,0,0)-dt.datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (dt.datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days
            BoundsArray[mm,0] = (dt.datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1] = (dt.datetime(TheYear,TheMonth+1,1,0,0,0)-StartDate).days
        else:
            DaysArray[mm] = (dt.datetime(TheYear+1,1,1,0,0,0)-dt.datetime(TheYear,TheMonth,1,0,0,0)).days/2. + (dt.datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days	
            BoundsArray[mm,0] = (dt.datetime(TheYear,TheMonth,1,0,0,0)-StartDate).days+1
            BoundsArray[mm,1] = (dt.datetime(TheYear+1,1,1,0,0,0)-StartDate).days
        TheMonth=TheMonth+1
        if (TheMonth == 13):
            TheMonth = 1
            TheYear = TheYear+1
	    
    return DaysArray,BoundsArray

#*********************************************************************************************
# WriteOut
def WriteOut(OutFile,OutputData,Latitudes,Longitudes,TheStyr,TheStmn,TheEdyr,TheEdmn,nmons,nlats,nlons,PassingOutDict): # styr,stmon,edyr,edmon to create 3D file
    # Write out
    if (nlats > 0):
        print('Writing out map of time average slices: ',OutFile)
        # We'll need LatBounds and LonBounds for this one
        LatBounds = np.empty((len(Latitudes),2),dtype='float')
        LonBounds = np.empty((len(Longitudes),2),dtype='float')

        LatBounds[:,0] = Latitudes - ((Latitudes[1]-Latitudes[0])/2.)
        LatBounds[:,1] = Latitudes + ((Latitudes[1]-Latitudes[0])/2.)

        LonBounds[:,0] = Longitudes - ((Longitudes[1]-Longitudes[0])/2.)
        LonBounds[:,1] = Longitudes + ((Longitudes[1]-Longitudes[0])/2.)    
    if (nmons > 0):
        print('Writing out time series of space average: ',OutFile)
        # We'll need time elements for this one
        TimPoints,TimBounds = MakeDaysSince(TheStyr,TheStmn,TheEdyr,TheEdmn)
        nTims = len(TimPoints)
    		    

    #pdb.set_trace()

    # No need to convert float data using given scale_factor and add_offset to integers - done within writing program (packV = (V-offset)/scale
    # Not sure what this does to float precision though...

    # Create a new netCDF file - have tried zlib=True,least_significant_digit=3 (and 1) - no difference
    ncfw = nc4.Dataset(OutFile,'w',format='NETCDF4_CLASSIC') # need to try NETCDF4 and also play with compression but test this first
    
    # Set up the dimension names and quantities
    # Go through each dimension and set up the variable and attributes for that dimension if needed
    if (nmons > 0):
        ncfw.createDimension('time',nmons)

        MyVarT = ncfw.createVariable('time','f4',('time',))
        MyVarT.standard_name = 'time'
        MyVarT.long_name = 'time'
        MyVarT.units = 'days since 1979-1-1 00:00:00'
        MyVarT.start_year = str(TheStyr)
        MyVarT.end_year = str(TheEdyr)
        MyVarT[:] = TimPoints

    if (nlats > 0):
        ncfw.createDimension('latitude',nlats)
        ncfw.createDimension('longitude',nlons)

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
    
    if (nmons > 0) & (nlats > 0):
        MyVarD = ncfw.createVariable(PassingOutDict['VarName'],'f4',('time','latitude','longitude',),fill_value = PassingOutDict['mdi'],zlib=True,least_significant_digit=4)
    elif (nmons == 0) & (nlats > 0):
        MyVarD = ncfw.createVariable(PassingOutDict['VarName'],'f4',('latitude','longitude',),fill_value = PassingOutDict['mdi'],zlib=True,least_significant_digit=4)
    if (nmons > 0) & (nlats == 0):
        MyVarD = ncfw.createVariable(PassingOutDict['VarName'],'f4',('time',),fill_value = PassingOutDict['mdi'],zlib=True,least_significant_digit=4)
       
    MyVarD.standard_name = PassingOutDict['StdName']
    MyVarD.long_name = PassingOutDict['LongName']
    MyVarD.units = PassingOutDict['Unit']
    MyVarD.valid_min = PassingOutDict['vmin']
    MyVarD.valid_max = PassingOutDict['vmax']
    MyVarD.missing_value = PassingOutDict['mdi']
    
    # Provide the data to the variable - depending on howmany dimensions there are
    if (nmons > 0) & (nlats > 0):
        MyVarD[:,:,:] = OutputData[:,:,:]
    elif (nmons == 0) & (nlats > 0):
        MyVarD[:,:] = OutputData[:,:]
    if (nmons > 0) & (nlats == 0):
        MyVarD[:] = OutputData[:]
    #pdb.set_trace()    	

    ncfw.close()
    #pdb.set_trace()
    
    return #WriteNetCDF
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************

if __name__ == '__main__':
    main(sys.argv[1:])
