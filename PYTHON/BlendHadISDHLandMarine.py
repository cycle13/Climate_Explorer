# Python 3
# 
# Author: Kate Willett
# Created: 25 February 2019
# Last update: 25 February 2019
# Location: /data/local/hadkw/HADCRUH2/UPDATE2018/PROGS/PYTHON/	
# GitHub: https://github.com/Kate-Willett/Climate_Explorer					
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code reads in the HadISDH-marine data in its HadISDH-land like format and
# outputs a blended with HadISDH-land to produce HadISDH-blend.
#
# Where land and sea present blend uses minimum of 25% of lowest fraction
# OBS, Sampling and Full uncertainties are combined in quadrature - using land fraction as a weight
# n_obs, n_grids, n_stations, n_pseudostations are just stored and not blended
# 
# <references to related published material, e.g. that describes data set>
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# IMPORTS:
#import os
#import datetime as dt
#import numpy as np
#import sys, getopt
#import math
#from math import sin, cos, sqrt, atan2, radians
#import scipy.stats
#import matplotlib
#matplotlib.use('Agg') 
#import matplotlib.pyplot as plt
#import calendar
#import gc
#import netCDF4 as ncdf
#import copy
#import pdb
#
# Kate:#
#from ReadNetCDF import GetGrid
#from ReadNetCDF import GetGrid4
#from GetNiceTimes import MakeDaysSince
# 
# -----------------------
# DATA
# -----------------------
# land data
# /data/local/hadkw/HADCRUH2/UPDATE2018/STATISTICS/GRIDS/HadISDH.land<var>.<version>_FLATgrid<homogtype>5by5_anoms8110_<thenmon><thenyear>_cf.nc
# marine data
# /data/local/hadkw/HADCRUH2/UPDATE2018/STATISTICS/GRIDS/HadISDH.marine<var>.<version>_BClocalSHIP5by5both_anoms8110_<thenmon><thenyear>_cf.nc
#
# Attributes - a list of things to use as global attributes
# /data/local/hadkw/HADCRUH2/UPDATE2018/PROGS/PYTHON/HadISDHBlend_attributes.dat
#
# Land (1)/ Sea (0) mask:
# /data/local/hadkw/HADCRUH2/UPDATE2018/OTHERDATA/HadCRUT.4.3.0.0.land_fraction.nc
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# modul load scitools/experimental-current
# python BlendHadISDHLandMarine.py 
# 
# -----------------------
# OUTPUT
# -----------------------
# /data/local/hadkw/HADCRUH2/UPDATE2018/STATISTICS/GRIDS/
# SHIP
# HadISDH-blendq.<version>_FLATgrid<homogtype>BClocalSHIPboth5by5_anoms8110_<nowmon><nowyear>_cf.nc
# SHIPS AND BUOYS
# HadISDH-blendq.<version>_FLATgrid<homogtype>BClocalboth5by5_anoms8110_<nowmon><nowyear>_cf.nc
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
#
# Version 1 (25 February 2019)
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
# Based original IDL code: Blend_HadISDHLandOcean_APR2016.pro
# -----------------------
#
#*******************************************************************
# IMPORTS:
import os
import datetime as dt
import numpy as np
import sys, getopt
import math
from math import sin, cos, sqrt, atan2, radians
import scipy.stats
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import calendar
import gc
import netCDF4 as ncdf
import copy
import pdb

# Kate:
from ReadNetCDF import GetGrid
from ReadNetCDF import GetGrid4
from GetNiceTimes import MakeDaysSince

# Set up variables
StYr = 1973 
EdYr = 2018
StMn = 1 
EdMn = 12
platform = 'ship' # 'all'

if (platform == 'ship'):
    PT = 'SHIP'
else:
    PT = ''
    
NowMon = 'FEB'
NowYear = '2019'
LandMonYear = 'JAN2019'
MarineMonYear = 'FEB2019'
ClimStart = 1981
ClimEnd = 2010
RefPeriod = str(ClimStart)+' to '+str(ClimEnd)
ClimPeriod = str(ClimStart)[2:4]+str(ClimEnd)[2:4]
LandVersion = '4.1.0.2018f'
MarineVersion = '1.0.0.2018f'
BlendVersion = '1.0.0.2018f'

################################################################################################################
# SUBROUTINES #
##################################################
# BlendIt #
##################################################
def BlendIt(TheLand, TheOcean, ThePctLandMask, TheMDI, IsUnc):
    '''This programs works over each gridbox to blend the land and marine
        data. It uses HadCRUT.4.3.0.0.land_fraction.nc to get a land/ocean coverage
	but if there is both a land and marine box present then there will be a min/max blend of 0.25/0.75.
	
	Apply minimum blend of 25% if both land and marine are present.
	
	This does a cursory uncertainty blend in quadrature of the obs, 
	sampling and full uncertainty.
	INPUT:
	TheLand - time, lat, lon np.array of land data
	TheOcean - time, lat, lon np.array of marine data
	ThePctLandMask - lat, lon np.array of percent land (1) and sea (0)
        TheMDI - scalar value of missing data indicator 
	IsUnc - Boolean True if uncertainty field, False if not
	
	OUTPUT:
	BlendedHadISDH - time, lat, lon np.array of blended land and marine data '''

    # Create empty array to fill
    BlendedHadISDH = np.empty((len(TheLand[:,0,0]),len(TheLand[0,:,0]),len(TheLand[0,0,:])),dtype = float)
    BlendedHadISDH.fill(TheMDI)

    # Fill with ocean where there is only ocean
    BlendedHadISDH[np.where((TheLand == TheMDI) & (TheOcean > TheMDI))] = TheOcean[np.where((TheLand == TheMDI) & (TheOcean > TheMDI))]

    # Fill with land where there is only land
    BlendedHadISDH[np.where((TheLand > TheMDI) & (TheOcean == TheMDI))] = TheLand[np.where((TheLand > TheMDI) & (TheOcean == TheMDI))]

    # Where there are both land and marine values fill with weighted mean or weighted quadrature sum
    if (IsUnc):
    
    # Its uncertainty!
    # obvs np.sqrt(1) is 1 so makes no difference but I'm being clear that this is a mean so should be /np.sqrt(N) which in this case N = 1!!!
        BlendedHadISDH[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))] = np.sqrt((PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]*(TheLand[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]**2)) 
            + ((1 - PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))])*(TheOcean[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]**2)))

    else:
    
    # Its just a mean of normal value - weights should sum to 1 so no need to /N as N = 1!!!
        BlendedHadISDH[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))] = (PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]*TheLand[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]) 
            + ((1 - PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))])*TheOcean[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))])

    print('Test the Blending of BlendedHadISDH from TheLand, TheOcean and PctLandArr')
    return BlendedHadISDH
 
##################################################
# Write_Netcdf_Variable_All #
##################################################
def Write_Netcdf_Variable_All(outfile, var, vlong, vunit, RefPeriod, TheMDI, data_arr):
    '''
    This is basically a tweak of the utils.py version but set up to work here
    and to be consistent with HadISDH netCDF files.
    
    Create the netcdf variable
    :param obj outfile: output file object
    :param string var: variable name
    :param list vlong: long variable name
    :param str vunit: unit of variable
    :param str RefPeriod: period of climatological reference   
    :param np array data_arr: times, lats, lons data to write
    
    '''
    # If there is no time dimension
    if (len(np.shape(data_arr)) == 2):
        nc_var = outfile.createVariable(var, np.dtype('float'), ('latitude','longitude',), zlib = True, fill_value = TheMDI) # with compression
    # If there is a time dimension
    else:
        if (len(data_arr[:,0,0]) > 12):
            nc_var = outfile.createVariable(var, np.dtype('float64'), ('time','latitude','longitude',), zlib = True, fill_value = TheMDI) # with compression
        else:
            nc_var = outfile.createVariable(var, np.dtype('float64'), ('month','latitude','longitude',), zlib = True, fill_value = TheMDI) # with compression

    nc_var.long_name = vlong
    nc_var.units = vunit
    nc_var.missing_value = TheMDI
    nc_var.reference_period = RefPeriod
    
    # We're not using masked arrays here - hope that's not a problem
    nc_var.valid_min = np.min(data_arr[np.where(data_arr != TheMDI)]) 
    nc_var.valid_max = np.max(data_arr[np.where(data_arr != TheMDI)]) 
    nc_var[:] = data_arr
        
    return # write_netcdf_variable_all
    
###################################################################
# Write_NetCDF_All #
###################
def Write_NetCDF_All(filename, land_var_list, marine_var_list, blend_var_list,
                     land_data_vals, marine_data_vals, blend_data_vals, 
		     land_var_long, marine_var_long, blend_var_long,
		     land_var_standard, marine_var_standard, blend_var_standard,
		     lats, lons, start_year, end_year, RefPeriod, unit, TheMDI): 

    '''
    This is basically a copy of utils.py version but tweaked to work here and now be consistent with HadISDH.land files
    
    2 sigma uncertainties are passed in here and written out!!!
    
    Write the relevant fields out to a netCDF file.
    
    :param str filename: output filename  
    :param list of land elements land_var_list: the land only elements
    :param list of marine elements marine_var_list: the marine only elements
    :param list of blend elements blend_var_list: the blend only elements
    :param list of np arrays land_data_vals: the land mean and actual station counts data array for [times?, lats, lons]
    :param list of np arrays marine_data_vals: the marine n_grids, n_obs, mean and actual pseudo station counts data array for [times?, lats, lons]
    :param list of np arrays blend_data_vals: the blended values and uncertainties data arrays for [times, lats, lons]
    :param list of land elements land_var_long: the land only elements
    :param list of marine elements marine_var_long: the marine only elements
    :param list of blend elements blend_var_long: the blend only elements
    :param list of land elements land_var_standard: the land only elements
    :param list of marine elements marine_var_standard: the marine only elements
    :param list of blend elements blend_var_standard: the blend only elements
    :param array lats: the latitudes
    :param array lons: the longitudes
    :param int start_year: the start year
    :param int end_year: the end year
    :param string RefPeriod: the climatological reference period.
    :param string unit: the variable units
    :param float TheMDI: the missing data indicator

    '''

    # remove file
    if os.path.exists(filename):
        os.remove(filename)

    outfile = ncdf.Dataset(filename,'w', format='NETCDF4')

    # Set up dimensions
    # Get times in terms of days since 1973-1-1 00:00:00
    DaysSince = MakeDaysSince(start_year,1,end_year,12,'month')
#    print('Test days since')
#    pdb.set_trace()
    time_dim = outfile.createDimension('time',len(DaysSince))
    month_dim = outfile.createDimension('month',12)
    lat_dim = outfile.createDimension('latitude',len(lats)) # as TRC of box edges given, size = # box centres to be written
    lon_dim = outfile.createDimension('longitude',len(lons))
    
    #***********
    # set up basic variables linked to dimensions
    # make time variable
    nc_var = outfile.createVariable('time', np.dtype('int'), ('time'), zlib = True) # with compression!!!
    nc_var.long_name = "time"
    nc_var.units = "days since 1973-1-1 00:00:00"
    nc_var.standard_name = "time"
    nc_var.axis = "T"
    nc_var.calendar = "gregorian"
    nc_var.start_year = str(start_year)
    nc_var.end_year = str(end_year)
    nc_var.start_month = '1'
    nc_var.end_month = '12'
    nc_var[:] = DaysSince

    # make month variable
    nc_var = outfile.createVariable('month', np.dtype('int'), ('month'), zlib = True) # with compression!!!
    nc_var.long_name = "month of year"
    nc_var.units = "month"
    nc_var.standard_name = "time in months"
    #nc_var.start_year = str(StartYear)
    #nc_var.end_year = str(EndYear)
    nc_var.start_month = '1'
    nc_var.end_month = '12'
    nc_var[:] = range(12)
    
    # make latitude variable
    nc_var = outfile.createVariable('latitude', np.dtype('float32'), ('latitude'), zlib = True) # with compression!!!
    nc_var.long_name = "latitude"
    nc_var.units = "degrees_north"
    nc_var.standard_name = "latitude"
    nc_var.point_spacing = "even"
    nc_var.axis = "X"
    nc_var[:] = lats
    
    # make longitude variable
    nc_var = outfile.createVariable('longitude', np.dtype('float32'), ('longitude'), zlib = True) # with compression!!!
    nc_var.long_name = "longitude"
    nc_var.units = "degrees_east"
    nc_var.standard_name = "longitude"
    nc_var.point_spacing = "even"
    nc_var.axis = "X"
    nc_var[:] = lons

    #***********
    # create variables for blended
    for b in range(len(blend_var_list)):
    
        Write_Netcdf_Variable_All(outfile, blend_var_list[b], blend_var_long[b], unit, RefPeriod, TheMDI, blend_data_vals[b])

    # create variables for land
    for l in range(len(land_var_list)):
    
        Write_Netcdf_Variable_All(outfile, land_var_list[l], land_var_long[l], 'standard', RefPeriod, -1, land_data_vals[l])

    # create variables for land
    for m in range(len(marine_var_list)):
    
        Write_Netcdf_Variable_All(outfile, marine_var_list[m], marine_var_long[m], 'standard', RefPeriod, -1, marine_data_vals[m])

    # Global Attributes 
    
    # Read these from file
    attr_file = os.path.join(os.getcwd(), "HadISDHBlend_attributes.dat")

    try:
        with open(attr_file,'r') as infile:        
            lines = infile.readlines()
        
    except IOError:
        print("Attributes file not found at " + attr_file)
        raise IOError
    
    attributes = {}
    
    for line in lines:
        split_line = line.split()
        
        attributes[split_line[0]] = " ".join(split_line[1:])    
    
    # Set the attributes
    for attr in attributes:
        
        outfile.__setattr__(attr, attributes[attr])
 
    outfile.file_created = dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d, %H:%M")
    outfile.Conventions = 'CF-1.5' 
    outfile.Metadata_Conventions = 'Unidata Dataset Discovery v1.0,CF Discrete Sampling Geometries Conventions'
    outfile.featureType = 'gridded'
    
    outfile.close()

    return # Write_NetCDF_All
    
################################################################################################################
# MAIN #
################################################################################################################
def main():
    # Set up this run files, directories and dates/clims: years, months, ship or all
    # Read in
    # Variables - short names
    # Quantities to be blended
    var_loop_lower = ['t','td','q','e','rh','tw','dpd']
    homogtype = ['IDPHA','PHADPD','IDPHA','IDPHA','IDPHA','IDPHA','PHA']
    var_loop = ['T','Td','q','e','RH','Tw','DPD']
    var_loop_full = ['Air Temperature','Dew Point Temperature','Specific Humidity','Vapor Pressure','Relative Humidity','Wet Bulb Temperature','Dew Point Pepression'] # This is the ReadInfo
    
    # INput land variables
    InLVarList = ['mean_n_stations','actual_n_stations','_anoms','_abs','_clims','_std','_obserr','_samplingerr','_combinederr']
    # start point and end point of Land vars with variable prefixes
    InLPointers = [2,9]
    
    # INput marine variables
    InMVarList = ['_n_grids','_n_obs','_clims_n_grids','_clims_n_obs','_clim_std_n_grids','_clim_std_n_obs','_uSAMP_n_grids','_usbarSQ_n_grids',
                  '_anoms','_abs','_clims','_clim_std','_anoms_uOBS','_abs_uOBS','_anoms_uSAMP','_abs_uSAMP','_anoms_uFULL','_abs_uFULL']
    # start point and end point of Marine vars with variable prefixes
    InMPointers = [0,18]
        
    # Shared output quantities to be blended, land only and marine only
    OutBVarList = ['_anoms','_abs','_clims','_climstd','_obserr','_abs_obserr','_samplingerr','_abs_samplingerr','_combinederr','_abs_combinederr']
    OutLVarList = ['Land_mean_n_stations','Land_actual_n_stations']
    OutMVarList = ['Marine_n_grids','Marine_n_obs','Marine_mean_pseudo_stations','Marine_actual_pseudo_stations',
                   'Marine_clim_n_obs','Marine_clim_n_grids','Marine_climstd_n_obs','Marine_climstd_n_grids']
    
    # Output variables list elements ????? have the var_loop_lower[v] prefixed
    OutBVarLong = ['Air Temperature','Dew Point Temperature','Specific Humidity','Vapor Pressure','Relative Humidity','Wet Bulb Temperature','Dew Point Pepression'] # This is the ReadInfo
    OutLVarLong = ['Air Temperature','Dew Point Temperature','Specific Humidity','Vapor Pressure','Relative Humidity','Wet Bulb Temperature','Dew Point Pepression'] # This is the ReadInfo
    OutMVarLong = ['Air Temperature','Dew Point Temperature','Specific Humidity','Vapor Pressure','Relative Humidity','Wet Bulb Temperature','Dew Point Pepression'] # This is the ReadInfo
    
    OutBVarStandard = ['air temperarature','dew point temperature','specific humidity','vapor pressure','relative humidity','wet bulb temperature','dew point depression'] # This is the ReadInfo
    OutLVarStandard = ['air temperarature','dew point temperature','specific humidity','vapor pressure','relative humidity','wet bulb temperature','dew point depression'] # This is the ReadInfo
    OutMVarStandard = ['air temperarature','dew point temperature','specific humidity','vapor pressure','relative humidity','wet bulb temperature','dew point depression'] # This is the ReadInfo

    units_loop = ['degrees C','degrees C','g/kg','hPa','%rh','degrees C','degrees C','standard']
        
    # Missing Data Indicator
    MDI = -1e30
    
    # Input and Output directory:
    WorkingDir = '/data/local/hadkw/HADCRUH2/UPDATE'+str(EdYr)
    DataDir = '/STATISTICS/GRIDS/'
    OtherDataDir = '/OTHERDATA/'
    
    # Input Files
    InFilLandBit1 = WorkingDir+DataDir+'HadISDH.land' # append to var+InFilLandBit2+homogtype+InFilBit3
    InFilLandBit2 = '.'+LandVersion+'_FLATgrid'
    InFilLandBit3 = '5by5_anoms8110_'+LandMonYear+'_cf.nc'

    InFilMarineBit1 = WorkingDir+DataDir+'HadISDH.marine' # append to var+InFilMarineBit2
    InFilMarineBit2 = '.'+MarineVersion+'_BClocal'+PT+'5by5both_anoms8110_'+MarineMonYear+'_cf.nc'

    # Output Files
    OutFilBit1 = WorkingDir+DataDir+'HadISDH.blend' # append to var+OutFilBit2+homogtype+InFilBit3
    OutFilBit2 = '.'+BlendVersion+'_FLATgrid'
    OutFilBit3 = 'BClocal'+PT+'both5by5_anoms8110_'+NowMon+NowYear+'_cf.nc'
        
    # Set up necessary dates - dates for output are just counts of months from 0 to 54?...
    Ntims = ((EdYr + 1) - StYr) * 12
    TimesArr = np.arange(Ntims) # months since January 1973
    
    # Read in the pct land data
    Filee = WorkingDir+OtherDir+'HadCRUT.4.3.0.0.land_fraction.nc' 
    LatInfo = ['latitude']
    LonInfo = ['longitude']
    PctLand, LatList, LonList = GetGrid(Filee,['land_area_fraction'],LatInfo,LonInfo)

    # Make PctLand a min and max of 0.25 or 0.75 respectively
    # We only use it when both land and marine are present so it doesn't matter in land only or marine only cases    
    PctLand[np.where(PctLand < 0.25] = 0.25
    PctLand[np.where(PctLand > 0.75] = 0.75
    # Now make PctLand have a time dimension the same length as the data
    3DPctLand = np.repeat(PctLand[np.newaxis,:,:],Ntims,axis=0)
    
#########
    # Loop through each variable to create the blended file
        
    # Which variable to loop through?
    for v,var in enumerate(var_loop):

        # Loop through quantities to be blended and build list
	BlendedList = []
	LandList = []
	LandDataList = []
	MarineList = []
	MarineDataList = []

        # For this var read in all land elements
        Filee = WorkingDir+InFilLandBit1+var+InFilLandBit2+homogtype[v]+InFilLandBit3
        LatInfo = ['latitude']
        LonInfo = ['longitude']
	TmpVarList = []
	TmpVarList = InLVarList
	TmpVarList[InLPointers[0]:InLPointers[1]] = [var+'_'+i for i in TmpVarList[InLPointers[0]:InLPointers[1]]]
        TmpDataList, LatList, LonList = GetGrid4(Filee,TmpVarList,LatInfo,LonInfo)
        # This comes out as:
        # TmpDataList: a list of np arrays (times, lats{87.5N to 87.5S], lons[-177.5W to 177.5E])
        # LatList: an NLats np array of lats centres (87.5N to 87.5S)
        # LonList: an NLons np array of lons centres (87.5N to 87.5S)
    
        # Get lat and lon counts
        NLats = len(LatList)
        NLons = len(LonList)
	
	# change MDI for counts to -1
	TmpDataList[0][np.where(TmpDataList[0] <= 0)] = -1
	TmpDataList[1][np.where(TmpDataList[1] <= 0)] = -1
    
        # Fill in lists
        LandList.append(TmpDataList[0]) # mean_n_stations 	
        LandList.append(TmpDataList[1]) # actual_n_stations
	LandDataList.append(TmpDataList[2]) # _anoms	
	LandDataList.append(TmpDataList[3]) # _abs	
	LandDataList.append(TmpDataList[4]) # _clims	
	LandDataList.append(TmpDataList[5]) # _climstd	
	LandDataList.append(TmpDataList[6]) # _obserr	
	LandDataList.append(TmpDataList[6]) # _abs_obserr	
	LandDataList.append(TmpDataList[7]) # _samplingerr	
	LandDataList.append(TmpDataList[7]) # _abs_samplingerr	
	LandDataList.append(TmpDataList[8]) # _combinederr	
	LandDataList.append(TmpDataList[8]) # _abs_combinederr	
	
        # For this var read in all marine elements
        Filee = WorkingDir+InFilMarineBit1+var+InFilMarineBit2
        LatInfo = ['latitude']
        LonInfo = ['longitude']
	TmpVarList = []
	TmpVarList = InMVarList
	TmpVarList[InMPointers[0]:InMPointers[1]] = [var+'_'+i for i in TmpVarList[InMPointers[0]:InMPointers[1]]]
        TmpDataList, LatList, LonList = GetGrid4(Filee,TmpVarList,LatInfo,LonInfo)
        # This comes out as:
        # TmpDataList: a list of np arrays (times, lats{87.5N to 87.5S], lons[-177.5W to 177.5E])
        # LatList: an NLats np array of lats centres (87.5N to 87.5S)
        # LonList: an NLons np array of lons centres (87.5N to 87.5S)
        
        # Fill in lists
        MarineList.append(TmpDataList[0]) # _n_grids 	
        MarineList.append(TmpDataList[1]) # _n_obs
        MarineList.append(TmpDataList[6]) # _mean_pseudo_stations
        MarineList.append(TmpDataList[7]) # _actual_pseudo_stations
        MarineList.append(TmpDataList[2]) # _clim_n_grids
        MarineList.append(TmpDataList[3]) # _clim_n_obs
        MarineList.append(TmpDataList[4]) # _climstd_n_grids
        MarineList.append(TmpDataList[5]) # _climstd_n_obs
	MarineDataList.append(TmpDataList[8]) # _anoms	
	MarineDataList.append(TmpDataList[9]) # _abs	
	MarineDataList.append(TmpDataList[10]) # _clims	
	MarineDataList.append(TmpDataList[11]) # _climstd	
	MarineDataList.append(TmpDataList[12]) # _obserr	
	MarineDataList.append(TmpDataList[13]) # _abs_obserr	
	MarineDataList.append(TmpDataList[14]) # _samplingerr	
	MarineDataList.append(TmpDataList[15]) # _abs_samplingerr	
	MarineDataList.append(TmpDataList[16]) # _combinederr	
	MarineDataList.append(TmpDataList[17]) # _abs_combinederr	
	
	# Now loop through the Blended quantities and blend
        for b,bvar in enumerate(OutBVarList):

	    # Is this an uncertainty quantity?
            if (b > 3): # Set IsUnc = True
	        IsUnc = True
            else:
                IsUnc = False
		
	    BlendedField = BlendIt(LandDataList[b],MarineDataList[b],3DPctLand,MDI,IsUnc)	    
	    BlendedList.append(BlendedField)
	    
##############
        # Write out combined file - uncertainties are all 2 sigma!!!.

        Write_NetCDF_All(OutFilBit1+var+OutFilBit2+homogtype[v]+OutFilBit3,
                         OutLVarList,OutMVarList,[var+'_'+i for i in OutBVarList],
			 LandList,MarineList,BlendList,
			 OutLVarLong,OutMVarLong,[var_loop_full[v]+' '+i for i in OutBVarLong],
			 OutLVarStandard,OutMVarStandard,[var_loop_full[v]+' '+i for i in OutBVarStandard],
                         LatList,LonList,StYr,EdYr,RefPeriod,units_loop[v],MDI)

#########

    print('And we are done!')

if __name__ == '__main__':
    
    main(sys.argv[1:])


pro Blend_HadISDHLandOcean_APR2016

# Modified to read in from already created NOCS grids:
# extract_NOCSvars_FEB2014.pro
#	> UPDATE2014/OTHERDATA/NOCSv2.0_ocean*

#*******************************************************************
# EDITABLES:
mdi = -1e30

# What year?
styr = 1973 
edyr = 2016

# What climatology period? (land is generally 1976-2005 until I sort out a new version)
clstyr = 1981
cledyr = 2010
clst   = clstyr - styr
cled   = cledyr - styr
climchoice = strmid(strcompress(clstyr,/remove_all),2,2)+strmid(strcompress(cledyr,/remove_all),2,2)

# What versions?
LandVersion   = '3.0.0.2016p'
MarineVersion = '1.0.0.2016p'
BlendVersion  = '1.0.0.2016p'

# What Types?
LandType   = 'standard' # this is actually 'FLATgrid*5by5' where IDPHA, PHADPD and PHA are used
#MarineType = 'OBSclim2BClocal' # this is the QC'd and bias corrected version
#BlendType  = 'FULL' # FULL=QC and bias corrected marine, homogenised land
MarineType = 'OBSclim2BClocalship' # this is the QC'd and bias corrected version
BlendType  = 'FULLship' # FULL=QC and bias corrected marine, homogenised land

# What working versions?
thenmon  = 'JAN'
thenyear = '2017'
nowmon   = 'JAN'
nowyear  = '2017'

# Work with anomalies or absolutes?
anomyes = 0	#0=anomalies, 1=absolutes

#*******************************************************************
# FILEPATHS AND VARIABLES 

updatedir  = 'UPDATE20'+strmid(strcompress(edyr,/remove_all),2,2)
workingdir = '/data/local/hadkw/HADCRUH2/'+updatedir

inMardir   = workingdir+'/OTHERDATA/'
inLanddir  = workingdir+'/STATISTICS/GRIDS/'
outdir     = workingdir+'/STATISTICS/GRIDS/'

inLandSeaMask = workingdir+'/OTHERDATA/new_coverpercentjul08.nc'

inMarfilANM   = inMardir+MarineType+'_5x5_monthly_renorm19812010_anomalies_from_daily_both_relax.nc'
inMarfilABS   = inMardir+MarineType+'_5x5_monthly_from_daily_both_relax.nc'
inMarfilCLM   = inMardir+MarineType+'_5x5_monthly_climatology_from_daily_both_relax.nc'
inMarfilSDV   = inMardir+MarineType+'_5x5_monthly_stdev_from_daily_both_relax.nc'

IF (LandType EQ 'standard') THEN BEGIN
  inLandq   = inLanddir+'HadISDH.landq.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandRH  = inLanddir+'HadISDH.landRH.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLande   = inLanddir+'HadISDH.lande.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandTw  = inLanddir+'HadISDH.landTw.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandT   = inLanddir+'HadISDH.landT.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandTd  = inLanddir+'HadISDH.landTd.'+LandVersion+'_FLATgridPHADPD5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandDPD = inLanddir+'HadISDH.landDPD.'+LandVersion+'_FLATgridPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
ENDIF

outMarfilq   = outdir+'HadISDH.marineq.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfile   = outdir+'HadISDH.marinee.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilTd  = outdir+'HadISDH.marineTd.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilT   = outdir+'HadISDH.marineT.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilTw  = outdir+'HadISDH.marineTw.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilRH  = outdir+'HadISDH.marineRH.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilDPD = outdir+'HadISDH.marineDPD.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'

outBlendfilq   = outdir+'HadISDH.blendq.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfile   = outdir+'HadISDH.blende.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilTd  = outdir+'HadISDH.blendTd.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilT   = outdir+'HadISDH.blendT.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilTw  = outdir+'HadISDH.blendTw.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilRH  = outdir+'HadISDH.blendRH.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilDPD = outdir+'HadISDH.blendDPD.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'

nyrs     = (edyr+1)-styr
nmons    = nyrs*12
int_mons = indgen(nmons)

nlons = 72
nlats = 36
lons  = (findgen(nlons)*5.)-177.5
lats  = (findgen(nlats)*(-5.))+87.5
lats  = REVERSE(lats)

#landarr   = make_array(nlons,nlats,nmons,/float,value=mdi)
#marinearr = make_array(nlons,nlats,nmons,/float,value=mdi)
blendarrANM  = make_array(nlons,nlats,nmons,/float,value=mdi)
blendarrABS  = make_array(nlons,nlats,nmons,/float,value=mdi)

#************************************************************************************
# read in land sea mask
inn   = NCDF_OPEN(inLandSeaMask)
varid = NCDF_VARID(inn,'pct_land')
NCDF_VARGET,inn,varid,pctl
NCDF_CLOSE,inn

#pctl = REVERSE(pctl,2)

# blend using coastal gridboxes as 25-75% weighting
# adjust pctl such that when value is not 0 or 100 it is at least 25 or max 75
# thus land and sea in coastal boxes will always be represented
lows = WHERE(pctl LT 25. AND pctl NE 0.,count)
IF (count GT 0) THEN pctl(lows) = 25.
highs = WHERE(pctl GT 75. AND pctl NE 100.,count)
IF (count GT 0) THEN pctl(highs) = 75.
#stop

# Loop through each variable
for loo = 0,6 DO BEGIN

  # read in marine ANOMALIES
  inn = NCDF_OPEN(inMarFilANM) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity_anomalies')		
    1: varid = NCDF_VARID(inn,'vapor_pressure_anomalies')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature_anomalies')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature_anomalies')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature_anomalies')		
    5: varid = NCDF_VARID(inn,'relative_humidity_anomalies')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature_anomalies')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrANM
  marinearrANM = reverse(marinearrANM,2) # CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrANM LT -100,count)
  if (count GT 0) THEN marinearrANM(bads) = mdi
  NCDF_CLOSE,inn
  # read in marine ACTUALS
  inn = NCDF_OPEN(inMarFilABS) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrABS
  marinearrABS = reverse(marinearrABS,2) # CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrABS LT -100,count)
  if (count GT 0) THEN marinearrABS(bads) = mdi
  NCDF_CLOSE,inn
  # read in marine CLIMATOLOGY (1981-2010)
  inn = NCDF_OPEN(inMarFilCLM) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrCLM
  marinearrCLM = reverse(marinearrCLM,2) # CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrCLM LT -100,count)
  if (count GT 0) THEN marinearrCLM(bads) = mdi
  NCDF_CLOSE,inn
  # read in marine STDEVS
  inn = NCDF_OPEN(inMarFilSDV) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrSDV
  marinearrSDV = reverse(marinearrSDV,2) # CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrSDV LT -100,count)
  if (count GT 0) THEN marinearrSDV(bads) = mdi
  NCDF_CLOSE,inn

  #stop

  # read in HadISDH data and renormalise to 1981-2010
  # could also mask by where uncertainties are larger than the standard deviation
  CASE loo OF
    0: BEGIN
      inn=NCDF_OPEN(inLandq)
      Avarid=NCDF_VARID(inn,'q_anoms') 
      ABvarid=NCDF_VARID(inn,'q_abs') 
    END
    1: BEGIN
      inn=NCDF_OPEN(inLande)
      Avarid=NCDF_VARID(inn,'e_anoms') 
      ABvarid=NCDF_VARID(inn,'e_abs') 
    END
    2: BEGIN
      inn=NCDF_OPEN(inLandTd)
      Avarid=NCDF_VARID(inn,'td_anoms') 
      ABvarid=NCDF_VARID(inn,'td_abs') 
    END
    3: BEGIN
      inn=NCDF_OPEN(inLandT)
      Avarid=NCDF_VARID(inn,'t_anoms') 
      ABvarid=NCDF_VARID(inn,'t_abs') 
    END
    4: BEGIN
      inn=NCDF_OPEN(inLandTw)
      Avarid=NCDF_VARID(inn,'tw_anoms') 
      ABvarid=NCDF_VARID(inn,'tw_abs') 
    END
    5: BEGIN
      inn=NCDF_OPEN(inLandRH)
      Avarid=NCDF_VARID(inn,'rh_anoms') 
      ABvarid=NCDF_VARID(inn,'rh_abs') 
    END
    6: BEGIN
      inn=NCDF_OPEN(inLandDPD)
      Avarid=NCDF_VARID(inn,'dpd_anoms') 
      ABvarid=NCDF_VARID(inn,'dpd_abs') 
    END
    
  ENDCASE
  #varid=NCDF_VARID(inn,'q_abs') 
  #NCDF_VARGET,inn,Avarid,Alandarr
  NCDF_VARGET,inn,Avarid,landarrANM
  NCDF_VARGET,inn,ABvarid,landarrABS
  bads = where(landarrANM LT -100.,count)
  if (count GT 0) THEN landarrANM(bads) = mdi
  bads = where(landarrABS LT -100.,count)
  if (count GT 0) THEN landarrABS(bads) = mdi
  timvarid = NCDF_VARID(inn,'time') 
  NCDF_VARGET,inn,timvarid,landtims
  NCDF_CLOSE,inn
  #stop

  # renormalise the land array anomalies to chosen clim if no 7605 or marine array if not 8110?
  FOR i = 0,nlons-1 DO BEGIN
    FOR j = 0,nlats-1 DO BEGIN
      # SORT OUT WHETHER TO RENORM LAND OR MARINE!!!
      IF (climchoice EQ '7605') THEN subarr = marinearrANM(i,j,*) ELSE subarr = landarrANM(i,j,*)
      gots = WHERE(subarr NE mdi,count)
      IF (count GT nmons/2) THEN BEGIN	# no point if there isn't much data
        monarr = REFORM(subarr,12,nyrs)
        FOR mm=0,11 DO BEGIN
          submon = monarr(mm,*)
	  gots = WHERE(submon NE mdi,count)
	  IF (count GT 20) THEN BEGIN
	    clims = submon(clst:cled)
	    gotclims = WHERE(clims NE mdi,count)
	    IF (count GT 15) THEN submon(gots) = submon(gots)-MEAN(clims(gotclims))
	  ENDIF ELSE submon(*,*) = mdi
	  monarr(mm,*) = submon
        ENDFOR
      IF (climchoice EQ '7605') THEN marinearrANM(i,j,*) = REFORM(monarr,nmons) ELSE landarrANM(i,j,*) = REFORM(monarr,nmons)
      ENDIF ELSE BEGIN
        IF (climchoice EQ '7605') THEN marinearrANM(i,j,*) = mdi ELSE landarrANM(i,j,*) = mdi
      ENDELSE
    ENDFOR
  ENDFOR
  #stop

  #Now do the blending
  FOR i=0,nlons-1 DO BEGIN
    FOR j=0,nlats-1 DO BEGIN
      FOR nn=0,nmons-1 DO BEGIN
        IF (pctl(i,j) EQ 100.) THEN BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*0.25)+(landarrANM(i,j,nn)*0.75) 
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*0.25)+(landarrABS(i,j,nn)*0.75) 
        ENDIF ELSE IF (pctl(i,j) EQ 0.) THEN BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*0.75)+(landarrANM(i,j,nn)*0.25) 
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*0.75)+(landarrABS(i,j,nn)*0.25) 
        ENDIF ELSE BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*((100.-pctl(i,j))/100.))+(landarrANM(i,j,nn)*(pctl(i,j)/100.))  
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*((100.-pctl(i,j))/100.))+(landarrABS(i,j,nn)*(pctl(i,j)/100.))  
        ENDELSE
      ENDFOR  
    ENDFOR
  ENDFOR

  # output gridded Blended product
  CASE loo OF
    0: wilma = NCDF_CREATE(outBlendfilq,/clobber)
    1: wilma = NCDF_CREATE(outBlendfile,/clobber)
    2: wilma = NCDF_CREATE(outBlendfilTd,/clobber)
    3: wilma = NCDF_CREATE(outBlendfilT,/clobber)
    4: wilma = NCDF_CREATE(outBlendfilTw,/clobber)
    5: wilma = NCDF_CREATE(outBlendfilRH,/clobber)
    6: wilma = NCDF_CREATE(outBlendfilDPD,/clobber)
  ENDCASE
  
  tid   = NCDF_DIMDEF(wilma,'time',nmons)
  clmid = NCDF_DIMDEF(wilma,'month',12)
  latid = NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid = NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar = NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar  = NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar  = NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  CASE loo OF
    0:   anomvar = NCDF_VARDEF(wilma,'q_anoms',[lonid,latid,tid],/FLOAT)
    1:   anomvar = NCDF_VARDEF(wilma,'e_anoms',[lonid,latid,tid],/FLOAT)
    2:   anomvar = NCDF_VARDEF(wilma,'td_anoms',[lonid,latid,tid],/FLOAT)
    3:   anomvar = NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
    4:   anomvar = NCDF_VARDEF(wilma,'tw_anoms',[lonid,latid,tid],/FLOAT)
    5:   anomvar = NCDF_VARDEF(wilma,'rh_anoms',[lonid,latid,tid],/FLOAT)
    6:   anomvar = NCDF_VARDEF(wilma,'dpd_anoms',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   absvar = NCDF_VARDEF(wilma,'q_abs',[lonid,latid,tid],/FLOAT)
    1:   absvar = NCDF_VARDEF(wilma,'e_abs',[lonid,latid,tid],/FLOAT)
    2:   absvar = NCDF_VARDEF(wilma,'td_abs',[lonid,latid,tid],/FLOAT)
    3:   absvar = NCDF_VARDEF(wilma,'t_abs',[lonid,latid,tid],/FLOAT)
    4:   absvar = NCDF_VARDEF(wilma,'tw_abs',[lonid,latid,tid],/FLOAT)
    5:   absvar = NCDF_VARDEF(wilma,'rh_abs',[lonid,latid,tid],/FLOAT)
    6:   absvar = NCDF_VARDEF(wilma,'dpd_abs',[lonid,latid,tid],/FLOAT)
  ENDCASE
  
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','months beginning Jan 1973'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','valid_min',0.
  NCDF_ATTPUT,wilma,'latitude','long_name','Latitude'
  NCDF_ATTPUT,wilma,'latitude','units','Degrees'
  NCDF_ATTPUT,wilma,'latitude','valid_min',-90.
  NCDF_ATTPUT,wilma,'latitude','valid_max',90.
  NCDF_ATTPUT,wilma,'longitude','long_name','Longitude'
  NCDF_ATTPUT,wilma,'longitude','units','Degrees'
  NCDF_ATTPUT,wilma,'longitude','valid_min',-180.
  NCDF_ATTPUT,wilma,'longitude','valid_max',180.

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_anoms','long_name','Blended Specific Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_anoms','units','g/kg'
      NCDF_ATTPUT,wilma,'q_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'q_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_anoms','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_anoms','long_name','Blended Vapour Pressure monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_anoms','units','hPa'
      NCDF_ATTPUT,wilma,'e_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'e_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_anoms','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_anoms','long_name','Blended Dew Point Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'td_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'td_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_anoms','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_anoms','long_name','Blended Air Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'t_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'t_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_anoms','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_anoms','long_name','Blended Wet Bulb Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'tw_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'tw_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_anoms','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_anoms','long_name','Blended Relative Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_anoms','units','%rh'
      NCDF_ATTPUT,wilma,'rh_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'rh_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_anoms','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_anoms','long_name','Blended Dew Point Depression monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_anoms','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_abs','long_name','Blended Specific Humidity monthly mean'
      NCDF_ATTPUT,wilma,'q_abs','units','g/kg'
      NCDF_ATTPUT,wilma,'q_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'q_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_abs','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_abs','long_name','Blended Vapour Pressure monthly mean'
      NCDF_ATTPUT,wilma,'e_abs','units','hPa'
      NCDF_ATTPUT,wilma,'e_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'e_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_abs','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_abs','long_name','Blended Dew Point Temperature monthly mean'
      NCDF_ATTPUT,wilma,'td_abs','units','deg C'
      NCDF_ATTPUT,wilma,'td_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'td_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_abs','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_abs','long_name','Blended Air Temperature monthly mean'
      NCDF_ATTPUT,wilma,'t_abs','units','deg C'
      NCDF_ATTPUT,wilma,'t_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'t_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_abs','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_abs','long_name','Blended Wet Bulb Temperature monthly mean'
      NCDF_ATTPUT,wilma,'tw_abs','units','deg C'
      NCDF_ATTPUT,wilma,'tw_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'tw_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_abs','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_abs','long_name','Blended Relative Humidity monthly mean'
      NCDF_ATTPUT,wilma,'rh_abs','units','%rh'
      NCDF_ATTPUT,wilma,'rh_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'rh_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_abs','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_abs','long_name','Blended Dew Point Depression monthly mean'
      NCDF_ATTPUT,wilma,'dpd_abs','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'dpd_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_abs','missing_value',mdi
    END

  ENDCASE

  current_time=SYSTIME()
  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"HadISDH monthly mean blended land and marine surface climate monitoring product from 1973 onwards. "+$
                                         "Quality control, homogenisation (land), bias correction (marine), uncertainty estimation, averaging over gridboxes (no smoothing "+$
					 "or interpolation). This combines HadISDH.land."+LandVersion+" and HadISDH.marine."+MarineVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"HadISDH monthly mean blended land and marine surface climate monitoring product from 1973 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"Met Office Hadley Centre (UK), National Climatic Data Centre (USA), Climatic Research Unit (UK), "+$
                                         "National Physical Laboratory (UK), Maynooth University (Ireland), National Oceanography Centre Southampton (NOCS)"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"HadISD.2.1.0.2016p (Dunn et al., 2016), www.metoffice.gov.uk/hadobs/hadisd/ and ICOADS.3.0.0 and ICOADS.3.0.1 (Freemand et al., 2016), icoads.noaa.gov"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"Willett, K. M., Dunn, R. J. H., Thorne, P. W., Bell, S., de Podesta, M., Parker, D. E., Jones, "+$
                                       "P. D., and Williams Jr., C. N.: HadISDH land surface multi-variable humidity and temperature "+$
				       "record for climate monitoring, Clim. Past, 10, 1983-2006, doi:10.5194/cp-10-1983-2014, 2014. and Willett et al. in prep."
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"HadISDH.blend."+BlendVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF

  NCDF_VARPUT, wilma,timesvar, landtims
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,anomvar,blendarrANM
  NCDF_VARPUT, wilma,absvar,blendarrABS

  NCDF_CLOSE,wilma

  # output gridded Marine product
  CASE loo OF
    0: wilma = NCDF_CREATE(outMarfilq,/clobber)
    1: wilma = NCDF_CREATE(outMarfile,/clobber)
    2: wilma = NCDF_CREATE(outMarfilTd,/clobber)
    3: wilma = NCDF_CREATE(outMarfilT,/clobber)
    4: wilma = NCDF_CREATE(outMarfilTw,/clobber)
    5: wilma = NCDF_CREATE(outMarfilRH,/clobber)
    6: wilma = NCDF_CREATE(outMarfilDPD,/clobber)
  ENDCASE
  
  tid   = NCDF_DIMDEF(wilma,'time',nmons)
  clmid = NCDF_DIMDEF(wilma,'month',12)
  latid = NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid = NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar = NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar  = NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar  = NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  CASE loo OF
    0:   anomvar = NCDF_VARDEF(wilma,'q_anoms',[lonid,latid,tid],/FLOAT)
    1:   anomvar = NCDF_VARDEF(wilma,'e_anoms',[lonid,latid,tid],/FLOAT)
    2:   anomvar = NCDF_VARDEF(wilma,'td_anoms',[lonid,latid,tid],/FLOAT)
    3:   anomvar = NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
    4:   anomvar = NCDF_VARDEF(wilma,'tw_anoms',[lonid,latid,tid],/FLOAT)
    5:   anomvar = NCDF_VARDEF(wilma,'rh_anoms',[lonid,latid,tid],/FLOAT)
    6:   anomvar = NCDF_VARDEF(wilma,'dpd_anoms',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   absvar = NCDF_VARDEF(wilma,'q_abs',[lonid,latid,tid],/FLOAT)
    1:   absvar = NCDF_VARDEF(wilma,'e_abs',[lonid,latid,tid],/FLOAT)
    2:   absvar = NCDF_VARDEF(wilma,'td_abs',[lonid,latid,tid],/FLOAT)
    3:   absvar = NCDF_VARDEF(wilma,'t_abs',[lonid,latid,tid],/FLOAT)
    4:   absvar = NCDF_VARDEF(wilma,'tw_abs',[lonid,latid,tid],/FLOAT)
    5:   absvar = NCDF_VARDEF(wilma,'rh_abs',[lonid,latid,tid],/FLOAT)
    6:   absvar = NCDF_VARDEF(wilma,'dpd_abs',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   climvar = NCDF_VARDEF(wilma,'q_clims',[lonid,latid,clmid],/FLOAT)
    1:   climvar = NCDF_VARDEF(wilma,'e_clims',[lonid,latid,clmid],/FLOAT)
    2:   climvar = NCDF_VARDEF(wilma,'td_clims',[lonid,latid,clmid],/FLOAT)
    3:   climvar = NCDF_VARDEF(wilma,'t_clims',[lonid,latid,clmid],/FLOAT)
    4:   climvar = NCDF_VARDEF(wilma,'tw_clims',[lonid,latid,clmid],/FLOAT)
    5:   climvar = NCDF_VARDEF(wilma,'rh_clims',[lonid,latid,clmid],/FLOAT)
    6:   climvar = NCDF_VARDEF(wilma,'dpd_clims',[lonid,latid,clmid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   stdvar = NCDF_VARDEF(wilma,'q_std',[lonid,latid,clmid],/FLOAT)
    1:   stdvar = NCDF_VARDEF(wilma,'e_std',[lonid,latid,clmid],/FLOAT)
    2:   stdvar = NCDF_VARDEF(wilma,'td_std',[lonid,latid,clmid],/FLOAT)
    3:   stdvar = NCDF_VARDEF(wilma,'t_std',[lonid,latid,clmid],/FLOAT)
    4:   stdvar = NCDF_VARDEF(wilma,'tw_std',[lonid,latid,clmid],/FLOAT)
    5:   stdvar = NCDF_VARDEF(wilma,'rh_std',[lonid,latid,clmid],/FLOAT)
    6:   stdvar = NCDF_VARDEF(wilma,'dpd_std',[lonid,latid,clmid],/FLOAT)
  ENDCASE
  
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','months beginning Jan 1973'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','valid_min',0.
  NCDF_ATTPUT,wilma,'latitude','long_name','Latitude'
  NCDF_ATTPUT,wilma,'latitude','units','Degrees'
  NCDF_ATTPUT,wilma,'latitude','valid_min',-90.
  NCDF_ATTPUT,wilma,'latitude','valid_max',90.
  NCDF_ATTPUT,wilma,'longitude','long_name','Longitude'
  NCDF_ATTPUT,wilma,'longitude','units','Degrees'
  NCDF_ATTPUT,wilma,'longitude','valid_min',-180.
  NCDF_ATTPUT,wilma,'longitude','valid_max',180.

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_anoms','long_name','Marine Specific Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_anoms','units','g/kg'
      NCDF_ATTPUT,wilma,'q_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'q_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_anoms','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_anoms','long_name','Marine Vapour Pressure monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_anoms','units','hPa'
      NCDF_ATTPUT,wilma,'e_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'e_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_anoms','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_anoms','long_name','Marine Dew Point Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'td_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'td_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_anoms','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_anoms','long_name','Marine Air Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'t_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'t_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_anoms','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_anoms','long_name','Marine Wet Bulb Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'tw_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'tw_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_anoms','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_anoms','long_name','Marine Relative Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_anoms','units','%rh'
      NCDF_ATTPUT,wilma,'rh_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'rh_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_anoms','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_anoms','long_name','Marine Dew Point Depression monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_anoms','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_abs','long_name','Marine Specific Humidity monthly mean'
      NCDF_ATTPUT,wilma,'q_abs','units','g/kg'
      NCDF_ATTPUT,wilma,'q_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'q_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_abs','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_abs','long_name','Marine Vapour Pressure monthly mean'
      NCDF_ATTPUT,wilma,'e_abs','units','hPa'
      NCDF_ATTPUT,wilma,'e_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'e_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_abs','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_abs','long_name','Marine Dew Point Temperature monthly mean'
      NCDF_ATTPUT,wilma,'td_abs','units','deg C'
      NCDF_ATTPUT,wilma,'td_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'td_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_abs','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_abs','long_name','Marine Air Temperature monthly mean'
      NCDF_ATTPUT,wilma,'t_abs','units','deg C'
      NCDF_ATTPUT,wilma,'t_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'t_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_abs','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_abs','long_name','Marine Wet Bulb Temperature monthly mean'
      NCDF_ATTPUT,wilma,'tw_abs','units','deg C'
      NCDF_ATTPUT,wilma,'tw_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'tw_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_abs','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_abs','long_name','Marine Relative Humidity monthly mean'
      NCDF_ATTPUT,wilma,'rh_abs','units','%rh'
      NCDF_ATTPUT,wilma,'rh_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'rh_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_abs','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_abs','long_name','Marine Dew Point Depression monthly mean'
      NCDF_ATTPUT,wilma,'dpd_abs','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'dpd_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_abs','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_clims','long_name','Marine Specific Humidity monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_clims','units','g/kg'
      NCDF_ATTPUT,wilma,'q_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'q_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_clims','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_clims','long_name','Marine Vapour Pressure monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_clims','units','hPa'
      NCDF_ATTPUT,wilma,'e_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'e_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_clims','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_clims','long_name','Marine Dew Point Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_clims','units','deg C'
      NCDF_ATTPUT,wilma,'td_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'td_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_clims','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_clims','long_name','Marine Air Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_clims','units','deg C'
      NCDF_ATTPUT,wilma,'t_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'t_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_clims','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_clims','long_name','Marine Wet Bulb Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_clims','units','deg C'
      NCDF_ATTPUT,wilma,'tw_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'tw_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_clims','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_clims','long_name','Marine Relative Humidity monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_clims','units','%rh'
      NCDF_ATTPUT,wilma,'rh_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'rh_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_clims','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_clims','long_name','Marine Dew Point Depression monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_clims','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'dpd_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_clims','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_std','long_name','Marine Specific Humidity monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_std','units','g/kg'
      NCDF_ATTPUT,wilma,'q_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'q_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_std','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_std','long_name','Marine Vapour Pressure monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_std','units','hPa'
      NCDF_ATTPUT,wilma,'e_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'e_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_std','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_std','long_name','Marine Dew Point Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_std','units','deg C'
      NCDF_ATTPUT,wilma,'td_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'td_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_std','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_std','long_name','Marine Air Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_std','units','deg C'
      NCDF_ATTPUT,wilma,'t_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'t_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_std','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_std','long_name','Marine Wet Bulb Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_std','units','deg C'
      NCDF_ATTPUT,wilma,'tw_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'tw_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_std','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_std','long_name','Marine Relative Humidity monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_std','units','%rh'
      NCDF_ATTPUT,wilma,'rh_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'rh_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_std','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_std','long_name','Marine Dew Point Depression monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_std','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'dpd_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_std','missing_value',mdi
    END

  ENDCASE

  current_time=SYSTIME()
  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"HadISDH monthly mean marine surface climate monitoring product from 1973 onwards. "+$
                                         "Quality control, bias correction, uncertainty estimation, averaging over gridboxes (no smoothing "+$
					 "or interpolation)."
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"HadISDH monthly mean marine surface climate monitoring product from 1973 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"Met Office Hadley Centre (UK), National Oceanography Centre Southampton (NOCS)"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"ICOADS.3.0.0 and ICOADS.3.0.1 (Freeman et al., 2016), http://icoads.noaa.gov/"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"Willett, K. M., Dunn, R. J. H., Parker, D. E., Kennedy, J. F. and Berry, D. I.: "+$
                                       "HadISDH marine surface multi-variable humidity and temperature "+$
				       "record for climate monitoring, in prep."
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"HadISDH.marine."+MarineVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF

  NCDF_VARPUT, wilma,timesvar, landtims
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,anomvar,marinearrANM
  NCDF_VARPUT, wilma,absvar,marinearrABS
  NCDF_VARPUT, wilma,climvar,marinearrCLM
  NCDF_VARPUT, wilma,stdvar,marinearrSDV

  NCDF_CLOSE,wilma


ENDFOR

stop
end
