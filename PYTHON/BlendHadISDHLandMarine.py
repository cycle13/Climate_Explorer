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
# modul load scitools/default-current
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
EdYr = 2019
StMn = 1 
EdMn = 12
platform = 'ship' # 'all'

if (platform == 'ship'):
    PT = 'SHIP'
else:
    PT = ''
    
NowMon = 'JAN'
NowYear = '2020'
LandMonYear = 'JAN2020'
MarineMonYear = 'JAN2020'
ClimStart = 1981
ClimEnd = 2010
RefPeriod = str(ClimStart)+' to '+str(ClimEnd)
ClimPeriod = str(ClimStart)[2:4]+str(ClimEnd)[2:4]
LandVersion = '4.2.0.2019f'
MarineVersion = '1.0.0.2019f'
BlendVersion = '1.0.0.2019f'

################################################################################################################
# SUBROUTINES #
##################################################
# BlendIt #
##################################################
def BlendIt(TheLand, TheOcean, PctLandArr, TheMDI, IsUnc):
    '''This programs works over each gridbox to blend the land and marine
        data. It uses HadCRUT.4.3.0.0.land_fraction.nc to get a land/ocean coverage
	but if there is both a land and marine box present then there will be a min/max blend of 0.25/0.75.
	
	Apply minimum blend of 25% if both land and marine are present.
	
	This does a cursory uncertainty blend in quadrature of the obs, 
	sampling and full uncertainty.
	INPUT:
	TheLand - time, lat, lon np.array of land data
	TheOcean - time, lat, lon np.array of marine data
	PctLandArr - lat, lon np.array of percent land (1) and sea (0)
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
        BlendedHadISDH[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))] = (PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]*TheLand[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))]) + ((1 - PctLandArr[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))])*TheOcean[np.where((TheLand > TheMDI) & (TheOcean > TheMDI))])

#    ## Test it
#    pdb.set_trace()

#    print('Test the Blending of BlendedHadISDH from TheLand, TheOcean and PctLandArr')
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
    OutBVarList = ['_anoms','_abs','_clims','_obserr','_abs_obserr','_samplingerr','_abs_samplingerr','_combinederr','_abs_combinederr']
    OutLVarList = ['Land_mean_n_stations','Land_actual_n_stations','_land_std']
    OutMVarList = ['Marine_n_grids','Marine_n_obs','Marine_mean_pseudo_stations','Marine_actual_pseudo_stations',
                   'Marine_clim_n_obs','Marine_clim_n_grids','Marine_climstd_n_obs','Marine_climstd_n_grids','_marine_clims_std']
    
    # Output variables list elements ????? have the var_loop_lower[v] prefixed
    OutBVarLong = ['monthly mean anomalies','monthly mean actuals','monthly mean climatologies','monthly mean anomaly observation uncertainty (2 sigma)','monthly mean actual observation uncertainty (2 sigma)',
                   'monthly mean anomaly sampling uncertainty','monthly mean actual sampling uncertainty','monthly mean anomaly full uncertainty','monthly mean actual full uncertainty']
    OutLVarLong = ['mean number of stations over land','actual number of stations over land','land monthly mean standard deviations']
    OutMVarLong = ['number of 1x1 daily grids with data over ocean','number of marine observations','mean number of pseudo stations over ocean','actual number of pseudo stations over ocean',
		   'number of marine observations in the climatology','number of 1x1 daily grids in the climatology over ocean','number of marine observations in the climatological standard deviation',
		   'number of 1x1 daily grids in the climatological standard deviation over ocean','marine monthly mean climatological standard deviation']
    
    OutBVarStandard = ['_anoms','_abs','_clims','_obserr','_abs_obserr','_samplingerr','_abs_samplingerr','_combinederr','_abs_combinederr']
    OutLVarStandard = ['Land_mean_n_stations','Land_actual_n_stations','_land_std']
    OutMVarStandard = ['Marine_n_grids','Marine_n_obs','Marine_mean_pseudo_stations','Marine_actual_pseudo_stations','Marine_clim_n_obs','Marine_clim_n_grids','Marine_climstd_n_obs','Marine_climstd_n_grids','_marine_clims_std']

    units_loop = ['degrees C','degrees C','g/kg','hPa','%rh','degrees C','degrees C','standard']
        
    # Missing Data Indicator
    MDI = -1e30
    
    # Input and Output directory:
    WorkingDir = '/data/users/hadkw/WORKING_HADISDH/UPDATE'+str(EdYr)
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
    Filee = WorkingDir+OtherDataDir+'HadCRUT.4.3.0.0.land_fraction.nc' 
    LatInfo = ['latitude']
    LonInfo = ['longitude']
    PctLand, LatList, LonList = GetGrid(Filee,['land_area_fraction'],LatInfo,LonInfo)

    # Make PctLand a min and max of 0.25 or 0.75 respectively
    # We only use it when both land and marine are present so it doesn't matter in land only or marine only cases    
    PctLand[np.where(PctLand < 0.25)] = 0.25
    PctLand[np.where(PctLand > 0.75)] = 0.75
    # Now make PctLand have a time dimension the same length as the data
    # PctLand is alread np.shape(1,26,72) so no need for np.newaxis in the first dimension
    FullPctLand = np.repeat(PctLand[:,:,:],Ntims,axis=0)
    ##Test It
    #pdb.set_trace()
    
#########
    # Loop through each variable to create the blended file
        
    # Which variable to loop through?
    for v,var in enumerate(var_loop_lower):
        
	# I haven't finished the land Td yet so skip this
	
        if (var == 'td'):
            continue
	    
        print('Working on: ',var,v)

        # Loop through quantities to be blended and build list
        BlendedList = []
        LandList = []
        LandDataList = []
        MarineList = []
        MarineDataList = []

        # For this var read in all land elements
        Filee = InFilLandBit1+var_loop[v]+InFilLandBit2+homogtype[v]+InFilLandBit3
        LatInfo = ['latitude']
        LonInfo = ['longitude']
        TmpVarList = []
#        pdb.set_trace()        
        TmpVarList = InLVarList.copy()
        TmpVarList[InLPointers[0]:InLPointers[1]] = [var+i for i in TmpVarList[InLPointers[0]:InLPointers[1]]]
        TmpDataList, LatList, LonList = GetGrid4(Filee,TmpVarList,LatInfo,LonInfo)
#        pdb.set_trace()
        # This comes out as:
        # TmpDataList: a list of np arrays (times, lats{87.5N to 87.5S], lons[-177.5W to 177.5E])
        # LatList: an NLats np array of lats centres (87.5N to 87.5S)
        # LonList: an NLons np array of lons centres (87.5N to 87.5S)
    
        # Get lat and lon counts
        NLats = len(LatList)
#        print(LatList)
        NLons = len(LonList)
	
	# change MDI for counts to -1
        TmpDataList[0][np.where(TmpDataList[0] <= 0)] = -1
        TmpDataList[1][np.where(TmpDataList[1] <= 0)] = -1
    
        # Fill in lists
        LandList.append(np.copy(TmpDataList[0])) # mean_n_stations 	
        LandList.append(np.copy(TmpDataList[1])) # actual_n_stations
        LandList.append(np.copy(TmpDataList[5])) # _std
        LandDataList.append(np.copy(TmpDataList[2])) # _anoms	
        LandDataList.append(np.copy(TmpDataList[3])) # _abs	
        LandDataList.append(np.copy(TmpDataList[4])) # _clims	
        LandDataList.append(np.copy(TmpDataList[6])) # _obserr	
        LandDataList.append(np.copy(TmpDataList[6])) # _abs_obserr	
        LandDataList.append(np.copy(TmpDataList[7])) # _samplingerr	
        LandDataList.append(np.copy(TmpDataList[7])) # _abs_samplingerr	
        LandDataList.append(np.copy(TmpDataList[8])) # _combinederr	
        LandDataList.append(np.copy(TmpDataList[8])) # _abs_combinederr	
	
        # For this var read in all marine elements
        Filee = InFilMarineBit1+var_loop[v]+InFilMarineBit2
        LatInfo = ['latitude']
        LonInfo = ['longitude']
        TmpVarList = []
        TmpVarList = InMVarList.copy()
#        pdb.set_trace()
        TmpVarList[InMPointers[0]:InMPointers[1]] = [var+i for i in TmpVarList[InMPointers[0]:InMPointers[1]]]
#        print(TmpVarList)
        TmpDataList, LatList, LonList = GetGrid4(Filee,TmpVarList,LatInfo,LonInfo)
#        pdb.set_trace()
        # This comes out as:
        # TmpDataList: a list of np arrays (times, lats{87.5N to 87.5S], lons[-177.5W to 177.5E])
        # LatList: an NLats np array of lats centres (87.5N to 87.5S)
        # LonList: an NLons np array of lons centres (87.5N to 87.5S)
        
        # Fill in lists

# NOT SURE THIS FLIPPING IS THE CASE ANYMORE
#	# NOTE: ALL MARINE ARRAYS NEED TO HAVE THEIR LATITUDES FLIPPED!!!
#        LatList = LatList[::-1]
#        pdb.set_trace()
#        TmpDataList[0] = TmpDataList[0][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[0])) # _n_grids 	
#        TmpDataList[1] = TmpDataList[1][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[1])) # _n_obs
#        TmpDataList[7] = TmpDataList[7][::-1,:]
        MarineList.append(np.copy(TmpDataList[7])) # _mean_pseudo_stations
#        TmpDataList[6] = TmpDataList[6][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[6])) # _actual_pseudo_stations
#        TmpDataList[2] = TmpDataList[2][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[2])) # _clim_n_grids
#        TmpDataList[3] = TmpDataList[3][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[3])) # _clim_n_obs
#        TmpDataList[4] = TmpDataList[4][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[4])) # _climstd_n_grids
#        TmpDataList[5] = TmpDataList[5][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[5])) # _climstd_n_obs
#        TmpDataList[11] = TmpDataList[11][:,::-1,:]
        MarineList.append(np.copy(TmpDataList[11])) # _climstd
#        TmpDataList[8] = TmpDataList[8][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[8])) # _anoms	
#        TmpDataList[9] = TmpDataList[9][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[9])) # _abs	
#        TmpDataList[10] = TmpDataList[10][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[10])) # _clims	
#        TmpDataList[12] = TmpDataList[12][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[12])) # _obserr	
#        TmpDataList[13] = TmpDataList[13][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[13])) # _abs_obserr	
#        TmpDataList[14] = TmpDataList[14][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[14])) # _samplingerr	
#        TmpDataList[15] = TmpDataList[15][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[15])) # _abs_samplingerr	
#        TmpDataList[16] = TmpDataList[16][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[16])) # _combinederr	
#        TmpDataList[17] = TmpDataList[17][:,::-1,:]
        MarineDataList.append(np.copy(TmpDataList[17])) # _abs_combinederr	
	
	# Now loop through the Blended quantities and blend
        for b,bvar in enumerate(OutBVarList):

            print('Working on: ',bvar)
	    
	    # Is this an uncertainty quantity?
            if (b > 3): # Set IsUnc = True
                IsUnc = True
            else:
                IsUnc = False
		
            BlendedField = BlendIt(LandDataList[b],MarineDataList[b],FullPctLand,MDI,IsUnc)	    
            BlendedList.append(BlendedField)
	    
##############
        # Write out combined file - uncertainties are all 2 sigma!!!.

        Write_NetCDF_All(OutFilBit1+var_loop[v]+OutFilBit2+homogtype[v]+OutFilBit3,
                         [OutLVarList[0],OutLVarList[1],var+OutLVarList[2]],[OutMVarList[0],OutMVarList[1],OutMVarList[2],OutMVarList[3],OutMVarList[4],OutMVarList[5],OutMVarList[6],OutMVarList[7],var+OutMVarList[8]],[var+i for i in OutBVarList],
			 LandList,MarineList,BlendedList,
			 [OutLVarLong[0],OutLVarLong[1],var_loop_full[v]+OutLVarLong[2]],[OutMVarLong[0],OutMVarLong[1],OutMVarLong[2],OutMVarLong[3],OutMVarLong[4],OutMVarLong[5],OutMVarLong[6],OutMVarLong[7],var_loop_full[v]+OutMVarLong[8]],[var_loop_full[v]+' '+i for i in OutBVarLong],
			 [OutLVarStandard[0],OutLVarStandard[1],var_loop_full[v]+OutLVarStandard[2]],[OutMVarStandard[0],OutMVarStandard[1],OutMVarStandard[2],OutMVarStandard[3],OutMVarStandard[4],OutMVarStandard[5],OutMVarStandard[6],OutMVarStandard[7],var_loop_full[v]+OutMVarStandard[8]],[var_loop_full[v]+' '+i for i in OutBVarStandard],
                         LatList,LonList,StYr,EdYr,RefPeriod,units_loop[v],MDI)

#########

    print('And we are done!')

if __name__ == '__main__':
    
    main()
